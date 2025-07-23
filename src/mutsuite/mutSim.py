###
# Modifies a certain proportion of reads at a specified location
# Only reads with a string of 'M' (matches) running across the modification site are considered for simulation. Other reads are not considered for mutations
#quality scores are left as is
# some notes:
# 1) Reads that have the deletion too close to one of the ends are not modified.
# 2) Nucleotides are added to the end of the sequence. If the original sequence had soft clipping a the end, the added nucleotides are added as clipped
# 3) All unmodified reads must be paired - otherwise they weren't included in bedtools bam2fastq
####
import argparse
import copy
import os
import pysam
import random
import re
import subprocess
import logging

from collections import defaultdict
from mutSimUtils import replaceRead, cleanIndelsFromReadEnds

logger = logging.getLogger(__name__)


#print("SETTING SEED")
#random.seed(123)

#cigar variables
MATCH = 0
INSERTION = 1
DELETION = 2
CLIPPING = 4

class ReadGroupInfo:
    def __init__(self, read_group_id, read_group_sample, read_group_platform, read_group_library):
        self.read_group_id = read_group_id
        self.read_group_sample = read_group_sample
        self.read_group_platform = read_group_platform
        self.read_group_library = read_group_library

    def __str__(self):
        return "ID: " + self.read_group_id + ", SM: " + self.read_group_sample + ", PL: " + self.read_group_platform + ", LB: " + self.read_group_library

    def __repr__(self):
        return self.__str__()

    def to_dict(self):
        return { 'ID':self.read_group_id,
                 'SM':self.read_group_sample,
                 'PL':self.read_group_platform,
                 'LB':self.read_group_library
                 }

cigar_explode_pattern = re.compile("(\d+)(\w)")
def explodeCigar(cigar_string):
    exploded_cigar = ""
    for (cigar_count,cigar_char) in re.findall(cigar_explode_pattern,cigar_string):
        exploded_cigar += cigar_char * int(cigar_count)
    return exploded_cigar

cigar_unexplode_pattern = re.compile(r'((\w)\2{0,})')
def unexplodeCigar(explode_cigar_string):
    cigar = ""
    for (cigar_str,cigar_char) in re.findall(cigar_unexplode_pattern,explode_cigar_string):
        cigar += str(len(cigar_str)) + cigar_char
    return cigar

def read_command_output(command):
    """
    Runs a shell command and returns an iter to read the output

    Args:
        command: shell command to run

    Returns:
        iter to read the output
    """

    p = subprocess.Popen(command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,shell=True,
#            encoding='utf-8',universal_newlines=True)
            universal_newlines=True,
            bufsize=-1) #bufsize system default
    return iter(p.stdout.readline, b'')


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--settings_file", help="Settings file containing [setting_name][tab][setting_value]", default=None)
    argparser.add_argument("--downsample_number", help="Depth after downsampling", type=int, default=None)

    argparser.add_argument("--mut_frequency", help="Frequency of simulated mutated reads from alternate into sample", type=float, default=None)
    argparser.add_argument("--mut_count", help="Count of simulated mutated reads from alternate into sample. If set, mut_frequency is ignored", type=float, default=None)

    argparser.add_argument("--qual_add", help="This number is added to all quality scores", type=int, default=0)
    argparser.add_argument("--mut_chr", help="Mutation chromosome", required=True)
    argparser.add_argument("--mut_loc", help="Mutation location", type=int, required=True)
    argparser.add_argument("--reference", help="Reference genome location", required=True)
    argparser.add_argument("--output_root", help="Output file root", required=True)
    argparser.add_argument("--unaltered_bam", help="Source aligned bam file with unaltered reads", required=True)
    argparser.add_argument("--unaltered_namesorted_bam", help="Source aligned bam file with unaltered reads, name sorted using samtools sort -n", default=None)
    argparser.add_argument("--altered_bam", help="Aligned bam file with altered reads to be inserted into simulated file. If none, reads will be simulated from the reference genome around the simulation location.", default=None, required=False)
    argparser.add_argument("--padding_region_size", help="Size of region around the mut_loc to extract reads from for analysis.", default=1000)
    argparser.add_argument("--mut_overlap_buffer", help="Reads starting or ending within this number of bp of the mut_loc will not be modified", default=5)
    argparser.add_argument("--only_include_altered_with_indel", help="Of reads from the altered file, only those with an indel are considered for simulating mutations",action='store_true')
    argparser.add_argument("--only_include_altered_crispresso_modified_indel", help="Of reads from the altered file, only those with an indel at the cut window as measured by CRISPResso are considered for simulating mutations", action='store_true')
    argparser.add_argument("--discard_indels_on_read_edge", help="Discard reads with indels at the edge of the altered reads", action='store_true')
    argparser.add_argument("--only_insert_first_indel", help="Only insert the indel from the first altered_bam read. This will produce samples that have two alleles and are easier for other indel callers to identify.", action='store_true')
    argparser.add_argument("--only_insert_specific_mutation", help="Only insert a specific mutation, specified by this parameter like 'SA' (substitute A), 'D1' (Delete 1bp) or 'IACCT' (Insert ACCT)", default=None)
    argparser.add_argument("--clearTags", help="If set, all sam tags are cleared from alignments.", action='store_true', default=False)
    argparser.add_argument("--keep_intermediate_files", help="If set, intermediate files are kept.", action='store_true', default=False)
    argparser.add_argument("--print_unmodified_read_info", help="If set, print information about unmodified reads - these are likely reads that were not modified because the alteredRead had no changes", action='store_true')
    argparser.add_argument("--new_read_group_id", help="Read group ID for simulated reads", default="simulated")
    argparser.add_argument("--new_read_group_sample", help="Read group sample ID for simulated reads", default="simulated")
    argparser.add_argument("--new_read_group_platform", help="Read group platform ID for simulated reads", default="ILLUMINA")
    argparser.add_argument("--new_read_group_library", help="Read group library ID for simulated reads", default="simulated")
    argparser.add_argument("--new_ctl_read_group_id", help="Read group ID for simulated reads in control sample", default="simulatedCtl")
    argparser.add_argument("--new_ctl_read_group_sample", help="Read group sample ID for simulated reads in control sample", default="simulatedCtl")
    argparser.add_argument("--new_ctl_read_group_platform", help="Read group platform ID for simulated reads in control sample", default="ILLUMINA")
    argparser.add_argument("--new_ctl_read_group_library", help="Read group library ID for simulated reads in control sample", default="simulatedCtl")
    argparser.add_argument("--debug", help="Print debug messages", action='store_true', default=False)
    
    args = argparser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # read in settings file if provided - note that these override command line arguments (and their defaults)
    if args.settings_file is not None:
        if not os.path.exists(args.settings_file):
            raise Exception("Settings file does not exist at " + args.settings_file)
        with open(args.settings_file) as settings_file:
            for line in settings_file:
                setting_name, setting_value = line.strip().split('\t')
                if hasattr(args, setting_name):
                    setattr(args, setting_name, setting_value)


    simulated_read_group_info = ReadGroupInfo(args.new_read_group_id, args.new_read_group_sample, args.new_read_group_platform, args.new_read_group_library)
    simulated_control_read_group_info = ReadGroupInfo(args.new_ctl_read_group_id, args.new_ctl_read_group_sample, args.new_ctl_read_group_platform, args.new_ctl_read_group_library)

    if args.mut_frequency is None and args.mut_count is None:
        raise Exception("--mut_frequency or --mut_count must be provided to specify the number of mutated reads to simulate")
    if args.mut_frequency is not None and args.mut_count is not None:
        raise Exception("--mut_frequency and --mut_count cannot both be provided")

    if args.altered_bam is None and args.only_insert_specific_mutation is None:
        raise Exception("either --altered_bam or --only_insert_specific_mutation must be specified")

    if not os.path.exists(args.reference):
        raise Exception("Reference file does not exist")

    if not os.path.exists(args.unaltered_bam):
        raise Exception("Unaltered bam file " + args.unaltered_bam + " does not exist")
    check_bam_indexed(args.unaltered_bam)

    input_bam_is_big = False
    if os.path.getsize(args.unaltered_bam) > 1000000000:  # 1GB
        input_bam_is_big = True
        logger.info("Input bam file is big - using big bam simulation method")

    if input_bam_is_big:
        sim_bam_muts_into_big_bam(args.unaltered_bam, args.altered_bam, args.reference, args.mut_chr, args.mut_loc, args.output_root, 
                                  padding_region_size=args.padding_region_size, mut_overlap_buffer=args.mut_overlap_buffer, downsample_number=args.downsample_number,
                                  mut_frequency=args.mut_frequency, mut_count=args.mut_count, qual_add=args.qual_add, only_include_altered_with_indel=args.only_include_altered_with_indel,
                                    only_include_altered_crispresso_modified_indel=args.only_include_altered_crispresso_modified_indel, discard_indels_on_read_edge=args.discard_indels_on_read_edge,
                                    only_insert_first_indel=args.only_insert_first_indel, only_insert_specific_mutation=args.only_insert_specific_mutation, clearTags=args.clearTags,
                                    keep_intermediate_files=args.keep_intermediate_files, print_unmodified_read_info=args.print_unmodified_read_info, simulated_read_group_info=simulated_read_group_info,
                                    simulated_control_read_group_info=simulated_control_read_group_info)
    else:
        sim_bam_muts_into_small_bam(args.unaltered_namesorted_bam, args.unaltered_bam, args.altered_bam, args.reference, args.mut_chr, args.mut_loc, args.output_root,
                                    padding_region_size=args.padding_region_size, mut_overlap_buffer=args.mut_overlap_buffer, downsample_number=args.downsample_number,
                                    mut_frequency=args.mut_frequency, mut_count=args.mut_count, qual_add=args.qual_add, only_include_altered_with_indel=args.only_include_altered_with_indel,
                                    only_include_altered_crispresso_modified_indel=args.only_include_altered_crispresso_modified_indel, discard_indels_on_read_edge=args.discard_indels_on_read_edge,
                                    only_insert_first_indel=args.only_insert_first_indel, only_insert_specific_mutation=args.only_insert_specific_mutation, clearTags=args.clearTags,
                                    keep_intermediate_files=args.keep_intermediate_files, print_unmodified_read_info=args.print_unmodified_read_info, simulated_read_group_info=simulated_read_group_info,
                                    simulated_control_read_group_info=simulated_control_read_group_info)

def sim_bam_muts_into_big_bam(unaltered_bam, altered_bam, reference, mut_chr, mut_loc, output_root, padding_region_size=1000, mut_overlap_buffer=5, downsample_number=None,
                            mut_frequency=None, mut_count=None, qual_add=0, only_include_altered_with_indel=False, only_include_altered_crispresso_modified_indel=False, 
                            discard_indels_on_read_edge=False, only_insert_first_indel=False, only_insert_specific_mutation=None, clearTags=False, keep_intermediate_files=False, print_unmodified_read_info=False, 
                            simulated_read_group_info=None, simulated_control_read_group_info=None):

    unaffected_reads_sam = output_root + ".unaffected_reads.sam"
    potential_reads_sam = output_root + ".potential_reads.sam"
    # first create bam for other chromosomes
    logger.debug('Extracting reads from unaltered bam ' + unaltered_bam + ' around the cut site ' + mut_chr + ':' + str(mut_loc))
    logger.debug('Extracting reads from other chromosomes to unaffected file ' + unaffected_reads_sam)

    unaltered_bam_handle = pysam.AlignmentFile(unaltered_bam, "rb")
    chroms = unaltered_bam_handle.references

    # write header to output bam
    subprocess.check_output('samtools view -H ' + unaltered_bam + ' > ' + unaffected_reads_sam, shell=True)
    # write aligned reads to output bam
    subprocess.check_output('samtools view -F 2 ' + unaltered_bam + ' ' + mut_chr + ' >> ' + unaffected_reads_sam, shell=True)
    # for other chrs just write them to this other bam
    for chrom in chroms:
        if chrom != mut_chr:
            subprocess.check_output('samtools view -f 2 ' + unaltered_bam + ' ' + chrom + ' >> ' + unaffected_reads_sam, shell=True)

    #here, we only include reads that start (or whose mates start) within the padding_region_size of the mutation site. All other reads are assumed to be unaffected.
    if padding_region_size < 300:
        logger.warning('Padding region size is less than 300 bp. This may result in reads that end near the mutation position being discarded.')

    logger.debug('Iterating through reads on chromosome ' + mut_chr + ' to separate reads that could be potentially affected by the simulation')
    total_chr_read_count = 0
    total_chr_printed_count = 0
    region_start = mut_loc - padding_region_size
    region_end = mut_loc + padding_region_size
    subprocess.check_output('samtools view -H ' + unaltered_bam + ' > ' + potential_reads_sam, shell=True)
    with open(unaffected_reads_sam, 'a') as unaffected_fout, open(potential_reads_sam, 'a') as potential_fout:
        for line in read_command_output('samtools view -f 2 ' + unaltered_bam + ' ' + mut_chr):
            if line == '':
                break
            total_chr_read_count += 1
            mate_is_in_region = False
            line_els = line.strip().split('\t')
            sam_flag = line_els[1]
            read_mate_is_unmapped = int(sam_flag) & 0x8
            if not read_mate_is_unmapped:
                read_mate_reference_name = line_els[6]
                if read_mate_reference_name == '=':
                    read_mate_reference_name = line_els[2]
                read_mate_aln_loc = int(line_els[7])
                if read_mate_reference_name == mut_chr:
                    if region_start <= read_mate_aln_loc <= region_end:
                        mate_is_in_region = True
            aln_start = int(line_els[3])
            if region_start <= aln_start <= region_end or mate_is_in_region: ## only check this read (and its mate) start positions - some reads may have hard-to-calculate read ends
                potential_fout.write(line)
                total_chr_printed_count += 1
            else:
                unaffected_fout.write(line)

    unaltered_bam_handle.close()

    logger.debug(f'Printed {total_chr_printed_count}/{total_chr_read_count} reads that overlap target on chromosome {mut_chr} to potential_reads_sam')
    
    #convert to bams
    logger.debug('Converting sam files to bam files')
    unaffected_reads_bam = output_root + ".unaffected_reads.bam"
    potential_reads_bam = output_root + ".potential_reads.bam"
    convert_command = f'samtools sort -o {unaffected_reads_bam} {unaffected_reads_sam} && samtools index {unaffected_reads_bam}'
    logger.debug('Running command: ' + convert_command)
    subprocess.check_output(convert_command, shell=True)
    os.remove(unaffected_reads_sam)

    convert_command = f'samtools sort -o {potential_reads_bam} {potential_reads_sam} && samtools index {potential_reads_bam}'
    logger.debug('Running command: ' + convert_command)
    subprocess.check_output(convert_command, shell=True)
    os.remove(potential_reads_sam)



    logger.debug('Simulating mutations into small bam file ' + potential_reads_bam)
    sub_output_root = output_root + ".sub"
    sim_bam_muts_into_small_bam(unaltered_namesorted_bam=None, unaltered_bam=potential_reads_bam, altered_bam=altered_bam, reference=reference, mut_chr=mut_chr, mut_loc=mut_loc, output_root=sub_output_root, 
                                padding_region_size=padding_region_size, mut_overlap_buffer=mut_overlap_buffer, downsample_number=downsample_number,
                                mut_frequency=mut_frequency, mut_count=mut_count, qual_add=qual_add, only_include_altered_with_indel=only_include_altered_with_indel, only_include_altered_crispresso_modified_indel=only_include_altered_crispresso_modified_indel,
                                discard_indels_on_read_edge=discard_indels_on_read_edge, only_insert_first_indel=only_insert_first_indel, only_insert_specific_mutation=only_insert_specific_mutation, clearTags=clearTags, keep_intermediate_files=keep_intermediate_files, print_unmodified_read_info=print_unmodified_read_info,
                                simulated_read_group_info=simulated_read_group_info, simulated_control_read_group_info=simulated_control_read_group_info)
    

    #downsample unmodified bam if necessary
    if downsample_number is not None:
        try:
            sourceAln = pysam.AlignmentFile(unaltered_bam,"rb")
        except Exception as e:
            raise Exception("Could not open unaltered bam file " + unaltered_bam + ": " + str(e))
        depth_at_target = 0
        for read in sourceAln.fetch(mut_chr,mut_loc,mut_loc+1):
            depth_at_target += 1
        sourceAln.close()
        downsamplePct = min(1,(float(downsample_number)/float(depth_at_target)))
        logger.debug('Downsampling unmodified reads to ' + str(downsample_number) + ' reads (' + str(downsamplePct) + '%)')
        downsampled_unaffected_reads_bam = output_root + ".unaffected_reads.downsampled.bam"
        downsample_command = 'samtools view -s ' + str(downsamplePct) + ' -b ' + unaffected_reads_bam + ' > ' + downsampled_unaffected_reads_bam
        logger.debug('Running command: ' + downsample_command)
        subprocess.check_output(downsample_command, shell=True)
        os.remove(unaffected_reads_bam)
        unaffected_reads_bam = downsampled_unaffected_reads_bam

    # merge into final bam
    logger.debug('Merging bams into final output bam')
    merge_command = f'samtools merge -f {output_root}.bam {unaffected_reads_bam} {sub_output_root}.bam && samtools index {output_root}.bam'
    logger.debug('Running command: ' + merge_command)
    subprocess.check_output(merge_command, shell=True)

    logger.debug('Merging bams into final control bam')
    merge_command = f'samtools merge -f {output_root}.ctl.bam {unaffected_reads_bam} {sub_output_root}.ctl.bam && samtools index {output_root}.ctl.bam'
    logger.debug('Running command: ' + merge_command)
    subprocess.check_output(merge_command, shell=True)

    copy_command = f'cp {sub_output_root}.mutations.txt {output_root}.mutations.txt'
    logger.debug('Running command: ' + copy_command)
    subprocess.check_output(copy_command, shell=True)

    copy_command = f'cp {sub_output_root}.mutatedReads.txt {output_root}.mutatedReads.txt'
    logger.debug('Running command: ' + copy_command)
    subprocess.check_output(copy_command, shell=True)

    copy_command = f'cp {sub_output_root}.info {output_root}.info'
    logger.debug('Running command: ' + copy_command)
    subprocess.check_output(copy_command, shell=True)

    if not keep_intermediate_files:
        logger.debug('Removing intermediate files')
        files = [unaffected_reads_bam, sub_output_root + ".bam", sub_output_root + ".bam.bai", sub_output_root + ".ctl.bam", sub_output_root + ".ctl.bam.bai",
                 sub_output_root + ".ctl.unsorted.bam", sub_output_root + ".unaltered_namesorted.bam", sub_output_root + ".unaltered_namesorted.bam.bai",
                 output_root + ".potential_reads.bam", output_root + ".potential_reads.bam.bai", output_root + ".unaffected_reads.bam", output_root + ".unaffected_reads.bam.bai",
                 sub_output_root + ".mutatedReads.txt", sub_output_root + ".mutations.txt", sub_output_root + ".info"]
        for f in files:
            if os.path.exists(f):
                os.remove(f)


def sim_bam_muts_into_small_bam(unaltered_namesorted_bam, unaltered_bam, altered_bam, reference, mut_chr, mut_loc, output_root, padding_region_size=1000, mut_overlap_buffer=5, downsample_number=None,
                            mut_frequency=None, mut_count=None, qual_add=0, only_include_altered_with_indel=False, only_include_altered_crispresso_modified_indel=False, 
                            discard_indels_on_read_edge=False, only_insert_first_indel=False, only_insert_specific_mutation=None, clearTags=False, keep_intermediate_files=False, print_unmodified_read_info=False, 
                            simulated_read_group_info=None, simulated_control_read_group_info=None):
    """
    Simulate mutations into a small bam file (e.g. not a genome-wide bam file)
    """

    if unaltered_namesorted_bam is None:
        new_unaltered_namesorted_bam = output_root + ".unaltered_namesorted.bam"
        if not os.path.exists(new_unaltered_namesorted_bam):
            logger.warning('Unaltered namesorted bam does not exist. Creating it now at ' + new_unaltered_namesorted_bam)
            subprocess.check_output('samtools sort -n -o ' + new_unaltered_namesorted_bam + ' ' + unaltered_bam, shell=True)
        unaltered_namesorted_bam = new_unaltered_namesorted_bam

    if not os.path.exists(unaltered_namesorted_bam):
        raise Exception("Unaltered namesorted bam file " + unaltered_namesorted_bam + " does not exist")

    downsample_number_unset = False
    if downsample_number is None:
        downsample_number_unset = True
        try:
            sourceAln = pysam.AlignmentFile(unaltered_bam,"rb")
        except Exception as e:
            raise Exception("Could not open unaltered bam file " + unaltered_bam + ": " + str(e))
        this_depth = 0
        for read in sourceAln.fetch(mut_chr,mut_loc,mut_loc+1):
            this_depth += 1
        sourceAln.close()
        print('Downsample number not provided. Setting to ' + str(this_depth))
        downsample_number = this_depth

    #get reference sequence
    refStart = mut_loc - padding_region_size
    refEnd = mut_loc + padding_region_size
    refFile = pysam.Fastafile(reference)
    refSeq = refFile.fetch(mut_chr,refStart,refEnd)

    try:
        sourceAln = pysam.AlignmentFile(unaltered_bam,"rb")
    except Exception as e:
        raise Exception("Could not open unaltered bam file " + unaltered_bam + ": " + str(e))

    sourceHeaderSim = sourceAln.header.to_dict()
    sourceHeaderCtl = sourceAln.header.to_dict()
    #make new read groups
    rg = { 'ID': simulated_read_group_info.read_group_id,
            'SM': simulated_read_group_info.read_group_sample,
            'PL': simulated_read_group_info.read_group_platform,
            'LB': simulated_read_group_info.read_group_library
            }
    sourceHeaderSim['RG'] = [rg]
    rg = { 'ID': simulated_control_read_group_info.read_group_id,
            'SM': simulated_control_read_group_info.read_group_sample,
            'PL': simulated_control_read_group_info.read_group_platform,
            'LB': simulated_control_read_group_info.read_group_library
            }
    sourceHeaderCtl['RG'] = [rg]

    unsortedOutName = "" # file to write final reads to (this will be sorted as a last step)
    mutationsFileName = "" # output text file name for list of mutations inserted
    unsortedOutName = output_root+".unsorted.sam"
    destAln = pysam.AlignmentFile(unsortedOutName,"w",header=sourceHeaderSim)
    unsortedCtlOutName = output_root+".ctl.unsorted.bam"
    destAlnCtl = pysam.AlignmentFile(unsortedCtlOutName,"wb",header=sourceHeaderCtl)

    mutatedReadsFileName = output_root+".mutatedReads.txt"
    mutationsFileName = output_root+".mutations.txt"

    #####
    #first, get reads that overlap with deletion site, make sure they are good quality, and then add them to a list
    allReadAtTargetCount = 0
    goodReadAtTargetCount = 0
    readsAtTarget = []
    readsAtTargetLookup = {} #dict for reads that overlap with target - these are the reads that we'll be precisely counting (e.g. only want 10 modified and 90 unmodified reads)
    readsAtTargetToDiscard = {} #dict for reads that overlap with target that we'll discard (e.g. reads that have indels at the edge of the read)
    logger.info("Getting WT reads that overlap with " + mut_chr + ":" + str(mut_loc))
    for read in sourceAln.fetch(mut_chr, mut_loc, mut_loc + 1):
        if read.query_name in readsAtTargetLookup:
            continue
        if read.query_length is None or read.query_alignment_start is None or read.query_alignment_end is None:
            continue
        #include soft-clipped bases in these rStart and rEnd locations
        rStart = read.reference_start - read.query_alignment_start
        rEnd = read.reference_end + (read.query_length - read.query_alignment_end)
        rSeq = read.query_sequence
        rQual = read.query_qualities

        allReadAtTargetCount += 1
        read_is_good = 0
        if read.is_unmapped:
            readsAtTargetToDiscard[read.query_name] = 1
            continue

        if rEnd < mut_loc + mut_overlap_buffer or rStart > mut_loc - mut_overlap_buffer:
            readsAtTargetToDiscard[read.query_name] = 1
            continue

        readInd = 0
    #	print str(read.cigar)
    #	print 'read starts' + str(rStart) + ' and ends '+ str(rEnd)
        for op,count in read.cigar:
            #if match, check for overlap with mutation site
            if op == MATCH:
                tStart = rStart + readInd
                tEnd = rStart + readInd + count
                # print("from " + str(tStart) + " to " + str(tEnd) + "are match")

                if rStart + readInd <= mut_loc and rStart + readInd + count > mut_loc:
                    read_is_good = 1
                readInd += count
            elif op == INSERTION:
    #			tStart = rStart + readInd
    #			tEnd = rStart + readInd + count
    #			print "from " + str(tStart) + " to " + str(tEnd) + "are insertion"
                readInd += count
            elif op == DELETION:
    #			tStart = rStart + readInd
    #			tEnd = rStart + readInd + count
    #			print "from " + str(tStart) + " to " + str(tEnd) + "are deletion"
                pass
            elif op == CLIPPING:
    #			tStart = rStart + readInd
    #			tEnd = rStart + readInd + count
    #			print "from " + str(tStart) + " to " + str(tEnd) + "are clipping"
                readInd += count
            else:
                raise Exception("got unrecognized op '"+str(op)+"' in cigar string " + str(read.cigar))

        if read_is_good:
            if read.query_name not in readsAtTargetLookup:  # only add read if mate is not already in list
                readsAtTarget.append(read)
                readsAtTargetLookup[read.query_name] = 1
        else:
            readsAtTargetToDiscard[read.query_name] = 1

    #check to make sure reads printed at on-target are paired (some variant callers may want the paired read to be in the window)
    #called with 'until_eof=True' because it doesn't have an index (name-sorted and can't be indexed)
    save = pysam.set_verbosity(0) # suppress message about not being able to find index
    try:
        sourceAlnNamesorted = pysam.AlignmentFile(unaltered_namesorted_bam,"rb")
    except Exception as e:
        raise Exception("Could not open unaltered namesorted bam file " + unaltered_namesorted_bam + ": " + str(e))
    pysam.set_verbosity(save)
    sourceAlnNamesortedIter = sourceAlnNamesorted.fetch(until_eof=True)
    readNameAppearanceCount = {} #keep track of how many times we see each read name
    for read in sourceAlnNamesortedIter:
        if read.query_name in readNameAppearanceCount:
            readNameAppearanceCount[read.query_name] += 1
        else:
            readNameAppearanceCount[read.query_name] = 1
    sourceAlnNamesorted.close()
    logger.info('got ' + str(len(readsAtTarget)) + ' reads at target')
    logger.info('got ' + str(len(readsAtTargetLookup)) + ' reads at target to precisely count')
    logger.info('got ' + str(len(readsAtTargetToDiscard)) + ' reads at target to discard')

    goodReadsAtTarget = []
    for read in readsAtTarget:
        if read.query_name in readNameAppearanceCount and readNameAppearanceCount[read.query_name] > 1 and read.query_name not in readsAtTargetToDiscard:
            goodReadAtTargetCount += 1
            goodReadsAtTarget.append(read)
        else:
            readsAtTargetToDiscard[read.query_name] = 1
    #done making sure reads are paired

    if goodReadAtTargetCount == 0:
        raise Exception('No good reads at target!')

    logger.info("read " + str(allReadAtTargetCount) + " reads, kept " + str(goodReadAtTargetCount) + " reads at target in unmodified bam")
    sourceAln.close()

    print('mut frequency: ' + str(mut_frequency))
    print('mut count: ' + str(mut_count))
    #set up mutation and downsample counts
    if mut_frequency is not None:
        numChanges = int(round(len(goodReadsAtTarget) * mut_frequency))
    elif mut_count is not None:
        numChanges = int(mut_count)

    numToPrint = downsample_number
    if downsample_number_unset:
        numToPrint = len(goodReadsAtTarget)

    if len(goodReadsAtTarget) < numToPrint:
        raise Exception('Not enough unmodified reads to print the requested sequencing depth. Need at least ' + str(numChanges) + ' unmodified reads, but only have ' + str(len(goodReadsAtTarget)) + ' unmodified reads.')

    ####
    #next, create altered_reads that will be inserted into the simulated file
    totR1Count = 0 #these are for counts from the altered_bam file
    alteredR1Count = 0
    totUnpairedCount = 0
    alteredUnpairedCount = 0
    totR2Count = 0
    alteredR2Count = 0
    #if we are only inserting a specific mutation, make a read with that mutation and copy it a bunch of times
    if only_insert_specific_mutation is not None:
        alteredReadsTmp = []
        new_read = copy.deepcopy(goodReadsAtTarget[0])
        if 'S' in only_insert_specific_mutation:
            sub_nuc = only_insert_specific_mutation[1:].upper()
            if sub_nuc not in ['A', 'C', 'G', 'T']:
                raise Exception('Cannot parse substitution nucleotide: ' + sub_nuc + ' from mutation string "' + only_insert_specific_mutation + '". Must be one of A, C, G, T. Please check the --only_insert_specific_mutation parameter.')
            old_base = refFile.fetch(mut_chr, mut_loc - 1, mut_loc)
            if old_base.upper() == sub_nuc:
                logger.warning('Substituting ubstitute base ' + old_base + ' with ' + sub_nuc + ' at position ' + str(mut_loc) + '. The base is already ' + sub_nuc + '. Please check the --only_insert_specific_mutation parameter.')
            sub_loc = mut_loc
            old_len = len(new_read.query_sequence)
            pre_bases = int(old_len/2) # how many bases pre-substitution
            new_start = sub_loc - pre_bases
            post_bases = int(old_len/2) - 1
            new_end = sub_loc + post_bases
            new_seq = refFile.fetch(mut_chr, new_start-1, sub_loc-1) + sub_nuc + refFile.fetch(mut_chr, sub_loc, new_end)
            tmp1 = refFile.fetch(mut_chr, new_start - 1, sub_loc-1)
            print('first part: ' + tmp1 + ' with length ' + str(len(tmp1)))
            tmp2 = refFile.fetch(mut_chr, sub_loc, new_end)
            print('second part: ' + tmp2 + ' with length ' + str(len(tmp2)))
            new_read.query_sequence = new_seq
            new_read.query_qualities = [30] * len(new_seq)
            new_read.reference_start = new_start - 1
            new_read.cigarstring = str(pre_bases + post_bases + 1) + "M"
        elif 'I' in only_insert_specific_mutation:
            insert_nucs = only_insert_specific_mutation[1:].upper()
            for nuc in insert_nucs:
                if nuc not in ['A', 'C', 'G', 'T']:
                    raise Exception('Cannot parse insertion nucleotide: ' + nuc + ' from mutation string "' + only_insert_specific_mutation + '". Must be one of A, C, G, T. Please check the --only_insert_specific_mutation parameter.')
            insert_loc = mut_loc
            old_len = len(new_read.query_sequence)
            pre_bases = int(old_len/2) # how many bases pre-insertion
            new_start = insert_loc - pre_bases
            post_bases = int(old_len/2) - len(insert_nucs) # how many bases post-insertion
            new_end = insert_loc + post_bases
            new_seq = refFile.fetch(mut_chr, new_start-1, insert_loc-1) + insert_nucs + refFile.fetch(mut_chr, insert_loc-1, new_end-1)
            new_read.query_sequence = new_seq
            new_read.query_qualities = [30] * len(new_seq)
            new_read.reference_start = new_start - 1
            new_read.cigarstring = str(pre_bases) + "M" + str(len(insert_nucs)) + "I" + str(post_bases) + "M"
        elif 'D' in only_insert_specific_mutation:
            del_len = only_insert_specific_mutation[1:]
            if not del_len.isdigit() or int(del_len) <= 0:
                raise Exception('Cannot parse deletion length: ' + del_len + ' from mutation string "' + only_insert_specific_mutation + '". Must be a number. Please check the --only_insert_specific_mutation parameter.')
            del_len = int(del_len)
            insert_loc = mut_loc
            old_len = len(new_read.query_sequence)  # create a new read with the same length as one in the alteredReads set
            pre_post_bases = int(old_len/2)
            new_start = insert_loc - pre_post_bases
            new_end = insert_loc + pre_post_bases + del_len
            new_seq = refFile.fetch(mut_chr, new_start-1, insert_loc-1) + refFile.fetch(mut_chr, insert_loc + del_len-1, new_end-1)
            new_read.query_sequence = new_seq
            new_read.query_qualities = [30] * len(new_seq)
            new_read.reference_start = new_start - 1
            new_read.cigarstring = str(pre_post_bases) + "M" + str(del_len) + "D" + str(pre_post_bases) + "M"
        else:
            raise Exception('Cannot parse indel string: ' + only_insert_specific_mutation + '. Must be of the form <location>I<insertedBP> or <location>D<deletionLen>')

        for i in range(numChanges):
            copy1 = copy.deepcopy(new_read)
            alteredReadsTmp.append(copy1)
        alteredReads = alteredReadsTmp
    #else, read reads from altered_bam that should be simulated into the output bam
    else:  
        alteredAln = pysam.AlignmentFile(altered_bam, "rb")
        alteredReads = []
        seenAlteredReads = {} #for synchronizing r1 and r2 reads
        for read in alteredAln.fetch(mut_chr, mut_loc, mut_loc + 1):
            if read.cigarstring is None:
                continue
            if not read.is_paired:
                totUnpairedCount += 1
            elif read.is_read1:
                totR1Count += 1
            else:
                totR2Count += 1
            if only_include_altered_crispresso_modified_indel is True:
                c2_value = read.get_tag("c2")
                m = re.search('MODS=D(\d+);I(\d+);S(\d+)', c2_value)
                if m is None:
                    raise Exception('Cannot find crispresso modifications for read ' + read.query_name + ' in tag c2. Please set --only_include_altered_crispresso_modified_indel to False.')
                deletion_count = int(m.group(1))
                insertion_count = int(m.group(2))
                indel_count = deletion_count + insertion_count
                if indel_count == 0:
                    continue
            if discard_indels_on_read_edge:
                read = cleanIndelsFromReadEnds(read)
            if only_include_altered_with_indel is True and not ("I" in read.cigarstring or "D" in read.cigarstring):
                continue
            if read.query_name not in seenAlteredReads: # only add read if mate is not already in list
                if not read.is_paired:
                    alteredReads.append(read)
                    alteredUnpairedCount += 1
                elif read.is_read1:
                    alteredReads.append(read)
                    alteredR1Count += 1
                else:
                    alteredReads.append(read)
                    alteredR2Count += 1
                seenAlteredReads[read.query_name] = True
        random.shuffle(alteredReads)
        alteredAln.close()

        filter_indel_str = ""
        if only_include_altered_crispresso_modified_indel is True:
            filter_indel_str = "that contain modified crispresso indels "
        elif only_include_altered_with_indel is True:
            filter_indel_str = "that contain indels "

        logger.info("keeping %d/%d modified Unpaired reads %sfor replacement"%(alteredUnpairedCount, totUnpairedCount, filter_indel_str))
        logger.info("keeping %d/%d modified R1 reads %sfor replacement"%(alteredR1Count, totR1Count, filter_indel_str))
        logger.info("keeping %d/%d modified R2 reads %sfor replacement"%(alteredR2Count, totR2Count, filter_indel_str))



        alteredCount = alteredR1Count + alteredR2Count + alteredUnpairedCount
        logger.info('Total altered count: ' + str(alteredCount))
        if alteredCount < numChanges:
            raise Exception('Not enough altered reads to make the requested number of changes. Need at least ' + str(numChanges) + ' altered reads, but only have ' + str(alteredCount) + ' altered reads.')

        #if we are only inserting a single indel, just use the first altered read and copy it a bunch of times
        if only_insert_first_indel:
            alteredReadsTmp = []
            first_altered_read = alteredReads[0]
            for i in range(numChanges):
                copy1 = copy.deepcopy(first_altered_read)
                alteredReadsTmp.append(copy1)
            alteredReads = alteredReadsTmp

    #####
    #next, randomize the list of good reads, and make edits to the first X reads
    readsToReplace = {} #dict holds reads that were changed
    readsAtTargetToPrint = {} #dict holds reads to print
    random.shuffle(goodReadsAtTarget)

    found_status_print_ids = {}
    found_status_mod_ids = {}
    logger.info("making changes in " + str(numChanges) + "/" + str(len(goodReadsAtTarget)) + " unmodified reads")


    currInd = 0
    currReplaceCount = 0
    currPrintCount = 0
    while currReplaceCount < numChanges:
    #	print 'replacement ' + str(i)
        if currInd >= len(goodReadsAtTarget):
            logger.info('Requested more modified reads than were read in the modified file.')
        read = goodReadsAtTarget[currInd]
        if read.query_name in readsAtTargetToDiscard or read.query_name in readsToReplace or read.query_name in readsAtTargetToPrint:
            currInd += 1
            continue
        readsToReplace[read.query_name] = read
        readsAtTargetToPrint[read.query_name] = 1
        found_status_mod_ids[read.query_name] = False
        found_status_print_ids[read.query_name] = False
        currReplaceCount += 1
        currPrintCount += 1
        currInd += 1

    while currPrintCount < numToPrint:
        if currInd >= len(goodReadsAtTarget):
            logger.info('Requested more modified reads than were read in the modified file.')
        read = goodReadsAtTarget[currInd]
        if read.query_name in readsAtTargetToDiscard or read.query_name in readsToReplace or read.query_name in readsAtTargetToPrint:
            currInd += 1
            continue
        readsAtTargetToPrint[read.query_name] = 1
        found_status_print_ids[read.query_name] = False
        currPrintCount += 1
        currInd += 1

    ####
    # finally, iterate over all the reads. Replace reads if they are in readToReplace
    # reset the qualities -- when the sequence is changed, the quality is also changed.
    readReadPairs = 0
    printedNotAtCutSiteReadPairs = 0
    printedAsIsReadPairs = 0
    printedChangedReadPairs = 0
    printedToControlReadPairs = 0

    # write mutations to file
    mutated_reads_file = open(mutatedReadsFileName, "w")

    def get_read_replacement_info_str(old_read, new_read, altered_read, replacement_info):
        """Creates a string with information about a modified read for printing to a file

        Args:
            old_read (pysam.AlignedSegment): wt unmodified read
            new_read (pysam.AlignedSegment): read containing source modificiation
            altered_read (pysam.AlignedSegment): wt read with modification applied
            replacement_info (list): list of strings describing the modifications
        """
        min_start = old_read.reference_start
        mut_chr = old_read.reference_name
        if new_read.reference_start < min_start:
            min_start = new_read.reference_start
        if altered_read.reference_start < min_start:
            min_start = altered_read.reference_start
        old_spaces = " "*(old_read.reference_start - min_start)
        new_spaces = " "*(new_read.reference_start - min_start)
        altered_spaces = " "*(altered_read.reference_start - min_start)
        this_ref_seq = refFile.fetch(mut_chr, min_start, min_start + 300)
        return(
            "ref: %10s %10s %20s %s"%(min_start,min_start+300," ",this_ref_seq) +
            "\nold: %10s %10s %20s%s %s %s"%(old_read.reference_start,old_read.reference_end,old_read.cigarstring, old_spaces,old_read.query_sequence,old_read.query_name) +
            "\nalt: %10s %10s %20s%s %s %s"%(altered_read.reference_start,altered_read.reference_end,altered_read.cigarstring, altered_spaces,altered_read.query_sequence,altered_read.query_name) +
            "\nnew: %10s %10s %20s%s %s %s"%(new_read.reference_start,new_read.reference_end,new_read.cigarstring, new_spaces,new_read.query_sequence,','.join(replacement_info)) +
            "\n\n")


    # we will be printing read pairs, so only print downsamplePct/2
    downsamplePct = (float(downsample_number)/float(allReadAtTargetCount))
    save = pysam.set_verbosity(0)  # suppress message about not being able to find index
    sourceAlnNamesorted = pysam.AlignmentFile(unaltered_namesorted_bam, "rb")
    pysam.set_verbosity(save)
    mutationDict = defaultdict(int)  # holds counts of all mutations
    sourceAlnNamesortedIter = sourceAlnNamesorted.fetch(until_eof=True)  # until_eof because it doesn't have an index (the bam is namesorted)
    read1 = next(sourceAlnNamesortedIter)
    readsDiscardedBecauseNoPairFound = 0
    readsUnmodifiedBecauseMutUnmodified = 0 # number of produced reads that are unmodified because the mutated read/altered read is unmodified
    for read2 in sourceAlnNamesortedIter:  # iterate through namesorted bam. R1 should be followed by R2, but this will allow us to skip unpaired reads
        if (read1.query_name != read2.query_name):  # if we found an unpaired read, skip it
            if read1.query_name in readsToReplace:  # well, first make sure it isn't in readsToReplace or readsAtTargetToPrint
                raise Exception('READ IS UNPAIRED BUT IN REPLACE: ' + str(read1))
            if read1.query_name in readsAtTargetToPrint:
                raise Exception('READ IS UNPAIRED BUT IN AT TARGET: ' + str(read1))
            read1 = read2
            readsDiscardedBecauseNoPairFound += 1
            # raise Exception("Namesorted bam does not contain pairs. Try increasing the size of the region used to create the unmodified bam.\nr1: " + read1.query_name + " is not " + read2.query_name)
            continue

        if read1.query_name in readsAtTargetToDiscard:
            read1 = next(sourceAlnNamesortedIter, None)
            continue

        readReadPairs += 1

        read1.query_qualities = [min(x+qual_add, 41) for x in read1.query_qualities]
        read2.query_qualities = [min(x+qual_add, 41) for x in read2.query_qualities]

        read1_ctl = copy.deepcopy(read1)
        read1_ctl.query_name = read1_ctl.query_name + "_ctl"
        if clearTags:
            read1_ctl.set_tags
            read1.set_tags([])
        read1_ctl.set_tag('RG', 'simulatedCtl')
        read1_ctl.set_tag('LB', 'simulatedCtl')
        read1.set_tag('RG', 'simulated')
        read1.set_tag('LB', 'simulated')

        read2_ctl = copy.deepcopy(read2)
        read2_ctl.query_name = read2_ctl.query_name + "_ctl"
        if clearTags:
            read2_ctl.set_tags
            read2.set_tags([])
        read2_ctl.set_tag('RG', 'simulatedCtl')
        read2_ctl.set_tag('LB', 'simulatedCtl')
        read2.set_tag('RG', 'simulated')
        read2.set_tag('LB', 'simulated')

        # if this read pair is at the target location,
        # read1.query_name is the same as read2.query_name
        if read1.query_name in readsAtTargetLookup:
            if read1.query_name in readsAtTargetToPrint:
                # write these reads to control
                destAlnCtl.write(read1_ctl)
                destAlnCtl.write(read2_ctl)
                printedToControlReadPairs += 1
                # replace and print (or just print below in else)
                this_read_mutation_dict = defaultdict(int)  # holds counts of mutations for this read pair (so we don't add mutations to the mutationDict twice)
                if read1.query_name in readsToReplace:
                    alteredRead = alteredReads.pop(0)
                    origRead1 = copy.deepcopy(read1)
                    origRead2 = copy.deepcopy(read2)
                    try:
                        newread1, replacement1Info, modification1Info = replaceRead(read1,alteredRead,refSeq,refStart, refEnd)
                        newread1.query_qualities = [min(x+qual_add,41) for x in newread1.query_qualities]
                        nonGenomicModArr1Info = []
                        for x in replacement1Info:
                            if x in modification1Info:
                                nonGenomicModArr1Info.append(x)
                    except Exception as e:
                        logger.critical("Error replacing read1: " + str(e))
                        logger.critical("Read1: " + str(read1))
                        logger.critical("Altered read: " + str(alteredRead))
                        logger.critical("Read2: " + str(read2))
                        raise e
                    try:
                        newread2, replacement2Info, modification2Info = replaceRead(read2,alteredRead,refSeq,refStart, refEnd)
                        newread2.query_qualities = [min(x+qual_add,41) for x in newread2.query_qualities]
                        nonGenomicModArr2Info = [] # create array of mutations introduced from alteredRead that are changes to genomic sequence
                        for x in replacement2Info:
                            if x in modification2Info:
                                nonGenomicModArr2Info.append(x)
                    except Exception as e:
                        logger.critical("Error replacing read2: " + str(e))
                        logger.critical("Read2: " + str(read2))
                        logger.critical("Altered read: " + str(alteredRead))
                        logger.critical("Read1: " + str(read1))
                        raise e
                    read1_is_mod = False
                    read2_is_mod = False
                    if read1.query_sequence != newread1.query_sequence:
                        for mut in nonGenomicModArr1Info:
                            this_read_mutation_dict[mut] += 1
                        mutated_reads_file.write(get_read_replacement_info_str(origRead1, newread1, alteredRead, nonGenomicModArr1Info))
                        read1_is_mod = True
                        newread1.set_tag('ZS', 'mod_simulated')
                        newread1.set_tag('ZR', ','.join(replacement1Info))
                        newread1.set_tag('ZI', ','.join(modification1Info))
                        newread2.set_tag('MC', newread1.cigarstring)  # set mate cigar to match
                    else:
                        newread1.set_tag('SS', 'mate_mod_simulated')
                    if read2.query_sequence != newread2.query_sequence:
                        for mut in nonGenomicModArr2Info:
                            this_read_mutation_dict[mut] += 1
                        mutated_reads_file.write(get_read_replacement_info_str(origRead2, newread2, alteredRead, nonGenomicModArr2Info))
                        read2_is_mod = True
                        newread2.set_tag('ZS', 'mod_simulated')
                        newread2.set_tag('ZR', ','.join(replacement2Info))
                        newread2.set_tag('ZI', ','.join(modification2Info))
                        newread1.set_tag('MC', newread2.cigarstring) # set mate cigar to match
                    else:
                        newread2.set_tag('SS', 'mate_mod_simulated')
                    if read1_is_mod or read2_is_mod:
                        printedChangedReadPairs += 1
                        for mut in this_read_mutation_dict:  # add mutations from this read pair to mutationDict
                            mutationDict[mut] += 1
                    if not read1_is_mod and not read2_is_mod:
                        readsUnmodifiedBecauseMutUnmodified += 1
                        if print_unmodified_read_info:
                            logger.warning('Error: Could not add mutation to read pair!')
                            logger.warning('Unmodified: ' + read1.query_name + ' ' + read2.query_name)
                            logger.warning('Pre:')
                            logger.warning(str(origRead1))
                            logger.warning(str(origRead2))
                            logger.warning('Altered:')
                            logger.warning(str(alteredRead))
                            logger.warning('Post:')
                            logger.warning(str(newread1))
                            logger.warning(str(newread2))
                    destAln.write(newread1)
                    destAln.write(newread2)
                else:
                    printedAsIsReadPairs += 1
                    destAln.write(read1)
                    destAln.write(read2)

        # if the read is not at the target location, print it out with downsampling
        else:
            if random.random() >= downsamplePct:
                read1 = next(sourceAlnNamesortedIter, None)
                continue
            else:
                destAln.write(read1)
                destAln.write(read2)
                printedNotAtCutSiteReadPairs += 1

                destAlnCtl.write(read1_ctl)
                destAlnCtl.write(read2_ctl)
                printedToControlReadPairs += 1



        # finally, read the next potential read pair (read 2 will be read at beginning of loop)
        read1 = next(sourceAlnNamesortedIter, None)

    destAln.close()
    destAlnCtl.close()
    sourceAlnNamesorted.close()

    logger.info('Discarding %s reads in unmodified sample with no pair found'%readsDiscardedBecauseNoPairFound)
    if readsUnmodifiedBecauseMutUnmodified > 0:
        logger.info('Note that ' + str(readsUnmodifiedBecauseMutUnmodified) + ' reads were unmodified because the altered/mutated read was unmodified (no mutations could be inserted). ' + \
            'To avoid this, ensure that the altered read has an indel at the target location or use the flag --only_include_altered_with_indel.')

    mutationsFile = open(mutationsFileName, "w")
    mutationsFile.write("LOC\tMUTATION\tINFO\tCOUNT\n")
    for mutationKey in mutationDict:
        if 'I' in mutationKey:
            mutation = 'I'
            mutation_els = mutationKey.split(mutation)
            loc = str(int(mutation_els[0]))
            info = mutation_els[1]
            mutationsFile.write(loc+"\t"+mutation+"\t"+info + "\t" + str(mutationDict[mutationKey]) + "\n")
        elif 'D' in mutationKey:
            mutation = 'D'
            mutation_els = mutationKey.split(mutation)
            loc = str(int(mutation_els[0]) - 1)
            info = mutation_els[1]
            mutationsFile.write(loc+"\t"+mutation+"\t"+info + "\t" + str(mutationDict[mutationKey]) + "\n")
        elif 'S' in mutationKey:
            mutation = 'S'
            mutation_els = mutationKey.split(mutation)
            loc = str(int(mutation_els[0]))
            info = mutation_els[1]
            mutationsFile.write(loc+"\t"+mutation+"\t"+info + "\t" + str(mutationDict[mutationKey]) + "\n")
    mutationsFile.close()

    logger.info("sorting...")
    subprocess.check_output("samtools sort -o " + output_root + ".bam " + unsortedOutName, shell=True)
    subprocess.check_output("samtools index " + output_root + ".bam", shell=True)
    if not keep_intermediate_files:
        logger.debug('removing ' + unsortedOutName)
        subprocess.check_output("rm " + unsortedOutName, shell=True)

    logger.info("sorting control...")
    subprocess.check_output("samtools sort -o " + output_root + ".ctl.bam " + unsortedCtlOutName, shell=True)
    subprocess.check_output("samtools index " + output_root + ".ctl.bam", shell=True)
    if not keep_intermediate_files:
        logger.debug('removing ' + unsortedCtlOutName)
        subprocess.check_output("rm " + unsortedCtlOutName, shell=True)

    totPrinted = printedNotAtCutSiteReadPairs + printedAsIsReadPairs + printedChangedReadPairs
    outstring = "Finished\n" 
    if only_insert_specific_mutation is not None:
        outstring += "Inserted mutation: " + only_insert_specific_mutation + "\n"
    else:
        if only_insert_first_indel:
            outstring += "Only inserted the first indel\n"
        outstring += "From Modified: keeping %d/%d modified Unpaired read pairs %sfor replacement\n"%(alteredUnpairedCount, totUnpairedCount, filter_indel_str)
        outstring += "From Modified: keeping %d/%d modified R1 read pairs %sfor replacement\n"%(alteredR1Count, totR1Count, filter_indel_str)
        outstring += "From Modified: keeping %d/%d modified R2 read pairs %sfor replacement\n"%(alteredR2Count, totR2Count, filter_indel_str)
    outstring += "From Unmodified: read " + str(readReadPairs) + " read pairs\n"
    outstring += "printed " + str(printedNotAtCutSiteReadPairs) + " read pairs not at the cut site (downsample pct was " + str(downsamplePct) + ")\n"
    outstring += "printed " + str(printedAsIsReadPairs) + " read pairs at cut site without modification\n"
    outstring += "printed " + str(printedChangedReadPairs) + " read pairs at cut site with modification\n"
    outstring += "printed " + str(totPrinted) + " read pairs to the treatment bam " + output_root + ".bam\n"
    outstring += "printed " + str(printedToControlReadPairs) + " read pairs to the control bam " + output_root + ".ctl.bam\n"


    logger.info(outstring)
    f1 = open(output_root+".info", "w")
    f1.write(outstring)
    f1.close()

def check_bam_indexed(bam_file_loc):
    """
    Check if a bam file is indexed
    """
    try:
        output = subprocess.check_output("samtools view " + bam_file_loc + " chr2345", stderr=subprocess.STDOUT, shell=True, encoding='utf-8')
        #should produce "[main_samview] region "chr2345" specifies an invalid region or unknown reference. Continue anyway."
        if 'Could not retrieve index file for' in output:
            raise Exception("Could not retrieve index file for " + bam_file_loc)
    except Exception as e:
        raise Exception("Could not open bam file " + bam_file_loc + ": " + str(e))

if __name__ == '__main__':
    main()