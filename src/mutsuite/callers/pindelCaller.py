from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil


class PindelCaller(Caller):

    def __init__(self, Config):
        swapChr = Config.get('Simulation', 'chr')  # location of simulated indels
        swapLoc = Config.getint('Simulation', 'loc')  # location of simulated indels
        fdrTolerances = [int(i) for i in Config.get('Simulation', 'FDRtolerance').split(",")]
        pindelRange = max(fdrTolerances) + 100  # pindel will look for indels within this range around the swap loc

        # get mean insert size
        unalteredBam = Config.get('Simulation', 'unalteredBam')
        mean_insert_size = -1
        with open(unalteredBam+".insertSizeMetrics", "r") as ub:
            for line in ub:
                if not line.startswith('## METRICS CLASS'):
                    continue
                metrics_line_els = ub.readline().split("\t")
                metrics_data_line_els = ub.readline().split("\t")
                metrics_dict = dict(zip(metrics_line_els, metrics_data_line_els))
                mean_insert_size = float(metrics_dict['MEDIAN_INSERT_SIZE'])

        if mean_insert_size < 0:
            raise Exception("Could not parse mean insert size from " + unalteredBam + ".insertSizeMetrics")

        self.call_chr = swapChr
        self.call_start = swapLoc - pindelRange
        self.call_end = swapLoc + pindelRange

        self.reference = Config.get('Simulation', 'reference')
        self.mean_insert_size = mean_insert_size

        # pindel_env_command is like: conda run -n conda_pindel bash -c "
        self.pindel_env_command = Config.get('Callers', 'pindel_env_command')

        self.reference_name = os.path.basename(self.reference).replace(".fa","")
        self.reference_date = "20101123"

        self.caller_threads = Config.getint('Callers', 'caller_threads', fallback=1)

    def get_name(self):
        return 'Pindel'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.pindel", sample_bam)
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)
        
        outfileName = os.path.join(outfolderName,"pindel")
        pindelConfigFile = outfileName + ".config"
        pFile = open(pindelConfigFile, "w")
        #pFile.write("/broad/hptmp/kendell/circle/simulations/alignNA12878/1kbAroundCutsite.BWA.bam " + str(meanInsertSize) + " NORMAL\n")
        pFile.write(control_bam + " " + str(self.mean_insert_size) + " NORMAL\n")
        pFile.write(sample_bam + " " + str(self.mean_insert_size) + " TUMOR\n")
        pFile.close()


#		r/--report_inversions           report inversions (default true)
#		-t/--report_duplications         report tandem duplications (default true)
#		-l/--report_long_insertions      report insertions of which the full sequence cannot be deduced because of their length (default true)
#		-k/--report_breakpoints          report breakpoints (default true)
#		-s/--report_close_mapped_reads   report reads of which only one end (the one closest to the mapped read of the paired-end read) could
#						    be mapped (default false)

        #pindel2vcf:
#       -G/--gatk_compatible  calls genotypes which could either be homozygous or heterozygous not as ./1 but as 0/1, to ensure compatibility with GATK
#       -mc/--min_coverage  The minimum number of reads to provide a genotype (default 10)

        ignoreFlags = "-r false -t false -l false -k false -s false"
        # command = "pindel -f /data/pinello/COMMON_DATA/REFERENCE_GENOMES/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa " + \
        tmp_prefix = outfileName + ".tmp."
        finished_file = outfileName + ".pindel.finished"
        clean_command = "rm -f " + tmp_prefix + "* && " 
        command = self.pindel_env_command + " pindel -f " + self.reference + " " + \
                ignoreFlags + " -i " + pindelConfigFile + " -o " + tmp_prefix + " -T " + str(self.caller_threads) + " -c " + self.call_chr + ":" + str(self.call_start) + "-" + str(self.call_end) + " && " + \
                    "pindel2vcf -P " + tmp_prefix + " -r " + self.reference + " -R " + self.reference_name + " -d " + self.reference_date + " -v " + outfileName + ".vcf -G -mc 1 && " + \
                    clean_command + \
                    " touch " + finished_file

        if self.pindel_env_command.endswith('"'):
            command += '"'

        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)

        vcf_file = outfileName + ".vcf"

        return (command, os.path.isfile(finished_file) and os.path.isfile(vcf_file))

    def get_results(self, sample_bam, control_bam, use_cached=True):
#From pindel documentation
#					The header line contains the following data:
#1) The index of the indel/SV (57 means that 57 insertions precede this insertion in the file)
#2) The type of indel/SV: I for insertion, D for deletion, INV for inversion, TD for tandem duplication
#3) The length of the SV
#4) "NT" (to indicate that the next number is the length of non-template sequences inserted; insertions are fully covered by the NT-fields, deletions can have NT bases if the deletion is not ‘pure’, meaning that while bases have been deleted, some bases have been inserted between the breakpoints)
#5) the length(s) of the NT fragment(s)
#6) the sequence(s) of the NT fragment(s)
#7-8) the identifier of the chromosome the read was found on
#9-10-11) BP: the start and end positions of the SV
#12-13-14) BP_range: if the exact position of the SV is unclear since bases at the edge of one read-half could equally well be appended to the other read-half. In the deletion example, ACA could be on any side of the gap, so the original deletion could have been between 1337143 and 1338815, between 1337144 and 1338816, or between 1337145 and 133817, or between 1337146 and 133818. BP-range is used to indicate this range.
#15) "Supports": announces that the total count of reads supporting the SV follow.
#16) The number of reads supporting the SV
#17) The number of unique reads supporting the SV (so not counting duplicate reads)
#18) +: supports from reads whose anchors are upstream of the SV
#19-20) total number of supporting reads and unique number of supporting reads whose anchors are upstream of the SV.
#21) -: supports from reads whose anchors are downstream of the SV
#22-23) total number of supporting reads and unique number of supporting reads whose anchors are downstream of the SV
#24-25) S1: a simple score, (“# +” + 1)* (“# -” + 1) ;
#26-27) SUM_MS: sum of mapping qualities of anchor reads, The reads with variants or unmapped are called split-read, whose mate is called anchor reads. We use anchor reads to narrow down the search space to speed up and increase sensitivity;
#28) the number of different samples scanned
#29-30-31) NumSupSamples?: the number of samples supporting the SV, as well as the number of samples having unique reads supporting the SV (in practice, these numbers are the same)
#32+) Per sample: the sample name, followed by the total number of supporting reads whose anchors are upstream, the total number of unique supporting reads whose anchors are upstream, the total number of supporting reads whose anchors are downstream, and finally the total number of unique supporting reads whose anchors are downstream.

        outfolderName = re.sub(r".bam$", ".caller.pindel", sample_bam)
        outfileName = os.path.join(outfolderName,"pindel")
        vcf_file = outfileName + ".vcf"

        cachedFile = outfolderName + ".results.txt"
        passing_indels = {}
        all_muts = {}

        if os.path.isfile(cachedFile) and use_cached:
            with open (cachedFile,'r') as f:
                header = f.readline()
                for line in f:
                    if line.strip() == "":
                        continue
                    mutation = Mutation.init_from_string(line)
                    if mutation.passing_filter and mutation.type in ['I', 'D']:
                        passing_indels[mutation.name] = mutation.count
                    all_muts[mutation.name] = mutation
            return passing_indels, all_muts


        #pindel's vcf doesn't include the normal sample
        vcf = cyvcf2.VCF(vcf_file)
        if len(vcf.samples) > 1:
            raise Exception('Expected only one sample in pindel vcf file ' + vcf_file + ' but found ' + str(len(vcf.samples)))

        sim_sample_index = 0

        for variant in vcf:
            is_passing = variant.FILTER is None
            ref_depths = variant.gt_ref_depths
            alt_depths = variant.gt_alt_depths
            sim_ref_depth = ref_depths[sim_sample_index]
            sim_alt_depth = alt_depths[sim_sample_index]

            # if multiple ALT alleles are called, just take the sum of them, and call the modification 'Multiple'
            if len(variant.ALT) > 1:
                thisKey = str(variant.start) + " M " + ",".join(variant.ALT)
                if is_passing:
                    passing_indels[thisKey] = int(sim_alt_depth)

            else:
                alt_seq = variant.ALT[0]
                seq_len_change = len(alt_seq) - len(variant.REF)

                if ((len(alt_seq) > 1 and len(variant.REF) > 1) # complex indel
                        or (seq_len_change != 0 and alt_seq[0] != variant.REF[0])):  # is not substitution and changed bases actually replace reference base
                    thisKey = str(variant.start) + " M " + variant.REF + "."+ alt_seq
                    if len(alt_seq) != len(variant.REF) and is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                elif seq_len_change > 0:  # insertion
                    # LOC INDEL INFO => COUNT
                    if alt_seq[0:len(variant.REF)] != variant.REF:
                        raise Exception("first base of insertion does not match reference: Ref:" + variant.REF + " Alt:" + alt_seq)
                    ins_seq = alt_seq[len(variant.REF):]
                    thisKey = str(variant.start) + " I " + ins_seq
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                elif seq_len_change < 0:   # deletion
                    # LOC INDEL INFO => COUNT
                    thisKey = str(variant.start) + " D " + str(-1 * seq_len_change)
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                else: # ..pindel doesn't find subs..
                    # LOC SUB INFO => COUNT
                    thisKey = str(variant.start) + " S " + alt_seq

            all_muts[thisKey] = Mutation(thisKey, variant.CHROM, variant.POS, 
                                variant.var_type, variant.REF, ",".join(variant.ALT),
                                variant.FILTER, is_passing,
                                None, None, int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))

        vcf.close()

        with open (cachedFile,'w') as f:
            f.write(Mutation.get_mutation_header()+"\n")
            for mutation in all_muts.values():
                f.write(mutation.get_mutation_string()+"\n")

        return passing_indels, all_muts

    def clean(self, sample_bam, control_bam):
        outfolderName = re.sub(r".bam$", ".caller.pindel", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)
