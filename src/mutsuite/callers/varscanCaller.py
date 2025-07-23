from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil

class VarscanCaller(Caller):

    def __init__(self, Config):

        self.reference = Config.get('Simulation','reference')
        self.varscan_env_command = Config.get('Callers', 'varscan_env_command')

    def get_name(self):
        return('Varscan')

    def run_caller(self,sample_bam,control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.varscan", sample_bam)
        
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)
        
        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)

        outfileRoot = os.path.join(outfolderName,"varscan")
        finishedFile = outfileRoot + ".finished"
        varscan_indel_vcf_file = outfileRoot + ".indel.vcf"
        varscan_snp_vcf_file = outfileRoot + ".snp.vcf"

        sim_mpileup = outfileRoot + ".sim.pileup"
        ctl_mpileup = outfileRoot + ".ctl.pileup"

        command = self.varscan_env_command + " samtools mpileup -f " + self.reference + " " + sample_bam + " > " + sim_mpileup + " " + \
                                          " && samtools mpileup -f " + self.reference + " " + control_bam + " > " + ctl_mpileup + " " + \
        " && varscan somatic " + ctl_mpileup + " " + sim_mpileup + " " + outfileRoot + " --output-vcf --min-coverage 1 --min-coverage-normal 1 --min-coverage-tumor 1 --somatic-p-value 1 " + \
        " && touch " + finishedFile
        #" && varscan somatic " + ctl_mpileup + " " + sim_mpileup + " " + outfileRoot + " --mpileup 1 --output-vcf --min-coverage 1 --min-coverage-normal 1 --min-coverage-tumor 1 --somatic-p-value 1 " + \
        
        if self.varscan_env_command.endswith('"'):
            command += '"'

        with open(outfileRoot + ".run.sh", 'w') as f:
            f.write(command)

        is_finished = os.path.isfile(finishedFile) and os.path.isfile(varscan_indel_vcf_file) and os.path.isfile(varscan_snp_vcf_file)
        return (command, is_finished)

    def get_results(self, sample_bam,control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.varscan", sample_bam)
        outfileRoot = os.path.join(outfolderName,"varscan")
        finishedFile = outfileRoot + ".finished"
        varscan_indel_vcf_file = outfileRoot + ".indel.vcf"
        varscan_snp_vcf_file = outfileRoot + ".snp.vcf"

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

        if not os.path.isfile(finishedFile):
            raise Exception('Command has not finished running yet (missing finished file ' + finishedFile + ')')  


        vcf = cyvcf2.VCF(varscan_indel_vcf_file)
        samples = vcf.samples
        ctl_sample_index = 0
        sim_sample_index = 1
        if samples[ctl_sample_index] != 'NORMAL':
            raise Exception("vcf " + varscan_indel_vcf_file + " does not have expected sample order (expecting NORMAL as sample " + str(ctl_sample_index) + ")")
        if samples[sim_sample_index] != 'TUMOR':
            raise Exception("vcf " + varscan_indel_vcf_file + " does not have expected sample order (expecting TUMOR as sample " + str(sim_sample_index) + ")")

        for variant in vcf:
            is_passing = variant.FILTER is None
            (sim_ref_fw_read_count, sim_ref_rv_read_count, sim_alt_fw_read_count, sim_alt_rv_read_count) = variant.format('DP4')[sim_sample_index].split(",")
            sim_alt_depth = int(sim_alt_fw_read_count) + int(sim_alt_rv_read_count)
            sim_ref_depth = int(sim_ref_fw_read_count) + int(sim_ref_rv_read_count)

            (ctl_ref_fw_read_count, ctl_ref_rv_read_count, ctl_alt_fw_read_count, ctl_alt_rv_read_count) = variant.format('DP4')[ctl_sample_index].split(",")
            ctl_alt_depth = int(ctl_alt_fw_read_count) + int(ctl_alt_rv_read_count)
            ctl_ref_depth = int(ctl_ref_fw_read_count) + int(ctl_ref_rv_read_count)

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
                    ins_seq = alt_seq[len(variant.REF):]
                    thisKey = str(variant.start) + " I " + ins_seq
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                elif seq_len_change < 0:   # deletion
                    # LOC INDEL INFO => COUNT
                    thisKey = str(variant.start) + " D " + str(-1 * seq_len_change)
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                else: 
                    # LOC SUB INFO => COUNT
                    thisKey = str(variant.start) + " S " + alt_seq

            all_muts[thisKey] = Mutation(thisKey, variant.CHROM, variant.POS, 
                                variant.var_type, variant.REF, ",".join(variant.ALT),
                                variant.FILTER, is_passing,
                                int(ctl_ref_depth), int(ctl_alt_depth), int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))
        vcf.close()

        vcf = cyvcf2.VCF(varscan_snp_vcf_file)
        samples = vcf.samples
        ctl_sample_index = 0
        sim_sample_index = 1
        if samples[ctl_sample_index] != 'NORMAL':
            raise Exception("vcf " + varscan_snp_vcf_file + " does not have expected sample order (expecting NORMAL as sample " + str(ctl_sample_index) + ")")
        if samples[sim_sample_index] != 'TUMOR':
            raise Exception("vcf " + varscan_snp_vcf_file + " does not have expected sample order (expecting TUMOR as sample " + str(sim_sample_index) + ")")

        for variant in vcf:
            is_passing = variant.FILTER is None
            (sim_ref_fw_read_count, sim_ref_rv_read_count, sim_alt_fw_read_count, sim_alt_rv_read_count) = variant.format('DP4')[sim_sample_index].split(",")
            sim_alt_depth = int(sim_alt_fw_read_count) + int(sim_alt_rv_read_count)
            sim_ref_depth = int(sim_ref_fw_read_count) + int(sim_ref_rv_read_count)

            (ctl_ref_fw_read_count, ctl_ref_rv_read_count, ctl_alt_fw_read_count, ctl_alt_rv_read_count) = variant.format('DP4')[ctl_sample_index].split(",")
            ctl_alt_depth = int(ctl_alt_fw_read_count) + int(ctl_alt_rv_read_count)
            ctl_ref_depth = int(ctl_ref_fw_read_count) + int(ctl_ref_rv_read_count)

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
                    ins_seq = alt_seq[len(variant.REF):]
                    thisKey = str(variant.start) + " I " + ins_seq
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                elif seq_len_change < 0:   # deletion
                    # LOC INDEL INFO => COUNT
                    thisKey = str(variant.start) + " D " + str(-1 * seq_len_change)
                    if is_passing:
                        passing_indels[thisKey] = int(sim_alt_depth)
                else: 
                    # LOC SUB INFO => COUNT
                    thisKey = str(variant.start) + " S " + alt_seq

            all_muts[thisKey] = Mutation(thisKey, variant.CHROM, variant.POS, 
                                variant.var_type, variant.REF, ",".join(variant.ALT),
                                variant.FILTER, is_passing,
                                int(ctl_ref_depth), int(ctl_alt_depth), int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))

        vcf.close()

        with open (cachedFile,'w') as f:
            f.write(Mutation.get_mutation_header()+"\n")
            for mutation in all_muts.values():
                f.write(mutation.get_mutation_string()+"\n")

        return passing_indels, all_muts

    def clean(self, sample_bam, control_bam):
        outfolderName = re.sub(r".bam$", ".caller.varscan", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)
