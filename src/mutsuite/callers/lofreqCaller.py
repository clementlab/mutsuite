from callers.caller import Caller
from mutation import Mutation
import configparser
import cyvcf2
import os
import re
import shutil

class LofreqCaller(Caller):
    def __init__(self, Config):
        self.reference = Config.get('Simulation', 'reference')
        # lofreq_env_command is like: conda run -n conda_lofreq bash -c "
        self.lofreq_env_command = Config.get('Callers', 'lofreq_env_command')



    def get_name(self):
        return 'Lofreq'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.lofreq", sample_bam)
        
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)

        
        outfileName = os.path.join(outfolderName,"lofreq")
        lofreqBam = outfileName + ".lofreq_indelqual.sim.bam"
        lofreqCtlBam = outfileName + ".lofreq_indelqual.ctl.bam"
        finishedFile = outfileName + ".finished"
        lofreqIndelFileName = os.path.join(outfolderName, "lofreq_somatic_final.indels.vcf.gz")
        lofreqSNVFileName = os.path.join(outfolderName, "lofreq_somatic_final.snvs.vcf.gz")
        
        command = self.lofreq_env_command + " rm -f " + outfileName + "* && " + \
            "lofreq indelqual --dindel -f " + self.reference + " -o " + lofreqBam + " " + sample_bam + " && samtools index " + lofreqBam + " && " + \
            "lofreq indelqual --dindel -f " + self.reference + " -o " + lofreqCtlBam + " " + control_bam + " && samtools index " + lofreqCtlBam + " && " + \
            "lofreq somatic -f " + self.reference + " " + "-t " + lofreqBam + " -n " + lofreqCtlBam + " --min-cov 0 --call-indels -o " + outfileName + "_" + \
            " && touch " + finishedFile

        if self.lofreq_env_command.endswith('"'):
            command += '"'

        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)

        is_finished = os.path.isfile(finishedFile) and os.path.isfile(lofreqIndelFileName) and os.path.isfile(lofreqSNVFileName)
        return (command, is_finished)

    def get_results(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.lofreq", sample_bam)
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

        lofreqIndelFileName = os.path.join(outfolderName, "lofreq_somatic_final.indels.vcf.gz")
        indel_vcf = cyvcf2.VCF(lofreqIndelFileName)
        for variant in indel_vcf:
            is_passing = variant.FILTER is None

            dp4 = variant.INFO.get('DP4')
            (ref_fw_read_count, ref_rv_read_count, alt_fw_read_count, alt_rv_read_count) = dp4

            sim_ref_depth = alt_fw_read_count + alt_rv_read_count
            sim_alt_depth = alt_fw_read_count + alt_rv_read_count
            ctl_ref_depth = None
            ctl_alt_depth = None

            # if multiple ALT alleles are called, just take the sum of them, and call the modification 'Multiple'
            # variant.ALT is an array of ALT alleles
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
                                None, None, int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))

        lofreqSNVFileName = os.path.join(outfolderName, "lofreq_somatic_final.snvs.vcf.gz")
        for variant in cyvcf2.VCF(lofreqSNVFileName):
            is_passing = variant.FILTER is None
            ref_depths = variant.gt_ref_depths
            alt_depths = variant.gt_alt_depths


            dp4 = variant.INFO.get('DP4')
            (ref_fw_read_count, ref_rv_read_count, alt_fw_read_count, alt_rv_read_count) = dp4

            sim_ref_depth = alt_fw_read_count + alt_rv_read_count
            sim_alt_depth = alt_fw_read_count + alt_rv_read_count
            ctl_ref_depth = None
            ctl_alt_depth = None

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
                                None, None, int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))

        with open (cachedFile,'w') as f:
            f.write(Mutation.get_mutation_header()+"\n")
            for mutation in all_muts.values():
                f.write(mutation.get_mutation_string()+"\n")

        return passing_indels, all_muts

    def clean(self, sample_bam, control_bam):
        outfolderName = re.sub(r".bam$", ".caller.lofreq", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)
