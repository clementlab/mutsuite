from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil

# ignoring this caller because it only calls SNPS? and includes the following restrictions:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268238/
# - Site is >10 bp from a predicted indel of quality ≥50.
# - Maximum mapping quality at the site is ≥40.
# - Fewer than three SNV calls in a 10 bp window around the site.
# - Site is covered by at least three reads.
# - Consensus quality ≥20.
# - Single Nucleotide Polymorphism (SNP) quality ≥20.

class SomaticSniperCaller(Caller):

    def __init__(self, Config):
        self.reference = Config.get('Simulation', 'reference')
        # somaticsniper_env_command is like: conda run -n conda_somaticsniper bash -c "
        self.somaticsniper_env_command = Config.get('Callers', 'somaticsniper_env_command')

    def get_name(self):
        return 'SomaticSniper'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.somaticSniper", sample_bam)
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)
            
        outfileName = os.path.join(outfolderName,"somaticSniper")
        vcfName = outfileName + ".vcf"
        finishedFile = outfileName + ".finished"

        # -N  INT number of haplotypes in the sample (for -c/-g) [2]
        # -F  STRING select output format (vcf or classic) [classic]
        command = self.somaticsniper_env_command + "bam-somaticsniper -N 100 -F vcf -f " + self.reference + " " + sample_bam + " " + control_bam + " " + vcfName + \
            " && touch " + finishedFile

        if self.somaticsniper_env_command.endswith('"'):
            command += '"'

        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)

        return (command, os.path.isfile(finishedFile) and os.path.isfile(vcfName))

    def get_results(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.somaticSniper", sample_bam)
        vcfFileName = os.path.join(outfolderName,"somaticSniper.vcf")

        passing_indels = {}
        all_muts = {}
        cachedFile = outfolderName + ".results.txt"

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



        vcf = cyvcf2.VCF(vcfFileName)
        samples = vcf.samples
        ctl_sample_index = 0
        sim_sample_index = 1
        if samples[ctl_sample_index] != 'NORMAL':
            raise Exception("vcf " + vcfFileName + " does not have expected sample order (expecting NORMAL as sample " + str(ctl_sample_index) + ")")
        if samples[sim_sample_index] != 'TUMOR':
            raise Exception("vcf " + vcfFileName + " does not have expected sample order (expecting TUMOR as sample " + str(sim_sample_index) + ")")

        for variant in vcf:
            is_passing = variant.FILTER is None
            (sim_ref_fw_read_count, sim_ref_rv_read_count, sim_alt_fw_read_count, sim_alt_rv_read_count) = variant.format('DP4')[sim_sample_index]
            sim_alt_depth = sim_alt_fw_read_count + sim_alt_rv_read_count
            sim_ref_depth = sim_ref_fw_read_count + sim_ref_rv_read_count

            (ctl_ref_fw_read_count, ctl_ref_rv_read_count, ctl_alt_fw_read_count, ctl_alt_rv_read_count) = variant.format('DP4')[ctl_sample_index]
            ctl_alt_depth = ctl_alt_fw_read_count + ctl_alt_rv_read_count
            ctl_ref_depth = ctl_ref_fw_read_count + ctl_ref_rv_read_count

            print('got depths: sim_ref_depth: ' + str(sim_ref_depth) + ' sim_alt_depth: ' + str(sim_alt_depth) + ' ctl_ref_depth: ' + str(ctl_ref_depth) + ' ctl_alt_depth: ' + str(ctl_alt_depth))

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
                                int(ctl_ref_depth), int(ctl_alt_depth), int(sim_ref_depth), int(sim_alt_depth),
                                str(variant))
        vcf.close()

        with open (cachedFile,'w') as f:
            f.write(Mutation.get_mutation_header()+"\n")
            for mutation in all_muts.values():
                f.write(mutation.get_mutation_string()+"\n")

        return passing_indels, all_muts

    def clean(self, sample_bam, control_bam):
        outfolderName = re.sub(r".bam$", ".caller.somaticSniper", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)

