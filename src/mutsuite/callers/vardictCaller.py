from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil


class VardictCaller(Caller):

    def __init__(self, Config):
        self.reference = Config.get('Simulation', 'reference')
        # vardict_env_command is like: conda run -n conda_vardict bash -c "
        self.vardict_env_command = Config.get('Callers', 'vardict_env_command')
        swapChr = Config.get('Simulation', 'chr')  # location of simulated indels
        swapLoc = Config.getint('Simulation', 'loc')  # location of simulated indels
        fdrTolerances = [int(i) for i in Config.get('Simulation', 'FDRtolerance').split(",")]
        callRange = max(fdrTolerances) + 100  # vardict will look for indels within this range around the swap loc

        self.call_chr = swapChr
        self.call_start = swapLoc - callRange
        self.call_end = swapLoc + callRange

    def get_name(self):
        return 'VarDict'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.vardict", sample_bam)
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)

        outfileName = os.path.join(outfolderName,"vardict")
        vcf_file = outfileName + ".vcf"
        finishedFile = outfileName + ".finished"

        obs_threshold = 0.0001
        obs_count_threshold = 1
        """
        -f double
            The threshold for allele frequency, default: 0.01 or 1%
        -R Region
            The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
        -r minimum reads
            The minimum # of variance reads, default 2
        """
        """
        var2vcf_paired.pl
        -A
            Indicate to output all variants at the same position. By default, only the variant with the highest allele frequency is converted to VCF.
        """
       
        # vardict requries piping which doesn't work well with conda run
        # vardict_env_command is like: conda run -n kc_vardict bash -c "
        command = self.vardict_env_command + " vardict -f " + str(obs_threshold) + \
            " -G " + self.reference + " -b '" + sample_bam + "|" + control_bam + "' " \
            " -R " + self.call_chr + ":" + str(self.call_start) +"-" + str(self.call_end) + \
            " -r " + str(obs_count_threshold) + " | testsomatic.R | var2vcf_paired.pl -A -N 'simulated|simulatedCtl' -f " + str(obs_threshold) + " > " + vcf_file + \
            " && touch " + finishedFile

        if self.vardict_env_command.endswith('"'):
            command += '"'

        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)

        return (command, os.path.isfile(finishedFile))

    def get_results(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.vardict", sample_bam)
        vardict_vcf_name = os.path.join(outfolderName,"vardict.vcf")
    
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


        vcf = cyvcf2.VCF(vardict_vcf_name)
        samples = vcf.samples
        sim_sample_index = 0
        ctl_sample_index = 1
        if samples[sim_sample_index] != 'simulated':
            raise Exception("vcf " + vardict_vcf_name + " does not have expected sample order (expecting simulated as sample 0)")

        for variant in vcf:
            is_passing = variant.FILTER is None
            ref_depths = variant.gt_ref_depths
            alt_depths = variant.gt_alt_depths
            sim_ref_depth = ref_depths[sim_sample_index]
            sim_alt_depth = alt_depths[sim_sample_index]
            ctl_ref_depth = ref_depths[ctl_sample_index]
            ctl_alt_depth = alt_depths[ctl_sample_index]

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
        outfolderName = re.sub(r".bam$", ".caller.vardict", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)

