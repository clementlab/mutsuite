from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil


class MutectCaller(Caller):
    def __init__(self, Config):
        self.reference = Config.get('Simulation', 'reference')
        # assert that the genome.dict file exists for Mutect
        if not os.path.isfile(self.reference + ".dict") and not os.path.isfile(re.sub(r".fa$",".dict",self.reference)):
            raise Exception("Mutect requires a .dict file for the reference genome. Please create one using Picard's CreateSequenceDictionary.jar")
        # gatk_env_command is like: conda run -n conda_gatk bash -c "
        self.gatk_env_command = Config.get('Callers', 'gatk_env_command')
        
        swapChr = Config.get('Simulation', 'chr')  # location of simulated indels
        swapLoc = Config.getint('Simulation', 'loc')  # location of simulated indels
        fdrTolerances = [int(i) for i in Config.get('Simulation', 'FDRtolerance').split(",")]
        callRange = max(fdrTolerances) + 100  # mutect will look for indels within this range around the swap loc
        self.call_chr = swapChr
        self.call_start = swapLoc - callRange
        self.call_end = swapLoc + callRange

    def get_name(self):
        return 'Mutect'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.mutect", sample_bam)
        
        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)

        outfileName = os.path.join(outfolderName,"mutect")
        unfilteredVCF = outfileName + ".unfiltered.vcf"
        filteredVCF = outfileName + ".filtered.vcf"
        finishedFile = outfileName + ".finished"

        command = self.gatk_env_command + "gatk Mutect2 -R " + self.reference + " -I " + sample_bam + " -I " + control_bam + \
            " -normal simulatedCtl -O " + unfilteredVCF + " -L " + self.call_chr + ":" + str(self.call_start) + "-" + str(self.call_end) + " && " + \
            " gatk FilterMutectCalls -R " + self.reference + " -V " + unfilteredVCF + " -O  " + filteredVCF + " && " +\
            "touch " + finishedFile
        
        if self.gatk_env_command.endswith('"'):
            command += '"'

        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)

        is_finished = os.path.exists(finishedFile) and os.path.isfile(unfilteredVCF) and os.path.isfile(filteredVCF)

        return (command, is_finished)


    def get_results(self, sample_bam, control_bam, use_cached=True):
        passing_indels = {}
        all_muts = {}
        outfolderName = re.sub(r".bam$", ".caller.mutect", sample_bam)
        outfileName = os.path.join(outfolderName,"mutect")
        filteredVCF = outfileName + ".filtered.vcf"

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

        vcf = cyvcf2.VCF(filteredVCF)

        if len(vcf.samples) != 2:
            raise Exception('Expected only one sample in mutect vcf file ' + filteredVCF + ' but found ' + str(len(vcf.samples)))
        if vcf.samples[0] != 'simulated':
            raise Exception('Expected first sample in mutect vcf file ' + filteredVCF + ' to be "simulated" but found ' + vcf.samples[0])
        if vcf.samples[1] != 'simulatedCtl':
            raise Exception('Expected second sample in mutect vcf file ' + filteredVCF + ' to be "simulatedCtl" but found ' + vcf.samples[0])

        sim_sample_index = 0
        ctl_sample_index = 1


        for variant in vcf:
            is_passing = variant.FILTER is None
            ref_depths = variant.gt_ref_depths
            alt_depths = variant.gt_alt_depths
            sim_ref_depth = ref_depths[sim_sample_index]
            sim_alt_depth = alt_depths[sim_sample_index]
            ctl_ref_depth = ref_depths[ctl_sample_index]
            ctl_alt_depth = alt_depths[ctl_sample_index]

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
        outfolderName = re.sub(r".bam$", ".caller.mutect", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)
