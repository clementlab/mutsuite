from callers.caller import Caller
from mutation import Mutation
import cyvcf2
import os
import re
import shutil


class StrelkaCaller(Caller):
    def __init__(self, Config):
        self.reference = Config.get('Simulation', 'reference')
        # assert that the genome.dict file exists for Mutect
        self.strelka_env_command = Config.get('Callers', 'strelka_env_command')
        self.caller_threads = Config.get('Callers', 'caller_threads', fallback=1)

    def get_name(self):
        return 'Strelka'

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.strelka", sample_bam)

        cachedFile = outfolderName + ".results.txt"
        if os.path.isfile(cachedFile) and use_cached:
            return ("", True)

        #create folder if it doesn't exist
        if not os.path.exists(outfolderName):
            os.makedirs(outfolderName)

        outfileName = os.path.join(outfolderName,"strelka")
        strelka_indel_vcf = os.path.join(outfolderName,"results/variants/somatic.indels.vcf.gz")
        strelka_snv_vcf = os.path.join(outfolderName,"results/variants/somatic.snvs.vcf.gz")
        finishedFile = outfileName + ".finished"

        manta_outfolder = os.path.join(outfolderName,"manta")
        if not os.path.exists(manta_outfolder):
            os.makedirs(manta_outfolder)

        mantaVCF = os.path.join(manta_outfolder, "results/variants/candidateSmallIndels.vcf.gz")

        command = self.strelka_env_command + " rm -rf " + outfolderName + " && rm -rf " + manta_outfolder + " && configManta.py " + \
            " --normalBam " + control_bam + \
            " --tumorBam " + sample_bam + \
            " --referenceFasta " + self.reference + \
            " --runDir " + manta_outfolder + \
            " && " + manta_outfolder + "/runWorkflow.py -m local -j " + str(self.caller_threads) + \
            " && configureStrelkaSomaticWorkflow.py " + \
            " --normalBam " + control_bam + \
            " --tumorBam " + sample_bam + \
            " --referenceFasta " + self.reference + \
            " --indelCandidates " + mantaVCF + \
            " --runDir " + outfolderName + \
            " && " + outfolderName + "/runWorkflow.py -m local -j " + str(self.caller_threads) + \
            " && touch " + finishedFile
        
        if self.strelka_env_command.endswith('"'):
            command += '"'
        
        with open(outfileName + ".run.sh", 'w') as f:
            f.write(command)
        
        is_finished = os.path.isfile(finishedFile) and os.path.isfile(strelka_indel_vcf) and os.path.isfile(strelka_snv_vcf)

        return (command, is_finished)


    def get_results(self, sample_bam, control_bam, use_cached=True):
        outfolderName = re.sub(r".bam$", ".caller.strelka", sample_bam)
        strelka_indel_vcf = os.path.join(outfolderName,"results/variants/somatic.indels.vcf.gz")
        strelka_snv_vcf = os.path.join(outfolderName,"results/variants/somatic.snvs.vcf.gz")

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


        # Somatic SNVs:
        # refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
        # altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
        # tier1RefCounts = First comma-delimited value from $refCounts
        # tier1AltCounts = First comma-delimited value from $altCounts
        # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
        # Somatic indels:
        # tier1RefCounts = First comma-delimited value from FORMAT/TAR
        # tier1AltCounts = First comma-delimited value from FORMAT/TIR
        # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
        vcf = cyvcf2.VCF(strelka_indel_vcf)

        if len(vcf.samples) != 2:
            raise Exception('Expected only one sample in strelka vcf file ' + strelka_indel_vcf + ' but found ' + str(len(vcf.samples)))
        if vcf.samples[0] != 'NORMAL':
            raise Exception('Expected first sample in strelka vcf file ' + strelka_indel_vcf + ' to be "simulated" but found ' + vcf.samples[0])
        if vcf.samples[1] != 'TUMOR':
            raise Exception('Expected second sample in strelka vcf file ' + strelka_indel_vcf + ' to be "simulatedCtl" but found ' + vcf.samples[1])

        ctl_sample_index = 0
        sim_sample_index = 1


        for variant in vcf:
            is_passing = variant.FILTER is None

            sim_alt_depth = int(variant.format('TIR')[sim_sample_index][0])
            sim_ref_depth = int(variant.format('TAR')[sim_sample_index][0])

            ctl_alt_depth = int(variant.format('TIR')[ctl_sample_index][0])
            ctl_ref_depth = int(variant.format('TAR')[ctl_sample_index][0])

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

        vcf = cyvcf2.VCF(strelka_snv_vcf)

        ctl_sample_index = 0
        sim_sample_index = 1

        if len(vcf.samples) != 2:
            raise Exception('Expected only one sample in strelka vcf file ' + strelka_snv_vcf + ' but found ' + str(len(vcf.samples)))
        if vcf.samples[0] != 'NORMAL':
            raise Exception('Expected first sample in strelka vcf file ' + strelka_snv_vcf + ' to be "simulated" but found ' + vcf.samples[ctl_sample_index])
        if vcf.samples[1] != 'TUMOR':
            raise Exception('Expected second sample in strelka vcf file ' + strelka_snv_vcf + ' to be "simulatedCtl" but found ' + vcf.samples[sim_sample_index])

        for variant in vcf:
            is_passing = variant.FILTER is None
            
            alt_base = variant.ALT[0]
            ref_base = variant.REF
            sim_alt_depth = int(variant.format(alt_base +"U")[sim_sample_index][0])
            sim_ref_depth = int(variant.format(ref_base +"U")[sim_sample_index][0])

            ctl_alt_depth = int(variant.format(alt_base +"U")[ctl_sample_index][0])
            ctl_ref_depth = int(variant.format(ref_base +"U")[ctl_sample_index][0])

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
        outfolderName = re.sub(r".bam$", ".caller.strelka", sample_bam)
        if os.path.isdir(outfolderName):
            shutil.rmtree(outfolderName)


