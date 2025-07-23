from callers.caller import Caller
from mutation import Mutation


class SimCaller(Caller):

    def __init__(self, Config):
        self.call_chr = Config.get('Simulation', 'chr')  # location of simulated indels

    def run_caller(self, sample_bam, control_bam):
        return []

    def get_results(self, sample_bam, control_bam):
        mutsIn = sample_bam.replace(".bam", ".mutations.txt")
        passing_indels = {}
        all_muts = {}
        with open(mutsIn, "r") as mutsFile:
            headLine = mutsFile.readline()
            #LOC MUTATIONTYPE INFO COUNT
            for line in mutsFile:
                line.strip()
                lineEls = line.split()

                mutation_loc = int(lineEls[0]) - 1

                # mutation callers return the base before deleted base
                if lineEls[1] in ['D']:
                    mutation_loc -= 1

                mutationKey = str(int(lineEls[0]) -1) + " " + lineEls[1] + " " + lineEls[2]

                if lineEls[1] in ['I', 'D']:
                    passing_indels[mutationKey] = int(lineEls[3])

                all_muts[mutationKey] = Mutation(mutationKey, self.call_chr, int(lineEls[0]) - 1,
                                    lineEls[1], None, lineEls[2],
                                    None, True,
                                    None, None, None, int(lineEls[3]),
                                    line)

        return passing_indels, all_muts
