class Caller:
    def __init__(self):
        """"Instantiate a new caller object.
        Any variables unique to this caller should be set here.
        For example, if this caller requires the path to the fastq reference, that should be set here.
        """
        pass

    def get_name(self):
        """Get the name of the caller.
        This name will be used in output that compares this caller to other callers
        Returns: A string with a name identifying the caller
        """
        pass

    def run_caller(self, sample_bam, control_bam, use_cached=True):
        """Extract indels from a pair of bams.
        sample_bam is the bam with simulated indels
        control_bam is the bam with no simulated indels
        use_cached is a boolean indicating whether to use cached results if they exist. Cached results are accessible even after 'clean' has been called.
        Returns: a tuple of:
            A string with a command to run this indel for this pair of bams
            A boolean for whether the command has already been completed successfully (and should not be run again)
        """
        pass

    def get_results(self, sample_bam, control_bam, use_cached=True):
        """Extract indel calling results.
        sample_bam is the bam with simulated indels
        control_bam is the bam with no simulated indels
        use_cached is a boolean indicating whether to use cached results if they exist. Cached results are accessible even after 'clean' has been called.
        run_caller has been run previously. This function collects the results for indel calling by this caller
        Returns: 
            a dictionary of mutations passing filters: 'LOC TYPE INFO' => count
                where TYPE is either 'I' for insertion, 'D' for deletion, 'S' for substitution or 'M' for mixed/multiple
                if TYPE is 'I' INFO will be the sequence of the insertion
                if TYPE is 'D' INFO will be the length of the deletion
                if TYPE is 'S' INFO will be the substituted base
                e.g. if a 25bp deletion was observed at position 100 supported by 5 reads, the key would be '100 D 25' and the value would be 5
                e.g. if a 2bp insertion of AT was observed at position 50 supported by 3 reads, the key would be '50 I AT' and the value would be 3
            a dictionary of all mutations: 'LOC TYPE INFO' => Mutation object
        """
        pass
    
    def clean(self, sample_bam, control_bam):
        """
        Remove all intermediate files associated with calling mutations for this sample
        """
        pass