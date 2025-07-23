class Mutation:
    def __init__(self, name, chrom, pos, type, ref_seq, alt_seq, filter_str, passing_filter, wt_ref_count, wt_mut_count, sim_ref_count, sim_mut_count, vcf_line):
        """"Instantiate a new mutation object.
        """
        self.name = name
        self.chrom = chrom
        self.pos = pos
        self.type = type
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.filter_str = filter_str
        self.passing_filter = passing_filter
        self.wt_ref_count = wt_ref_count
        self.wt_mut_count = wt_mut_count
        self.sim_ref_count = sim_ref_count
        self.sim_mut_count = sim_mut_count
        self.vcf_line = vcf_line

    def print_mutation(self):
        print("Mutation: name: " + self.name +
                "\nchrom: " + self.chrom +
                "\npos: " + str(self.pos) +
                "\ntype: " + str(self.type) +
                "\nref_seq: " + str(self.ref_seq) +
                "\nalt_seq: " + str(self.alt_seq) +
                "\nfilter_str: " + str(self.filter_str) +
                "\npassing_filter: " + str(self.passing_filter) +
                "\nwt_ref_count: " + str(self.wt_ref_count) +
                "\nwt_mut_count: " + str(self.wt_mut_count) +
                "\nsim_ref_count: " + str(self.sim_ref_count) +
                "\nsim_mut_count: " + str(self.sim_mut_count) +
                "\nvcf_line: " + str(self.vcf_line))

    def __str__(self):
        return "Mutation: %s %s %s %s ref:%s alt:%s filter:%s %s %s %s %s"%(self.name, self.chrom, self.pos, self.type, 
                    self.ref_seq, self.alt_seq, 
                    'PASS' if self.passing_filter else 'NOT PASS', 
                    str(self.wt_ref_count), str(self.wt_mut_count), str(self.sim_ref_count), str(self.sim_mut_count))

    # two functions for storing mutations in a flat file
    def get_mutation_string(self):
        return "\t".join(str(x) for x in [self.name, self.chrom, self.pos, self.type, self.ref_seq, self.alt_seq, self.filter_str, self.passing_filter, self.wt_ref_count, self.wt_mut_count, self.sim_ref_count, self.sim_mut_count, self.vcf_line])

    @classmethod
    def get_mutation_header(cls):
        return "\t".join(str(x) for x in ['Name', 'Chrom', 'Pos', 'Type', 'Ref_seq', 'Alt_seq', 'Filter_str', 'Passing_filter', 'Wt_ref_count', 'Wt_mut_count', 'Sim_ref_count', 'Sim_mut_count', 'Vcf_line'])
    
    @classmethod
    def init_from_string(cls, mutation_str):
        str_els = mutation_str.strip().split("\t")
        name = str_els.pop(0)
        chrom = str_els.pop(0)
        pos = int(str_els.pop(0))
        type = str_els.pop(0)
        ref_seq = str_els.pop(0)
        alt_seq = str_els.pop(0)
        filter_str = str_els.pop(0)
        passing_filter_str = str_els.pop(0)
        passing_filter = False
        if 'true' in passing_filter_str.lower():
            passing_filter = True

        wt_ref_count_val = str_els.pop(0)
        if 'none' in wt_ref_count_val.lower():
            wt_ref_count = None
        else:
            wt_ref_count = int(wt_ref_count_val)

        wt_mut_count_val = str_els.pop(0)
        if 'none' in wt_mut_count_val.lower():
            wt_mut_count = None
        else:
            wt_mut_count = int(wt_mut_count_val)

        sim_ref_count_val = str_els.pop(0)
        if 'none' in sim_ref_count_val.lower():
            sim_ref_count = None
        else:
            sim_ref_count = int(sim_ref_count_val)

        sim_mut_count_val = str_els.pop(0)
        if 'none' in sim_mut_count_val.lower():
            sim_mut_count = None
        else:
            sim_mut_count = int(sim_mut_count_val)

        vcf_line = '\t'.join(str_els)
        return cls(name, chrom, pos, type, ref_seq, alt_seq, filter_str, passing_filter, wt_ref_count, wt_mut_count, sim_ref_count, sim_mut_count, vcf_line)