from collections import defaultdict

def aggregate_indels(caller_mutations):
    """Aggregate counts of each mutations type across simulated samples"""
    aggregated_mutations = {}
    for name in caller_mutations:
        for key in caller_mutations[name]:
            if key not in aggregated_mutations:
                aggregated_mutations[key] = 0
            aggregated_mutations[key] += caller_mutations[name][key]
    return aggregated_mutations

def aggregate_mutations(caller_mutations, passing_filter=False):
    """Aggregate counts of each mutations type across simulated samples
        if passing_filter == True, only mutations passing caller filters are aggregated
    """
    aggregated_mutations = defaultdict(int)
    for sim_name in caller_mutations:
        for key in caller_mutations[sim_name]:
            if passing_filter and not caller_mutations[sim_name][key].passing_filter:
                continue
            mut_count = caller_mutations[sim_name][key].sim_mut_count
            aggregated_mutations[key] += mut_count
    return aggregated_mutations

def print_mutations(mutations, fileName=None):
    """Prints counts of mutations for a single caller across simulated samples
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    """
    out_string = "Simulated Sample\tBP\tMUTATION\tLOC\tCOUNT\n"
    for name in sorted(mutations):
        for key in sorted(mutations[name]):
            out_string += name + "\t" + key.replace(" ", "\t")

    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)


def print_aggregate_mutations(mutations, fileName=None):
    """Prints counts of aggregated mutations across simulated samples
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    """
    out_string = "BP\tMUTATION\tLOC\tCOUNT\n"
    for key in sorted(mutations):
        out_string += key.replace(" ", "\t")

    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)


def print_all_aggregate_mutations(all_aggregated_mutations, aggregated_simulated_mutations, callers, fileName=None, mutation_types_to_report=None):
    """For all callers, prints counts of aggregated mutations
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    """
    all_mutations = {}  # keep track of all possible keys
    for mutation in aggregated_simulated_mutations:
        if mutation not in all_mutations:
            all_mutations[mutation] = 0

    head_string = "BP\tMUTATION\tLOC\tSIMULATED"
    caller_names = []
    for caller in callers:
        caller_name = caller.get_name()
        caller_names.append(caller_name)
        head_string += "\t"+caller_name
        for mutation in all_aggregated_mutations[caller_name]:
            if mutation not in all_mutations:
                all_mutations[mutation] = 0

    out_string = ""
    sorted_mutations = sorted(all_mutations, key=lambda x: int(x.split(" ")[0]), reverse=False)
    for mutation in sorted_mutations:
        if mutation_types_to_report is not None and mutation.split(" ")[1] not in mutation_types_to_report:
            continue
        out_string += mutation.replace(" ", "\t")
        sim_count = 0
        if mutation in aggregated_simulated_mutations:
            sim_count = aggregated_simulated_mutations[mutation]
        out_string += "\t" + str(sim_count)

        for caller_name in caller_names:
            caller_count = 0
            if caller_name in all_aggregated_mutations and mutation in all_aggregated_mutations[caller_name]:
                caller_count = all_aggregated_mutations[caller_name][mutation]
            out_string += "\t" + str(caller_count)
        out_string += "\n"

    out_string = head_string + "\n" + out_string
    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)
    return out_string

def print_all_mutations_by_sample(called_mutations,simulated_mutations,callers,fileName=None, mutation_types_to_report=None, collapse_replicates=False, require_passing=False):
    """For all callers, prints counts of aggregated mutations
    called_mutations: mutations called by callers in form of called_mutations[caller_name][sample_name][mutation] = Mutation object
    simulated_mutations: mutations simulated in form of simulated_mutations[sample_name][mutation] = Mutation object
    callers: list of callers
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    mutation_types_to_report: list of mutation types to include (e.g. ['I', 'D'])
    collapse_replicates: if true, mutation counts will be collapsed by replicate
    require_passing: if true, only mutations passing caller filters will be included
    """
    sim_mutations_by_sample = defaultdict(lambda:defaultdict(int))
    #first, count all simulated mutations for each sample
    for sample_name in simulated_mutations:
        sample_name_for_output = sample_name
        if collapse_replicates:
            sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
        for mutation in simulated_mutations[sample_name]:
            if mutation not in sim_mutations_by_sample[sample_name]:
                sim_mutations_by_sample[sample_name_for_output][mutation] += simulated_mutations[sample_name][mutation].sim_mut_count

    #next iterate through mutations 
    called_mutations_by_sample = defaultdict(lambda:defaultdict(int)) # called_mutation_by_sample[sample_name][mutation] = count (collapsed by sample_name_for_output)
    called_mutations_by_sample_and_caller = defaultdict(lambda:defaultdict(lambda:defaultdict(int))) # called_mutation_by_sample[sample_name][caller][mutation] = count (collapsed by sample_name_for_output)
    caller_names = []
    sample_names = {}
    sample_output_names = {}
    for caller in callers:
        caller_name = caller.get_name()
        caller_names.append(caller_name)
        for sample_name in called_mutations[caller_name]:
            sample_names[sample_name] = 1
            sample_name_for_output = sample_name
            if collapse_replicates:
                sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
            sample_output_names[sample_name_for_output] = 1
            for mutation in called_mutations[caller_name][sample_name]:
                this_mutation_loc, this_mutation_type, this_mutation_info = mutation.split(" ")
                if mutation_types_to_report is not None and this_mutation_type not in mutation_types_to_report:
                    continue
                if require_passing and not called_mutations[caller_name][sample_name][mutation].passing_filter:
                    continue
                called_mutations_by_sample[sample_name_for_output][mutation] += 1
                called_mutations_by_sample_and_caller[sample_name_for_output][caller_name][mutation] += called_mutations[caller_name][sample_name][mutation].sim_mut_count

    head_string = "SAMPLE\tBP\tMUTATION\tLOC\tSIMULATED"
    head_string += "\t"+"\t".join(caller_names)
    out_string = ""
    for sample_name in sample_output_names:
        sorted_mutations = sorted(called_mutations_by_sample[sample_name], key=lambda x: int(x.split(" ")[0]), reverse=False)
        for mutation in sorted_mutations:
            out_string += sample_name + "\t" + mutation.replace(" ", "\t")
            out_string += "\t" + str(sim_mutations_by_sample[sample_name][mutation])
            for caller_name in caller_names:
                out_string += "\t" + str(called_mutations_by_sample_and_caller[sample_name][caller_name][mutation])
            out_string += "\n"

    out_string = head_string + "\n" + out_string
    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)
    return out_string

def print_recovery_reads_by_sample(called_mutations,simulated_mutations,callers,fileName=None,mutation_types_to_report=None,nearby_distance_cutoff=10,false_positive_distance_cutoffs=[1000],collapse_replicates=False, require_passing=False):
    """For all callers, prints read counts for simulated mutations and read counts for recovered mutations for each sample
    if mutation is within distance_cutoff of a simulated mutation, it is considered recovered_nearby
    called_mutations: mutations called by callers in form of called_mutations[caller_name][sample_name][mutation] = Mutation object
    simulated_mutations: mutations simulated in form of simulated_mutations[sample_name][mutation] = Mutation object
    callers: list of callers
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    mutation_types_to_report: list of mutation types to include (e.g. ['I', 'D'])
    nearby_distance_cutoff: distance from simulated mutation for a mutation to be considered recovered_nearby
    false_positive_distance_cutoffs: list of distances from simulated mutation for a mutation to be considered a false positive (e.g. reporting false-positives (subject to require_passing and mutation_types_to_report) called within 100bp, 1000bp, etc.)
    collapse_replicates: if true, mutation counts will be collapsed by replicate
    require_passing: if true, only mutations passing caller filters will be included
    """

    max_false_positive_cutoff = max(false_positive_distance_cutoffs)

    sim_mutations_list_by_sample = defaultdict(list)
    #first, aggregate all simulated mutations
    sim_count_by_sample = defaultdict(int)
    for sample_name in simulated_mutations:
        sample_name_for_output = sample_name
        if collapse_replicates:
            sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
        for mutation in simulated_mutations[sample_name]:
            if mutation not in sim_mutations_list_by_sample[sample_name]:
                sim_mutations_list_by_sample[sample_name].append(mutation)
                sim_count_by_sample[sample_name_for_output] += simulated_mutations[sample_name][mutation].sim_mut_count


    #next iterate through mutations 
    head_string = "SAMPLE\tCLASS\tSIMULATED"
    caller_names = []
    sample_names = {}
    sample_output_names = {}
    sample_caller_exact_count = defaultdict(lambda: defaultdict(int))
    sample_caller_nearby_count = defaultdict(lambda: defaultdict(int))
    sample_caller_false_positive_count = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for caller in callers:
        caller_name = caller.get_name()
        caller_names.append(caller_name)
        head_string += "\t"+caller_name
        for sample_name in called_mutations[caller_name]:
            sample_names[sample_name] = 1
            sample_name_for_output = sample_name
            if collapse_replicates:
                sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
            sample_output_names[sample_name_for_output] = 1
            for mutation in called_mutations[caller_name][sample_name]:
                this_mutation_loc, this_mutation_type, this_mutation_info = mutation.split(" ")
                if mutation_types_to_report is not None and this_mutation_type not in mutation_types_to_report:
                    continue
                if require_passing and not called_mutations[caller_name][sample_name][mutation].passing_filter:
                    continue
                found_mutation_as_simulated = False
                min_distance_to_sim_mutation = max_false_positive_cutoff * 2
                for sim_mutation in sim_mutations_list_by_sample[sample_name]:
                    sim_mutation_loc, sim_mutation_type, sim_mutation_info = sim_mutation.split(" ")
                    this_dist_to_sim_mutation = abs(int(this_mutation_loc) - int(sim_mutation_loc))
                    if this_dist_to_sim_mutation < min_distance_to_sim_mutation:
                        min_distance_to_sim_mutation = this_dist_to_sim_mutation
                    if sim_mutation == mutation:
                        sample_caller_exact_count[sample_name_for_output][caller_name] += called_mutations[caller_name][sample_name][mutation].sim_mut_count
                        found_mutation_as_simulated = True
                    elif sim_mutation_type == this_mutation_type:
                        if this_dist_to_sim_mutation <= nearby_distance_cutoff:
                            sample_caller_nearby_count[sample_name_for_output][caller_name] += called_mutations[caller_name][sample_name][mutation].sim_mut_count
                            found_mutation_as_simulated = True
                        elif this_mutation_type == "D":
                            if int(this_mutation_loc) <= int(sim_mutation_loc) and int(this_mutation_loc) + int(this_mutation_info) >= int(sim_mutation_loc): #if deletion overlaps simulated deletion
                                sample_caller_nearby_count[sample_name_for_output][caller_name] += called_mutations[caller_name][sample_name][mutation].sim_mut_count
                                found_mutation_as_simulated = True

                    if not found_mutation_as_simulated:
                        for false_positive_cutoff in false_positive_distance_cutoffs:
                            if min_distance_to_sim_mutation < false_positive_cutoff:
                                sample_caller_false_positive_count[false_positive_cutoff][sample_name_for_output][caller_name] += called_mutations[caller_name][sample_name][mutation].sim_mut_count


    out_string = ""

    sample_output_names_arr = sorted(sample_output_names.keys())
    
    for sample_name in sample_output_names_arr:
        out_string += sample_name + "\tsim\t" + str(sim_count_by_sample[sample_name]) + "\t"
        for caller_name in caller_names:
            out_string += str(sample_caller_exact_count[sample_name][caller_name]) + "\t"
        out_string += "\n"

        out_string += sample_name + "\tnearby\t" + str(sim_count_by_sample[sample_name]) + "\t"
        for caller_name in caller_names:
            out_string += str(sample_caller_nearby_count[sample_name][caller_name]) + "\t"
        out_string += "\n"

        for false_positive_cutoff in false_positive_distance_cutoffs:
            out_string += sample_name + "\tfalse-positive " + str(false_positive_cutoff) + "\t0\t"
            for caller_name in caller_names:
                out_string += str(sample_caller_false_positive_count[false_positive_cutoff][sample_name][caller_name]) + "\t"
            out_string += "\n"

    out_string = head_string + "\n" + out_string
    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)
    return out_string

def print_recovery_mutations_by_sample(called_mutations,simulated_mutations,callers,fileName=None,mutation_types_to_report=None,nearby_distance_cutoff=10,false_positive_distance_cutoffs=[1000],collapse_replicates=False, require_passing=False):
    """For all callers, prints read counts for simulated mutations and read counts for recovered mutations for each sample
    if mutation is within distance_cutoff of a simulated mutation, it is considered recovered_nearby
    called_mutations: mutations called by callers in form of called_mutations[caller_name][sample_name][mutation] = Mutation object
    simulated_mutations: mutations simulated in form of simulated_mutations[sample_name][mutation] = Mutation object
    callers: list of callers
    fileName: If given, results are written to the file, otherwise they are printed to stdOut"
    mutation_types_to_report: list of mutation types to include (e.g. ['I', 'D'])
    nearby_distance_cutoff: distance from simulated mutation for a mutation to be considered recovered_nearby
    false_positive_distance_cutoffs: list of distances from simulated mutation for a mutation to be considered a false positive (e.g. reporting false-positives (subject to require_passing and mutation_types_to_report) called within 100bp, 1000bp, etc.)
    collapse_replicates: if true, mutation counts will be collapsed by replicate
    require_passing: if true, only mutations passing caller filters will be included
    """

    max_false_positive_cutoff = max(false_positive_distance_cutoffs)

    sim_mutations_list_by_sample = defaultdict(list)
    #first, aggregate all simulated mutations
    sim_count_by_sample = defaultdict(int)
    for sample_name in simulated_mutations:
        sample_name_for_output = sample_name
        if collapse_replicates:
            sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
        for mutation in simulated_mutations[sample_name]:
            if mutation not in sim_mutations_list_by_sample[sample_name]:
                sim_mutations_list_by_sample[sample_name].append(mutation)
                sim_count_by_sample[sample_name_for_output] += 1


    #next iterate through mutations 
    head_string = "SAMPLE\tCLASS\tSIMULATED"
    caller_names = []
    sample_names = {}
    sample_output_names = {}
    sample_caller_exact_count = defaultdict(lambda: defaultdict(int))
    sample_caller_nearby_count = defaultdict(lambda: defaultdict(int))
    sample_caller_false_positive_count = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for caller in callers:
        caller_name = caller.get_name()
        caller_names.append(caller_name)
        head_string += "\t"+caller_name
        for sample_name in called_mutations[caller_name]:
            sample_names[sample_name] = 1
            sample_name_for_output = sample_name
            if collapse_replicates:
                sample_name_for_output = '_'.join(sample_name.split("_")[:-1])
            sample_output_names[sample_name_for_output] = 1
            for mutation in called_mutations[caller_name][sample_name]:
                this_mutation_loc, this_mutation_type, this_mutation_info = mutation.split(" ")
                if mutation_types_to_report is not None and this_mutation_type not in mutation_types_to_report:
                    continue
                if require_passing and not called_mutations[caller_name][sample_name][mutation].passing_filter:
                    continue
                found_mutation_as_simulated = False
                min_distance_to_sim_mutation = max_false_positive_cutoff * 2
                for sim_mutation in sim_mutations_list_by_sample[sample_name]:
                    sim_mutation_loc, sim_mutation_type, sim_mutation_info = sim_mutation.split(" ")
                    this_dist_to_sim_mutation = abs(int(this_mutation_loc) - int(sim_mutation_loc))
                    if this_dist_to_sim_mutation < min_distance_to_sim_mutation:
                        min_distance_to_sim_mutation = this_dist_to_sim_mutation
                    if sim_mutation == mutation:
                        sample_caller_exact_count[sample_name_for_output][caller_name] += 1
                        found_mutation_as_simulated = True
                    elif sim_mutation_type == this_mutation_type:
                        if this_dist_to_sim_mutation <= nearby_distance_cutoff:
                            sample_caller_nearby_count[sample_name_for_output][caller_name] += 1
                            found_mutation_as_simulated = True
                        elif this_mutation_type == "D":
                            if int(this_mutation_loc) <= int(sim_mutation_loc) and int(this_mutation_loc) + int(this_mutation_info) >= int(sim_mutation_loc): #if deletion overlaps simulated deletion
                                sample_caller_nearby_count[sample_name_for_output][caller_name] += 1
                                found_mutation_as_simulated = True

                    if not found_mutation_as_simulated:
                        for false_positive_cutoff in false_positive_distance_cutoffs:
                            if min_distance_to_sim_mutation < false_positive_cutoff:
                                sample_caller_false_positive_count[false_positive_cutoff][sample_name_for_output][caller_name] += 1


    out_string = ""

    sample_output_names_arr = sorted(sample_output_names.keys())
    
    for sample_name in sample_output_names_arr:
        out_string += sample_name + "\tsim\t" + str(sim_count_by_sample[sample_name]) + "\t"
        for caller_name in caller_names:
            out_string += str(sample_caller_exact_count[sample_name][caller_name]) + "\t"
        out_string += "\n"

        out_string += sample_name + "\tnearby\t" + str(sim_count_by_sample[sample_name]) + "\t"
        for caller_name in caller_names:
            out_string += str(sample_caller_nearby_count[sample_name][caller_name]) + "\t"
        out_string += "\n"

        for false_positive_cutoff in false_positive_distance_cutoffs:
            out_string += sample_name + "\tfalse-positive " + str(false_positive_cutoff) + "\t0\t"
            for caller_name in caller_names:
                out_string += str(sample_caller_false_positive_count[false_positive_cutoff][sample_name][caller_name]) + "\t"
            out_string += "\n"

    out_string = head_string + "\n" + out_string
    if fileName is not None:
        fh = open(fileName, "w")
        fh.write(out_string)
        fh.close
    else:
        print(out_string)
    return out_string
