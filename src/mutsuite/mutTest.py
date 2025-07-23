from runReplaceReads import runReplaceReads
from runCallersHelpers import aggregate_mutations, print_mutations, print_aggregate_mutations, print_all_aggregate_mutations, print_all_mutations_by_sample, print_recovery_reads_by_sample, print_recovery_mutations_by_sample


from callers.simCaller import SimCaller
import concurrent.futures
import configparser
import glob
import os
import sys
from subprocess import call

def simulateSamples(config):
    maxReps = config.getint('Simulation','reps')
    reps = range(0,maxReps)

    # first simulate samples by replacing reads
    workFolder = config.get('runningSettings','workFolder')
    samples, simulateReadsCommands = runReplaceReads(config, workFolder)
    if len(simulateReadsCommands) > 0:
        print("got " + str(len(simulateReadsCommands)) + " commands to simulate samples")
        for command in simulateReadsCommands:
            call(command, shell=True)
    else:
        print('Finished simulating samples')

    return samples

def call_command(command):
    print("Running command " + command)
    return call(command,shell=True)

def runCallers(config, callers, samples):
    # read simulated indels
    simCaller = SimCaller(config)
    simulated_mutations = {}
    for sample_name, simulated_bam, simulated_control in samples:
        sample_simulated_indels, sample_simulated_mutations = simCaller.get_results(simulated_bam, simulated_control)
        simulated_mutations[sample_name] = sample_simulated_mutations
    aggregated_simulated_mutations = aggregate_mutations(simulated_mutations, passing_filter=True)


    n_processes_for_callers = config.getint('Callers','synchronousCallers')

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_processes_for_callers) as executor:
        checked_times = 0
        while(True):
            commands = []
            for caller in callers:
                caller_name = caller.get_name()
                for name, simulated_bam, simulated_control in samples:
                    this_command, is_finished = caller.run_caller(simulated_bam,simulated_control)
                    if not is_finished:
                        commands.append(this_command)

            if len(commands) > 0:
                print ("Got " + str(len(commands)) + " commands to run callers")

            # run commands in parallel
            caller_results = executor.map(call_command, commands)
            for result in caller_results:
                pass

            #check to see if they are finished
            all_are_finished = True
            for caller in callers:
                for name, simulated_bam, simulated_control in samples:
                    this_command, is_finished = caller.run_caller(simulated_bam,simulated_control)
                    if not is_finished:
                        all_are_finished = False
            if all_are_finished:
                break

            checked_times += 1

            if checked_times > 3:
                raise Exception('Tried to run caller commands 3 times and they still did not finish')

            if len(commands) == 0:
                break

            print('Finished running ' + caller_name + ' commands')




    all_mutations = {} #contains information for each simulated sample
    all_aggregated_indels = {} #aggregates information for each indel type across each sample

    for caller in callers:
        caller_name = caller.get_name()
        print('Running caller ' + str(caller_name))
        caller_mutations = {}
        for name, simulated_bam, simulated_control in samples:
            this_command, is_finished = caller.run_caller(simulated_bam,simulated_control)
            if not is_finished:
                caller_mutations[name] = None
            else:
                caller_passing_indels, caller_all_mutations = caller.get_results(simulated_bam,simulated_control)
                caller_mutations[name] = caller_all_mutations
        all_mutations[caller_name] = caller_mutations
        caller_aggregated_indels = aggregate_mutations(caller_mutations, passing_filter=True)
        all_aggregated_indels[caller_name] = caller_aggregated_indels
        print("\tFinished reading results from calling indels with " + caller_name)


    return all_aggregated_indels, aggregated_simulated_mutations, all_mutations, simulated_mutations


if __name__ == "__main__":
    #parse settings file
    if len(sys.argv) < 2:
        raise Exception("Couldn't find settings file")
    settingsFile = sys.argv[1]
    if settingsFile == "":
        settingsFile = "simSettings.txt"
    if not os.path.isfile(settingsFile):
        raise Exception("Couldn't find settings file")
    workFolder = os.path.abspath(settingsFile) + ".output/"
    if not os.path.isdir(workFolder):
        os.makedirs(workFolder)
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    config.read(settingsFile)
    
    config['runningSettings'] = {}
    config['runningSettings']['workFolder'] = workFolder

    samples = simulateSamples(config)
    samples_file = settingsFile+".samples.txt" # what percent of mutations were recovered
    with open(samples_file,'w') as s_out:
        for sample_name, simulated_bam, simulated_control in samples:
            s_out.write('\t'.join([sample_name, simulated_bam, simulated_control]) + '\n')

    useAllCallers = False
    callerNames = config.get('Callers','CallerList',fallback=None)
    if callerNames is None or 'all' in callerNames.lower():
        useAllCallers = True
    else:
        callerNames = [x.strip() for x in config.get('Callers','CallerList').lower().split(",")]

    callers = []
    # create each caller and append it to to the list of callers

    if useAllCallers or 'pindel' in callerNames:
        from callers.pindelCaller import PindelCaller #additional callers can be implemented and imported here
        pindelCaller = PindelCaller(config)
        callers.append(pindelCaller)
    if useAllCallers or 'lofreq' in callerNames:
        from callers.lofreqCaller import LofreqCaller
        lofreqCaller = LofreqCaller(config)
        callers.append(lofreqCaller)
    if useAllCallers or 'vardict' in callerNames:
        from callers.vardictCaller import VardictCaller
        vardictCaller = VardictCaller(config)
        callers.append(vardictCaller)
    if useAllCallers or 'somaticsniper' in callerNames:
        from callers.somaticSniperCaller import SomaticSniperCaller
        somaticsniperCaller = SomaticSniperCaller(config)
        callers.append(somaticsniperCaller)
    if useAllCallers or 'varscan' in callerNames:
        from callers.varscanCaller import VarscanCaller
        varscanCaller = VarscanCaller(config)
        callers.append(varscanCaller)
    if useAllCallers or 'mutect' in callerNames:
        from callers.mutectCaller import MutectCaller
        mutectCaller = MutectCaller(config)
        callers.append(mutectCaller)
    if useAllCallers or 'strelka' in callerNames:    
        from callers.strelkaCaller import StrelkaCaller
        strelkaCaller = StrelkaCaller(config)
        callers.append(strelkaCaller)
    if useAllCallers or 'strelkaTier2' in callerNames:
        from callers.strelkaTier2Caller import StrelkaTier2Caller
        strelkaTier2Caller = StrelkaTier2Caller(config)
        callers.append(strelkaTier2Caller)

    print('Running ' + str(len(callers)) + ' callers')

    all_aggregated_indels, aggregated_simulated_mutations, called_mutations, simulated_mutations = runCallers(config, callers, samples)

    FDR_distance_cutoffs = config.get('Simulation','FDRtolerance',fallback=None)
    if FDR_distance_cutoffs is not None:
        FDR_distance_cutoffs = [int(x.strip()) for x in FDR_distance_cutoffs.split(",")]
    else:
        FDR_distance_cutoffs = [1000]

    mutation_types_to_report = config.get('Reporting','MutationTypesToReport',fallback=None)
    if 'all' in mutation_types_to_report.lower():
        mutation_types_to_report = None
    if mutation_types_to_report is not None:
        mutation_types_to_report = [x.strip() for x in mutation_types_to_report.split(",")]

    collapse_replicates = config.getboolean('Reporting','CollapseReplicates',fallback=False)
    
    report_file = settingsFile+".report.txt"
    print_all_aggregate_mutations(all_aggregated_indels,aggregated_simulated_mutations,callers,fileName=report_file, mutation_types_to_report=mutation_types_to_report)
    report_by_name_unfiltered_file = settingsFile+".report_by_sample_unfiltered.txt"
    print_all_mutations_by_sample(called_mutations,simulated_mutations,callers,fileName=report_by_name_unfiltered_file,
                                  collapse_replicates=collapse_replicates, mutation_types_to_report=mutation_types_to_report, require_passing=False)
    report_by_name_file = settingsFile+".report_by_sample_passing_filters.txt"
    print_all_mutations_by_sample(called_mutations,simulated_mutations,callers,fileName=report_by_name_file,
                                  collapse_replicates=collapse_replicates, mutation_types_to_report=mutation_types_to_report, require_passing=True)
    recovery_report_reads_file = settingsFile+".read_recovery_report_by_sample.txt" # what number of reads were recovered
    print_recovery_reads_by_sample(called_mutations,simulated_mutations,callers,fileName=recovery_report_reads_file, false_positive_distance_cutoffs=FDR_distance_cutoffs, 
                                   collapse_replicates=collapse_replicates, mutation_types_to_report=mutation_types_to_report, require_passing=True)
    recovery_report_mutations_file = settingsFile+".mutation_recovery_report_by_sample.txt" # what percent of mutations were recovered
    print_recovery_mutations_by_sample(called_mutations,simulated_mutations,callers,fileName=recovery_report_mutations_file, false_positive_distance_cutoffs=FDR_distance_cutoffs, 
                                    collapse_replicates=collapse_replicates, mutation_types_to_report=mutation_types_to_report, require_passing=True)
    print('wrote report to ' + report_file)

    keep_caller_intermediate = config.getboolean('Callers','keepIntermediateFiles', fallback=False)
    if not keep_caller_intermediate:
        print('Cleaning callers')
        for caller in callers:
            print('Cleaning ' + str(caller.get_name()))
            for name, simulated_bam, simulated_control in samples:
                caller.clean(simulated_bam,simulated_control)

    print("Finished")

