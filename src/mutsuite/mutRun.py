from subprocess import check_call
import configparser
import os
import stat
import sys


def mutRun(config, work_folder, debug=False):
    """Runs a simulation replacing reads with edits

    Args:
        config (Configparser): configuration settings
        work_folder (str): path to work folder
    """

    depths = [int(i) for i in config.get('Simulation', 'depths').split(",")]
    pct_mut = None
    if config.has_option('Simulation', 'pctMut'):
        pct_mut = [float(i) for i in config.get('Simulation', 'pctMut').split(",")]
    count_mut = None
    if config.has_option('Simulation', 'countMut'):
        count_mut = [float(i) for i in config.get('Simulation', 'countMut').split(",")]
    qual_add_arr = [int(i) for i in config.get('Simulation', 'addQual').split(",")]
    max_reps = config.getint('Simulation', 'reps', fallback=1)
    swap_chr = config.get('Simulation', 'chr')
    swap_loc = config.getint('Simulation', 'loc')
    use_only_indel_flag = config.getboolean('Simulation', 'useOnlyIndels', fallback=False)
    use_only_crispresso_modified_flag = config.getboolean('Simulation', 'useOnlyCRISPRessoModifiedReads', fallback=False)
    only_insert_first_indel_flag = config.getboolean('Simulation', 'onlyInsertFirstIndel', fallback=False)
    reference = config.get('Simulation', 'reference')
    use_bowtie2 = config.getboolean('Simulation', 'useBowtie2', fallback=False)
    source_bam = config.get('Simulation', 'sourceBam')
    altered_bam = config.get('Simulation', 'alteredBam', fallback=None)
    unaltered_bam = config.get('Simulation', 'unalteredBam') #downsampled from sourceBam, does not have to exist
    unaltered_namesorted_bam = config.get('Simulation', 'unalteredNamesortedBam') #downsampled from sourceBam, does not have to exist
    simulate_range = config.getint('Simulation', 'simulateRange', fallback=5000)  # bp around cut site to simulate reads for
    mut_sim_script = config.get('Simulation', 'mutSimScript', fallback="mutSim.py")
    keep_intermediate_files = config.getboolean('Simulation', 'keepIntermediateFiles', fallback=False)
    n_processes = config.getint('Simulation', 'n_processes', fallback=20) #processes/threads for alignment for each teask

    reps = range(0, max_reps)

    specific_mutation_arr = [None]
    if config.has_option('Simulation', 'onlyInsertSpecificMutation'):
        only_insert_specific_mutation_str = config.get('Simulation', 'onlyInsertSpecificMutation')
        if only_insert_specific_mutation_str.lower().strip() != "none":
            specific_mutation_arr = [i.strip() for i in only_insert_specific_mutation_str.split(",") if i.strip() != ""]
            for s in specific_mutation_arr:
                if 'I' not in s and 'D' not in s and 'S' not in s:
                    raise Exception("Specific mutation '" + s + "' cannot be parsed")

    use_only_indel_string = ""
    if use_only_indel_flag:
        use_only_indel_string = " --only_include_altered_with_indel "

    use_only_crispresso_modified_string = ""
    if use_only_crispresso_modified_flag:
        use_only_crispresso_modified_string = " --only_include_altered_crispresso_modified_indel "

    only_insert_first_indel_string = ""
    if only_insert_first_indel_flag:
        only_insert_first_indel_string = " --only_insert_first_indel "


    keep_intermediate_files_string = ""
    if keep_intermediate_files:
        keep_intermediate_files_string = " --keep_intermediate_files "

    if specific_mutation_arr == [None] and not os.path.isfile(altered_bam):
        raise Exception("Cannot find altered bam '" + altered_bam + "''")

    if not os.path.isfile(reference):
        raise Exception("Cannot find reference '" + reference + "''")

    processed_bams_created = False
    if os.path.isfile(unaltered_bam) and os.path.isfile(unaltered_namesorted_bam):
        processed_bams_created = True

    simulate_distance = simulate_range * 2
    simulate_distance_suffix = "b"
    if simulate_distance > 1000:
        simulate_distance = int(simulate_distance / 1000)
        simulate_distance_suffix = "kb"
    if simulate_distance > 1000:
        simulate_distance = int(simulate_distance / 1000)
        simulate_distance_suffix = "mb"

    if not processed_bams_created:
        if debug:
            print('Creating bam files around cut site for simulation')
        files_to_delete = []
        # check whether source bam exists (this is really big, so sometimes it is deleted)
        if not os.path.isfile(source_bam):
            raise Exception("Cannot find source bam '" + source_bam + "''")
        sim_start = swap_loc - simulate_range
        sim_end = swap_loc + simulate_range
        temp_unaltered = os.path.join(work_folder, str(simulate_distance) + simulate_distance_suffix + "AroundCutsite.fromUnaltered.bam")
        # samtools view -Sb ../../NA12878/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam chr2:72928870-72938870 > 1kbAroundCutsite.novoalign.bam
        cmd = f"samtools view -Sb {source_bam} {swap_chr}:{sim_start}-{sim_end} > {temp_unaltered}"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        files_to_delete.append(temp_unaltered)

        temp_sort = temp_unaltered + ".qsort.bam"
        # samtools sort -n 1kbAroundCutsite.novoalign.bam 1kbAroundCutsite.novoalign.bam.qsort
        cmd = f"samtools sort -n {temp_unaltered} -o {temp_sort}"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        files_to_delete.append(temp_sort)

        # bedtools bamtofastq -i 1kbAroundCutsite.novoalign.bam.qsort.bam -fq 1kbAroundCutsite.novoalign.bam.R1.fastq -fq2 1kbAroundCutsite.novoalign.bam.R2.fastq
        cmd = f"bedtools bamtofastq -i {temp_sort} -fq {temp_sort}.R1.fastq -fq2 {temp_sort}.R2.fastq"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        files_to_delete.append(temp_sort+".R1.fastq")
        files_to_delete.append(temp_sort+".R2.fastq")


        if use_bowtie2:
            bowtie_reference = reference.replace(".fasta", "").replace(".fa", "")
            aligned_sam = os.path.join(work_folder, str(simulate_distance) + simulate_distance_suffix + "AroundCutsite.Bowtie2.sam")
            cmd = f"bowtie2 -p {n_processes} -x {bowtie_reference} -1 {temp_sort}.R1.fastq -2 {temp_sort}.R2.fastq > {aligned_sam}"
        else:
            aligned_sam = os.path.join(work_folder, str(simulate_distance) + simulate_distance_suffix + "AroundCutsite.BWA.sam")
            cmd = f"bwa mem -t {n_processes} {reference} {temp_sort}.R1.fastq {temp_sort}.R2.fastq > {aligned_sam}"
        # bwa mem /seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta 1kbAroundCutsite.novoalign.bam.R1.fastq 1kbAroundCutsite.novoalign.bam.R2.fastq > 1kbAroundCutsite.BWA.sam
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)
        files_to_delete.append(aligned_sam)

        # filter out secondary alignments (-F 2048)
        # samtools view -F 2048 -Sb 1kbAroundCutsite.BWA.sam | samtools sort - -o 1kbAroundCutsite.BWA.bam && samtools index 1kbAroundCutsite.BWA.bam
        cmd = f"samtools view -F 2048 -Sb {aligned_sam} | samtools sort - -o {unaltered_bam} && samtools index {unaltered_bam}"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        # create name-sorted bam for operating on read pairs
        # samtools view -F 2048 -Sb 1kbAroundCutsite.BWA.sam | samtools sort -n -o 1kbAroundCutsite.BWA.nameSort.bam
        cmd = f"samtools view -F 2048 -Sb {aligned_sam} | samtools sort -n -o {unaltered_namesorted_bam}"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        # java -Xmx8g -jar /seq/software/picard/current.20170216/bin/picard.jar CollectInsertSizeMetrics I=1kbAroundCutsite.BWA.bam O=1kbAroundCutsite.BWA.bam.insertSizeMetrics H=1kbAroundCutsite.BWA.bam.insertSizeMetrics.pdf
        cmd = f"picard CollectInsertSizeMetrics -I {unaltered_bam} -O {unaltered_bam}.insertSizeMetrics -H {unaltered_bam}.insertSizeMetrics.pdf"
        if debug:
            print("Calling: "+cmd)
        check_call(cmd, shell=True)

        if not keep_intermediate_files:
            for file in files_to_delete:
                check_call(["rm",file])


    if count_mut is not None:
        mut_string = " --mut_count "
        mut_pct_arr = count_mut
    else:
        mut_string = " --mut_freq "
        mut_pct_arr = pct_mut

    command_array = []
    samples = []
    for d in depths:
        for p in mut_pct_arr:
            for q in qual_add_arr:
                for s in specific_mutation_arr:
                    for r in reps:
                        s_string = ""
                        if s is not None:
                            only_insert_specific_mutation = " --only_insert_specific_mutation " + s + " "
                            s_string = "_s" + s

                        name = f"sim_d{d}_p{p}_q{q}{s_string}_r{r}"
                        outfile_root = os.path.join(work_folder, name)
                        simulated_bam = outfile_root + ".bam"
                        simulated_control = outfile_root + ".ctl.bam"
                        command = f"python {mut_sim_script} --downsample_number {d} {mut_string}{p} --qual_add {q} " + \
                                f"--mut_chr {swap_chr} --mut_loc {swap_loc} --reference {reference} " + \
                                f"--unaltered_bam {unaltered_bam} --unaltered_namesorted_bam {unaltered_namesorted_bam} --altered_bam {altered_bam} " + \
                                f"--output_root {outfile_root} {use_only_indel_string} {use_only_crispresso_modified_string} {only_insert_first_indel_string} " + \
                                f"{only_insert_specific_mutation} {keep_intermediate_files_string} ; touch {outfile_root}.finished"
                        if (not os.path.isfile(outfile_root+".finished") or not os.path.isfile(simulated_bam)):
                            command_array.append(command)

                        samples.append((name, simulated_bam, simulated_control))

    if debug:
        print('Got ' + str(len(command_array)) + ' commands to simulate samples')
    return samples, command_array


if __name__ == "__main__":

    if len(sys.argv) < 2:
        raise Exception("Couldn't find settings file")
    settings_file = sys.argv[1]
    if settings_file == "":
        settings_file = "simSettings.txt"
    if not os.path.isfile(settings_file):
        raise Exception("Couldn't find settings file")

    work_folder = os.path.abspath(settings_file) + ".output/"
    if not os.path.isdir(work_folder):
        os.makedirs(work_folder)

    config = configparser.ConfigParser(inline_comment_prefixes=['#'])
    config.read(settings_file)
    print("Read setting file " + settings_file)

    samples, command_array = mutRun(config, work_folder)
    print('Got ' + str(len(command_array)) + ' commands to run')

    if len(command_array) > 0:
        script_name = settings_file + ".run.sh"
        q_commands_file = script_name + ".commands.txt"

        print('command array: ' + str(command_array))
        with open(q_commands_file,'w') as f:
            f.write("\n".join(command_array))

        print (f'{len(command_array)} commands printed to {q_commands_file}')
    else:
        print ("No Commands")