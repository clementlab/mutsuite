####
# Modifies a certain proportion of reads at a specified location
# Only reads with a string of 'M' (matches) running across the modification site are considered for swapping. Other reads are deleted
# quality scores are left as is
# some notes:
# 1) Reads that have the deletion too close to one of the ends are discarded
# 2) nucleotides are added to the end of the sequence. If the original sequence had soft clipping a the end, the added nucleotides are added as clipped
# 3) All reads are paired - otherwise they weren't included in bedtools bam2fastq
####
import pysam
from Bio import SeqIO
import re
import copy

cigar_explode_pattern = re.compile("(\d+)(\w)")


def explodeCigar(cigar_string):
    exploded_cigar = ""
    for (cigar_count, cigar_char) in re.findall(cigar_explode_pattern, cigar_string):
        exploded_cigar += cigar_char * int(cigar_count)
    return exploded_cigar


cigar_unexplode_pattern = re.compile(r'((\w)\2{0,})')


def unexplodeCigar(exploded_cigar_string):
    cigar = ""
    for (cigar_str, cigar_char) in re.findall(cigar_unexplode_pattern, exploded_cigar_string):
        cigar += str(len(cigar_str)) + cigar_char
    return cigar

def unexplodeMDZ(exploded_mdz_arr):
    """
        :param exploded_mdz_arr: array of mdz operations (= for matched bases, "A" for sub A, "^A" e.g. for deleted A)
        :return: mdz string
    """
    mdz = ""
    match_count_so_far = 0
    deleted_char_str = ""

    for el in exploded_mdz_arr:
        if el == "=":  # match base
            if deleted_char_str != "":
                mdz += "^" + deleted_char_str
                deleted_char_str = ""
            match_count_so_far += 1
        elif el.startswith('^'):  # deletion base
            if match_count_so_far != 0:
                mdz += str(match_count_so_far)
                match_count_so_far = 0
            deleted_char_str += el[1]
        else:  # substitution base
            if deleted_char_str != "":
                mdz += "^" + deleted_char_str
                deleted_char_str = ""
            if match_count_so_far != 0:
                mdz += str(match_count_so_far)
                match_count_so_far = 0
            else:  # if transitioning from a previous deletion or substitution
                if len(mdz) > 0:
                    mdz += "0"
            mdz += el

    if deleted_char_str != "":
        mdz += "^" + deleted_char_str
        deleted_char_str = ""
    if match_count_so_far != 0:
        mdz += str(match_count_so_far)
        match_count_so_far = 0
        
    return mdz


def unexplodeReplacementInfo(replacement_info_els):
    """
    Collapse adjacent replacements into a single string (e.g. multiple insertions and multiple deletions)
    e.g. [(10,'D',1), (11,'D',1), (12,'I','A'), (12,'I','C'), (15,'I','C')] -> ['10D2', '12IAC', '15IC']
    :param replacement_info_els: array of tuples (loc, op, info) where 
        loc is the location of the operation, 
        op is one of S, D, I, R (substitution, deletion, insertion, restoration (of a deleted base)) and 
        info is the information for the operation (e.g. the base for substitution, the number of bases for deletion, and the string of bases for insertion)
    :return: array of strings of the form loc+op+info
    """

    replacement_info = []

    first_deleted_loc = 0
    deleted_count_so_far = 0

    first_inserted_loc = 0
    inserted_char_str = ""

    last_replacement_info_loc = -1 # merge infos that are adjacent

    for (loc, op, info) in replacement_info_els:
        # check for unrecorded insertions and deletions if the next operation is not the same or not adjacent
        if deleted_count_so_far > 0 and (op != "D" or loc != last_replacement_info_loc+1):
            replacement_info.append(str(first_deleted_loc)+"D"+str(deleted_count_so_far))
            deleted_count_so_far = 0
        if inserted_char_str != "" and (op != "I" or loc != last_replacement_info_loc):
            replacement_info.append(str(first_inserted_loc)+"I"+inserted_char_str)
            inserted_char_str = ""

        last_replacement_info_loc = loc

        if op == "S":  # substitution
            replacement_info.append(str(loc)+op+info)
        elif op == "D":
            if deleted_count_so_far == 0:
                first_deleted_loc = loc
                deleted_count_so_far = 1
            else:
                deleted_count_so_far += 1
        elif op == "I":
            if inserted_char_str == "":
                inserted_char_str = info
                first_inserted_loc = loc
            else:
                inserted_char_str += info
        elif op == "R":  # restore/undelete / revert/uninsert
            replacement_info.append(str(loc)+op+info)
        else:
            raise Exception("Invalid replacement operation: " + str(op) + ' in ' + str(replacement_info_els))

    #finally check one last time for unrecorded insertions and deletions
    if deleted_count_so_far > 0:
        replacement_info.append(str(first_deleted_loc)+"D"+str(deleted_count_so_far))
    if inserted_char_str != "":
        replacement_info.append(str(first_inserted_loc)+"I"+inserted_char_str)
        
    return replacement_info

def replaceRead(_wt_read, _new_read, _genome_seq, _genome_start, _genome_end=None, debug=False):
    """ Replaces bases from _new_read into _wt_read
        If _new_read and _wt_read do not overlap, wt_read is returned
        If _wt_read is unmapped, wt_read is returned

        Returns:
            pysam.AlignedSegment: _wt_read with modifications from _new_read
            replacement_info: array of replacement positions and details (these replacements were made from _new_read to _wt_read, e.g. if the _wt_read has a substitution which is restored to the genomic sequence by the _new_read that will appear here)
            modification_info: array of modification positions and details (these are modifications as compared to the genome, e.g. if the wt read has a substitution, that will appear here)

    """
    wt_read = copy.deepcopy(_wt_read)
    new_read = copy.deepcopy(_new_read)

    if _genome_end is None:
        _genome_end = _genome_start + len(_genome_seq)

    if wt_read.is_unmapped or new_read.is_unmapped:
        return wt_read, [], []

    # print("wt: " + str(wt_read))
    # print("new: " + str(new_read))

    # wt is unaltered read
    wt_start = wt_read.reference_start
    wt_start_with_clipping = wt_read.reference_start - wt_read.query_alignment_start  # includes soft-clipped bases
    wt_end = wt_read.reference_end - 1
    wt_end_with_clipping = wt_read.reference_end + (wt_read.query_length - wt_read.query_alignment_end)

    # new is altered read from which we'll pull mutations
    new_start = new_read.reference_start
    #new_start_with_clipping = new_read.reference_start - new_read.query_alignment_start
    new_end = new_read.reference_end - 1
    #new_end_with_clipping = new_read.reference_end + (new_read.query_length - new_read.query_alignment_end)

    #print('wt_start: ' + str(wt_start) + ' wt_end: ' + str(wt_end) + ' new_start: ' + str(new_start) + ' new_end: ' + str(new_end))


    if wt_start > new_end:
        return _wt_read, [], get_read_modifications(wt_read, _genome_seq, _genome_start, _genome_end)
    if wt_end < new_start:
        return _wt_read, [], get_read_modifications(wt_read, _genome_seq, _genome_start, _genome_end)

    # print ("wt_start: " + str(wt_start) + " wt_end " + str(wt_end))
    # print ("new_start: " + str(new_start) + " new_end " + str(new_end))

    explode_wt_cigar = explodeCigar(str(wt_read.cigarstring))
    explode_new_cigar = explodeCigar(str(new_read.cigarstring))

    final_read = ""
    final_cigar = ""
    final_mdz_els = []  # keep this as an arr because some operations are multiple characters (e.g. ^G)
    final_nm_count = 0  # keep track of NM for number of mismatches
    replacement_info_els = []  # array of replacement els to be collapsed - some operations are multiple bases
    modification_info_els = []  # array of changes from the genome (not just those inserted by the new read) 

    if debug:
        print ("wt: " + wt_read.query_sequence)
        print ("wt: " + wt_read.cigarstring)
        print ("wt: " + explode_wt_cigar)
        print ("wt: " + str(wt_read.query_qualities))
        print ("new: " + new_read.query_sequence)
        print ("new: " + new_read.cigarstring)
        print ("new: " + explode_new_cigar)
        print ("new: " + str(new_read.query_qualities))

    wt_read_ind = 0
    wt_cigar_ind = 0
    final_genomic_loc = wt_start
    # add soft/hard clipped bases to final read
    while explode_wt_cigar[wt_cigar_ind] == "H":
        final_cigar += "H"
        wt_cigar_ind += 1

    while explode_wt_cigar[wt_cigar_ind] == "S":
        final_cigar += "S"
        final_read += wt_read.query_sequence[wt_read_ind]
        wt_cigar_ind += 1
        wt_read_ind += 1

    if debug:
        print('216 wt_read_ind: ' + str(wt_read_ind) + ' wt_cigar_ind: ' + str(wt_cigar_ind) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' final_genomic_loc: ' + str(final_genomic_loc))

    new_read_ind = 0
    new_cigar_ind = 0
    # skip soft/hard clipped bases in new
    while explode_new_cigar[new_cigar_ind] == "H":
        new_cigar_ind += 1

    while explode_new_cigar[new_cigar_ind] == "S":
        new_cigar_ind += 1
        new_read_ind += 1

    # if WT read starts to the left of the NEW read, just take bases from the WT read
    while final_genomic_loc < new_start:
        if debug:
            print('231 wt_read_ind: ' + str(wt_read_ind) + ' wt_cigar_ind: ' + str(wt_cigar_ind) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' final_genomic_loc: ' + str(final_genomic_loc) + " new_start: " + str(new_start))
            print('final_genomic_loc: ' + str(final_genomic_loc) + '_genome_start: ' + str(_genome_start))
            print('genome at pos: ' + str(_genome_seq[final_genomic_loc - _genome_start]))
        if explode_wt_cigar[wt_cigar_ind] == "D":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_mdz_els.append("^" + _genome_seq[final_genomic_loc - _genome_start])
            modification_info_els.append((final_genomic_loc+1, "D", 1))
            wt_cigar_ind += 1
            final_genomic_loc += 1
        elif explode_wt_cigar[wt_cigar_ind] == "I":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            modification_info_els.append((final_genomic_loc, "I", wt_read.query_sequence[wt_read_ind]))
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "S":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "H":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            wt_cigar_ind += 1
        else:
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            base_from_wt = wt_read.query_sequence[wt_read_ind]
            base_in_genome = _genome_seq[final_genomic_loc - _genome_start]
            final_read += base_from_wt
            if base_from_wt.upper() != base_in_genome.upper():
                final_mdz_els.append(base_from_wt)
                final_nm_count += 1
                modification_info_els.append((final_genomic_loc + 1, "S", base_from_wt))
            else:
                final_mdz_els.append("=")
            wt_cigar_ind += 1
            wt_read_ind += 1
            final_genomic_loc += 1

    # print('added wt bases that didn\'t overlap with new read')
    # print('final_read: ' + str(final_read) + ' final cigar: ' + str(final_cigar))
    # print('New start should be >= final_genomic_loc: ' + str(new_start) + ' >= ' + str(final_genomic_loc))
    # print('235 final_genomic_loc: ' + str(final_genomic_loc) + ' new_start ' + str(new_start))

    # print("236 wt len is : " + str(len(wt_read.query_sequence)) + " final_read is: " + final_read + " wt_start " + str(wt_start) + " wt_end " + str(wt_end) + " new start: " + str(new_start) + " new_end: " + str(new_end))

    new_temp_genomic_loc = new_start
    # print('new_temp_genomic_loc: ' + str(new_temp_genomic_loc) + ' final_genomic_loc: ' + str(final_genomic_loc))
    # if wt read starts in middle of new, catch the new_read pointers up to where we start pulling from the new_read sequence
    # in simple cases, final_genomic_loc will be the wt_start
    while new_temp_genomic_loc < final_genomic_loc:
        #print('245 catching up new_temp_genomic_loc: ' + str(new_temp_genomic_loc) + ' final_genomic_loc: ' + str(final_genomic_loc) + ' explode_new_cigar[' + str(new_cigar_ind) + ']:' + explode_new_cigar[new_cigar_ind])
        if str(explode_new_cigar[new_cigar_ind]) == "H":
            new_cigar_ind += 1
        elif str(explode_new_cigar[new_cigar_ind]) == "S":
            new_cigar_ind += 1
            new_read_ind += 1
        elif str(explode_new_cigar[new_cigar_ind]) == "D":
            new_cigar_ind += 1
            new_temp_genomic_loc += 1
        elif str(explode_new_cigar[new_cigar_ind]) == "I":
            new_cigar_ind += 1
            new_read_ind += 1
        else:  # M
            new_cigar_ind += 1
            new_read_ind += 1
            new_temp_genomic_loc += 1

    # at this point the final_genomic_loc should be where new and wt overlap
    # print('Ready to start adding matching bases')
    # print('286 final_genomic_loc: ' + str(final_genomic_loc))

    # print("checking if first cigar is deletion (explode_new_cigar["+str(new_cigar_ind)+"]:" + explode_new_cigar[new_cigar_ind]+")")
    # print("finalCigar: " + final_cigar)
    # print("finalRead: " + final_read)
    # print('orig wt start: ' + str(wt_start))

    # print("newSeq: " + new_read.query_sequence)
    # print("explodenewcigar: " + explode_new_cigar)
    # print('checking first overlapping position is deletion in wt: ' + str(explode_new_cigar[new_cigar_ind]) + ' ' + str(new_cigar_ind))
    # if first position in wt is a deletion in mut, add the bases from wt to final
    if explode_new_cigar[new_cigar_ind] == "D":
        # print('first overlapping position in wt is a deletion in mut')
        # first determine how long the deletion string was
        new_cigar_del_start = new_cigar_ind
        new_cigar_del_end = new_cigar_ind
        genomic_loc_del_end = final_genomic_loc
        while new_cigar_del_end < len(explode_new_cigar):
            if explode_new_cigar[new_cigar_del_end] != "D":
                break
            if explode_new_cigar[new_cigar_del_end] not in ['I','S','H','P']:  # these cigar actions don't consume the reference
                genomic_loc_del_end += 1
            new_cigar_del_end += 1

        # find start of del by backtracking through mut cigar
        genomic_loc_del_start = final_genomic_loc
        while new_cigar_del_start > 0:
            if explode_new_cigar[new_cigar_del_start] != "D":
                break
            if explode_new_cigar[new_cigar_del_start] not in ['I','S','H','P']:  # these cigar actions don't consume the reference
                genomic_loc_del_start -= 1
            new_cigar_del_start -= 1
        genomic_loc_del_end -= 1 # deletion ends the base before non-D cigar
        new_cigar_del_end -= 1
        genomic_loc_del_start += 1  # deletion starts the next base after non-D cigar
        new_cigar_del_start += 1
        deletion_length = new_cigar_del_end - new_cigar_del_start + 1
        # print('del start was ' + str(deletion_length) + 'bp from ' + str(new_cigar_del_start) + ' ( ' + str(genomic_loc_del_start) + 'bp) and end was ' + str(new_cigar_del_end) + ' at genomic loc ' + str(genomic_loc_del_end))

        # insert bases from genome into final before deletion
        #   find cigar loc where deletion ends in wt, then insert bases/changes from wt before that into final
        #   There is some ambiguity about how this should happen. Here, if a wt read starts in the middle of a deletion in mut/new
        #     if the wt cigar has M we insert the bases from genome into the final read before the deletion, then continue with the deletion
        #     This could potentially ignore indels in mut/new that occur before this deletion, but technically those indels don't overlap this wt read...
        temp_wt_cigar_ind = 0  # index in cigar
        temp_wt_read_ind = 0  # index in read
        temp_genome_loc = final_genomic_loc
        num_bases_to_add_from_genome = 0
        while temp_genome_loc <= genomic_loc_del_end:
            if explode_wt_cigar[temp_wt_cigar_ind] in ['S']:  # these cigar actions consume the read but bases shouldn't be added from genome for S
                temp_wt_read_ind += 1
            if explode_wt_cigar[temp_wt_cigar_ind] in ['M', 'I', '=', 'X']:  # these cigar actions consume the read (and S)
                num_bases_to_add_from_genome += 1
                temp_wt_read_ind += 1
            if explode_wt_cigar[temp_wt_cigar_ind] in ['M', 'D', 'N', '=', 'X']:  # these cigar actions consume the reference
                temp_genome_loc += 1
            temp_wt_cigar_ind += 1
        # print('now temp_genome_loc ' + str(temp_genome_loc) + ' is > genomic_loc_del_end ' + str(genomic_loc_del_end))
        
        wt_read_ind = temp_wt_read_ind 
        wt_cigar_ind = temp_wt_cigar_ind
        new_cigar_ind = new_cigar_del_end + 1
        final_genomic_loc = genomic_loc_del_end + 1

        # print('364 final_cigar: ' + str(final_cigar))
        # print('Location in WT of end of deletion: temp_wt_cigar_ind: ' + str(temp_wt_cigar_ind) + \
        #     ' temp_wt_read_ind: ' + str(temp_wt_read_ind) + ' temp_genome_loc: ' + str(temp_genome_loc))
        
        # print("deletion length: " + str(deletion_length))

        bases_to_add_from_genome = _genome_seq[genomic_loc_del_start - num_bases_to_add_from_genome - _genome_start:genomic_loc_del_start - _genome_start]
        # print('num bases to add from genome: ' + str(num_bases_to_add_from_genome) +
        #       ' (' + str(genomic_loc_del_start) + ' - ' + str(num_bases_to_add_from_genome) + ' - ' + str(_genome_start) + ' :' +
        #       str(genomic_loc_del_start) + ' - ' + str(_genome_start) + ') : ' + str(bases_to_add_from_genome))

        final_cigar += "M"*num_bases_to_add_from_genome + "D"*deletion_length
        final_mdz_els = ["="]*num_bases_to_add_from_genome + ['^'+_genome_seq[x - _genome_start] for x in range(genomic_loc_del_start, genomic_loc_del_end+1)]  # this is an arr because we can't reverse multiple characters (e.g. ^A)
        # update start of read
        wt_read.reference_start = genomic_loc_del_start - num_bases_to_add_from_genome
        final_read += bases_to_add_from_genome
        for i in range(genomic_loc_del_start, genomic_loc_del_end+1):
            replacement_info_els.append((i+1, 'D', 1))
            modification_info_els.append((i+1, 'D', 1))

    # print('At 419')
    # print("finalCigar: " + final_cigar)
    # print("finalRead: " + final_read)
    # print('wt_read reference start: ' + str(wt_read.reference_start))
    # print('replacement info els: ' + str(replacement_info_els))

    # print("new_read_ind: " + str(new_read_ind) + " new_cigar_ind " + str(new_cigar_ind) + " final_genomic_loc: " + str(final_genomic_loc))
    # print("at 1: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind, new_cigar_ind) + \
    #  "wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind, wt_cigar_ind) + \
    #  "final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc, final_read, final_cigar))

    # print('wt read: ' + str(wt_read.query_sequence) + ' length ' + str(len(wt_read.query_sequence)))
    # print('getting overlapping positions')
    # overlapping locations, or where we still have bases in new_read
    tmp_wt_genomic_loc = final_genomic_loc  # this can sometimes be different than final_genomic_loc if wt has indels
    check_wt_for_replacement_info = True  # if we're in the region where wt and new overlap, check for sequence changes

    while len(final_read) < len(wt_read.query_sequence):
        # first, if we've reached the end of the new read, break
        if new_read_ind >= len(new_read.query_sequence):
            break
        # then, if the new read has additional bases, add them as appropriate
        if explode_new_cigar[new_cigar_ind] == "H":
            break
        if explode_new_cigar[new_cigar_ind] == "S":
            break
        elif explode_new_cigar[new_cigar_ind] == "D":
            final_cigar += "D"
            final_mdz_els.append('^'+_genome_seq[final_genomic_loc - _genome_start])
            replacement_info_els.append((final_genomic_loc+1, 'D', 1))
            modification_info_els.append((final_genomic_loc+1, 'D', 1))
            new_cigar_ind += 1
            final_genomic_loc += 1
        elif explode_new_cigar[new_cigar_ind] == "I":
            final_cigar += "I"
            replacement_info_els.append((final_genomic_loc, 'I', new_read.query_sequence[new_read_ind]))
            modification_info_els.append((final_genomic_loc, 'I', new_read.query_sequence[new_read_ind]))
            final_read += new_read.query_sequence[new_read_ind]
            new_cigar_ind += 1
            new_read_ind += 1
        else:  # match
            final_cigar += explode_new_cigar[new_cigar_ind]
            base_to_add = new_read.query_sequence[new_read_ind]
            final_read += base_to_add
            # print('378 adding base to read: ' + str(base_to_add) + '(now ' + str(final_read) + ')  at loc: ' + str(final_genomic_loc) + ' len: ' + str(len(final_read)) + ' len: ' + str(len(wt_read.query_sequence)))
            if base_to_add.upper() != _genome_seq[final_genomic_loc - _genome_start].upper():
                final_mdz_els.append(base_to_add)
                final_nm_count += 1
                modification_info_els.append((final_genomic_loc+1, 'S', base_to_add))
                # print('adding modification info el: ' + str((final_genomic_loc+1, 'S', base_to_add)))
            else:
                final_mdz_els.append("=")
            # print('new read loc: ' + str(new_read_ind) + ' len: ' + str(len(new_read.query_sequence)) + ' wt_read_ind: ' + str(wt_read_ind) + ' len: ' + str(len(wt_read.query_sequence)))
            if check_wt_for_replacement_info:
                # print('398 checking wt for replacement info: explode_wt_cigar[' + str(wt_cigar_ind) + ']: ' + str(explode_wt_cigar[wt_cigar_ind]) + ' base_to_add: ' + str(base_to_add) + ' wt_read.query_sequence['+str(wt_read_ind) + ']: ' + str(wt_read.query_sequence[wt_read_ind]) + ' wt_read.query_sequence: ' + str(wt_read.query_sequence))
                if explode_wt_cigar[wt_cigar_ind] == "S":
                    pass
                elif explode_wt_cigar[wt_cigar_ind] == "I":  # insertions will be restored below in the wt catchup loop
                    pass
                elif explode_wt_cigar[wt_cigar_ind] == "D":
                    replacement_info_els.append((final_genomic_loc+1, 'R', _genome_seq[final_genomic_loc - _genome_start]))
                elif base_to_add.upper() != wt_read.query_sequence[wt_read_ind].upper():
                    replacement_info_els.append((final_genomic_loc + 1, 'S', base_to_add))
            new_cigar_ind += 1
            new_read_ind += 1
            final_genomic_loc += 1

        while tmp_wt_genomic_loc < final_genomic_loc:
            if wt_read_ind >= len(wt_read.query_sequence):
                check_wt_for_replacement_info = False
                break
            if explode_wt_cigar[wt_cigar_ind] in ['S', 'H']: # at end of a read
                check_wt_for_replacement_info = False
                break
            if explode_wt_cigar[wt_cigar_ind] == 'I':
                replacement_info_els.append((tmp_wt_genomic_loc - 1, 'R', wt_read.query_sequence[wt_read_ind]))
            if explode_wt_cigar[wt_cigar_ind] in ['M', 'I', '=', 'X']:  # these cigar actions consume the read
                wt_read_ind += 1
            if explode_wt_cigar[wt_cigar_ind] in ['M', 'D', 'N', '=', 'X']:  # these cigar actions consume the reference
                tmp_wt_genomic_loc += 1
            wt_cigar_ind += 1
        # if we've reached the end of the wt sequence
        if wt_read_ind >= len(wt_read.query_sequence):
            check_wt_for_replacement_info = False

    # print("469 len is : " + str(len(wt_read.query_sequence)) + " final is: " + final_read + " " + final_cigar + " wt_start " + str(wt_start) + " new start: " + str(new_start))
    # print("at 2: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind,new_cigar_ind)+\
    # 	"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
    # 	"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

    # pprint.pprint(str(wt_read))
    # refPos = wt_read.get_reference_positions(full_length=True)
    # print("pos: " + str(refPos))
    # pprint.pprint(str(new_read))
    # refPos = new_read.get_reference_positions(full_length=True)
    # print("pos: " + str(refPos))
    # print('Adding remaining bases from WT (if any)')
    # print('final read: ' + str(final_read))
    # print('final cigar: ' + str(final_cigar))
    # print('replacement info els: ' + str(replacement_info_els))
    # print("len final read: " + str(len(final_read)) + " wt len: " + str(len(wt_read.query_sequence)))
    # print("at 3: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind, new_cigar_ind)+\
    # 		"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind, wt_cigar_ind)+\
    # 		"final_genomic_loc: %s final_read: %s finalCigar: %s replacement_els %s\n" % (final_genomic_loc,final_read,final_cigar,str(replacement_info_els)))

    #print('490 final genomic_loc: ' + str(final_genomic_loc) + ' wt_end_with_clipping: ' + str(wt_end_with_clipping) + ' final read: ' + str(final_read) + ' len(final_read): ' + str(len(final_read)) + ' len(wt_read.query_sequence): ' + str(len(wt_read.query_sequence)))
    # finished adding from indices where bases overlap
    # add from wt 
    while final_genomic_loc < wt_end_with_clipping and len(final_read) < len(wt_read.query_sequence) and wt_cigar_ind < len(explode_wt_cigar):
        #print('493 final_genomic_loc: ' + str(final_genomic_loc) + ' wt_end: ' + str(wt_end) + ' len(final_read): ' + str(len(final_read)) + ' len(wt_read.query_sequence): ' + str(len(wt_read.query_sequence)))
        #print('wt_cigar_ind: ' + str(wt_cigar_ind) + " (len " + str(len(explode_wt_cigar)) + ') wt_read_ind: ' + str(wt_read_ind) + \
        #       ' final_genomic_loc: ' + str(final_genomic_loc) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' replacement_info_els: ' + str(replacement_info_els) + ' mdz: ' + str(final_mdz_els))
        #print('wt_cigar: ' + str(explode_wt_cigar))
#            print('cigar val: ' + str(explode_wt_cigar[wt_cigar_ind]))
        #print('wt read: ' + str(wt_read))
        #print('new_read: ' + str(new_read))
        if explode_wt_cigar[wt_cigar_ind] == "H":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            wt_cigar_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "S":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "D":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_mdz_els.append("^"+_genome_seq[final_genomic_loc - _genome_start])
            modification_info_els.append((final_genomic_loc+1, 'D', 1))
            wt_cigar_ind += 1
            final_genomic_loc += 1
        elif explode_wt_cigar[wt_cigar_ind] == "I":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            modification_info_els.append((final_genomic_loc, 'I', wt_read.query_sequence[wt_read_ind]))
            wt_cigar_ind += 1
            wt_read_ind += 1
        else:
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            base_to_add = wt_read.query_sequence[wt_read_ind]
            final_read += base_to_add
            if base_to_add.upper() != _genome_seq[final_genomic_loc - _genome_start].upper():
                final_mdz_els.append(base_to_add)
                final_nm_count += 1
                modification_info_els.append((final_genomic_loc+1, 'S', base_to_add))
            else:
                final_mdz_els.append("=")
            wt_cigar_ind += 1
            wt_read_ind += 1
            final_genomic_loc += 1

    # print("at 4: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind,new_cigar_ind)+\
    # 		"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
    # 		"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

    # finally, add from genome
    if final_cigar[-1] == "S": # if read ends with soft clipped bases, add soft-clipped bases (M can't come after S at end of cigar)
        while len(final_read) < len(wt_read.query_sequence):
            final_cigar += 'S'
            final_read += _genome_seq[final_genomic_loc - _genome_start]
            final_genomic_loc += 1

    while len(final_read) < len(wt_read.query_sequence):
        final_cigar += 'M'
        final_read += _genome_seq[final_genomic_loc - _genome_start]
        final_mdz_els.append("=")
        final_genomic_loc += 1
        
    # add trailing hardclipped bases (if any)
    temp_h_wt_pos = len(explode_wt_cigar) - 1
    while (explode_wt_cigar[temp_h_wt_pos] == "H"):
        final_cigar += "H"
        temp_h_wt_pos -= 1
    # print("at 5: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind,new_cigar_ind)+\
    # 	"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
    # 	"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

    # if re.search("14D",new_read.cigarstring):
    # 	print ("wt_start: " + str(wt_start) + " wt_end " + str(wt_end))
    # 	print ("new_start: " + str(new_start) + " new_end " + str(new_end))
    # 	print ("wt: " + wt_read.query_sequence)
    # 	print ("wt: " + wt_read.cigarstring)
    # 	print ("wt: " + explode_wt_cigar)
    # 	print ("wt: " + str(wt_read.query_qualities))
    # 	print ("new: " + new_read.query_sequence)
    # 	print ("new: " + new_read.cigarstring)
    # 	print ("new: " + explode_new_cigar)
    # 	print ("new: " + str(new_read.query_qualities))
    #
    # 	print ("final: " + final_read)
    # 	print ("final: " + final_cigar)
    #   print ("final mdz: " + str(final_mdz_els))

    wt_quals = wt_read.query_qualities
    wt_read.query_sequence = final_read
    wt_read.cigarstring = unexplodeCigar(final_cigar)
    wt_read.query_qualities = wt_quals
    if wt_read.has_tag("NM"):
        wt_read.set_tag("NM", final_nm_count)
    if wt_read.has_tag("MD"):
        wt_read.set_tag("MD", unexplodeMDZ(final_mdz_els))

    replacement_info = unexplodeReplacementInfo(replacement_info_els)
    modification_info = unexplodeReplacementInfo(modification_info_els)

    return wt_read, replacement_info, modification_info

def modifyRead(_wt_read, _genome_seq, _genome_start, mod_type, mod_chrom, mod_pos, substituted_bp=None, inserted_bp=None, deletion_length=None, debug=False):
    """ Modifies _wt_read by adding a modification to the read
        If _new_read and _wt_read do not overlap, wt_read is returned
        If _wt_read is unmapped, wt_read is returned

        Returns:
            pysam.AlignedSegment: _wt_read with modifications from _new_read
            replacement_info: array of replacement positions and details (these replacements were made to _wt_read, e.g. if the _wt_read has a substitution which is restored to the genomic sequence, that will appear here)
            modification_info: array of modification positions and details (these are modifications as compared to the genome, e.g. if the wt read has a substitution, that will appear here)

    """
    wt_read = copy.deepcopy(_wt_read)

    _genome_end = _genome_start + len(_genome_seq)

    if wt_read.is_unmapped:
        return wt_read, [], []

    # print("wt: " + str(wt_read))
    # print("new: " + str(new_read))

    # wt is unaltered read
    wt_start = wt_read.reference_start
    wt_start_with_clipping = wt_read.reference_start - wt_read.query_alignment_start  # includes soft-clipped bases
    wt_end = wt_read.reference_end - 1
    wt_end_with_clipping = wt_read.reference_end + (wt_read.query_length - wt_read.query_alignment_end)


    new_start = mod_pos
    new_end = mod_pos
    if mod_type == "I":
        if inserted_bp == None:
            raise Exception('mod_type is "I" but inserted_bp is ' + str(inserted_bp))
        new_end += 1
    elif mod_type == "D":
        if deletion_length == None:
            raise Exception('mod_type is "D" but deletion_length is ' + str(deletion_length))
        new_end += deletion_length
    elif mod_type == "S":
        if substituted_bp == None:
            raise Exception('mod_type is "S" but substituted_bp is ' + str(substituted_bp))
        new_end += 1
    else:
        raise Exception('mod_type is not "I", "D", or "S"')

    print('wt_start: ' + str(wt_start) + ' new end: ' + str(new_end))
    if wt_start > new_end:
        return _wt_read, [], get_read_modifications(wt_read, _genome_seq, _genome_start, _genome_end)
    if wt_end < new_start:
        return _wt_read, [], get_read_modifications(wt_read, _genome_seq, _genome_start, _genome_end)

    explode_wt_cigar = explodeCigar(str(wt_read.cigarstring))

    final_read = ""
    final_cigar = ""
    final_mdz_els = []  # keep this as an arr because some operations are multiple characters (e.g. ^G)
    final_nm_count = 0  # keep track of NM for number of mismatches
    replacement_info_els = []  # array of replacement els to be collapsed - some operations are multiple bases
    modification_info_els = []  # array of changes from the genome (not just those inserted by the new read) 

    if debug:
        print ("wt: " + wt_read.query_sequence)
        print ("wt: " + wt_read.cigarstring)
        print ("wt: " + explode_wt_cigar)
        print ("wt: " + str(wt_read.query_qualities))
        print ("mod_type: " + mod_type)
        print ("mod_chrom: " + mod_chrom)
        print ("mod_pos: " + str(mod_pos)) 
        print ("substituted_bp: " + str(substituted_bp))
        print ("inserted_bp: " + str(inserted_bp))
        print ("deletion_length: " + str(deletion_length))

    wt_read_ind = 0
    wt_cigar_ind = 0
    final_genomic_loc = wt_start
    # add soft/hard clipped bases to final read
    while explode_wt_cigar[wt_cigar_ind] == "H":
        final_cigar += "H"
        wt_cigar_ind += 1

    while explode_wt_cigar[wt_cigar_ind] == "S":
        final_cigar += "S"
        final_read += wt_read.query_sequence[wt_read_ind]
        wt_cigar_ind += 1
        wt_read_ind += 1

    if debug:
        print('669 wt_read_ind: ' + str(wt_read_ind) + ' wt_cigar_ind: ' + str(wt_cigar_ind) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' final_genomic_loc: ' + str(final_genomic_loc))

    new_read_ind = 0
    new_cigar_ind = 0

    # if WT read starts to the left of the NEW read, just take bases from the WT read
    while final_genomic_loc < new_start:
        if debug:
            print('231 wt_read_ind: ' + str(wt_read_ind) + ' wt_cigar_ind: ' + str(wt_cigar_ind) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' final_genomic_loc: ' + str(final_genomic_loc) + " new_start: " + str(new_start))
            print('final_genomic_loc: ' + str(final_genomic_loc) + '_genome_start: ' + str(_genome_start))
            print('genome at pos: ' + str(_genome_seq[final_genomic_loc - _genome_start]))
        if explode_wt_cigar[wt_cigar_ind] == "D":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_mdz_els.append("^" + _genome_seq[final_genomic_loc - _genome_start])
            modification_info_els.append((final_genomic_loc+1, "D", 1))
            wt_cigar_ind += 1
            final_genomic_loc += 1
        elif explode_wt_cigar[wt_cigar_ind] == "I":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            modification_info_els.append((final_genomic_loc, "I", wt_read.query_sequence[wt_read_ind]))
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "S":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            final_read += wt_read.query_sequence[wt_read_ind]
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "H":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            wt_cigar_ind += 1
        elif explode_wt_cigar[wt_cigar_ind] == "M":
            final_cigar += explode_wt_cigar[wt_cigar_ind]
            base_from_wt = wt_read.query_sequence[wt_read_ind]
            base_in_genome = _genome_seq[final_genomic_loc - _genome_start]
            final_read += base_from_wt
            if base_from_wt.upper() != base_in_genome.upper():
                final_mdz_els.append(base_from_wt)
                final_nm_count += 1
                modification_info_els.append((final_genomic_loc + 1, "S", base_from_wt))
            else:
                final_mdz_els.append("=")
            wt_cigar_ind += 1
            wt_read_ind += 1
            final_genomic_loc += 1
        else:
            raise Exception("Unexpected CIGAR operation (" + str(explode_wt_cigar[wt_cigar_ind] + ") in WT cigar: " + str(wt_read.cigarstring)))

    if debug:
        print('added wt bases that didn\'t overlap with new read')
        print('final_read: ' + str(final_read) + ' final cigar: ' + str(final_cigar))
        print('New start should be >= final_genomic_loc: ' + str(new_start) + ' >= ' + str(final_genomic_loc))
        print('721 final_genomic_loc: ' + str(final_genomic_loc) + ' new_start ' + str(new_start))

        print("723 wt len is : " + str(len(wt_read.query_sequence)) + " final_read is: " + final_read + " wt_start " + str(wt_start) + " wt_end " + str(wt_end) + " new start: " + str(new_start) + " new_end: " + str(new_end))

    # at this point the final_genomic_loc should be where new and wt overlap
    tmp_wt_genomic_loc = final_genomic_loc  # this can sometimes be different than final_genomic_loc if wt has indels
    check_wt_for_replacement_info = True  # if we're in the region where wt and new overlap, check for sequence changes

    if final_genomic_loc > new_start:
        raise Exception('Final genomic loc is greater than new start: ' + str(final_genomic_loc) + ' > ' + str(new_start))

    if mod_type == 'I':
        for i in range(len(inserted_bp)):
            replacement_info_els.append((final_genomic_loc, 'I', inserted_bp[i]))
            modification_info_els.append((final_genomic_loc, 'I', inserted_bp[i]))
            final_cigar += 'I'
            final_read += inserted_bp[i]
    elif mod_type == 'D':
        for i in range(deletion_length):
            final_cigar += "D"
            final_mdz_els.append('^'+_genome_seq[final_genomic_loc - _genome_start])
            replacement_info_els.append((final_genomic_loc+1, 'D', 1))
            modification_info_els.append((final_genomic_loc+1, 'D', 1))
            final_genomic_loc += 1
    elif mod_type == 'S':
        print('mod type here 743')
        print('final genomic loc: ' + str(final_genomic_loc))
        base_to_add = substituted_bp
        final_read += base_to_add
        if base_to_add.upper() != _genome_seq[final_genomic_loc - _genome_start].upper():
            final_mdz_els.append(base_to_add)
            final_nm_count += 1
            modification_info_els.append((final_genomic_loc+1, 'S', base_to_add))
        else:
            final_mdz_els.append("=")
        if check_wt_for_replacement_info:
            if explode_wt_cigar[wt_cigar_ind] == "S":
                pass
            elif explode_wt_cigar[wt_cigar_ind] == "I":  # insertions will be restored below in the wt catchup loop
                pass
            elif explode_wt_cigar[wt_cigar_ind] == "D":
                replacement_info_els.append((final_genomic_loc+1, 'R', _genome_seq[final_genomic_loc - _genome_start]))
            elif base_to_add.upper() != wt_read.query_sequence[wt_read_ind].upper():
                replacement_info_els.append((final_genomic_loc + 1, 'S', base_to_add))
        final_genomic_loc += 1

        #move pointer to correct position in wt
        while tmp_wt_genomic_loc < final_genomic_loc:
            if wt_read_ind >= len(wt_read.query_sequence):
                check_wt_for_replacement_info = False
                break
            if explode_wt_cigar[wt_cigar_ind] in ['S', 'H']: # at end of a read
                check_wt_for_replacement_info = False
                break
            if explode_wt_cigar[wt_cigar_ind] == 'I':
                replacement_info_els.append((tmp_wt_genomic_loc - 1, 'R', wt_read.query_sequence[wt_read_ind]))
            if explode_wt_cigar[wt_cigar_ind] in ['M', 'I', '=', 'X']:  # these cigar actions consume the read
                wt_read_ind += 1
            if explode_wt_cigar[wt_cigar_ind] in ['M', 'D', 'N', '=', 'X']:  # these cigar actions consume the reference
                tmp_wt_genomic_loc += 1
            wt_cigar_ind += 1

        # if we've reached the end of the wt sequence
        if wt_read_ind >= len(wt_read.query_sequence):
            check_wt_for_replacement_info = False

    if debug:
        print("len is : " + str(len(wt_read.query_sequence)) + " final is: " + final_read + " " + final_cigar + " wt_start " + str(wt_start) + " new start: " + str(new_start))
        print("at 2: wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
        	"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

        print('wt_read:' + str(wt_read))
        print('Adding remaining bases from WT (if any)')
        print('final read: ' + str(final_read))
        print('final cigar: ' + str(final_cigar))
        print('replacement info els: ' + str(replacement_info_els))
        print("len final read: " + str(len(final_read)) + " wt len: " + str(len(wt_read.query_sequence)))
    if len(final_read) < len(wt_read.query_sequence):
        if debug:
            print("at 3: wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind, wt_cigar_ind)+\
        		"final_genomic_loc: %s final_read: %s finalCigar: %s replacement_els %s\n" % (final_genomic_loc,final_read,final_cigar,str(replacement_info_els)))

        # add from wt 
        while final_genomic_loc < wt_end_with_clipping and len(final_read) < len(wt_read.query_sequence):
            # print('469 final_genomic_loc: ' + str(final_genomic_loc) + ' wt_end: ' + str(wt_end) + ' len(final_read): ' + str(len(final_read)) + ' len(wt_read.query_sequence): ' + str(len(wt_read.query_sequence)))
            # print('wt_cigar_ind: ' + str(wt_cigar_ind) + " (len " + str(len(explode_wt_cigar)) + ') wt_read_ind: ' + str(wt_read_ind) + ' final_genomic_loc: ' + str(final_genomic_loc) + ' final_read: ' + final_read + ' final_cigar: ' + final_cigar + ' replacement_info_els: ' + str(replacement_info_els))
            # print('wt_cigar: ' + str(explode_wt_cigar))
            # print('cigar val: ' + str(explode_wt_cigar[wt_cigar_ind]))
            if explode_wt_cigar[wt_cigar_ind] == "H":
                final_cigar += explode_wt_cigar[wt_cigar_ind]
                wt_cigar_ind += 1
            elif explode_wt_cigar[wt_cigar_ind] == "S":
                final_cigar += explode_wt_cigar[wt_cigar_ind]
                final_read += wt_read.query_sequence[wt_read_ind]
                wt_cigar_ind += 1
                wt_read_ind += 1
            elif explode_wt_cigar[wt_cigar_ind] == "D":
                final_cigar += explode_wt_cigar[wt_cigar_ind]
                final_mdz_els.append("^"+_genome_seq[final_genomic_loc - _genome_start])
                modification_info_els.append((final_genomic_loc+1, 'D', 1))
                wt_cigar_ind += 1
                final_genomic_loc += 1
            elif explode_wt_cigar[wt_cigar_ind] == "I":
                final_cigar += explode_wt_cigar[wt_cigar_ind]
                final_read += wt_read.query_sequence[wt_read_ind]
                modification_info_els.append((final_genomic_loc, 'I', wt_read.query_sequence[wt_read_ind]))
                wt_cigar_ind += 1
                wt_read_ind += 1
            else:
                final_cigar += explode_wt_cigar[wt_cigar_ind]
                base_to_add = wt_read.query_sequence[wt_read_ind]
                final_read += base_to_add
                if base_to_add.upper() != _genome_seq[final_genomic_loc - _genome_start].upper():
                    final_mdz_els.append(base_to_add)
                    final_nm_count += 1
                    modification_info_els.append((final_genomic_loc+1, 'S', base_to_add))
                else:
                    final_mdz_els.append("=")
                wt_cigar_ind += 1
                wt_read_ind += 1
                final_genomic_loc += 1

        # print("at 4: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind,new_cigar_ind)+\
        # 		"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
        # 		"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

        # finally, add from genome
        last_cigar_letter = "M"
        if final_cigar[-1] == "S": # if read ends with soft clipped bases, add soft-clipped bases (M can't come after S at end of cigar)
            last_cigar_letter = "S"
        while len(final_read) < len(wt_read.query_sequence):
            final_cigar += last_cigar_letter
            final_read += _genome_seq[final_genomic_loc - _genome_start]
            final_mdz_els.append("=")
            final_genomic_loc += 1
        
    # add trailing hardclipped bases (if any)
    temp_h_wt_pos = len(explode_wt_cigar) - 1
    while (explode_wt_cigar[temp_h_wt_pos] == "H"):
        final_cigar += "H"
        temp_h_wt_pos -= 1
    # print("at 5: new_read_ind: %s new_cigar_ind: %s\n" % (new_read_ind,new_cigar_ind)+\
    # 	"wt_read_ind: %s wt_cigar_ind: %s\n" % (wt_read_ind,wt_cigar_ind)+\
    # 	"final_genomic_loc: %s final_read: %s finalCigar: %s\n" % (final_genomic_loc,final_read,final_cigar))

    # if re.search("14D",new_read.cigarstring):
    # 	print ("wt_start: " + str(wt_start) + " wt_end " + str(wt_end))
    # 	print ("new_start: " + str(new_start) + " new_end " + str(new_end))
    # 	print ("wt: " + wt_read.query_sequence)
    # 	print ("wt: " + wt_read.cigarstring)
    # 	print ("wt: " + explode_wt_cigar)
    # 	print ("wt: " + str(wt_read.query_qualities))
    # 	print ("new: " + new_read.query_sequence)
    # 	print ("new: " + new_read.cigarstring)
    # 	print ("new: " + explode_new_cigar)
    # 	print ("new: " + str(new_read.query_qualities))
    #
    # 	print ("final: " + final_read)
    # 	print ("final: " + final_cigar)
    #   print ("final mdz: " + str(final_mdz_els))

    wt_quals = wt_read.query_qualities
    wt_read.query_sequence = final_read
    wt_read.cigarstring = unexplodeCigar(final_cigar)
    wt_read.query_qualities = wt_quals
    if wt_read.has_tag("NM"):
        wt_read.set_tag("NM", final_nm_count)
    if wt_read.has_tag("MD"):
        wt_read.set_tag("MD", unexplodeMDZ(final_mdz_els))

    replacement_info = unexplodeReplacementInfo(replacement_info_els)
    modification_info = unexplodeReplacementInfo(modification_info_els)

    return wt_read, replacement_info, modification_info


def cleanIndelsFromReadEnds(_oldRead,  num_match_bases_to_see=5):
    """Cleans indels from the beginning and end of reads until at least num_match_bases_to_see bases are seen.

    Args:
        _oldRead (pysam.AlignedSegment): Read to clean indels from
        num_match_bases_to_see (int, optional): Indels are cleaned until this number of matched bases are observed. Defaults to 5.

    Returns:
        _type_: _description_
    """
    tmp_read = cleanStartIndelsFromRead(_oldRead, num_match_bases_to_see)
    new_read = cleanEndIndelsFromRead(tmp_read, num_match_bases_to_see)
    return new_read


def cleanStartIndelsFromRead(_oldRead, num_match_bases_to_see=5):
    """
    This function takes a read and cleans up the start of the read
    by removing any indels that are at the start of the read.
    """

    cleanRead = copy.deepcopy(_oldRead)
    explodeOldCigar = explodeCigar(str(_oldRead.cigarstring))
    old_qualities = pysam.qualities_to_qualitystring(_oldRead.query_qualities)

    final_cigar = ""
    final_read_seq = ""
    final_read_qual = ""

    match_count = 0
    cigar_ind = 0
    read_ind = 0
    ref_ind = _oldRead.reference_start
    first_match_ref_ind = 0

    for cigar_ind in range(len(explodeOldCigar)):
        if explodeOldCigar[cigar_ind] == "M":
            # match: add seq and quality and increment alignment pos
            final_cigar += "M"
            final_read_seq += _oldRead.query_sequence[read_ind]
            final_read_qual += old_qualities[read_ind]
            if first_match_ref_ind == 0:
                first_match_ref_ind = ref_ind
            read_ind += 1
            ref_ind += 1
            match_count += 1
        else:
            final_cigar = ""
            final_read_seq = ""
            final_read_qual = ""
            match_count = 0
            first_match_ref_ind = 0

            if explodeOldCigar[cigar_ind] == "I":
                # insertions don't affect alignment position
                # don't add these bases to the quality or seq strings
                read_ind += 1
            elif explodeOldCigar[cigar_ind] == "D":
                # deletions affect alignment position
                # don't add to the quality or seq strings
                ref_ind += 1
            elif explodeOldCigar[cigar_ind] == "H":
                # hard-clipped bases aren't included in the read sequence
                pass
            elif explodeOldCigar[cigar_ind] == "S":
                # soft-clipped bases are included in the read sequence but don't affect alignment position
                read_ind += 1
            else:
                raise Exception("Unknown cigar operation: %s" % explodeOldCigar[cigar_ind])

        if match_count >= num_match_bases_to_see:
            break

    # add the rest of the sequence, qual, and cigar string
    cleanRead.reference_start = first_match_ref_ind
    cleanRead.query_sequence = final_read_seq + _oldRead.query_sequence[read_ind:]
    cleanRead.query_qualities = pysam.qualitystring_to_array(final_read_qual + old_qualities[read_ind:])
    cleanRead.cigarstring = unexplodeCigar(final_cigar + explodeOldCigar[cigar_ind + 1:])

    return (cleanRead)


def cleanEndIndelsFromRead(_oldRead, num_match_bases_to_see=5):
    """
    This function takes a read and cleans up the end of the read
    by removing any indels that are at the end of the read.
    """

    cleanRead = copy.deepcopy(_oldRead)
    explodeOldCigar = explodeCigar(str(_oldRead.cigarstring))
    old_qualities = pysam.qualities_to_qualitystring(_oldRead.query_qualities)

    final_cigar = ""
    final_read_seq = ""
    final_read_qual = ""

    match_count = 0
    cigar_ind = 0
    read_ind = len(_oldRead.query_sequence) - 1
    ref_ind = _oldRead.reference_start

    for cigar_ind in reversed(range(len(explodeOldCigar))):
        if explodeOldCigar[cigar_ind] == "M":
            # match: add seq and quality and increment alignment pos
            final_cigar += "M"
            final_read_seq += _oldRead.query_sequence[read_ind]
            final_read_qual += old_qualities[read_ind]
            read_ind -= 1
            ref_ind -= 1
            match_count += 1
        else:
            final_cigar = ""
            final_read_seq = ""
            final_read_qual = ""
            match_count = 0

            if explodeOldCigar[cigar_ind] == "I":
                # insertions don't affect alignment position
                # don't add these bases to the quality or seq strings
                read_ind -= 1
            elif explodeOldCigar[cigar_ind] == "D":
                # deletions affect alignment position
                # don't add to the quality or seq strings
                ref_ind -= 1
            elif explodeOldCigar[cigar_ind] == "H":
                # hard-clipped bases aren't included in the read sequence
                pass
            elif explodeOldCigar[cigar_ind] == "S":
                # soft-clipped bases are included in the read sequence but don't affect alignment position
                read_ind -= 1
            else:
                raise Exception("Unknown cigar operation: %s" % explodeOldCigar[cigar_ind])

        if match_count >= num_match_bases_to_see:
            break

    # add the rest of the sequence, qual, and cigar string
    cleanRead.query_sequence = _oldRead.query_sequence[:read_ind+1] + final_read_seq[::-1]
    cleanRead.query_qualities = pysam.qualitystring_to_array(''.join(old_qualities[:read_ind+1]+final_read_qual[::-1]))
    cleanRead.cigarstring = unexplodeCigar(explodeOldCigar[:cigar_ind] + final_cigar[::-1])

    return (cleanRead)


def get_read_modifications(wt_read, _genome_seq, _genome_start, _genome_end):
    """Gets modifications from a read (compared to genome)
    This function is used when returning from replace_reads early and not inserting changes

    Args:
        wt_read (pysam.AlignedSegment): read to get modifications from
        _genome_seq (str): genome sequence (may be a subsequence of the whole chromosome)
        _genome_start (int): start position of the genome_seq provided
        _genome_end (int): end position of the genome_seq provided
    """
    modification_info_els = []

    explode_wt_cigar = explodeCigar(wt_read.cigarstring)
    wt_start = wt_read.reference_start
    final_genomic_loc = wt_start
    wt_read_ind = 0
    wt_cigar_ind = 0
    while wt_cigar_ind < len(explode_wt_cigar):
        if wt_read_ind >= len(wt_read.query_sequence):
            raise Exception('Cigar string is incorrect for read: ' + str(wt_read) + ' sequence: ' + str(wt_read_ind) + ' wt_read_ind: ' + str(wt_read_ind) + ' wt_cigar_ind: ' + str(wt_cigar_ind) + ' explode_wt_cigar: ' + str(explode_wt_cigar))
        if str(explode_wt_cigar[wt_cigar_ind]) == "H":
            wt_cigar_ind += 1
        if str(explode_wt_cigar[wt_cigar_ind]) == "S":
            wt_cigar_ind += 1
            wt_read_ind += 1
        elif str(explode_wt_cigar[wt_cigar_ind]) == "D":
            modification_info_els.append((final_genomic_loc + 1, "D", 1))
            wt_cigar_ind += 1
            final_genomic_loc += 1
        elif str(explode_wt_cigar[wt_cigar_ind]) == "I":
            modification_info_els.append((final_genomic_loc, "I", wt_read.query_sequence[wt_read_ind]))
            wt_cigar_ind += 1
            wt_read_ind += 1
        else:
            base_from_wt = wt_read.query_sequence[wt_read_ind]
            if final_genomic_loc > _genome_start and final_genomic_loc < _genome_end:
                base_in_genome = _genome_seq[final_genomic_loc - _genome_start]
                if base_from_wt.upper() != base_in_genome.upper():
                    modification_info_els.append((final_genomic_loc + 1, "S", base_from_wt))
            wt_cigar_ind += 1
            wt_read_ind += 1
            final_genomic_loc += 1
    
    modification_info = unexplodeReplacementInfo(modification_info_els)
    return (modification_info)

def make_new_read_from_genome(refFile, chrom, start, end, readName, readQuals=None):
    """
    Make a new read from the reference genome
    """
    new_read = pysam.AlignedSegment()
    new_read.query_name = readName
    new_read.query_sequence = refFile.fetch(chrom, start, end)
    new_read.flag = 0
    new_read.reference_id = 0
    new_read.reference_start = start
    new_read.mapping_quality = 20
    new_read.cigarstring = str(end - start) + "M"

    new_read.query_qualities = [30] * len(new_read.query_sequence)
    if readQuals is not None:
        new_read.query_qualities = readQuals
    return new_read


def make_read_with_indel(refFile, mod_type, mod_chrom, mod_pos, mod_bases, read_to_alter=None):
    """
    Make a new read with a specific indel
    
    Args:
        refFile (pysam.FastaFile): Reference genome file
        mod_type (str): Type of modification (I, D, S)
        mod_chrom (str): Chromosome of modification
        mod_pos (int): Position of modification
        mod_bases (str): Bases to insert or delete
        read_to_alter (pysam.AlignedSegment, optional): Read to alter. Defaults to None.

    """
    if read_to_alter is None:
        new_read = make_new_read_from_genome(refFile, mod_chrom, mod_pos-1, mod_pos + 1, "newRead")
    else:
        new_read = copy.deepcopy(read_to_alter)
    if mod_type == 'S':
        sub_nuc = mod_bases
        sub_loc = int(mod_pos)
        old_len = len(new_read.query_sequence)
        pre_bases = int(old_len/2) # how many bases pre-substitution
        new_start = sub_loc - pre_bases
        post_bases = int(old_len/2) - 1
        new_end = sub_loc + post_bases
        new_seq = refFile.fetch(mod_chrom, new_start, sub_loc) + sub_nuc + refFile.fetch(mod_chrom, sub_loc+1, new_end)
        new_read.query_sequence = new_seq
        new_read.query_qualities = [30] * len(new_seq)
        new_read.reference_start = new_start
        new_read.cigarstring = str(pre_bases + post_bases + 1) + "M"
    elif mod_type == 'I':
        insert_loc = int(mod_pos)
        insert_nucs = mod_bases
        old_len = len(new_read.query_sequence)
        pre_bases = int(old_len/2) # how many bases pre-insertion
        new_start = insert_loc - pre_bases
        post_bases = int(old_len/2) - len(insert_nucs) # how many bases post-insertion
        new_end = insert_loc + post_bases
        new_seq = refFile.fetch(mod_chrom, new_start, insert_loc) + insert_nucs + refFile.fetch(mod_chrom, insert_loc, new_end)
        new_read.query_sequence = new_seq
        new_read.query_qualities = [30] * len(new_seq)
        new_read.reference_start = new_start
        new_read.cigarstring = str(pre_bases) + "M" + str(len(insert_nucs)) + "I" + str(post_bases) + "M"
    elif mod_type == 'D':
        del_loc = int(mod_pos)
        del_len = int(mod_bases)
        del_len = int(del_len)
        old_len = len(new_read.query_sequence)  # create a new read with the same length as one in the alteredReads set
        pre_post_bases = int(old_len/2)
        new_start = del_loc - pre_post_bases
        new_end = del_loc + pre_post_bases + del_len
        new_seq = refFile.fetch(mod_chrom, new_start, del_loc) + refFile.fetch(mod_chrom, del_loc + del_len, new_end)
        new_read.query_sequence = new_seq
        new_read.query_qualities = [30] * len(new_seq)
        new_read.reference_start = new_start
        new_read.cigarstring = str(pre_post_bases) + "M" + str(del_len) + "D" + str(pre_post_bases) + "M"
    else:
        raise Exception('Cannot parse indel type: ' + mod_type + '. Must be one of "I", "D", or "S"')

def prep_ref_file(ref_file_path):
    """
    Initializes a reference file for use in the script
    """
    ref_file = pysam.Fastafile(ref_file_path)
    return ref_file

if __name__ == "__main__":

    genome = "G"*100
    genome_start = 3

    a = pysam.AlignedSegment()
    qseq = "A"*20
    a.query_name = "read1"
    a.query_sequence = qseq
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 10
    a.mapping_quality = 20
    a.cigarstring = str(len(qseq)) + "M"
    a.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
    a.tags = (("NM", 1), ("RG", "L1"))

    b = pysam.AlignedSegment()
    qseq = "T"*20
    b.query_name = "read2"
    b.query_sequence = qseq
    b.flag = 0
    b.reference_id = 0
    b.reference_start = 15
    b.mapping_quality = 20
    b.cigarstring = str(len(qseq)) + "M"
    b.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
    b.tags = (("NM", 1), ("RG", "L1"))

    print("readA: " + str(a))
    print("readB: " + str(b))
    # b is old
    if (0):
        new = replaceRead(b, a)
        if new.query_sequence != "AAAAAAAAAAAAAAATTTTT":
            raise Exception("Unexpected result!")

    #    10
    # a: AAAAAAAAAAAAAAAAAAA
    # b:      TTTTTTTTTTTTTTTTTTTT
    # m: AAAAAAAAAA-----AAAAACCCCC
    # m2:       AAA---CCC
    m = pysam.AlignedSegment()
    qseq = "AAAAAAAAAAAAAAACCCCC"
    m.query_name = "mod"
    m.query_sequence = qseq
    m.flag = 0
    m.reference_id = 0
    m.reference_start = 10
    m.mapping_quality = 20
    m.cigarstring = "10M5D10M"
    m.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
    m.tags = (("NM", 1), ("RG", "L1"))

    m2 = pysam.AlignedSegment()
    qseq = "AAACCC"
    m2.query_name = "mod2"
    m2.query_sequence = qseq
    m2.flag = 0
    m2.reference_id = 0
    m2.reference_start = 18
    m2.mapping_quality = 20
    m2.cigarstring = "3M3D3M"
    m2.query_qualities = pysam.qualitystring_to_array("<"*len(qseq))
    m2.tags = (("NM", 1), ("RG", "L1"))

    new = replaceRead(b, m2, genome, genome_start)
    print("new:" + str(new))
    print("readA: " + str(a))
    print("readB: " + str(b))
