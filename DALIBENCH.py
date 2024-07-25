#DALIBENCH
#Author: Marcel Gabor
#Institute for Mathematics and Computer Science, University of Greifswald
#11/04/2024

import argparse
import sys
from datetime import datetime
from distutils.command import clean
import re
import os
import time
import pandas
import math
import random
import requests
import subprocess
import gzip
import warnings
import traceback
import contextlib
from rcsbsearchapi.search import SequenceQuery
from Bio import AlignIO
from Bio import SeqIO
import matplotlib

#defining path to script
def set_script_path():
    
    global SCRIPT_DIR
    SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
    os.chdir(SCRIPT_DIR)


#writing error messages to log
#void (str)
def errorFct(errorMsg):
    try:
        print(errorMsg)
        log_path = os.path.join(SCRIPT_DIR, "errorLog.txt")
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errorMsg)
            errFile.write(logEntry)
    except:
        print("Error log file could not be opened.")
        exit(1)

#safely closing files, writing to error log if it fails
#void (file handle, str, str)
def close_file_safely(file, file_path, errMsgForward):
    try:
        file.close()
    except:
        errorMsg = "%s\nCould not close %s." % (errMsgForward, file_path)
        errorFct(errorMsg)
        exit(1)
        
#create dir in specified path
#void str        
def create_dir_safely(path): 
    
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    except:
        errMsg = "Failed to create \"%s\" directory. Insufficient writing priviliges in specified output path?" % os.path.split(path)[1]
        errorFct(errMsg)
        exit(1)

#writing list of dictionaries to custom CSV file importable by pandas into dataframe
#void ([{},{},...], str, str)
def write_list_of_dicts_to_csv(mydict_list, output_full_path):
    
    try:
        csv_file = open(output_full_path, 'w')
    except:
        errMsg = "CSV file %s could not be opened for writing." % output_full_path
        errorFct(errMsg)
        return
    #writing the csv, likely to be optimizable
    #remove trailing delimiters, consider using more robust delimiter
    try:
        #usually all dicts should have same size. if function is to be generalized, a check needs to be implemented
        for key in mydict_list[0].keys():
            csv_file.write(key + ",")
        csv_file.write("\n")
        for mydict in mydict_list:
            for value in mydict.values():
                csv_file.write(str(value)+",")
            csv_file.write("\n")
    except:
        errMsg = "An error occured while writing to %s." % output_full_path
        close_file_safely(csv_file, output_full_path, errMsg)
        errorFct(errMsg)

    close_file_safely(csv_file, output_full_path, '')

#removing gaps and other symbols not in "valid_symbols" from given alignment sequence data
#(str, str) -> str
def clean_seq(seq, regex_str):

    cleaned_seq = re.sub(regex_str,'',seq)
    cleaned_seq = cleaned_seq.upper()
    
    return cleaned_seq

#function for custom selection of relevant metadata contained in the response body to be exported into CSV for further evaluation
#([str,str,...], {}) -> [{},{},...]
def parse_relevant_search_metadata(sample_seq_list, header_pdb_dict):
   
    relevant_item_dict = {}
    relevant_item_dict_list = []
    header_id_list = [row[1:3] for row in sample_seq_list]
    
    for header,internal_id in header_id_list:
        for entry in header_pdb_dict[header]:
            relevant_item_dict.clear()
            relevant_item_dict['identifier'] = entry['identifier']
            relevant_item_dict['original_fasta_header'] = header
            relevant_item_dict['score'] = entry['score']
            relevant_item_dict['sequence_identity'] = entry['services'][0]['nodes'][0]['match_context'][0]['sequence_identity']
            relevant_item_dict['evalue'] = entry['services'][0]['nodes'][0]['match_context'][0]['evalue']
            relevant_item_dict['bitscore'] = entry['services'][0]['nodes'][0]['match_context'][0]['bitscore']
            relevant_item_dict['alignment_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['alignment_length']
            relevant_item_dict['mismatches'] = entry['services'][0]['nodes'][0]['match_context'][0]['mismatches']
            relevant_item_dict['gaps_opened'] = entry['services'][0]['nodes'][0]['match_context'][0]['gaps_opened']
            relevant_item_dict['query_beg'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_beg']
            relevant_item_dict['query_end'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_end']
            relevant_item_dict['subject_beg'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_beg']
            relevant_item_dict['subject_end'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_end']
            relevant_item_dict['query_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['query_length']
            relevant_item_dict['subject_length'] = entry['services'][0]['nodes'][0]['match_context'][0]['subject_length']
            relevant_item_dict['internal_id'] = internal_id
            relevant_item_dict_list.append(relevant_item_dict.copy())

    return relevant_item_dict_list
 
#building sequence search query with set parameters and begin search. cast results to list.
#([str,str,int], str) -> []
def search_pdb_for_sequences(seq, aln_file_path, identity_cutoff):

    raw_response_list = []
    result_list = []
    e_val_cutoff = 1
    try:
        results = SequenceQuery(seq[0],e_val_cutoff,identity_cutoff)
        #request_options is currently not accessible via exec(), e.g. request_options["sort_by"]="sequence_identity"...
        result_list = list(results(results_verbosity="verbose", return_type="polymer_entity"))
    except:
        errMsg = "Failed to get results for sequence with header %s in file %s with %s and %s in sequence search." %(seq[1], aln_file_path, e_val_cutoff, identity_cutoff)
        errorFct(errMsg)

    for entry in result_list:
        raw_response_list.append(entry)
        
    if len(raw_response_list) == 0:
        return []

    return raw_response_list

#extracting sequences and headers from alignment file in FASTA format
#(str) -> [[str,str,int],[str,str,int],...]
def import_seq_list_from_fasta_aln(aln_file_path):

    seq_list = []

    try:
        aln_file = open(aln_file_path)
    except:
        errMsg = "Could not open alignment file %s." % aln_file_path
        errorFct(errMsg)
        return seq_list
    try:
        alignment = AlignIO.read(aln_file, "fasta")
    except:
        errMsg = "Could not read fasta file %s." %aln_file_path
        close_file_safely(aln_file, aln_file_path, errMsg)
        errorFct(errMsg)
        return seq_list

    close_file_safely(aln_file, aln_file_path, '')

    internal_id = 0

    for record in alignment:
        internal_id += 1
        seq_list.append([str(record.seq), str(record.id), internal_id])

    return seq_list

#reading file names in specified path and returning list of full paths
#function needs rework with os.path, and rename to get_file_paths_in_path
#(str) -> [str,str,...]
def get_file_names_in_path(path):

    path_file_list = []
 
    for element in os.listdir(path):
        full_path = os.path.join(path, element)
        if os.path.isfile(full_path):
            path_file_list.append(full_path)


    return path_file_list

#function that returns a list of PDB-IDs, one for each fasta header in CSV. default is the one with the highest SEQUENCE IDENTITY
#(str) -> ([str,str,...],{})
def get_ids_from_csv(csv_file_path, sort_by):
    
    id_list = []
    internal_id_pdb_name_dict = {}
    df = pandas.read_csv(csv_file_path)
    #remove next line, once trailing delimiters have been removed
    df = df.drop(['Unnamed: 16'], axis=1)
    #split df into one dataframe per fasta_header
    df_list = [x for _,x in df.groupby('original_fasta_header')]
    for split_df in df_list:
        if sort_by == "identity":
            split_df = split_df.sort_values(by=['sequence_identity'], ascending=False)
        elif sort_by == "score":
            split_df = split_df.sort_values(by=['score'], ascending=False)
        else:
            errMsg = "Can't sort by %s. Exiting." % sort_by
            errorFct(errMsg)
            exit(1)
        pdb_identifier = (split_df['identifier'].iloc[0])[0:4]
        id_list.append(pdb_identifier)
        internal_id_pdb_name_dict[pdb_identifier] = split_df['internal_id'].iloc[0]
    
    return id_list, internal_id_pdb_name_dict

#function to automatically send download requests to PDB according to a list of PDB-IDs,
#gathered from the generated CSV files
#inspired by https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script
#(str,str,int) -> {}
def download_pdb_files_by_CSV(csv_file_path, out_path, sort_by, verbosity):

    counter = 1
    id_list, internal_id_pdb_name_dict = get_ids_from_csv(csv_file_path, sort_by)
    id_list_len = len(id_list)
    
    for identifier in id_list:
        if verbosity > 0:
            msg = "Downloading \"%s\"... [%s/%s]" % (identifier, counter, id_list_len)
            print(msg)
            counter += 1
            
        ext_id = identifier + ".pdb"
        url = "https://files.rcsb.org/download/%s" % ext_id
        try:
            r = requests.get(url)
        except:
            errMsg = "Failed to get download request for ID %s." % identifier
            errorFct(errMsg)
            continue
        try:
            with open(os.path.join(out_path, ext_id), 'wb') as f:
                f.write(r.content)
        except:
            errMsg = "Failed to write to %s." % ext_id
            errorFct(errMsg)
    
    #maybe refactor so that id_list and internal_id_name_dict is passed to download_pdb_files_by_CSV directly
    return internal_id_pdb_name_dict

#function to select alignment sequence alphabet (e.g. protein, RNA, DNA)
#currently without choice    
#(str) -> [str,str,...]
def sel_alphabet(choice):
 
    if choice == "AA":
        #depending on the length of the protein sequence, different orders might be more efficient
        #see "Amino acid composition and protein dimension", Protein Sci. 17(12): 2187-2191 (2008)
        #maybe make into set. then rework of build_regex is required
        valid_symbols = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","O","S","U","T","W","Y","V","B","Z","X","J"]
        #IUPAC-IUB Joint Commission on Biochemical Nomenclature.Nomenclature and Symbolism for Amino Acids and Peptides. Eur. J. Biochem. 138: 9-37 (1984)
    elif choice == "DNA/RNA":
        valid_symbols = ["A","G","C","T"]
    
    return valid_symbols

#takes string of invalid symbols and returns a list
#this function mostly exists only for factoring consistency
#str -> [str,str,...]
def set_non_alphabet(non_alphabet):
    
    return list(non_alphabet)

#calculates the reference "sum-of-pairs-set" for structural alignments (lowercase chars mean not equivalent)
#which is a set containing double-nested 2-tuples ( i.e. ((x,y),(a,b)) )
#where x is the row index of the query structure, a is the row index of the sbjct structure
#y is the column index of the query structure, b is the column index of the sbjct structure
#if ((x,y),(a,b)) is in ref_ind_set, then this means that according to the structure alignment
#the sequence with index x at residue number y is equivalent to the sequence with index a at residue number b
#((int,int),(str,str)) -> set( ((),()) )
def calc_ref_SP_set(sequence_ids, aligned_sequences):

    query = aligned_sequences[0]
    sbjct = aligned_sequences[1]
    query_gaps = 0
    sbjct_gaps = 0
    ref_ind_set = set()

    #raw query and sbjct strings should have the same length in alignment (incl. gaps)
    if len(query) == len(sbjct):
        max_len = len(query)
    elif len(query) < len(sbjct):
        errMsg = "Warning: Aligned subject sequence is longer than query sequence."
        errorFct(errMsg)
        max_len = len(sbjct)
    else:
        errMsg = "Warning: Aligned query sequence is longer than subject sequence."
        errorFct(errMsg)
        max_len = len(query)
    for i in range(0, max_len):
        
        try:
            equiv_sbjct = True
            equiv_query = True
            #need to define a set of expected symbols or set of gap/unknown symbols
            if sbjct[i] == "-" or sbjct[i] == ".":
                sbjct_gaps += 1
                equiv_sbjct = False
            if query[i] == "-"  or sbjct[i] == ".":
                query_gaps += 1
                equiv_query = False
            #if one of the sequences has a gap, it's a mismatch
            if equiv_sbjct is False or equiv_query is False:
                continue
            if sbjct[i].islower():
                continue
            if query[i].islower():
                continue
            #is the following is an inefficient approach, when we know what symbols to expect?
            #if sbjct[i] not in valid_symbols:
             #   continue
            #if query[i] not in valid_symbols:
             #   continue
            ref_ind_set.add(((sequence_ids[0],i-query_gaps),(sequence_ids[1],i-sbjct_gaps)))
        except IndexError:
            break
    
    return ref_ind_set

#function to calculate MSA sum-of-pairs-set (different equivalence relation)
#its partly duplicate to "calc_ref_SP_set" and will most likely be refactored together with latter function
#
def calc_MSA_SP_set(sequence_ids, aligned_sequences):
    
    query = aligned_sequences[0]
    sbjct = aligned_sequences[1]
    query_gaps = 0
    sbjct_gaps = 0
    MSA_ind_set = set()
    #raw query and sbjct strings should have the same length in alignment (incl. gaps)
    if len(query) == len(sbjct):
        max_len = len(query)
    elif len(query) < len(sbjct):
        errMsg = "Warning: Aligned subject sequence is longer than query sequence."
        errorFct(errMsg)
        max_len = len(sbjct)
    else:
        errMsg = "Warning: Aligned query sequence is longer than subject sequence."
        errorFct(errMsg)
        max_len = len(query)
    for i in range(0, max_len):
        try:
            equiv_sbjct = True
            equiv_query = True
            #need to define a set of expected symbols or set of gap/unknown symbols
            if sbjct[i] == "-" or sbjct[i] == ".":
                sbjct_gaps += 1
                equiv_sbjct = False
            if query[i] == "-"  or sbjct[i] == ".":
                query_gaps += 1
                equiv_query = False
            #if one of the sequences has a gap, it's a mismatch
            if equiv_sbjct is False or equiv_query is False:
                continue
            if sbjct[i].upper() != query[i].upper():
                continue
            #is the following is an inefficient approach, when we know what symbols to expect?
            #if sbjct[i] not in valid_symbols:
            #   continue
            #if query[i] not in valid_symbols:
            #   continue
            MSA_ind_set.add(((sequence_ids[0],i-query_gaps),(sequence_ids[1],i-sbjct_gaps)))
        except IndexError:
            break
    
    return MSA_ind_set
        

#function for naive approach to maximizing intersection between two sets by shifting all element values in set_2 by some constant
#sets need to contain mutable elements
def argmax_intersection(set_1, set_2):
    
    #make copies of original sets
    orig_set_1 = set_1.copy()
    orig_set_2 = set_2.copy()

    #initial shift bounds +-1000
    shift_bounds = range(1000,-1001,-1)
    max_intersection = orig_set_1.intersection(orig_set_2)
    cur_max = len(max_intersection)
    cur_argmax = orig_set_2
    max_shift = 0
    for shift in shift_bounds:  
        new_set_2 = set()
        set_intersection = set()
        for element in set_2:
            y = element[0][1] + shift
            b = element[1][1] + shift
            x = element[0][0]
            a = element[1][0]
            new_set_2.add(((x,y),(a,b)))
        set_intersection = set_1.intersection(new_set_2)
        if len(set_intersection) > cur_max:
            max_intersection = set_intersection.copy()
            cur_max = len(max_intersection)
            cur_argmax = new_set_2.copy()
            max_shift = shift
        set_1 = orig_set_1
        set_2 = orig_set_2

    return cur_argmax, max_shift

#writes SP sets in job and resulting SPS to CSV for subsequent statistical evaluation
# 
def write_job_SPS_CSV(ref_set_list, aln_set_list, job_SPS):
    
    pass

#calculates SP scores for job and returns list of dictionaries containing information about each alignment
#highly specific signature
#([[{},{},...],[{},{},...],...], {}, {}) -> ([{},{},...], float)
def calc_SPS_scores_in_job(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict):
    current_max_ref_SP_set = set()
    save_dict = {}
    SPS_dict_list = []
    max_ref_SP_set_list = []
    argmax_aln_SP_set_list = []
    job_SPS = 0
    #for each alignment folder in job, i.e. for each sub job
    #for-else statement to continue outer loop in case of breaking of inner loop
    for dict_list in job_list_of_dict_lists:
        aln_SP_set = set()
        pair_SPS = 0
        pdb_name_1_const_check = dict_list[0]["query"][0:4]
        pdb_name_2_const_check = dict_list[0]["sbjct"][0:4]
        pdb_internal_id_1_const = internal_id_pdb_name_dict[pdb_name_1_const_check]
        pdb_internal_id_2_const = internal_id_pdb_name_dict[pdb_name_2_const_check]
        sub_job_counter = 0
        #for each chain alignment in PDB alignment 
        #calculating reference SP set for DALI alignment
        for dic in dict_list:
            ref_SP_set = set()
            pdb_name_1 = dic["query"][0:4]
            pdb_name_2 = dic["sbjct"][0:4]
            if (pdb_name_1 != pdb_name_1_const_check) or (pdb_name_2 != pdb_name_2_const_check):
                errMsg = "Differing PDB names detected within sub job. Skipping sub job with data:\n %s" % str(dict_list)
                errorFct(errMsg)
                break
            ref_SP_set = calc_ref_SP_set((pdb_internal_id_1_const, pdb_internal_id_2_const), (dic["DALI_query_seq"], dic["DALI_sbjct_seq"]))
            #DALI generates an unknown number of TXT files or alignments per TXT file depending on the number of chains per PDB file and what file acts as query or subject
            #for example if PDB1 contains 2 chains and PDB2 contains 3 chains, if PDB1 is query and PDB2 is subject, DALI will generate
            #2 TXT files, namely "PDB1A.txt" and "PDB1B.txt". each will contain 3 alignments. Thus we are left with 6 alignments.
            #namely: PDB1A vs. PDB2A, PDB1A vs. PDB2B, PDB1A vs. PDB2C, PDB1B vs. PDB2A, ... and so forth
            #the following check chooses the set with maximum length out of all possible chain alignments, so that we dont accidentally align
            #dissimilar chains with each other 
            if (len(ref_SP_set) > len(current_max_ref_SP_set)) or sub_job_counter == 0:
                current_max_ref_SP_set.clear()
                save_dict.clear()
                current_max_ref_SP_set = ref_SP_set.copy()
                save_dict = dic.copy()
            sub_job_counter += 1
        #calculate SP set for original alignment
        #outer loop    
        else:
            max_ref_SP_set_list.append(current_max_ref_SP_set.copy())
            aln_seq_1 = internal_id_raw_seq_dict[pdb_internal_id_1_const]
            aln_seq_2 = internal_id_raw_seq_dict[pdb_internal_id_2_const]
            aln_SP_set = calc_MSA_SP_set((pdb_internal_id_1_const, pdb_internal_id_2_const), (aln_seq_1, aln_seq_2))
            #extend saved aln dict with additional data
            save_dict["SP_DALI"] = len(current_max_ref_SP_set)
            save_dict["SP_orig_aln"] = len(aln_SP_set)
            save_dict["query_internal_id"] = pdb_internal_id_1_const
            save_dict["sbjct_internal_id"] = pdb_internal_id_2_const
            save_dict["orig_query_seq"] = aln_seq_1
            save_dict["orig_sbjct_seq"] = aln_seq_2
            #calculate SPS
            
            #try to maximize intersection size by shifting y and b in ((x,y),(a,b)) in aln_SP_set by some constants
            #try to find shift pattern to minimize maximization effort
            #this could actually be viable. but since we know from the search meta data where in the PDB the particular sequence match starts
            #we should use that information first
            #in any way, maximizing intersection increases sensitivity
            argmax_aln_SP_set, max_shift = argmax_intersection(current_max_ref_SP_set, aln_SP_set)
            argmax_aln_SP_set_list.append(argmax_aln_SP_set.copy())
            save_dict["max_shift"] = max_shift
            try:
                pair_SPS = len(argmax_aln_SP_set.intersection(current_max_ref_SP_set))/len(current_max_ref_SP_set)
            except ZeroDivisionError:
                pair_SPS = 0
            save_dict["pair_SPS"] = pair_SPS
            SPS_dict_list.append(save_dict.copy())
            
    #first unite all current_max_ref_SP_sets. unite aln_SP_sets. then intersect unions and divide by union of cur_max_ref_SP_sets.
    overall_ref_SP_set = set()
    overall_MSA_SP_set = set()
    for s in max_ref_SP_set_list:
        overall_ref_SP_set = overall_ref_SP_set.union(s)
    for s in argmax_aln_SP_set_list:
        overall_MSA_SP_set = overall_MSA_SP_set.union(s)
    try:
        job_SPS = len(overall_ref_SP_set.intersection(overall_MSA_SP_set))/len(overall_ref_SP_set)
    except ZeroDivisionError:
        #if length of denominator set is 0, it means that the intersection must be empty as well
        job_SPS = 0
        
    return SPS_dict_list, job_SPS

#import function for DALI outputs. one pass, line per line. might need rework, because it's hard to read and looks ugly. hopefully it's fast, at least.
#(str) -> ([{},{},...])
def import_DALI_aln(aln_file_path):

    try:
        alignment = open(aln_file_path, 'r')
    except:
        errorMsg = "Could not open %s for reading." % aln_file_path
        errorFct(errorMsg)
        exit(1)

    #read file line by line
    lineCount = 0
    queryFound = False
    sbjctFound = False
    query_list = []
    sbjct_list = []
    aln_dict_list = []
    aln_dict = {}
    aln_counter = 0
    while True:
        lineCount += 1
        line = alignment.readline()
        if not line:
            break

        #import jobname
        job_match = re.search(r"^#[ ]*Job:[ ]*(.*)", line)
        if job_match is not None:
            job_name = job_match.group(1)
            job_match = None

        #import header lines
        header_match = re.search(r"^[ ]*[0-9]+:[ ]*[a-zA-Z0-9-]*[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([A-Za-z0-9.: -]*)", line)
        if header_match is not None:
            aln_dict.clear()
            aln_dict['job_name'] = job_name
            aln_dict['Z-score'] = header_match.group(1)
            aln_dict['rmsd'] = header_match.group(2)
            aln_dict['lali'] = header_match.group(3)
            aln_dict['nres'] = header_match.group(4)
            aln_dict['perc_id'] = header_match.group(5)
            aln_dict['pdb_descr'] = header_match.group(6)
            aln_dict_list.append(aln_dict.copy())
            header_match = None
            continue
        
        if aln_counter > 0 and (line.startswith("# S") or line.startswith("# T")):
            overall_query = ''.join(query_list)
            overall_sbjct = ''.join(sbjct_list)
            aln_dict_list[aln_counter-1]['DALI_query_seq'] = overall_query
            aln_dict_list[aln_counter-1]['DALI_sbjct_seq'] = overall_sbjct
            #we're not using structural equivalences at the moment
            close_file_safely(alignment, aln_file_path, "")
            return aln_dict_list


        #set aln_counter
        aln_count_match = re.search(r"^No[ ]+([0-9]+):[ ]*Query=([a-zA-Z0-9]+)[ ]*Sbjct=([a-zA-Z0-9]+)", line)
        if aln_count_match is not None:
            aln_counter = int(aln_count_match.group(1))
            aln_dict_list[aln_counter-1]['query'] = aln_count_match.group(2)
            aln_dict_list[aln_counter-1]['sbjct'] = aln_count_match.group(3)
            #check if there's previous alignment. if we're here, it means that we're at the start of a new alignment
            if aln_counter > 1:
                overall_query = ''.join(query_list)
                overall_sbjct = ''.join(sbjct_list)
                aln_dict_list[aln_counter-2]['DALI_query_seq'] = overall_query
                aln_dict_list[aln_counter-2]['DALI_sbjct_seq'] = overall_sbjct
                query_list = []
                sbjct_list = []
            aln_count_match = None
            continue

        #search for Query sequence
        if queryFound is False:
            query_match = re.search(r"^Query ([a-zA-Z-]+)", line)
            if query_match is None:
                continue
            else:
                queryFound = True
                continue
        #search for Sbjct sequence
        if sbjctFound is False:
            sbjct_match = re.search(r"^Sbjct ([a-zA-Z-]+)", line)
            if sbjct_match is None:
                continue
            else:
                sbjctFound = True
                continue
            
        if (queryFound and sbjctFound) is True:
            query_list.append(query_match.group(1))
            queryFound = False
            sbjct_list.append(sbjct_match.group(1))
            sbjctFound = False
            continue

    #still assemble sequences, if end of file is reached
    if aln_counter <= 0:
        aln_dict.clear()
        #get query and sbjct name from dir name
        try:
            head_tail = os.path.split(aln_file_path)
            head_tail2 = os.path.split(head_tail[0])
            query_name = head_tail2[1][0:4]
            sbjct_name = head_tail2[1][5:9]
        except:
            query_name = "NA"
            sbjct_name = "NA"
            errMsg = "Path \"%s\" does not contain any information about query and subject name at expected location." % aln_file_path
            errorFct(errMsg)
        aln_dict['query'] = query_name
        aln_dict['sbjct'] = sbjct_name
        #set remaining dict entries
        aln_dict['job_name'] = job_name
        aln_dict['Z-score'] = "NA"
        aln_dict['rmsd'] = "NA"
        aln_dict['lali'] = "NA"
        aln_dict['nres'] = "NA"
        aln_dict['perc_id'] = "NA"
        aln_dict['pdb_descr'] = "NA"
        aln_dict['DALI_query_seq'] = "0"
        aln_dict['DALI_sbjct_seq'] = "0"
        aln_dict_list.append(aln_dict.copy())
        errMsg = "File \'%s\' contains no alignments." % aln_file_path
        errorFct(errMsg)
        close_file_safely(alignment, aln_file_path, errMsg)
        return aln_dict_list
    
    if (len(query_list) == len(sbjct_list)) and (queryFound == sbjctFound):
        overall_query = ''.join(query_list)
        overall_sbjct = ''.join(sbjct_list)
        aln_dict_list[aln_counter-1]['DALI_query_seq'] = overall_query
        aln_dict_list[aln_counter-1]['DALI_sbjct_seq'] = overall_sbjct           
    else:
        errorMsg = "Differing amounts of regex matches for Query and Sbjct sequences in %s" % aln_file_path
        close_file_safely(alignment, aln_file_path, errorMsg)
        errorFct(errorMsg)
        exit(1)

    close_file_safely(alignment, aln_file_path, '')
    return aln_dict_list

#import all alignment files in a dir, appends job_list_of_dict_list exactly ONCE
#(str) -> [{},{},...]
def import_aln_files_in_dir(aln_path):
    
    #dict_list for alignment dir
    aln_dict_list = []

    #get all TXT files
    txt_path_list = []
    for file in os.listdir(aln_path):
        if file.endswith(".txt"):
            txt_path_list.append(os.path.join(aln_path, file))
    #import all TXT files in dir and extend list
    for txt_path in txt_path_list:
        aln_dict_list.extend(import_DALI_aln(txt_path))

    return aln_dict_list

#import all alignments in one job
#(str, str) -> [[{},{},...], [{},{},...],...]
def import_aln_files_in_job(job_title, out_path):

    job_list_of_dict_lists = []
    job_file_list = []
    job_aln_path = os.path.join(out_path, job_title, "ALN")
    #list all subdirs in job dir
    job_file_list = os.listdir(job_aln_path)
    for file in job_file_list:
        if not os.path.isdir(os.path.join(job_aln_path, file)):
            job_file_list.remove(file)
    for aln_dir in job_file_list:
        aln_path = os.path.join(job_aln_path, aln_dir)
        job_list_of_dict_lists.append(import_aln_files_in_dir(aln_path))
        
    return job_list_of_dict_lists

#build regex from list of valid symbols for clean_seq() function
#([str,str,...]) -> str
def build_regex_for_seq_cleaning_whitelist(valid_symbols):

    regex_str_list = []
    regex_str_list = valid_symbols.copy()
    valid_symbols_lower = list(map(str.lower, valid_symbols.copy()))
    regex_str_list.extend(valid_symbols_lower)
    regex_str_list.append("]")
    regex_str_list.insert(0,"^")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

def build_regex_for_seq_cleaning_blacklist(invalid_symbols):
    
    regex_str_list = []
    regex_str_list = invalid_symbols.copy()
    invalid_symbols_lower = list(map(str.lower, invalid_symbols.copy()))
    regex_str_list.append("]")
    regex_str_list.insert(0,"[")
    regex_str = ''.join(regex_str_list)
    
    return regex_str

#takes a list of sequence-header pairs and returns a random sample of given size
#([[str,str,int],[str,str,int],...],int) -> [[str,str,int],[str,str,int],...]
def choose_random_sample_from_aln_file(seq_list, sample_size):

    sample_seq_list = []
    random_int_list = []
    counter = 0

    try:
        sample_size = int(sample_size)
    except:
        errMsg = "Sample size \"%s\" must be given as an integer." % sample_size
        errorFct(errMsg)
        exit(0)

    if sample_size >= len(seq_list):
        sample_seq_list = seq_list
        return sample_seq_list

    while counter < sample_size:

        random_int = random.randint(0, len(seq_list)-1)
        if random_int in random_int_list:
            continue
        random_int_list.append(random_int)
        sample_seq_list.append(seq_list[random_int].copy())
        counter += 1

    return sample_seq_list

#counts number of "real" proteins in alignment
#(str,str) -> (str,int,int,float) 
def count_homologues(aln_file_path, hom_aln_file_path):
    
    seq_list = import_seq_list_from_fasta_aln(aln_file_path)
    hom_seq_list = import_seq_list_from_fasta_aln(hom_aln_file_path)
    
    no_of_proteins = len(seq_list)
    no_of_proteins_and_homologues = len(hom_seq_list)
    prot_frac = no_of_proteins/no_of_proteins_and_homologues
    head, tail = os.path.split(hom_aln_file_path)
    
    #might want to write this to CSV
    msg = "There are a total of %s sequences in file %s. %s of them are proteins. Fraction: %s" % (no_of_proteins_and_homologues, tail, no_of_proteins, prot_frac)
    print(msg)
    
    return (tail, no_of_proteins_and_homologues, no_of_proteins, prot_frac)

#wrapper function to parse fasta alignment sequences and headers, query a PDB search and write relevant metadata to CSV in one function call
#verbosity option allows progress tracking in terminal
#([str,str,...],str,str,int,int) -> {}
def write_aln_file_search_hits_to_csv(valid_symbols, aln_file_path, output_path, verbosity, sample_size):

    seq_list = []
    regex_str = ""
    header_pdb_dict = {}

    if(verbosity>0):
        msg = "Parsing fasta sequences from %s." % aln_file_path
        print(msg)

    #rename seq_list, because it's not only a list of sequences 
    seq_list = import_seq_list_from_fasta_aln(aln_file_path)

    if sample_size != 0:
        sample_seq_list = choose_random_sample_from_aln_file(seq_list, sample_size)
    else:
        sample_seq_list = seq_list
    
    #double check whitelist, i.e. valid_symbols
    regex_str = build_regex_for_seq_cleaning_whitelist(valid_symbols)
    
    #make a copy of seq list before sequences are cleaned. returned by function.
    raw_sample_seq_list = sample_seq_list.copy()
    internal_id_raw_seq_dict = {}
    #fill dictionary to get sequence by internal sequence id or header later, returned by function. not elegant.
    for raw_sample in raw_sample_seq_list:
        internal_id_raw_seq_dict[raw_sample[2]] = raw_sample[0]
    
    if(verbosity>0):
        msg = "Searching for sequences in alignment file %s." % (os.path.split(aln_file_path)[1])
        seq_no = len(sample_seq_list)
        counter = 0
        print(msg)

    for seq in sample_seq_list:
        
        identity_cutoff = 1.00
        
        seq[0] = clean_seq(str(seq[0]), regex_str)
        if seq[1] in header_pdb_dict.keys():
            errMsg = "Duplicate header \"%s\" in alignment file %s." %(seq[1], os.path.split(aln_file_path)[1])
            errorFct(errMsg)
        
        header_pdb_dict[seq[1]] = search_pdb_for_sequences(seq, aln_file_path, identity_cutoff)

        #this loop keeps lowering the identity cutoff until proteins have been found
        #currently theres no good and easy way to sort search hits via the python interface
        while not header_pdb_dict[seq[1]]:
            if (verbosity>0):
                notification = "No matches for header \'%s\' with identity cutoff %s." %(seq[1], identity_cutoff)
                print(notification)

            header_pdb_dict[seq[1]] = search_pdb_for_sequences(seq, aln_file_path, identity_cutoff)
            if identity_cutoff < 0 or math.isclose(identity_cutoff, 0.00,abs_tol=1e-03):
                header_pdb_dict[seq[1]] = []
                break
            if identity_cutoff < 0.30:
                identity_cutoff = 0.00
            if math.isclose(identity_cutoff,1.00,abs_tol=1e-03): 
                identity_cutoff -= 0.02
            if identity_cutoff > 0.90:
                identity_cutoff -= 0.05
            if identity_cutoff < 0.90:
                identity_cutoff -= 0.2
        
        if(verbosity>0):
            counter += 1
            msg = "%s/%s. Header: \"%s\". Done!" % (counter, seq_no, seq[1])
            print(msg)
    
    #extract list of headers from list of sample and headers
    #header_list = [row[1] for row in sample_seq_list]

    if(verbosity>0):
            print("Parsing search metadata.")

    relevant_search_results_list = parse_relevant_search_metadata(sample_seq_list, header_pdb_dict)

    if(verbosity>0):
            msg = "Writing CSV to %s" % (output_path)
            print(msg)
    
        
    aln_file_name = os.path.basename(aln_file_path)
    csv_output_file_name = "%s.csv" % aln_file_name
    csv_output_full_path = os.path.join(output_path, csv_output_file_name)
    write_list_of_dicts_to_csv(relevant_search_results_list, csv_output_full_path)

    if(verbosity>0):
            print("Done!")
            
    return internal_id_raw_seq_dict

#function to write error log for errors encountered during subprocesses involving the shell
#void str        
def shell_err_fct(errMsg):
    try:
        log_path = os.path.join(SCRIPT_DIR, "shell_err_log.txt")
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errMsg)
            errFile.write(logEntry)
    except:
        print("Error log file could not be opened.")
        exit(1)

#returns list of PDB names from given dir of PDB files
#str -> [str,str,...]        
def get_pdb_name_list(pdb_path):
    
    pdb_path_list = get_file_names_in_path(pdb_path)
    pdb_name_list = []
    for pdb_path in pdb_path_list:
        head, tail = os.path.split(pdb_path)
        tail = os.path.splitext(tail)[0]
        pdb_name_list.append(tail)
        
    return pdb_name_list

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run import.pl of DaliLite.v5 to import PDBs in a given dir to DAT file format used by DALI
#void (str, str, str, int)
def DALI_import_PDBs(pl_bin_path, pdb_path, DAT_path, verbosity):

    pdb_name_list = get_pdb_name_list(pdb_path)
        
    import_pl_full_path = os.path.join(pl_bin_path, "import.pl")

    #create dir for output and err
    create_dir_safely(os.path.join(DAT_path, "err_logs"))
    
    create_dir_safely(os.path.join(DAT_path, "out_logs"))
    
    #this loop is necessary for clean error logging, since DaliLite.v5 seems to have bad logging
    for pdb_name in pdb_name_list:
        
        shell_input = []
        pdb_full_path = os.path.join(pdb_path, pdb_name+".pdb")
        err_file_path = os.path.join(DAT_path, "err_logs/", pdb_name)
        out_file_path = os.path.join(DAT_path, "out_logs/", pdb_name)
        
        #/.../<dali_dir>/import.pl --pdbfile <path> --pdbid <id> --dat <path> --verbose --clean
        shell_input = [import_pl_full_path, '--pdbfile', pdb_full_path, '--pdbid', pdb_name, '--dat', DAT_path, '--clean', '1', '>', out_file_path, '2', '>', err_file_path]
        
        #create file handles for stdout and stderr
        try:
            err_file = open(err_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % err_file_path
            #no need to close file when failed to open
            try:
                close_file_safely(err_file, err_file_path, errMsg)
            except UnboundLocalError:
                pass
            errorFct(errMsg)
            exit(1)
        try:
            out_file = open(out_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % out_file_path
            try:
                close_file_safely(out_file, out_file_path, errMsg)
            except UnboundLocalError:
                pass
            errorFct(errMsg)
            exit(1)

        #subprocess for importing PDBs
        process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
        
        try:
            out, err = process.communicate(timeout=15)
            if out is not None and verbosity>1:
                print(out.decode('ascii'))
        except subprocess.TimeoutExpired:
            process.kill()
            out, err = process.communicate()
            #print(out.decode('ascii'))
            errMsg = "Communication with import.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt."
            close_file_safely(err_file, err_file_path, errMsg)
            close_file_safely(out_file, out_file_path, errMsg)
            #shell_err_fct(err.decode('ascii'))
            errorFct(errMsg)
            exit(1)
            
        #close files
        close_file_safely(err_file, err_file_path, "")
        close_file_safely(out_file, out_file_path, "")
    
    #rewrite a sufficient check whether import has been successful
    for file_name in os.listdir(DAT_path):
        if file_name.endswith(".d"):
            errMsg = "Failed to convert PDBs to DAT format. Exiting."
            errorFct(errMsg)
            exit(1)

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run dali.pl of DaliLite.v5. Uses "pairwise" functionality of DALI for better control over what alignments are being made.
#writes stderr and stdout to files for each job and alignment. improved output dir structure, since DaliLite.v5 doesn't seem to support specifying an output dir
#dir structure: <jobtitle>/ALN/<queryPDBname_sbjctPDBname>/        
#void (str,str,str,str,str,int)
def DALI_all_vs_all_query(pl_bin_path, pdb_path, out_path, DAT_path, job_title, verbosity):

    #for some weird implementation reasons (probably Fortran though) DALI doesn't generate alignment files for job titles longer than 12 chars
    #and generates only empty alignments for job titles exactly 12 chars long
    #this is a DALI problem and can't be changed at the moment
    #i could, however probably find a workaround, but this isn't really worth the effort right now
    if len(job_title) > 11:
        errMsg = "Job title is longer than 11 characters. Please enter a shorter job title."
        errorFct(errMsg)
        exit(1)
    
    #saving cwd for later file system navigation
    wd = os.getcwd()
    
    #changing relative paths to absolute
    if os.path.isabs(pl_bin_path):
        dali_pl_full_path = os.path.join(pl_bin_path, "dali.pl")
    else:
        dali_pl_full_path = os.path.join(wd, pl_bin_path, "dali.pl")
        
    if not os.path.isabs(pdb_path):
        pdb_path = os.path.join(wd, pdb_path)
        
    if not os.path.isabs(out_path):
        out_path = os.path.join(wd, out_path)
    
    if not os.path.isabs(DAT_path):
        DAT_path = os.path.join(wd, DAT_path)

    pdb_name_list = get_pdb_name_list(pdb_path)
        
    
    #creating job dir
    job_path = os.path.join(out_path, job_title)
    create_dir_safely(job_path)
    
    
    #Dali's all-against-all feature behaves in unexpected ways. That's why I'm manually looping through the PDB list for generation of pairwise alignments.
    n = len(pdb_name_list)
    no_of_alns = sum(range(1,n))
    if n<2:
        errMsg = "Please provide more than one PDB file for alignment."
        errorFct(errMsg)
        exit(0)
    
    #we dont want redundant or trivial alignments: e.g. (1,1) is trivial, (1,2) = (2,1) are (effectively?) redundant
    
    counter = 0    
    for i in range(0, n):
        k = i+1
        for j in range(k, n):
            
            shell_input = []
            sub_job = "%s_%s" % (pdb_name_list[i], pdb_name_list[j])
            
            
            counter += 1 
            if verbosity>0:
                msg = "Calculating alignment %s. [%s/%s] " % (sub_job, counter, no_of_alns)
                print(msg)
            
            #creating ALN dir
            aln_path = os.path.join(job_path, "ALN")
            create_dir_safely(aln_path)
            
            #creating sub job dir
            sub_job_path = os.path.join(aln_path, sub_job)    
            create_dir_safely(sub_job_path)
                          
            pdb_full_path_query = os.path.join(pdb_path, pdb_name_list[i]+".pdb")
            pdb_full_path_sbjct = os.path.join(pdb_path, pdb_name_list[j]+".pdb")
            err_file_path = os.path.join(sub_job_path, "err_log")
            out_file_path = os.path.join(sub_job_path, "output_log")
            #/.../<dali_dir>/dali.pl --pdbfile1 <path> --pdbid1 <id> --pdbfile2 <path> --pdbid2 <id> --dat1 <path> --dat2 <path> --title <string> \
            # --outfmt "summary,alignments,equivalences,transrot" --clean
            shell_input = [dali_pl_full_path, '--pdbfile1', pdb_full_path_query, '--pdbid1', pdb_name_list[i], 
                           '--pdbfile2', pdb_full_path_sbjct, '--pdbid2', pdb_name_list[j], '--dat1', DAT_path,
                           '--dat2', DAT_path, '--title', job_title, '--outfmt', "summary,alignments,equivalences,transrot",
                           '--clean', '1', '>', out_file_path, '2', '>', err_file_path]

            #navigate to output dir in loop
            os.chdir(sub_job_path)
            
            #create file handles for stdout and stderr, only needed for Popen. rework maybe
            try:
                err_file = open(err_file_path, 'w+')
            except:
                errMsg = "Failed to open file %s." % err_file_path
                try:
                    close_file_safely(err_file, err_file_path, errMsg)
                except UnboundLocalError:
                    pass
                errorFct(errMsg)
                exit(1)
            try:
                out_file = open(out_file_path, 'w+')
            except:
                errMsg = "Failed to open file %s." % out_file_path
                try:
                    close_file_safely(out_file, out_file_path, errMsg)
                except UnboundLocalError:
                    pass
                errorFct(errMsg)
                exit(1)

            #subprocess for alignment generation
            #need to PIPE stdout and stderr for output forwarding    
            process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
    
            try:
                out, err = process.communicate(timeout=500) #timeout may need to be longer (or shorter), depending on size of alignments
                if out is not None and verbosity>1:
                    print(out.decode('ascii'))
            except subprocess.TimeoutExpired:
                process.kill()
                out, err = process.communicate()
                #print(out.decode('ascii'))
                errMsg = "Communication with dali.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt"
                close_file_safely(err_file, err_file_path, errMsg)
                close_file_safely(out_file, out_file_path, errMsg)
                #shell_err_fct(err.decode('ascii'))
                errorFct(errMsg)
                exit(1)
            
            #close files
            close_file_safely(err_file, err_file_path, "")
            close_file_safely(out_file, out_file_path, "")
            
            #change dir back to previous wd
            os.chdir(wd)
            
#wrapper function to calculate SPS of a (sub-)set of aligned sequences 
#(str, str, str, str, int, [str,str,...], int) -> (float, [{},{},...])
def calc_SPS_from_aln_sample(output_path, aln_file_path, job_title, dali_path, sample_size, valid_symbols, verbosity):
    
    job_out_path = os.path.join(output_path, job_title)
    create_dir_safely(job_out_path)
    
    job_data_out_path = os.path.join(job_out_path, "DATA")
    create_dir_safely(job_data_out_path)
    
    CSV_out_path = os.path.join(job_data_out_path, "CSV")
    create_dir_safely(CSV_out_path)
    
    #refactor write_aln_file_search_hits_to_csv
    internal_id_raw_seq_dict = write_aln_file_search_hits_to_csv(valid_symbols, aln_file_path, CSV_out_path, verbosity, sample_size)
    
    pdb_out_path = os.path.join(job_data_out_path, "PDB_lib")
    create_dir_safely(pdb_out_path)
    
    #refactor download_pdb_files_by_CSV
    aln_file_name = os.path.split(aln_file_path)[1]
    CSV_file_path = os.path.join(CSV_out_path, aln_file_name+".csv")
    internal_id_pdb_name_dict = download_pdb_files_by_CSV(CSV_file_path, pdb_out_path, "identity", verbosity)
    
    dat_out_path = os.path.join(job_data_out_path, "DAT_lib")
    create_dir_safely(dat_out_path)

    DALI_import_PDBs(dali_path, pdb_out_path, dat_out_path, verbosity)
    
    DALI_all_vs_all_query(dali_path, pdb_out_path, output_path, dat_out_path, job_title, verbosity)
    
    job_list_of_dict_lists = import_aln_files_in_job(job_title, output_path)
    
    SPS_dict_list, job_SPS = calc_SPS_scores_in_job(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict)
    
    print(SPS_dict_list)
    print(job_SPS)

    SPS_out_path = os.path.join(job_data_out_path, "SPS")
    create_dir_safely(SPS_out_path)
    SPS_csv_full_out_path = os.path.join(SPS_out_path, job_title+"_SPS.csv")
    #write SPS_dict_list to csv
    write_list_of_dicts_to_csv(SPS_dict_list, SPS_csv_full_out_path)
    
    return job_SPS, SPS_dict_list

#evaluate and create overview of search data
#void (str, str, str)    
def evaluate_search_data(csv_path, output_path, sort_by):
    
    ex_match_counter = 0
    total_counter = 0
    csv_file_list = get_file_names_in_path(csv_path)
    
    for csv_file_path in csv_file_list:
        df = pandas.read_csv(csv_file_path)
        new_df = pandas.DataFrame(columns=['original_fasta_header','sequence_identity','score','mismatches','gaps_opened'])
        row_list = []
        #remove next line, once trailing delimiters have been removed
        df = df.drop(['Unnamed: 16'], axis=1)
        #split df into one dataframe per fasta_header
        df_list = [x for _,x in df.groupby('original_fasta_header')]
        for split_df in df_list:
            if sort_by == "identity":
                split_df = split_df.sort_values(by=['sequence_identity'], ascending=False)
            elif sort_by == "score":
                split_df = split_df.sort_values(by=['score'], ascending=False)
            else:
                errMsg = "Can't sort by %s. Exiting." % sort_by
                errorFct(errMsg)
                exit(1)
            total_counter += 1
            #append highest score/identity row to df list
            if math.isclose(split_df.iloc[[0]]['sequence_identity'],1,abs_tol=10e-5):
                ex_match_counter += 1
            row_list.append(split_df.iloc[[0]].copy())
        #combine list into one df
        new_df = pandas.concat(row_list, ignore_index=True)
        plot_df = pandas.DataFrame(new_df[['sequence_identity','score']])
        plot_df['mean_identity'] = new_df['sequence_identity'].mean()
        csv_file_name = os.path.splitext(os.path.split(csv_file_path)[1])[0]
        #pandas.plot returns matplotlib.axes.Axes object
        axes_obj = plot_df.plot(title=csv_file_name)
        #we can get figure object from axes_obj using get_figure() method and save it to vector image
        axes_obj.get_figure().savefig(os.path.join(output_path,csv_file_name+".svg"))
    print("Exact matches: %s.Total number of sequences: %s." % (ex_match_counter, total_counter))

def update_local_pdb_index(raw_line, tmp_index_file_handle):

    tic_index_file = time.perf_counter()
    #write_pdb_index_file(pdb_path_list, index_file_full_path)
    toc_index_file = time.perf_counter()
    print(f"Indexing done in {toc_index_file - tic_index_file:0.4f} seconds.")
    
def repair_index(index_file_full_path):
    pass
   
#synchronizes copy of remote PDB archive with local copy
#void (str, int)  
def sync_pdb_copy(loc_db_out_path, index_file_full_path, tmp_index_file_full_path, verbosity):
    
    create_dir_safely(loc_db_out_path)
    #https://www.wwpdb.org/ftp/pdb-ftp-sites
    #rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb
    shell_input = []
    shell_input = ["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444", "rsync.rcsb.org::ftp_data/structures/divided/pdb/", loc_db_out_path]
    
    if verbosity>0:
        print("Synchronizing local PDB database with remote. This may take a while.")
    
    process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    if verbosity>0:
        print("Updating PDB index file...")

    try:
        tmp_index_file_handle = open(tmp_index_file_full_path, 'w')
    except:
        errMsg = "Cannot open file %s." % tmp_index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:    
        for line in process.stdout:
            if verbosity>0:
                print(line, end='\r')
            update_local_pdb_index(line, tmp_index_file_handle)
    except:
        errMsg1 = ""
        try:
            close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
        except:
            errMsg1 = "Error while closing temporary index file. Skipping."
        finally:
            repair_code = repair_index(index_file_full_path)
            errorFct("Error while updating index. Rewriting")
            if repair_code == 0:
                errorFct("Successfully repaired PDB index and updated FASTA accordingly. Try runnung sync again. Requires internet connection.")
            else:
                errMsg2 = "Repair code %s. We got a problem." % repair_code
                errorFct(errMsg1)
                errorFct(errMsg2)
    try:
        out, err = process.communicate()
    except:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with rsync subprocess failed. Killed it."
        shell_err_fct(err)
        errorFct(errMsg)
        exit(1)

#recursively get all pdb file paths in a specified path and return them as list
# str -> [str, str, ...]        
def get_pdb_path_list_recursively(loc_pdb_db_list):
    
    pdb_path_list = []
    
    for root, dirs, files in os.walk(loc_pdb_db_list):
        for file in files:
            if ".pdb" in file.lower() or ".ent" in file.lower():
                pdb_path_list.append(os.path.join(root, file))
                
    return pdb_path_list

#extracts chain numbers, pdbid and sequences from pdb file and writes them to FASTA file index
#also writes biopython output to dedicated log file
#consider parallelization
#void (str, file, file)
def write_loc_pdb_to_fasta_file(pdb_file_path, fasta_file_handle, biopython_log_file_handle):
    
    #make sure that this also works with uncompressed files
    try:
        pdb_file_handle = gzip.open(pdb_file_path, 'rt')
    except:
        errMsg = "Failed to open %s." % pdb_file_path
        errorFct(errMsg)
        exit(1)
    else:
        with gzip.open(pdb_file_path, 'rt') as pdb_file_handle:
            #supressing warnings for this line, specifically Bio.PDB.PDBExceptions.PDBConstructionWarning does NOT work in ANY way
            #not with BiopythonWarnings or BiopythonExperimentalWarnings either
            #since I am fed up with trying something that should work, but doesnt, I will resort to
            #redirecting all stdout and stderr streams to a different file
            with contextlib.redirect_stdout(biopython_log_file_handle):
                with contextlib.redirect_stderr(sys.stdout):
                    try:
                        chains = SeqIO.PdbIO.PdbSeqresIterator(pdb_file_handle)
                        #atom iterator, speed killer. but safer? implement option to choose, if slow or fast.
                        #also make sure to switch to atomiterator, if seqres not readable. find out exception
                        #chains = SeqIO.PdbIO.PdbAtomIterator(pdb_file_handle)
                    except:
                        errMsg = "Failed to parse chain data from %s." % pdb_file_path
                        errorFct(errMsg)
                        exit(1)
                    try:
                        SeqIO.write(chains, fasta_file_handle, "fasta")
                    except:
                        errMsg = "Error while writing %s to FASTA file." % pdb_file_path
                        errorFct(errMsg)
                        exit(1)

#creating fasta file from a local directory with PDB files recursively
#extends fasta file, if it already exists. now make this make sense
#void (str, str, str, str, int)
def create_fasta_and_index_file_from_pdb_collection(loc_pdb_db_path, fasta_file_full_path, index_file_full_path, biopython_log_full_path, verbosity):
    
    pdb_path_list = get_pdb_path_list_recursively(loc_pdb_db_path)
    #pdb_path_list = get_pdb_path_list_by_index_file

    seq = ""
    header = ""
    counter = 1
    number_of_pdbs = len(pdb_path_list)
    header_list = []

    if verbosity>0:
        print("Writing Fasta file...")
    try:
        #check integrity of/rewrite last entry
        fasta_file_handle = open(fasta_file_full_path, 'a')
    except:
        errMsg = "File %s could not be opened for writing." % fasta_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        biopython_log_file_handle = open(biopython_log_full_path)
    except:
        errMsg = "File %s could not be opened for writing." % biopython_log_full_path
        errorFct(errMsg)
        exit(1)

    try:
        for pdb_path in pdb_path_list:
            if counter>100:
                break
            write_loc_pdb_to_fasta_file(pdb_path, fasta_file_handle, biopython_log_file_handle)
            if verbosity>0:
                msg = "[%s\/%s]" % (str(counter), str(number_of_pdbs))
                print(msg, end="\r")
            counter += 1
    except:
        errMsg = "Error while writing to %s." % fasta_file_full_path
        print(traceback.format_exc())
        close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg)
        close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)
        
    
    close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg)
    close_file_safely(fasta_file_handle, fasta_file_full_path, '')


#creates diamond database from collection of pdb files in directory or index file
#void (str, str, str, int)        
def create_diamond_database(diamond_file_path, loc_pdb_db_path, db_out_path, verbosity):
    
    fasta_file_full_path = os.path.join(db_out_path, "pdb_database.fasta")
    index_file_full_path = os.path.join(db_out_path, "pdb_database.index")
    biopython_log_full_path = os.path.join(db_out_path, "biopython_log.txt")

    try:
        create_fasta_and_index_file_from_pdb_collection(loc_pdb_db_path, fasta_file_full_path, index_file_full_path, biopython_log_full_path, verbosity)
    except:
        errMsg = "Failed to create FASTA or INDEX file."
        errorFct(errMsg)
        exit(1)
    
    #<path_to_diamond_file> makedb --in <path_to_FASTA> --db <db_out_path>
    shell_input = []
    shell_input = [diamond_file_path, "makedb", "--in", fasta_file_full_path, "--db", db_out_path]
    
    if verbosity>0:
        print("Creating local database...")
    if verbosity>1:
        process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            print(line, end='')
    else:
        process = subprocess.Popen(shell_input, stdout=subprocess.DEVNULL)

    try:
        out, err = process.communicate()
        shell_err_fct(err)
    except:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with \"diamond makedb\" subprocess failed. Killed it."
        shell_err_fct(err)
        errorFct(errMsg)
        exit(1)


def main():

    tic = time.perf_counter()
    
    set_script_path()

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-p", "--createsearchplots", default=0, action="count", help="Creates search plots from CSV search metadata.", required=False)
    argParser.add_argument("-csv", "--csvdir", help="Path to CSV files with search metadata.", required=False)
    argParser.add_argument("-e", "--doitall", default=0, action="count", help="Do it all.", required=False)
    argParser.add_argument("-sync", "--syncdb", default=0, action="count", help="Sync local PDB database with remote.", required=False)
    argParser.add_argument("-dm", "--datamode", default=0, action="count",
                           help="data mode for gathering PDB files and other data from given FASTA alignment", required=False)
    argParser.add_argument("-db", "--locdb", default="PDB_db", help="Path to local PDB database.", required=False)
    argParser.add_argument("-af", "--alignmentfile", help="MSA file in FASTA format", required=False)
    argParser.add_argument("-bat", "--batchsearch", help="Path to MSA files in FASTA format", required=False)
    argParser.add_argument("-ss", "--samplesize", default=10, type = int,
                           help="Give sample size of random sequences taken from alignment file as integer. Default:10. 0 means full alignment.", required=False)
    argParser.add_argument("-sn", "--samplenumber",type = int, default = 1, help="Provide integer for how many times you want to sample. Default:1", required=False)
    argParser.add_argument("-dow", "--download", help="Path to CSV files to download PDB files from", required=False)
    argParser.add_argument("-bnch", "--benchmarkmode", default=0, action="count",
                           help="benchmark mode for generating alignments from PDB file database and measuring alignment scores", required=False)
    argParser.add_argument("-pdb", "--pdbdir", help="Path to directory with PDB files corresponding to proteins in MSA file", required=False)
    argParser.add_argument("-dat", "--datdir", help="Path to directory with DAT files", required=False)
    argParser.add_argument("-dali", "--dalidir", help="Path to DaliLite.v5/bin directory", required=False)
    argParser.add_argument("-jd", "--jobdir", default=0, action="count", help="Specify job dir. Currently no purpose.", required=False)
    argParser.add_argument("-al", "--align", default=0, action="count", help="Generate pairwise alignments for all PDB files specified under --pdbdir", required=False)
    argParser.add_argument("-im", "--importdat", default=0, action="count", help="Import all PDB files specified under --pdbdir into DAT format", required=False)
    argParser.add_argument("-ti", "--title", help="Job title. At most 11 characters.", required=False)
    argParser.add_argument("-o", "--output", help="Output path", required=False)
    argParser.add_argument("-s", "--AAscoring", choices=["BLOSUM62", "PAM"], help="Select type of AA-scoring: (BLOSUM, PAM).", required=False)
    argParser.add_argument("-a", "--alphabet", choices=["AA", "DNA/RNA"], default="AA", help="Select alphabet: (AA, DNA/RNA). Default: AA", required=False)
    argParser.add_argument("-na", "--nonalphabet", default="-.", help="Select non-alphabet. Default: -.", required=False)
    argParser.add_argument("-v", "--verbose", default=0, action="count",
                           help="Print progress to terminal. No effect on error logging.1: Basic output. 2: Doesn't work atm. Redirect DALI output.",
                           required=False)
    argParser.add_argument("-tdb", "--testdb", action="count", default=0, help="Testing diamond database creation.", required=False)
    argParser.add_argument("-dia", "--diamondfile", help="Path to diamond program file.", required=False)

    args = argParser.parse_args()
    
    if args.testdb>0:
        create_diamond_database(args.diamondfile, args.locdb, args.output, args.verbose)

    if args.syncdb>0:
        sync_pdb_copy(args.locdb, args.verbose)

    valid_symbols = sel_alphabet(args.alphabet)
    invalid_symbols = set_non_alphabet(args.nonalphabet)

    if args.doitall > 0:
        calc_SPS_from_aln_sample(args.output, args.alignmentfile, args.title, args.dalidir, args.samplesize, valid_symbols, args.verbose)

    if args.datamode > 0:
        if args.alignmentfile is not None:
            _ = write_aln_file_search_hits_to_csv(valid_symbols, args.alignmentfile, args.output, args.verbose, args.samplesize)
        if args.createsearchplots > 0:
            evaluate_search_data(args.csvdir, args.output, "identity")
        #wrap in function    
        elif args.batchsearch is not None:
            aln_file_list = get_file_names_in_path(args.batchsearch)
            for file_path in aln_file_list:
                csv_path = os.path.join(args.output,os.path.split(file_path)[1]+".csv")
                if os.path.exists(csv_path):
                    errMsg = "File %s already exists. Skipping." % csv_path
                    errorFct(errMsg)
                    continue
                _ = write_aln_file_search_hits_to_csv(valid_symbols, file_path, args.output, args.verbose, args.samplesize)
        if args.download is not None:        
            download_pdb_files_by_CSV(args.download, args.output)
    if args.benchmarkmode > 0:
        if args.importdat > 0:
            DALI_import_PDBs(args.dalidir, args.pdbdir, args.datdir, args.verbose)
        if args.align > 0:
            DALI_all_vs_all_query(args.dalidir, args.pdbdir, args.output, args.datdir, args.title, args.verbose)
        if args.jobdir > 0:
            job_list_of_dict_lists = import_aln_files_in_job(args.title, args.output)
            print(len(job_list_of_dict_lists))
            print(job_list_of_dict_lists)
    

    toc = time.perf_counter()
    print(f"Done in {toc - tic:0.4f} seconds.")
    exit(0)

if __name__ == "__main__":
    main()