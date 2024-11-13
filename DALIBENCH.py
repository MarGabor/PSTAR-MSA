#DALIBENCH
#Author: Marcel Gabor
#Institute for Mathematics and Computer Science, University of Greifswald
#11/04/2024

import argparse
import sys
from datetime import datetime
import re
import os
import time
import pandas
import math
import random
import requests
import subprocess
import gzip
import shutil
#import svg_stack
from itertools import tee
import lxml
import traceback
import contextlib
from rcsbsearchapi.search import Query, SequenceQuery
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
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
        valid_symbols = ["A","G","C","T","U"]
    
    return valid_symbols

#takes string of invalid symbols and returns a list
#this function mostly exists only for factoring consistency
#str -> [str,str,...]
def set_non_alphabet(non_alphabet):
    
    return list(non_alphabet)

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

#removing gaps and other symbols not in "valid_symbols" from given alignment sequence data
#preserves case
#(str, str) -> str
def clean_seq(seq, valid_symbols):

    regex_str = build_regex_for_seq_cleaning_whitelist(valid_symbols)

    cleaned_seq = re.sub(regex_str,'',seq)
    
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
# [sequence, header, number]
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
#(str) -> [str,str,...]
def get_file_paths_in_dir(path):

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
#((int,int),(str,str)) -> set( ((),()) )
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
            #define a set of expected symbols or set of gap/unknown symbols
            sbjct_gap = False
            query_gap = False
            if sbjct[i] == "-" or sbjct[i] == ".":
                sbjct_gaps += 1
                sbjct_gap = True
            if query[i] == "-"  or sbjct[i] == ".":
                query_gaps += 1
                query_gap = True
            if sbjct_gap or query_gap:
                continue
            MSA_ind_set.add(((sequence_ids[0],i-query_gaps),(sequence_ids[1],i-sbjct_gaps)))
        except IndexError:
            break
    
    return MSA_ind_set
        

#function for naive approach to maximizing intersection between two sets by shifting all element values in set_2 by some constant
#sets need to contain mutable elements
# (set(),set()) -> (set(),int)
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
def calc_job_SPS(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict):
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

#import function for DALI outputs. one pass, line per line. might need rework,
#because it's hard to read and looks ugly. hopefully it's fast(-ish), at least.
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
        
        
        #we're not using structural equivalences or transrot at the moment,
        #thus we can end parsing once we reached this point
        if aln_counter > 0 and (line.startswith("# S") or line.startswith("# T")):
            overall_query = ''.join(query_list)
            overall_sbjct = ''.join(sbjct_list)
            aln_dict_list[aln_counter-1]['DALI_query_seq'] = overall_query
            aln_dict_list[aln_counter-1]['DALI_sbjct_seq'] = overall_sbjct
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
            #if chain identifiers are used for naming of dirs
            if len(head_tail2[1]) == 11:
                query_name = head_tail2[1][0:5]
                sbjct_name = head_tail2[1][6:11]
            elif len(head_tail2[1]) == 9:
                query_name = head_tail2[1][0:4]
                sbjct_name = head_tail2[1][5:9]
            else:
                errMsg = "Sub job dir name %s is not of expected format." % head_tail[0]
                errorFct(errMsg)
                raise NameError
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
    #all upper case alphabet letters
    regex_str_list = valid_symbols.copy()
    valid_symbols_lower = list(map(str.lower, valid_symbols.copy()))
    #and all lower case alphabet letters
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
#irrelevant function
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
    
    #delete next line, if everything runs fine
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
        
        #pass regex_str to clean_seq() again, if something goes wrong
        seq[0] = clean_seq(str(seq[0]), valid_symbols)
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

    return internal_id_raw_seq_dict

#function to write error log. for errors encountered during subprocesses involving the shell
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
    
    ent_gz = False
    pdb_path_list = get_file_paths_in_dir(pdb_path)
    pdb_name_list = []
    for pdb_path in pdb_path_list:
        head, file_name = os.path.split(pdb_path)
        #this loop works for .pdb and for .ent.gz formats
        while "." in file_name:
            file_name = os.path.splitext(file_name)[0]
        #standard rsync pdb database has file names of format "pdbxxxx.ent.gz"
        if len(file_name) == 7:
            file_name = file_name[3:8]
            ent_gz = True
        pdb_name_list.append(file_name)
        
    return pdb_name_list, ent_gz

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run import.pl of DaliLite.v5 to import PDBs in a given dir to DAT file format used by DALI
#void (str, str, str, int)
def DALI_import_PDBs(pl_bin_path, pdb_path, DAT_path, verbosity):

    pdb_name_list, ent_gz = get_pdb_name_list(pdb_path)
        
    import_pl_full_path = os.path.join(pl_bin_path, "import.pl")

    #create dir for output and err
    create_dir_safely(os.path.join(DAT_path, "err_logs"))
    
    create_dir_safely(os.path.join(DAT_path, "out_logs"))
    
    #this loop is necessary for clean error logging, since DaliLite.v5 seems to have bad logging
    for pdb_name in pdb_name_list:
        
        shell_input = []
        if ent_gz:
            pdb_full_path = os.path.join(pdb_path, "pdb"+pdb_name+".ent.gz")
        else:
            pdb_full_path = os.path.join(pdb_path, pdb_name+".pdb")
        err_file_path = os.path.join(DAT_path, "err_logs/", pdb_name)
        out_file_path = os.path.join(DAT_path, "out_logs/", pdb_name)
        
        #/.../<dali_dir>/import.pl --pdbfile <path> --pdbid <id> --dat <path> --verbose --clean
        shell_input = [import_pl_full_path, '--pdbfile', pdb_full_path, '--pdbid', pdb_name.upper(), '--dat', DAT_path, '--clean', '1', '>', out_file_path, '2', '>', err_file_path]
        
        #create file handles for stdout and stderr
        try:
            err_file = open(err_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % err_file_path
            errorFct(errMsg)
            exit(1)
        try:
            out_file = open(out_file_path, 'w+')
        except:
            errMsg = "Failed to open file %s." % out_file_path
            errorFct(errMsg)
            exit(1)

        #subprocess for importing PDBs
        #shell input already sends stdout to out_file_path and stderr to err_file_path
        #probably just pipe output in the next line    
        process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
        
        try:
            out, err = process.communicate()
            #its really just a lot to print...
            #if out is not None and verbosity>1:
                #print(out)
        except:
            process.kill()
            out, err = process.communicate()
            errMsg = "Communication with import.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt."
            close_file_safely(err_file, err_file_path, errMsg)
            close_file_safely(out_file, out_file_path, errMsg)
            shell_err_fct(err)
            errorFct(errMsg)
            if out is not None:
                print(out)
            exit(1)
            
        #close files
        close_file_safely(err_file, err_file_path, "")
        close_file_safely(out_file, out_file_path, "")
    
    #rewrite a sufficient check of whether import has been successful
    for file_name in os.listdir(DAT_path):
        if file_name.endswith(".d"):
            errMsg = "Failed to convert PDB %s to DAT format. Exiting." % file_name
            errorFct(errMsg)
            exit(1)

#http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html
#creates subprocess to run dali.pl of DaliLite.v5. Uses "pairwise" functionality of DALI for better control over what alignments are being made.
#writes stderr and stdout to files for each job and alignment. improved output dir structure, since DaliLite.v5 doesn't seem to support specifying an output dir
#dir structure: <jobtitle>/ALN/<queryPDBname_sbjctPDBname>/
#if chain_list is passed, only the specified chains of the PDBs will be used in alignment       
#void (str,str,str,str,str,int,[str,str,...])
def DALI_all_vs_all_query(pl_bin_path, pdb_path, out_path, DAT_path, job_title, verbosity, chain_list=[]):

    #for some weird implementation reasons (probably Fortran though) DALI doesn't generate alignment files for job titles longer than 12 chars
    #and generates only empty alignments for job titles exactly 12 chars long
    #this is a DALI problem and can't be changed at the moment
    #i could, however probably find a workaround, but this isn't really worth the effort right now
    #maybe it's also just because job_title needs to be passed with "" around it.
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
    
    chains_used = False
    #if optional argument chain_list is passed, use that instead    
    if len(chain_list) != 0:
        pdb_name_list = list(chain_list).copy()
        chains_used = True
    else:
        pdb_name_list = get_pdb_name_list(pdb_path)
        
    
    #creating job dir
    job_path = os.path.join(out_path, job_title)
    create_dir_safely(job_path)
    
    #Dali's all-against-all feature behaves in unexpected ways. 
    #That's why manual looping through the PDB list for generation of pairwise alignments.
    #this is also beneficial, since it allows full control over what alignments are to be generated
    #and individual sub job error logging
    n = len(pdb_name_list)
    no_of_alns = sum(range(1,n))
    if n<2:
        errMsg = "Please provide more than one PDB file for alignment."
        errorFct(errMsg)
        exit(0)
    
    #we dont want redundant or trivial alignments:
    # e.g. (1,1) is trivial, (1,2) = (2,1) are (effectively?) redundant
    
    counter = 0    
    for i in range(0, n):
        k = i+1
        for j in range(k, n):
            
            if chains_used:
                try:
                    chain_id1 = str(pdb_name_list[i][0:4])+str(pdb_name_list[i][5])
                    chain_id2 = str(pdb_name_list[j][0:4])+str(pdb_name_list[j][5])
                except IndexError:
                    errMsg = "Chain identifier %s or %s is shorter than 6 chars." % (pdb_name_list[i], pdb_name_list[j])
                    errorFct(errMsg)
                    raise

            shell_input = []
            sub_job = "%s_%s" % (chain_id1, chain_id2)
            
            counter += 1 
            if verbosity>0:
                msg = "Calculating alignment %s. [%s/%s] " % (sub_job, counter, no_of_alns)
                print(msg, end="\r")
            
            #creating ALN dir
            aln_path = os.path.join(job_path, "ALN")
            create_dir_safely(aln_path)
            
            #creating sub job dir
            sub_job_path = os.path.join(aln_path, sub_job)    
            create_dir_safely(sub_job_path)
                          
            err_file_path = os.path.join(sub_job_path, "err_log")
            out_file_path = os.path.join(sub_job_path, "output_log")
            if chains_used:
                #/.../<dali_dir>/dali.pl --cd1 <chain_id1> --cd2 <chain_id2> --dat1 <path> --dat2 <path> --title <string> \
                # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
                shell_input = [dali_pl_full_path, '--cd1', chain_id1, '--cd2', chain_id2, 
                           '--dat1', DAT_path, '--dat2', DAT_path, '--title', job_title,
                           '--outfmt', "summary,alignments,equivalences,transrot",
                           '--clean', '1', '>', out_file_path, '2', '>', err_file_path]
            else:
                #/.../<dali_dir>/dali.pl --pdbfile1 <path> --pdbid1 <id> --pdbfile2 <path> --pdbid2 <id> --dat1 <path> --dat2 <path> --title <string> \
                # --outfmt "summary,alignments,equivalences,transrot" --clean 1> <out_log_file> 2> <err_log_file>
                #technically --dat argument makes no sense here
                pdb_full_path_query = os.path.join(pdb_path, pdb_name_list[i]+".pdb")
                pdb_full_path_sbjct = os.path.join(pdb_path, pdb_name_list[j]+".pdb")    
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
                errorFct(errMsg)
                exit(1)
            try:
                out_file = open(out_file_path, 'w+')
            except:
                errMsg = "Failed to open file %s." % out_file_path
                errorFct(errMsg)
                exit(1)

            #subprocess for alignment generation
            #need to PIPE stdout and stderr for output forwarding    
            process = subprocess.Popen(shell_input, stdout=out_file, stderr=err_file)
    
            try:
                #timeout may need to be longer (or shorter), depending on size of alignments
                #for the calculation of pairwise alignments of up to a few chains
                #500 seconds seems to be reasonable though
                out, err = process.communicate()
                #if out is not None and verbosity>1:
                    #print(out)
            except subprocess.TimeoutExpired:
                process.kill()
                out, err = process.communicate()
                errMsg = "Communication with dali.pl subprocess timed out. Killed it. Forwarding error to shell_err_log.txt"
                close_file_safely(err_file, err_file_path, errMsg)
                close_file_safely(out_file, out_file_path, errMsg)
                shell_err_fct(err)
                errorFct(errMsg)
                if out is not None:
                    print(out)
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
    
    SPS_dict_list, job_SPS = calc_job_SPS(job_list_of_dict_lists, internal_id_pdb_name_dict, internal_id_raw_seq_dict)
    
    print(SPS_dict_list)
    print(job_SPS)

    SPS_out_path = os.path.join(job_data_out_path, "SPS")
    create_dir_safely(SPS_out_path)
    SPS_csv_full_out_path = os.path.join(SPS_out_path, job_title+"_SPS.csv")
    #write SPS_dict_list to csv
    write_list_of_dicts_to_csv(SPS_dict_list, SPS_csv_full_out_path)
    
    return job_SPS, SPS_dict_list

#function to take a list of svg images and concatenates them to a large svg image
#void ([str,str,...],str)    
def append_svg_images(svg_path, out_img_full_path):
    
    svg_full_path_list = []
    svg_file_name_list = os.listdir(svg_path)
    for svg_file_name in svg_file_name_list:
        svg_full_path = os.path.join(svg_path, svg_file_name)
        if os.path.isfile(svg_full_path) is False:
            continue
        svg_full_path_list.append(svg_full_path)

    doc = svg_stack.Document()
    
    layout1 = svg_stack.HBoxLayout()
    for svg_full_path in svg_full_path_list:
        layout1.addSVG(svg_full_path, alignment=svg_stack.AlignCenter)
    doc.setLayout(layout1)
    doc.save(out_img_full_path)

#evaluate and create overview of search data
#void (str, str, str)    
def evaluate_rcsb_search_data(csv_path, output_path, sort_by):
    
    ex_match_counter = 0
    total_counter = 0
    csv_file_list = get_file_paths_in_dir(csv_path)
    
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

#evaluate diamond blastp TSV files and create overview plots for structure availability
#optimizable (?) iloc is massively slowing things down.
#try with df.set_value(), almost 3 times faster probably, but its deprecated. annoying  
# void (str,str,str,int)    
def evaluate_diamond_search_data(tsv_path, output_path, sort_by, verbosity):
    
    tsv_file_list = get_file_paths_in_dir(tsv_path)
    exact_match_dict = {}
    file_counter = 1
    no_of_files = len(tsv_file_list)

    for tsv_file_path in tsv_file_list:
        if verbosity>0:
            msg = "[%s/%s] Plotting %s..." % (file_counter, no_of_files, tsv_file_path)
            print(msg)
        tsv_ex_match_counter = 0
        tsv_total_counter = 0
        df = pandas.read_csv(tsv_file_path, sep='\t')
        new_df = pandas.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen'])
        row_list = []
        df['is_exact_match'] = 0
        ex_match_col_num = df.columns.get_loc('is_exact_match')
        mismatch_col_num = df.columns.get_loc('mismatch')
        len_col_num = df.columns.get_loc('length')
        qlen_col_num = df.columns.get_loc('qlen')
        #split df into one dataframe per fasta_header
        df_list = [x for _,x in df.groupby('qseqid')]
        for split_df in df_list:
            if sort_by == "identity":
                split_df = split_df.sort_values(by=['pident'], ascending=False)
            elif sort_by == "score":
                split_df = split_df.sort_values(by=['bitscore'], ascending=False)
            else:
                errMsg = "Can't sort by %s. Exiting." % sort_by
                errorFct(errMsg)
                exit(1)
            tsv_total_counter += 1
            #append highest score/identity row to df list
            #if math.isclose(split_df.iloc[[0]]['pident'],1,abs_tol=10e-1):
             #   ex_match_counter += 1
            if int(split_df.iloc[0,mismatch_col_num]) == 0 and int(split_df.iloc[0,len_col_num]) == int(split_df.iloc[0,qlen_col_num]):
                tsv_ex_match_counter += 1
                split_df.iloc[0,ex_match_col_num] = 100
            row_list.append(split_df.iloc[[0]].copy())
        tsv_file_name = os.path.splitext(os.path.split(tsv_file_path)[1])[0]
        exact_match_dict[tsv_file_name] = (tsv_ex_match_counter/tsv_total_counter)*100
        #combine list into one df
        new_df = pandas.concat(row_list, ignore_index=True)
        plot_df = pandas.DataFrame(new_df[['pident','is_exact_match']])
        plot_df['mean_identity'] = new_df['pident'].mean()
        plot_df['exact_match_perc'] = (tsv_ex_match_counter/tsv_total_counter)*100
        #pandas.plot returns matplotlib.axes.Axes object
        plot_df_len = len(plot_df)
        xtick_list = list(range(0,plot_df_len,math.floor(plot_df_len/5)))
        axes_obj = plot_df[['pident','is_exact_match',]].plot(title=tsv_file_name, kind="bar", yticks=[0,2.5,5,7.5,10,15,20,25,30,40,50,60,70,80,90,100],
                                                             xticks=xtick_list)
        plot_df[['mean_identity','exact_match_perc']].plot(kind="line", ax=axes_obj)
        #we can get figure object from axes_obj using get_figure() method and save it to vector image
        axes_obj.get_figure().savefig(os.path.join(output_path, tsv_file_name+".svg"))
        #matplot.pyplot.figure objects need to be explicitly closed. otherwise they are all open at the same time
        matplotlib.pyplot.close()
        file_counter += 1
    exact_perc_df = pandas.DataFrame.from_dict(data=exact_match_dict, orient='index', columns=['ex_match_perc'])
    exact_perc_df['mean_exact_perc'] = exact_perc_df['ex_match_perc'].mean()
    axes_obj = exact_perc_df[['ex_match_perc']].plot(title="Exact match distribution among test alignments", kind="bar")
    exact_perc_df[['mean_exact_perc']].plot(kind="line", ax=axes_obj, color='orange')
    axes_obj.set_xticklabels(exact_perc_df.index.values, rotation = 45)
    axes_obj.tick_params(axis='x', which='major', labelsize=2)
    axes_obj.get_figure().savefig(os.path.join(output_path, "mean_ex_match_perc.svg"))
    matplotlib.pyplot.close()

#recursively get all pdb file paths in a specified path and return them as list
# str -> [str, str, ...]        
def get_pdb_path_list_recursively(loc_pdb_db_path):
    
    pdb_path_list = []
    
    for root, dirs, files in os.walk(loc_pdb_db_path):
        for file in files:
            if ".pdb" in file.lower() or ".ent" in file.lower():
                pdb_path_list.append(os.path.join(root, file))
                
    return pdb_path_list

#extract pdb ID from pdb file path and write it to index file as line
# (str,file_handle) -> (str)
def write_pdb_index_entry_from_file_path(raw_line, index_file_handle):

    #this might need to be generalized more for other structure formats downloaded by rsync
    pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', raw_line)
    pdb_id = pdb_id_match.group(2)
    index_line = pdb_id+"\n"
    index_file_handle.write(index_line)
    
    return pdb_id

#wrapper function to write pdb index entries according to list of paths to pdb files
# void (str, [str,str,...])
def write_index_file_from_file_path_list(old_index_file_full_path, pdb_path_list):
    
    try:
        old_index_file_handle = open(old_index_file_full_path, "w")
    except:
        errMsg = "Cannot open file %s." % old_index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        old_index_file_handle.truncate()
        for pdb_path in pdb_path_list:
            write_pdb_index_entry_from_file_path(pdb_path, old_index_file_handle)
    except:
        errMsg = "Error while writing PDB IDs to index file %s." % old_index_file_full_path
        close_file_safely(old_index_file_handle, old_index_file_full_path, errMsg)
        exit(1)
            
    close_file_safely(old_index_file_handle, old_index_file_full_path, "")

#writes an index file from the contents of a set
# void (str, set())   
def write_index_file_from_set(index_file_full_path, pdb_id_set):
    
    try:
        index_file_handle = open(index_file_full_path, "w")
    except:
        errMsg = "Cannot open file %s." % index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        index_file_handle.truncate()
        for pdb_id in pdb_id_set:
            index_line = str(pdb_id).lower()+"\n"
            index_file_handle.write(index_line)
    except:
        errMsg = "Error while writing PDB IDs to index file %s." % index_file_full_path
        close_file_safely(index_file_handle, index_file_full_path, errMsg)
        exit(1)
            
    close_file_safely(index_file_handle, index_file_full_path, "")
        
    

#read index file and return set with pdb IDs
# (file_handle) -> set()    
def read_index_file(old_index_file_handle):
    
    old_index_set = set()

    pdb_ids = old_index_file_handle.readlines()
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.rstrip().lower()
        old_index_set.add(pdb_id)
        
    return old_index_set

#function to be called before editing the database fasta file
#creates a copy of the fasta file in the same location
# void (str,str)
def backup_fasta_file(fasta_file_full_path):
     
    dest_backup = os.path.join(os.path.split(fasta_file_full_path)[0], "pdb_db_backup.fasta")
    cp_process = subprocess.Popen(['cp', fasta_file_full_path, dest_backup], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    try:
        out, err = cp_process.communicate()
    except:
        cp_process.kill()
        out, err = cp_process.communicate()
        shell_err_fct(err)
        raise
    
    return dest_backup

#called when editing fasta file has been completed successfully
#double checks to see, if original fasta file is present in the same location
# void (str,str)    
def remove_fasta_backup(dest_backup, fasta_file_full_path):
    
    #just to be extra safe
    if os.path.exists(fasta_file_full_path):
        os.remove(dest_backup)

#called when exception is raised during editing of fasta file
#simply overwrites the (potentially) broken, edited fasta file with the backup and renames it
# void (str,str)        
def restore_fasta_file_from_backup(dest_backup, fasta_file_full_path):
    
    mv_process = subprocess.Popen(['mv', dest_backup, fasta_file_full_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    try:
        out, err = mv_process.communicate()
    except:
        mv_process.kill()
        out, err = mv_process.communicate()
        shell_err_fct(err)
        raise

#function for removal of specific entries in fasta file
#unfortunately theres no perfect way to do it at the moment
#lines need to be rewritten when we want to delete some
#but maybe the biopython parser is what makes it slow
#just try it later and see how performance does
# void (set(),str)    
def remove_fasta_db_entries_by_header(removal_set, diamond_db_path):

    fasta_file_full_path = os.path.join(diamond_db_path, "pdb_db.fasta")
    with open(fasta_file_full_path, 'r+') as fasta_file:
        for line in fasta_file:
            if line[0] == ">":
                if line[1:5].lower() in removal_set:
                    pass

#takes a SeqRecord iterator, if chain.id does not have expected format and rewrites the whole iterator
#fixing those entries that don't adhere to XXXX:X format in the header by getting pdb id from file name
#if not even the chain identifier can be read (accessing first char of chain.id throws IndexError)
#the header will simply contain the pdb id without chain information
# (iterator, str) -> (iterator)                
def update_chains_iterator(iterator, pdb_file_path):
         
    new_chains_list = []

    for chain in iterator:
        try:
            chain_id_str = str(chain.id)
            if chain_id_str[4] != ":":
                pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', pdb_file_path)
                pdb_id = pdb_id_match.group(2)
                new_chain_id = pdb_id+":"+chain_id_str[0].upper()
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
            else:
                #new_chain_id = str(chain.id)
                #new_chain_seq = Seq.Seq(str(chain.seq))
                #new_chain = SeqIO.SeqRecord(new_chain_seq, id=new_chain_id)
                new_chains_list.append(chain)
        except IndexError:
            pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', pdb_file_path)
            pdb_id = pdb_id_match.group(2)
            try:
                new_chain_id = (pdb_id+":"+chain_id_str[0]).upper()
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
            except IndexError:
                new_chain_id = pdb_id
                new_chain_seq = str(chain.seq)
                new_chain = SeqIO.SeqRecord(Seq.Seq(new_chain_seq), id=new_chain_id)
                new_chains_list.append(new_chain)
                
    return iter(new_chains_list)
                
#extracts chain numbers, pdbid and sequences from pdb file and writes them to FASTA file
#also writes biopython output to dedicated log file
#consider parallelization of pdb file parser, it's very slow like this
#i suspect the biopython pdb parser     
#void (str, file_handle, file_handle)
def write_loc_pdb_to_fasta_file(pdb_file_path, fasta_file_handle, biopython_log_file_handle):
    
    #make sure that this also works with uncompressed files
    try:
        pdb_file_handle = gzip.open(pdb_file_path, 'rt')
    except:
        errMsg = "Failed to open %s." % pdb_file_path
        errorFct(errMsg)
        exit(1)
    #supressing warnings for this line, specifically Bio.PDB.PDBExceptions.PDBConstructionWarning does NOT work in ANY way
    #not with BiopythonWarnings nor with BiopythonExperimentalWarnings
    #redirecting all stdout and stderr streams to a different file instead
    with contextlib.redirect_stdout(biopython_log_file_handle):
        with contextlib.redirect_stderr(sys.stdout):
            try:
                chains = SeqIO.PdbIO.PdbSeqresIterator(pdb_file_handle)
                #atom iterator is speed killer.
                #IMPORTANT REMARK: AtomIterator ends up producing FASTA files down the line that
                #confuse DIAMOND. The current FTP path of rsync also downloads DNA files.
                #something about the handling of HETATOMS in AtomIterator causes problems (wrong Biopython parser).
                #stick to using SeqResIterator for now 
                #chains = SeqIO.PdbIO.PdbAtomIterator(pdb_file_handle)  
            except:
                errMsg = "Failed to parse chain data from %s." % pdb_file_path
                close_file_safely(pdb_file_handle, pdb_file_path, errMsg)
                errorFct(errMsg)
                return
    
    update = False
    #copies of iterators are required for a) counting length, b) looping to check for incomplete headers
    #c) passing to update_chains_iterator() function 
    chains, chains_copy, chains_copy_2, chains_copy_3 = tee(chains, 4)
    chains_choice = chains

    #if parser fails to get valid chain ID, try to get PDB ID from file name
    #this is not a particularly elegant implementation, but at least it requires only
    #a minimal amount of rewriting iterators  
    for chain in chains_copy_2:
        try:
            chain_id_str = str(chain.id)
            if chain_id_str[4] != ":":
                new_chains = update_chains_iterator(chains_copy_3, pdb_file_path)
                update = True
                break
        except IndexError:
            new_chains = update_chains_iterator(chains_copy_3, pdb_file_path)
            update = True
            break
        
    if update:
        #check len of old iterator
        old_iter_len = 0
        for chain in chains_copy:
            old_iter_len +=1
        new_chains, new_chains_copy = tee(new_chains)
        #check len of new iterator
        new_iter_len = 0
        for chain in new_chains_copy:
            new_iter_len += 1
        if old_iter_len == new_iter_len:
            chains_choice = new_chains
        else:
            errMsg = "Warning: Failed to update chains for %s." % (pdb_file_path)
            errorFct(errMsg)

    try:
        SeqIO.write(chains_choice, fasta_file_handle, "fasta")
    except:
        #try:
            #   repair_recent_fasta_entries(fasta_file_handle, chains)
        #except:
            #   errMsg = "Error while attempting to repair most recent FASTA entries."
            #  errorFct(errMsg)
            # raise Exception
        print(traceback.format_exc())
        errMsg = "Error while writing %s to FASTA file. Skipping." % pdb_file_path
        close_file_safely(pdb_file_handle, pdb_file_path, errMsg)
        errorFct(errMsg)
    close_file_safely(pdb_file_handle, pdb_file_path, "")

#takes a set of pdb IDs and appends existing fasta file with the elements of this set
#throws FileNotFoundException, if indexed PDB ID is not in local PDB database
# void (set(),str,str,str,int)                
def add_fasta_db_entries(added_set, fasta_file_full_path, loc_pdb_db_path, biopython_log_full_path, verbosity):
    
    if verbosity>0:
        print("Appending FASTA file...")
    try:
        dest_backup = backup_fasta_file(fasta_file_full_path)
    except:
        head, tail = os.path.split(fasta_file_full_path)
        errMsg = "Failed to back up %s. Insufficient writing priviliges in %s?" % (tail, head)
        errorFct(errMsg)
        exit(1)
    
    try:
        biopython_log_file_handle = open(biopython_log_full_path)
    except:
        errMsg = "File %s could not be opened for writing." % biopython_log_full_path
        errorFct(errMsg)
        exit(1)

    try:
        counter = 1
        number_of_pdbs = len(added_set)
        with open(fasta_file_full_path, 'a') as fasta_file_handle:
            for pdb_id in added_set:
                if verbosity>0:
                    msg = "[%s/%s]              " % (str(counter), str(number_of_pdbs))
                    print(msg, end="\r")
                pdb_file_path = os.path.join(loc_pdb_db_path, pdb_id[1:3], "pdb"+pdb_id+".ent.gz")
                write_loc_pdb_to_fasta_file(pdb_file_path, fasta_file_handle, biopython_log_file_handle)
                counter += 1
    except:
        try:
            restore_fasta_file_from_backup(dest_backup, fasta_file_full_path)
        except:
            print(traceback.format_exc())
            errMsg = "Error while restoring FASTA file from backup %s." % dest_backup
            errorFct(errMsg)
        finally:
            print(traceback.format_exc())
            errMsg = "Error while writing to %s." % fasta_file_full_path
            close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg)
            close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg)
            errorFct(errMsg)
            exit(1)

    remove_fasta_backup(dest_backup, fasta_file_full_path)
    close_file_safely(biopython_log_file_handle, biopython_log_full_path, "")

#updates index file by reading the old index file and a new index file.
#the new index file is generated by recursively listing all pdb files that are in the local database.
#this ensures that the up-to-date index really only contains IDs that exist locally.    
#the index entries are then added to two different sets. the resulting difference set is subsequently used
#to write only the new entries to the fasta file (s. function "sync_pdb_copy()")    
#this funciton gets called for every rsync call
# (str,str,str,str) -> set()    
def update_pdb_index(loc_pdb_db_path, diamond_db_path, tmp_index_file_full_path, old_index_file_full_path):

    tmp_index_set = set()
    old_index_set = set()
    pdb_id = ""

    try:
        tmp_index_file_handle = open(tmp_index_file_full_path, 'w')
    except:
        errMsg = "Cannot open file %s." % tmp_index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:      
        pdb_path_list = get_pdb_path_list_recursively(loc_pdb_db_path)
        for pdb_path in pdb_path_list:
            pdb_id = write_pdb_index_entry_from_file_path(pdb_path, tmp_index_file_handle)
            tmp_index_set.add(pdb_id.lower())
    except:
        errMsg = "Error while writing index file %s." % tmp_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
    
    #if old index file does not exist, then write it  
    if os.path.exists(old_index_file_full_path) is False:
        write_index_file_from_file_path_list(old_index_file_full_path, pdb_path_list)
     
    try:
        old_index_file_handle = open(old_index_file_full_path, "r")
    except:
        errMsg = "Cannot open file %s." % old_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)
            
    try:
        old_index_set = read_index_file(old_index_file_handle)
    except:
        errMsg = "Error while reading index file %s." % old_index_file_full_path
        close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, errMsg)
        close_file_safely(old_index_file_handle, old_index_file_full_path, errMsg)
        exit(1)
        
    close_file_safely(tmp_index_file_handle, tmp_index_file_full_path, "")
    close_file_safely(old_index_file_handle, old_index_file_full_path, "")
  
    #for all entries that have been removed in recent sync
    #this will probably be a rare occurence
    #removal_set = old_index_set.difference(tmp_index_set)
    #if len(removal_set) != 0:
     #   remove_fasta_db_entries_by_header(removal_set)
        
    #set of all entries that got added in recent sync
    added_set = tmp_index_set.difference(old_index_set)
    
    return added_set


#replaces old index file with new one
#called after appending of fasta file finished successfully
# void (str,str)
def replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path):
    
    mv_process = subprocess.Popen(['mv', tmp_index_file_full_path, old_index_file_full_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        out, err = mv_process.communicate()
    except:
        mv_process.kill()
        out, err = mv_process.communicate()
        shell_err_fct(err)
        errMsg = "Communcation with mv_process to replace old index file failed. Killed it."
        errorFct(errMsg)
        exit(1)

#extra. biopython should handle fasta write errors just fine.
#the idea would be to traverse the file backwards with the cursor and delete everything
#between the first occurence of the PDB ID contained within "chains" and the end of the file
# void (file_handle,Bio.SeqRecord)
def repair_recent_fasta_entries(fasta_file_handle, chains):
    
    offset = 0
    while True:
        fasta_file_handle.seek(offset, 2)
        if fasta_file_handle.read(1) == ">":
            pass

#using this function probably does not have a big advantage over using "get_pdb_path_list_recursively()"
#though tying everything to the index file instead of to the file system structure might be useful later
# (str,str) -> [str,str,...]                      
def get_pdb_path_list_by_index_file(index_file_full_path, loc_pdb_db_path):
    
    pdb_path = ""
    pdb_path_list = []

    try:
        index_file_handle = open(index_file_full_path, 'r')
    except:
        errMsg = "Could not open %s." % index_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        for pdb_id in index_file_handle:
            pdb_id_stripped = pdb_id.rstrip()
            pdb_path = os.path.join(loc_pdb_db_path, pdb_id_stripped[1:3], "pdb"+pdb_id_stripped+".ent.gz")
            pdb_path_list.append(pdb_path)
    except:
        errMsg = "Error while reading %s." % index_file_full_path
        close_file_safely(index_file_handle, index_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)

    close_file_safely(index_file_handle, index_file_full_path, "")
    
    return pdb_path_list

#if there is an exception during the writing of the FASTA file in create_fasta_file_from_index_file(),
#then this function is called to remove the most recently added entries belonging to the PDB that was read
#while the exception occured.
#starts with the file handle at the end of file, works its way up to the last occurence of most_recent_pdb_id
#and deletes everything from there to end of file
# void (str, str)
def rm_most_recent_pdb_from_fasta_file(fasta_file_full_path, most_recent_pdb_path):
    
    #get PDB ID from file path
    pdb_id_match = re.search(r'(.+\/)*pdb(....)\..*', most_recent_pdb_path)
    pdb_id = pdb_id_match.group(2).upper()

    try:
        fasta_file_handle = open(fasta_file_full_path, 'rb+')
    except:
        errMsg = "Failed to open %s." % fasta_file_full_path
        errorFct(errMsg)
        raise
    try:
        offset = 0
        first_occurence_offset = 0
        while True:
            fasta_file_handle.seek(offset, 2)
            char_byte = fasta_file_handle.read(1)
            #every char in UTF-8 that is only decoded by exactly one byte has a 0 as a first bit
            if int.from_bytes(char_byte, byteorder='little') & (1 << 7):
                errMsg = "FASTA file contains char that is not encoded by exactly 1 byte in UTF-8. Exiting."
                errorFct(errMsg)
                raise UnicodeDecodeError('utf-8',char_byte,offset,offset+1,'Char is encoded using more than 1 byte.')
            if char_byte == ">".encode('utf-8'):
                header = fasta_file_handle.read(4)
                if header == str(pdb_id).encode('utf-8'):
                    first_occurence_offset = offset
                    offset -= 1
                    continue
                else:
                    break
            else:
                offset -= 1
                continue
        #sets the file handle end of line before first occurence of matching pdb header
        fasta_file_handle.seek(first_occurence_offset-1, 2)
        #deletes everything from there to eof
        fasta_file_handle.truncate()
            
    except:
        errMsg = "Error during repairing of %s." % fasta_file_full_path
        close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg)
        errorFct(errMsg)
        raise
        
    close_file_safely(fasta_file_handle, fasta_file_full_path, "")

#reads a fasta file and returns all pdb ids (not chain names) in the headers
# (str) -> set() 
def read_fasta_pdb_ids(fasta_file_full_path):
    
    try:
        fasta_file_handle = open(fasta_file_full_path, 'r')
    except:
        errMsg = "Failed to open %s." % fasta_file_full_path
        errorFct(errMsg)
        raise

    pdb_id_set = set()

    for line in fasta_file_handle:
        if line.startswith(">"):
            pdb_id_set.add(line[1:5])
    
    return pdb_id_set     

#creates fasta file from index file
#void (str, str, str, str, int)
def create_fasta_file_from_index_file(loc_pdb_db_path, fasta_file_full_path, index_file_full_path, biopython_log_full_path, verbosity):
    
    pdb_path_list = get_pdb_path_list_by_index_file(index_file_full_path, loc_pdb_db_path)

    counter = 1
    number_of_pdbs = len(pdb_path_list)

    if verbosity>0:
        print("Writing Fasta file...")
    try:
        fasta_file_handle = open(fasta_file_full_path, 'w')
    except:
        errMsg = "File %s could not be opened for writing." % fasta_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        biopython_log_file_handle = open(biopython_log_full_path, 'a+')
    except:
        errMsg = "File %s could not be opened for writing." % biopython_log_full_path
        errorFct(errMsg)
        exit(1)
    try:
        for pdb_path in pdb_path_list:
            most_recent_pdb_path = pdb_path
            write_loc_pdb_to_fasta_file(pdb_path, fasta_file_handle, biopython_log_file_handle)
            if verbosity>0:
                msg = "[%s/%s]" % (str(counter), str(number_of_pdbs))
                print(msg, end="\r")
            counter += 1
    except:
        print(traceback.format_exc())
        errMsg1 = "Error while writing to %s." % fasta_file_full_path
        close_file_safely(biopython_log_file_handle, biopython_log_full_path, errMsg1)
        #if an exception is thrown during writing of fasta file this will flush fasta file to memory
        close_file_safely(fasta_file_handle, fasta_file_full_path, errMsg1)
        try:
            dest_backup = backup_fasta_file(fasta_file_full_path)
        except:
            errMsg = "Error while backing up FASTA file. Exiting."
            errorFct(errMsg)
            exit(1)
        #removes most recent pdb chains from fasta file
        try:    
            rm_most_recent_pdb_from_fasta_file(fasta_file_full_path, most_recent_pdb_path)
            fasta_pdb_id_set = read_fasta_pdb_ids(fasta_file_full_path)
            #overwrite index file
            write_index_file_from_set(index_file_full_path, fasta_pdb_id_set)
        except:
            try:
                restore_fasta_file_from_backup(dest_backup, fasta_file_full_path)
            except:
                print(traceback.format_exc())
                errMsg3 = "Error while restoring FASTA file backup %s." % dest_backup
                errorFct(errMsg3)
            finally:
                print(traceback.format_exc())
                errMsg2 = "Error while attempting to repair FASTA file %s." % fasta_file_full_path
                errorFct(errMsg2)
        else:
            remove_fasta_backup(dest_backup, fasta_file_full_path)
        finally:
            errorFct(errMsg1)
            exit(1)
          
    close_file_safely(biopython_log_file_handle, biopython_log_full_path, '')
    close_file_safely(fasta_file_handle, fasta_file_full_path, '')


#creates diamond database from fasta file
#void (str, str, str, int)        
def create_diamond_database(diamond_file_path, diamond_db_path, fasta_file_full_path, verbosity):
    
    diamond_db_full_path = os.path.join(diamond_db_path, "diamond_db")

    #<path_to_diamond_file> makedb --in <path_to_FASTA> --db <db_out_path>
    shell_input = []
    shell_input = [diamond_file_path, "makedb", "--in", fasta_file_full_path, "--db", diamond_db_full_path]
    
    if verbosity>0:
        print("Creating local diamond database...")
    if verbosity>1:
        process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            print(line, end='\n')
    else:
        process = subprocess.Popen(shell_input, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
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

#think about when to call this.
#when the rsync process fails, e.g. internet connection is lost,
#broken pdb files might remain that might not be cleaned up by calling rsync again, since
#they're already indexed as present in the directory tree.
#the assumption here is that the pdb files are being written from top to bottom.
#when the writing process fails halfway though the file will be truncated, but might still seem intact.
#thus, in truncated files the trailing "END" statement will be missing. check if its there.
# void (str)        
def rm_incomplete_pdbs(loc_pdb_db_path):
    
    pdb_path_list = get_pdb_path_list_recursively(loc_pdb_db_path)
        
#synchronizes copy of remote PDB archive with local copy
#void (str,str,str,int)  
def sync_pdb_copy(diamond_file_path, loc_pdb_db_path, diamond_db_path, verbosity):
    
    tmp_index_file_full_path = os.path.join(diamond_db_path, "tmp_pdb_db.index")
    old_index_file_full_path = os.path.join(diamond_db_path, "pdb_db.index")
    fasta_file_full_path = os.path.join(diamond_db_path, "pdb_db.fasta")
    biopython_log_full_path = os.path.join(diamond_db_path, "biopython_log.txt")

    create_dir_safely(loc_pdb_db_path)
    create_dir_safely(diamond_db_path)
    #https://www.wwpdb.org/ftp/pdb-ftp-sites
    #rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb
    #there's also DNA/RNA in there, look into it to maybe only sync protein structures
    shell_input = []
    shell_input = ["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444", "rsync.rcsb.org::ftp_data/structures/divided/pdb/", loc_pdb_db_path]
    
    if verbosity>0:
        print("Synchronizing local PDB database with remote. This may take a while.")
    
    process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    if verbosity>0:
        counter = 1
        for line in process.stdout:
            if counter>1000:
                process.kill()
                break
            line = line.rstrip()+"                   "
            print(line, end="\r")
            counter += 1

    try:
        out, err = process.communicate()
    except:
        process.kill()
        #rm_incomplete_pdbs(loc_pdb_db_path)
        out, err = process.communicate()
        errMsg = "Communication with rsync subprocess failed. Killed it."
        shell_err_fct(err)
        errorFct(errMsg)
        exit(1)
                    
    if verbosity>0:
        print("Updating PDB index file...")
    
    added_set = update_pdb_index(loc_pdb_db_path, diamond_db_path, tmp_index_file_full_path, old_index_file_full_path)
    if len(added_set)>0 and os.path.exists(fasta_file_full_path):
        add_fasta_db_entries(added_set, fasta_file_full_path, loc_pdb_db_path, biopython_log_full_path, verbosity)
        replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path)
    elif os.path.exists(fasta_file_full_path) is False:
        create_fasta_file_from_index_file(loc_pdb_db_path, fasta_file_full_path, old_index_file_full_path, biopython_log_full_path, verbosity)
        replace_old_index_file(tmp_index_file_full_path, old_index_file_full_path)
    else:
        if verbosity>0:
            print("No new PDB entries for update.")
       
    create_diamond_database(diamond_file_path, diamond_db_path, fasta_file_full_path, verbosity)

#removes gaps and other symbols not in valid_symbols and writes new fasta file with cleaned sequences
# void (str,str,[str,str,...], int)    
def clean_fasta_file(aln_file_full_path, out_file_full_path, valid_symbols, verbosity):
    
    if verbosity>0:
        msg = "Cleaning %s..." % os.path.split(aln_file_full_path)[1]
        print(msg)

    try:
        aln_file_handle = open(aln_file_full_path)
    except:
        errMsg = "Could not open alignment file %s." % aln_file_full_path
        errorFct(errMsg)
        exit(1)
    try:
        records = AlignIO.read(aln_file_handle, "fasta")
    except:
        errMsg = "Could not read fasta file %s." % aln_file_full_path
        close_file_safely(aln_file_handle, aln_file_full_path, errMsg)
        errorFct(errMsg)
        exit(1)

    close_file_safely(aln_file_handle, aln_file_full_path, "")
    
    counter = 1
    no_of_records = len(records)
    for record in records:
        if verbosity>0:
            msg = "[%s/%s]" % (str(counter), str(no_of_records))
            print(msg, end='\r')
        cleaned_seq = clean_seq(str(record.seq), valid_symbols)
        #if this fails use seqrecord constructor
        record.seq = Seq.Seq(cleaned_seq)
        counter += 1

    try:
        SeqIO.write(records, out_file_full_path, "fasta")
    except:
        errMsg = "Error while writing %s." % out_file_full_path
        errorFct(errMsg)
        exit(1)

#clean sequences of all fasta files in given dir. writes new ones to out_path with the same name
# void (str,str,[str,str,...],int)        
def clean_fasta_files_in_dir(path_to_fasta_files, out_path, valid_symbols, verbosity):
    
    for fasta_file_name in os.listdir(path_to_fasta_files):
        fasta_file_full_path = os.path.join(path_to_fasta_files, fasta_file_name)
        out_file_full_path = os.path.join(out_path, os.path.split(fasta_file_full_path)[1])
        clean_fasta_file(fasta_file_full_path, out_file_full_path, valid_symbols, verbosity)

# ./diamond blastp -d <diamond_db_file_full_path> -q <query_fasta_file_full_path> -o <out_path> --iterate
#launches a subprocess for a diamond blastp query without restriction of sequence identity
# void (str,str,str,str,int)        
def diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path, verbosity):

    shell_input = []
    shell_input = [diamond_exe_full_path, 'blastp', '-d', diamond_db_file_full_path, '-q', query_fasta_file_full_path, '-o', out_file_full_path,
                   '--iterate', '--log', '--header', 'simple', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                   'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen']
    if verbosity>0:
        msg = "Performing blastp query for file %s." % os.path.split(query_fasta_file_full_path)[1]
        print(msg)
    
    blastp_process = subprocess.Popen(shell_input, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    try:
        out, err = blastp_process.communicate()
    except:
        blastp_process.kill()
        out, err = blastp_process.communicate()
        shell_err_fct(err)
        errMsg = "Failed to communicate with diamond blastp subprocess. Killed it."
        errorFct(errMsg)
        exit(1)

#starts a diamond blastp query, but with an entire dir of fasta files, writes them with the same file name to out_path
#void (str,str,str,str,int)        
def batch_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_dir_path, out_path, verbosity):
    
    create_dir_safely(out_path)

    for query_fasta_file_name in os.listdir(query_fasta_file_dir_path):
        query_fasta_file_full_path = os.path.join(query_fasta_file_dir_path, query_fasta_file_name)
        if os.path.isfile(query_fasta_file_full_path) is False:
            continue
        out_file_full_path = os.path.join(out_path, os.path.splitext(query_fasta_file_name)[0]+".tsv")
        diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, query_fasta_file_full_path, out_file_full_path, verbosity)

#checks if TSV file generated by diamond is properly tab-separated
# (str) -> ([int,int,...])        
def check_TSV_integrity(TSV_file_full_path):

    broken_line_list = []
    with open(TSV_file_full_path, 'r') as TSV_file_handle:
        header_match_number = 0
        line_count = 1
        for header in TSV_file_handle:
            header_match = re.findall(r'([^\t\n]+)(\t|\n)', header)
            header_match_number = len(header_match)
            break
        for line in TSV_file_handle:
            line_match = re.findall(r'([^\t\n]+)(\t|\n)', line)
            line_match_number = len(line_match)
            if line_match_number != header_match_number:
                errMsg = "Warning: Inconsistent data in file %s, line %s:\n%s" % (TSV_file_full_path, str(line_count), line_match)
                errorFct(errMsg)
            broken_line_list.append(line_count)
            line_count += 1
            
    return broken_line_list

#this function defines what an exact match is
#currently it's 0 mismatches and alignment over full length of query
# (str) -> (pandas.DataFrame)
def filter_exact_matches(TSV_file_full_path):
    
    df = pandas.read_csv(TSV_file_full_path, sep='\t')
    new_df = pandas.DataFrame(columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qlen'])
    row_list = []
    #remove next line, once trailing delimiters have been removed
    #df = df.drop(['Unnamed: 16'], axis=1)
    df['is_exact_match'] = 0
    ex_match_col_num = df.columns.get_loc('is_exact_match')
    mismatch_col_num = df.columns.get_loc('mismatch')
    len_col_num = df.columns.get_loc('length')
    qlen_col_num = df.columns.get_loc('qlen')
    #split df into one dataframe per fasta_header
    df_list = [x for _,x in df.groupby('qseqid')]
    for split_df in df_list:
        #next lines needs to be changed, if definition of exact match requires it
        split_df = split_df.sort_values(by=['pident'], ascending=False)
        if int(split_df.iloc[0,mismatch_col_num]) == 0 and int(split_df.iloc[0,len_col_num]) == int(split_df.iloc[0,qlen_col_num]):
            split_df.iloc[0,ex_match_col_num] = 100
            row_list.append(split_df.iloc[[0]].copy())
    new_df = pandas.concat(row_list, ignore_index=True)

    return new_df

#definition of exact match needs to be applied here
# (str) -> { {},{},... },{ {},{},... }
def get_exact_matches(TSV_file_full_path):
    
    #optional check, if I trust diamond to generate good data
    #broken_line_list = check_TSV_integrity(TSV_file_full_path)
    #print(broken_line_list)
    #if len(broken_line_list) != 0:
     #   attempt_TSV_repair()
    #apply filter according to exact match definition to get at most one chain id for each query header
    exact_matches_df = filter_exact_matches(TSV_file_full_path)
    #for each row in dataframe write dict and append it to list
    exact_match_dict = {}
    #iterating over dataframe rows is slow, look for better options
    for index, row in exact_matches_df.iterrows():
        row_dict = {}
        for col_name in exact_matches_df.columns:
            row_dict[col_name] = row[col_name]
        exact_match_dict[row["sseqid"]] = row_dict.copy()
        
    return exact_match_dict
  
#copies PDB files from local PDB database to <job_dir>/DATA/PDB_lib
#mainly implemented for consistency, no real advantage
# void ([str,str,...],str,str)
def cp_PDB_files_to_job_dir(chain_list, pdb_db_path, pdb_out_path):
    
    for chain in chain_list:
        try:
            pdb_id = chain[0:4].lower()
        except IndexError:
            errMsg = "Chain name has less than 4 chars."
            errorFct(errMsg)
            raise
        pdb_file_name = "pdb"+str(pdb_id)+".ent.gz"
        src_pdb_full_path = os.path.join(pdb_db_path, pdb_id[1:3], pdb_file_name)
        dest_pdb_full_path = os.path.join(pdb_out_path, pdb_file_name)
        try:
            shutil.copy(src_pdb_full_path, dest_pdb_full_path)
        except:
            errMsg = "Error while copying file %s to %s." % (src_pdb_full_path, dest_pdb_full_path)
            errorFct(errMsg)
            raise

#assigns numbers (temporary ids) to headers and sequences in alignment file, return dict
#this function is very similar to import_seq_list_from_fasta_aln(), but i dont want to change the latter at the moment
# (str) -> { (str,int),(str,int),... }        
def enumerate_aln_headers(aln_file_full_path):
    
    aln_dict = {}
    
    with open(aln_file_full_path, 'r') as aln_file_handle:
        tmp_id = 0
        try:
            alignment = AlignIO.read(aln_file_handle, "fasta")
        except:
            errMsg = "Failed to read alignment file %s." % aln_file_handle
            errorFct(errMsg)
            raise
        for record in alignment:
            aln_dict[record.id] = (str(record.seq), tmp_id)
            tmp_id += 1
            
    return aln_dict
                
def diamond_calc_job_SPS(job_list_of_dict_lists, exact_match_dict, aln_dict, job_title, verbosity):
    
    if verbosity>0:
        msg = "Calculating SPS for job %s." % job_title
        print(msg)

    ref_DALI_set = set()
    MSA_set = set()
    
    #each list in job_list_of_dict_lists represents a sub job, i.e. a dir of format xxxxX_xxxxX
    #for each sub job in job
    for sub_job_dict_list in job_list_of_dict_lists:
        #for each alignment in sub job
        for dali_aln in sub_job_dict_list:
            DALI_aln_index_set = set()
            query_range = set()
            sbjct_range = set()
            DALI_aln_index_set.clear()
            query_range.clear()
            sbjct_range.clear()
            try:
                query_chain_id = dali_aln["query"][0:4]+":"+dali_aln["query"][4]
                sbjct_chain_id = dali_aln["sbjct"][0:4]+":"+dali_aln["sbjct"][4]
            except IndexError:
                errMsg = "Query or subject name in DALI alignment %s contain less than 5 characters."
                errorFct(errMsg)
                exit(1)
            DALI_query_seq = dali_aln["DALI_query_seq"]
            DALI_sbjct_seq = dali_aln["DALI_sbjct_seq"]
            orig_query_header = exact_match_dict[query_chain_id]["qseqid"]
            orig_sbjct_header = exact_match_dict[sbjct_chain_id]["qseqid"]
            tmp_query_id = aln_dict[orig_query_header][1]
            tmp_sbjct_id = aln_dict[orig_sbjct_header][1]
            orig_query_seq = aln_dict[orig_query_header][0]
            orig_sbjct_seq = aln_dict[orig_sbjct_header][0]
            #calculate interval intersection
            query_sstart = exact_match_dict[query_chain_id]["sstart"]
            query_send = exact_match_dict[query_chain_id]["send"]
            query_range = set(range(query_sstart,query_send+1))
            sbjct_sstart = exact_match_dict[sbjct_chain_id]["sstart"]
            sbjct_send = exact_match_dict[sbjct_chain_id]["send"]
            sbjct_range = set(range(sbjct_sstart,sbjct_send+1))
            overlap_interval = query_range.intersection(sbjct_range)
            DALI_tmp_aln_index_set = calc_ref_SP_set((tmp_query_id, tmp_sbjct_id), (DALI_query_seq, DALI_sbjct_seq))
            #add only those elements to ref_DALI_set that overlap in the original chain, since non-overlapping regions
            #in the original MSA cannot be measured
            for element in DALI_tmp_aln_index_set:
                if element[0][1] in overlap_interval and element[1][1] in overlap_interval:
                    DALI_aln_index_set.add(element)
            ref_DALI_set = ref_DALI_set.union(DALI_aln_index_set.copy())
            orig_aln_index_set = calc_MSA_SP_set((tmp_query_id, tmp_sbjct_id), (orig_query_seq, orig_sbjct_seq))
            MSA_set = MSA_set.union(orig_aln_index_set.copy())
    print("DALI_set: ")
    print(ref_DALI_set)
    print("MSA_set: ")
    print(MSA_set)
    intersection_set = ref_DALI_set.intersection(MSA_set)
    if len(ref_DALI_set) == 0:
        job_SPS = 0
    else:
        job_SPS = len(intersection_set)/len(ref_DALI_set)
    
    return job_SPS

#takes alignment and calculates SPS score from pairwise structural alignments
#        
def calc_aln_score_with_diamond(aln_file_full_path, diamond_exe_full_path, diamond_db_file_full_path, pdb_db_path, out_path,
                                dali_path, sample_size, job_title, valid_symbols, verbosity):
    
    job_out_path = os.path.join(out_path, job_title)
    create_dir_safely(job_out_path)
    
    job_data_out_path = os.path.join(job_out_path, "DATA")
    create_dir_safely(job_data_out_path)
    
    aln_file_name = os.path.splitext(os.path.split(aln_file_full_path)[1])[0]
    TSV_out_path = os.path.join(job_data_out_path, "TSV")
    TSV_file_full_path = os.path.join(TSV_out_path, aln_file_name+".tsv")
    create_dir_safely(TSV_out_path)

    clean_fasta_out_path = os.path.join(job_data_out_path, "CLEAN_ALN")
    create_dir_safely(clean_fasta_out_path)
    clean_fasta_full_path = os.path.join(clean_fasta_out_path, aln_file_name)
    clean_fasta_file(aln_file_full_path, clean_fasta_full_path, valid_symbols, verbosity)

    #it's a bit of a waste to read the alignment file yet again and not import this information while the
    #file is read the first time, but it's definitely clearer this way
    #aln_dict has original fasta headers as keys, check for duplicate keys
    aln_dict = enumerate_aln_headers(aln_file_full_path)

    diamond_blastp_query(diamond_exe_full_path, diamond_db_file_full_path, clean_fasta_full_path, TSV_file_full_path, verbosity)
    #exact_match_dict has chains as keys, i.e. XXXX:X
    #be extra careful and watch out for leading spaces and wrong or missing \t separations
    #check for duplicate keys
    exact_match_dict = get_exact_matches(TSV_file_full_path)
    
    if sample_size != 0:
        #randomize key popping and write function
        key_list = list(exact_match_dict.keys())
        key_list_len = len(key_list)
        no_to_remove = key_list_len - sample_size
        if no_to_remove < 0:
            no_to_remove = 0
        i = 1
        while i <= no_to_remove:
            #generate random number between 0 and len(key_list)
            random_int = random.randint(0, len(key_list)-1)
            key = key_list[random_int]
            try:
                exact_match_dict.pop(key)
            except KeyError:
                continue
            i += 1
    
    pdb_out_path = os.path.join(job_data_out_path, "PDB_lib")
    create_dir_safely(pdb_out_path)

    cp_PDB_files_to_job_dir(exact_match_dict.keys(), pdb_db_path, pdb_out_path)

    dat_out_path = os.path.join(job_data_out_path, "DAT_lib")
    create_dir_safely(dat_out_path)

    DALI_import_PDBs(dali_path, pdb_out_path, dat_out_path, verbosity)
    
    #use xxxxX format for all dali related things, instead of xxxx:X
    DALI_all_vs_all_query(dali_path, pdb_out_path, out_path, dat_out_path, job_title, verbosity, chain_list=exact_match_dict.keys())

    job_list_of_dict_lists = import_aln_files_in_job(job_title, out_path)

    job_SPS = diamond_calc_job_SPS(job_list_of_dict_lists, exact_match_dict, aln_dict, job_title, verbosity)
    
    print(job_SPS)
        

#calculation rework
#__________________________________________________________________________________________________________________________________________________
def extract_exact_query_matches(aln_file_path_list, out_path, valid_symbols):
    for aln_file_path in aln_file_path_list:
        path, aln_file_name_ext = os.path.split(aln_file_path)
        aln_file_name, ext = os.path.splitext(aln_file_name_ext)
        #aln_list is a list of lists each of length 3, containing [sequence, header, entry_number]
        aln_list = import_seq_list_from_fasta_aln(aln_file_path)
        #VERY CAREFUL HERE, DOUBLE CHECK, TRIPLE CHECK, IF IT MAKES SENSE TO DO IT LIKE THIS
        #IF WE FORGOT A SINGLE LETTER, THIS RUINS ALL RESULTS
        #HOW DO WE HANDLE UNCOMMON AAs OR 'X' AND THE LIKE?
        regex_str = build_regex_for_seq_cleaning_whitelist(valid_symbols)
        for aln in aln_list:
            cleaned_seq = clean_seq(aln[0], regex_str)

#creates job dir from a list of alignments
def create_job(aln_path_list, out_path, job_name):
    job_path = os.path.join(out_path, job_name)
    create_dir_safely(job_path)
    extract_exact_query_matches(aln_path_list, job_path)

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
    argParser.add_argument("-db", "--locpdbdb", default="PDB_db", help="Path to local PDB database.", required=False)
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
    argParser.add_argument("-clf", "--cleanfasta", action="count", default=0, help="testing fasta cleaning", required=False)
    argParser.add_argument("-babp", "--batchblastp", action="count", default=0, help="testing fasta cleaning", required=False)
    argParser.add_argument("-ddb", "--diamonddbfile", help="Path to diamond database file.", required=False)
    argParser.add_argument("-svg", "--svgpath", help="Path to SVG files to be concatenated.", required=False)
    argParser.add_argument("-cr", "--calcrework", default=0, action="count", help="Placeholder for new calculation of job", required=False)
    args = argParser.parse_args()
    
    valid_symbols = sel_alphabet(args.alphabet)
    invalid_symbols = set_non_alphabet(args.nonalphabet)

    if args.cleanfasta>0:
        clean_fasta_files_in_dir(args.batchsearch, args.output, valid_symbols, args.verbose)
    
    if args.batchblastp>0:
        batch_blastp_query(args.diamondfile, args.diamonddbfile, args.batchsearch, args.output, args.verbose)

    if args.syncdb>0:
        sync_pdb_copy(args.diamondfile, args.locpdbdb, args.output, args.verbose)

    if args.doitall > 0:
        #calc_SPS_from_aln_sample(args.output, args.alignmentfile, args.title, args.dalidir, args.samplesize, valid_symbols, args.verbose)
        calc_aln_score_with_diamond(args.alignmentfile, args.diamondfile, args.diamonddbfile, args.locpdbdb, args.output,
                                args.dalidir, args.samplesize, args.title, valid_symbols, args.verbose)

    if args.datamode > 0:
        if args.alignmentfile is not None:
            _ = write_aln_file_search_hits_to_csv(valid_symbols, args.alignmentfile, args.output, args.verbose, args.samplesize)
        if args.createsearchplots > 0:
            evaluate_diamond_search_data(args.csvdir, args.output, "identity", args.verbose)
            #with import svg_stack:
                #append_svg_images(args.svgpath, args.output)
        #wrap in function    
        elif args.batchsearch is not None:
            aln_file_list = get_file_paths_in_dir(args.batchsearch)
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
    if args.calcrework > 0:
        create_job()
    

    toc = time.perf_counter()
    print(f"Done in {toc - tic:0.4f} seconds.")
    exit(0)

if __name__ == "__main__":
    main()