#!/usr/bin/env python
# _*_ coding: utf-8 _*_


#### Editting Log (27th Dec 2022)####
## 1, Correcting 'uniq' checking for samples from the same individual
## 2, Correcting 'Endo' checking using carriers only
## 3, Adding one column contaning WGS position coverage
## 4, Unifying column classes for latter filtering
## 5, Adding row and column names for output tables

import LT_basic as LT
import sys
import os
import itertools
from collections import defaultdict

#################
##### Usage #####
#################
# python  summary_inherited_ERV.py meta.tab constitute_IS.list samplesheet known_IS output_path IS_depth bams_list WGS_cov
#   
#   meta.tab: file including all individual: pos IS SS direction ID prime 
#   constitute_IS.list: (pos prime) This file come from filtering meta table with more than 50 IS and more than 50 SS and recovering one half fail in meta table
#   samplesheet: one column with sample name only
#   known_IS: known IS from locaTER in bed format
#   output_path: folder to write the output tables
#   IS_depth: IS bam coverage against all PCIP Bams
#   bams_list: PCIP Bam list
#   WGS_cov: IS bam coverage against one WGS
#   Reps: this is a file containing replicates PCIP from the same subject.


#####################
## Functions
#####################

def write_list(clipped_list, prefix):
    fo = LT.fun_open_file(prefix, 'w')
    for i in clipped_list:
        for j in range(len(i)-1):
            fo.write("{}\t".format(i[j]))
        fo.write("{}\n".format(i[-1]))
    fo.close()

def print_help():
    if len(sys.argv) < 9:
        LT.fun_print_help("meta.tab", "Endo.list", "samplesheet", "known_IS", "path_out", "IS_depth", "bams_list", "WGS_cov", "Reps")
        sys.exit(1)

def load_Endo_li(Endo):
    IS_li = LT.fun_read(Endo)
    l = []
    for i in IS_li:
        l.append(i[0])
    return l

def load_meta_uniq_series(meta):
    with open(meta,"r") as f_in:
        
        uniq_carrier = defaultdict(list)
        uniq_direction = defaultdict(list)
        uniq_ss_counts = defaultdict(list)
        for line in f_in:
            line = line.rstrip().split()
            try:
                uniq_carrier[line[0]].append(str(line[4]))
            except KeyError:
                uniq_carrier[line[0]] = [str(line[4])]
            try:
                uniq_direction[line[0]].append(str(line[3]))
            except KeyError:
                uniq_direction[line[0]] = [str(line[3])]
            try:
                uniq_ss_counts[line[0]].append(int(line[2]))
            except KeyError:
                uniq_ss_counts[line[0]] = [int(line[2])]
    return  uniq_carrier, uniq_direction, uniq_ss_counts

def load_replicates(replicates):
    '''
        this is a file containing replicates PCIP from the same subject. Format as follow:
        
        A-rep-1 A
        A-rep-2 A
        A-rep-3 A
        B-rep-1 B
        B-rep-2 B

    '''
    with open(replicates, "r") as f_in:
        Reps = {}
        for line in f_in:
            fields = line.strip().split("\t")
            try:
                Reps[fields[0]] = fields[1]
            except KeyError:
                Reps[fields[0]] = fields[1]
    return Reps

def one_individual(carrier_li, info_carrier_dic):
    names = [info_carrier_dic[i] if i in info_carrier_dic.keys() else "single" for i in carrier_li]
    if ("single" not in names)  and names.count(names[0]) == len(names):
        return True
    else:
        return False

def is_IS_from_one_individual(li, Replicates_dic):
    '''
        when one IS was carried by multiple PCIP samples, it's import to know if the samples correspoding to one individual or not.
        chech a list of carriers, return True if they from one individual
        ["BE090909O", "BE090909Y"],["BE0909-sperm", "BE0909-muscle"]
    '''
    if len(li) == 1:
        return True
    elif one_individual(li, Replicates_dic):
        return True
    else:
        return False

def is_in_Endo_li(IS, Endo_li):
        if IS in Endo_li:
            return True
        else:
            return False
    
def make_Endo_search_di(Endo_li):
    search_di = {}
    for i in Endo_li:
        chr = i.split(":")[0]
        pos = int(i.split(":")[1])
        try:
            search_di[chr].append(pos)
        except KeyError:
            search_di[chr] = [pos]
    return search_di

def is_close_to_IS(IS, search_di):
    query_chrom = str(IS.split(":")[0])
    query_pos = int(IS.split(":")[1])
    # Be careful of that chrm don't have constitutive ERV
    try:
        query_li = search_di[query_chrom]
        distance_li = [query_pos - i if i < query_pos else i - query_pos for i in query_li]
        min_distance = min(distance_li)
        if min_distance <= 10000:
            closest_pos = query_li[distance_li.index(min_distance)]
            closest_IS = "{}:{}".format(query_chrom, closest_pos)
            return min_distance, closest_IS
        else:
            return "no", "no"
    except KeyError:
        return "no", "no"

def make_known_search_di(known_bed):
    with open(known_bed, "r") as f_in:
        search_di = {}
        for line in f_in:
            line = line.rstrip().split()
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            try:
                search_di[chr].append(start)
                search_di[chr].append(end)
            except KeyError:
                search_di[chr] = [start]
                search_di[chr] = [end]
        return search_di

def meta_dict(meta):
    """
        write the meta as a dictionary in order to searching based on key 'chr27:23296039_BE963859523_5LTR'
    """
    out_dict = {}
    for record in meta:
        key = "_".join([record[0], record[4], record[5]])
        value = int(record[2])
        out_dict[key] = value
    return out_dict

def count_IS(IS, samp_li, meta_di):
    out_li = []
    for samp in samp_li:
        key = "_".join([IS[0], samp, IS[1]])
        try:
            num = meta_di[key]
        except KeyError:
            num = 0
        out_li.append(num)
    return out_li

def IS_table(IS_li, samp_li, meta_di):
    out_li = [samp_li]
    for IS in IS_li:
        li = count_IS(IS, samp_li, meta_di)
        out_li.append(li)
    return out_li

def load_IS_depth(IS_depth):
    with open(IS_depth) as f_in:
        f_out = {}
        for line in f_in:
            line = line.strip().split()
            key = ":".join(line[0:2])
            value = line[2:]
            f_out[key] = value 
    return f_out

def make_IS_depth_bam_index(bam_list):
    with open(bam_list) as f_in:
        bam_li = []
        for line in f_in:
            line = line.rstrip().split()
            bam_li.append(line)
    bam_3LTR = []
    bam_5LTR = []
    for i in bam_li:
        if i[0].find('3LTR') != -1:
            bam_3LTR.append(i)
        else:
            bam_5LTR.append(i) 
    bam_3LTR_index = [bam_li.index(i) for i in bam_3LTR]
    bam_5LTR_index = [bam_li.index(i) for i in bam_5LTR]
    return bam_3LTR_index, bam_5LTR_index 

def make_NumofLibs_NumofReads(IS_depth_di, IS, index):
    count_li = []
    for i in index:
        try:
            all_count_li = IS_depth_di[IS]
            count_li.append(int(all_count_li[i]))
        except KeyError:
            count_li = [0]
    NumofLibs = sum([1 if i > 0 else 0 for i in count_li])
    NumofReads = sum(count_li)
    return NumofLibs, NumofReads
    
def main():
    print_help()
    # read meta_data.tab
    meta = LT.fun_read(sys.argv[1])
    
    # read IS list need to organize
    # This list come from filtering meta table with more than 50 IS and more than 100 SS 
    # and recover one half fail or not exit in meta table
    IS_li = LT.fun_read(sys.argv[2])

    ##################################################
    ## summarize Endo ERVs molecular in one table
    ##################################################
    # read sample list
    fi = LT.fun_open_file(sys.argv[3])
    samp_li = []
    for line in fi:
        samp_li.append(line.rstrip('\n'))
    # convert meta to a dictionary for searching
    meta_di = meta_dict(meta)
    # write the output
    Endo_sum = IS_table(IS_li, samp_li, meta_di)
    if os.path.exists(sys.argv[5]):
        path = os.path.join(sys.argv[5], "Endo_summary.txt")
        row_names = [["chr:pos", "prime"]] + IS_li
        Endo = []
        for (a,b) in zip(row_names, Endo_sum):
            Endo.append((a+b))
        write_list(Endo, path)
    else:
        print("Output file dir does not exist!")
        sys.exit(1)

    ##################################################
    ## add filters for meta-fill table
    ##################################################
    uniq_carrier, uniq_direction, uniq_ss_counts = load_meta_uniq_series(sys.argv[1])
    Endo_li = load_Endo_li(sys.argv[2])
    Endo_search_di = make_Endo_search_di(Endo_li)
    known_search_di = make_known_search_di(sys.argv[4])
    IS_depth_di = load_IS_depth(sys.argv[6])
    IS_depth_index_3, IS_depth_index_5 = make_IS_depth_bam_index(sys.argv[7])
    IS_depth_WGS = load_IS_depth(sys.argv[8])
    Reps = load_replicates(sys.argv[9])

    meta_fill = []
    for entry in meta:
        IS, carrier, direction = entry[0], entry[4], entry[5]
        # annotate is it uniq or not?
        if is_IS_from_one_individual(uniq_carrier[IS], Reps):
            entry.append("{}".format(len(uniq_carrier[IS])))
            entry.append("one")
            entry.append("/".join(uniq_carrier[IS]))
        else:
            entry.append("{}".format(len(uniq_carrier[IS])))
            entry.append("multiple")
            entry.append("#")
        # annotate is it close to Endo or not?
        min_dis_Endo, close_IS_Endo = is_close_to_IS(IS, Endo_search_di)
        if min_dis_Endo != "no" and (carrier in uniq_carrier[close_IS_Endo]):
            entry.append(min_dis_Endo)
            entry.append(close_IS_Endo)
        else:
            entry.append("10000")
            entry.append("no")
        # annotate is it close to known or not?
        min_dis_known, close_IS_known = is_close_to_IS(IS, known_search_di)
        if min_dis_known != "no":
            entry.append(min_dis_known)
            entry.append(close_IS_known)
        else:
            entry.append("10000")
            entry.append("no")
        
        # annotate IS depth for PCIP Bams
        if direction == "3LTR":
            numoflibs, numofreads = make_NumofLibs_NumofReads(IS_depth_di, IS, IS_depth_index_3)
        else:
            numoflibs, numofreads = make_NumofLibs_NumofReads(IS_depth_di, IS, IS_depth_index_5)
        entry.append(numoflibs)
        entry.append(numofreads)     
        # annotate IS depth for WGS
        try:
            entry.append(str(IS_depth_WGS[IS][0]))
        except KeyError:
            entry.append("NA")
        meta_fill.append(entry)

    if os.path.exists(sys.argv[5]):
        path = os.path.join(sys.argv[5], "Meta_fill.txt")
        header = [["pos", "IS", "SS", "direction", "ID", "prime", "occurrence", "uniq", "carrier", "close-to-Endo", "which-Endo","close-to-Known", "which-known", "Numoflib", "Numofreads", "WGS-coverage" ]]
        write_list((header + meta_fill), path)
    else:
        print("Output file dir does not exist!")
        sys.exit(1)
    


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        LT.fun_print_error("user interrupted, abort!!!")
        sys.exit(0)

