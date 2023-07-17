#!/usr/bin/env python
# _*_ coding: utf-8 _*_

###### load modules ######
import sys
import os
import os.path
import LT_basic as LT
import pprint as pp

###### define functions ######
def filter_alignment_score(li, threshold=0.6):
    """
        filter out entries with alignments score less than 0.6
    """
    out_li = []
    for entry in li:
        if float(entry[2]) >= threshold:
            out_li.append(entry)
    return out_li

def de_multientries(li):
    """
        write a dict with each uniq locus as one key
        { chr:pos:[[alignment score], [flag of moving], [length of split_reads]], chr:pos:[[],[],[]], chr:pos[[],[],[]] ...}
    """
    out_dict = {}
    for entry in li:
        key = ":".join(entry[0:2])
        try:
            out_dict[key][0].append(float(entry[2]))
            out_dict[key][1].append(int(entry[3]))
            out_dict[key][2].append(int(entry[4]))
        except KeyError:
            out_dict[key]=[[float(entry[2])], [int(entry[3])], [int(entry[4])]]
    return out_dict

def number_format(di):
    """
        caculate the number of split reads for each uniq locus
        { chr:pos: int, chr:pos: int, ... }
    """
    out_dict = {}
    for locus in di.keys():
        out_dict[locus] = len(di[locus][1])
    return out_dict

def ltr_jc_clustering(DICT, threshold = 50):
    """
        cluster the LTR_JC within 50bp left and right windows
        [[{}, {}, {}, ...], [], [], ...]
    """
    clusters = []
    cluster = []
    IS = 0
    chr = ''
    # for key, value in DICT.items(): #####
    for i in sorted(DICT.items()):
        key, value = i[0], i[1]
        if key.split(":")[0] == chr and IS - threshold < int(key.split(":")[1]) < IS + threshold:
            cluster.append({key:value})
        else:
            clusters.append(cluster)
            cluster = []
            cluster.append({key:value})
            chr = key.split(":")[0] 
            IS = int(key.split(":")[1])
    clusters.append(cluster)
    return clusters[1:]

def is_pos_by_moving(li):
    if len([i for i in li if i > 0 ]) == len(li):
        return True
    else:
        return False

def proviral_jc_refine(DICT, fraction=.9, minimalMinorSupport=5, direction=''):
    drop_di = {}
    keep_di = {}
    if direction == 'left':
        for i in sorted(DICT.items()):
            key, value = i[0], i[1]
            if is_pos_by_moving(value[1]):
                chr = key.split(":")[0] 
                IS = int(key.split(":")[1])
                Nminorpeak = len(value[1])
                # moving 1 base
                if value[1][0] == 1:
                    if (":".join([chr,str(IS-1)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport: 
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 2 bases
                elif value[1][0] == 2:
                    if (":".join([chr,str(IS-2)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-1)]) in DICT.keys()) and DICT[":".join([chr,str(IS-1)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:    
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 3 bases
                elif value[1][0] == 3:
                    if (":".join([chr,str(IS-3)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-3)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-2)]) in DICT.keys()) and DICT[":".join([chr,str(IS-2)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-1)]) in DICT.keys()) and DICT[":".join([chr,str(IS-1)])][1][0] == 2:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 4 bases
                elif value[1][0] == 4:
                    if (":".join([chr,str(IS-4)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-4)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-3)]) in DICT.keys()) and DICT[":".join([chr,str(IS-3)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-3)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-2)]) in DICT.keys()) and DICT[":".join([chr,str(IS-2)])][1][0] == 2:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS-1)]) in DICT.keys()) and DICT[":".join([chr,str(IS-1)])][1][0] == 3:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS-1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
            else:
                keep_di[key] = value 
    else:
        for i in sorted(DICT.items()):
            key, value = i[0], i[1]
            if is_pos_by_moving(value[1]):
                chr = key.split(":")[0] 
                IS = int(key.split(":")[1])
                Nminorpeak = len(value[1])
                # moving 1 base
                if value[1][0] == 1:
                    if (":".join([chr,str(IS+1)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport: 
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 2 bases
                elif value[1][0] == 2:
                    if (":".join([chr,str(IS+2)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+1)]) in DICT.keys()) and DICT[":".join([chr,str(IS+1)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:    
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 3 bases
                elif value[1][0] == 3:
                    if (":".join([chr,str(IS+3)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+3)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+2)]) in DICT.keys()) and DICT[":".join([chr,str(IS+2)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+1)]) in DICT.keys()) and DICT[":".join([chr,str(IS+1)])][1][0] == 2:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
                # moving 4 bases
                elif value[1][0] == 4:
                    if (":".join([chr,str(IS+4)]) in DICT.keys()):
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+4)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+3)]) in DICT.keys()) and DICT[":".join([chr,str(IS+3)])][1][0] == 1:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+3)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+2)]) in DICT.keys()) and DICT[":".join([chr,str(IS+2)])][1][0] == 2:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+2)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    elif (":".join([chr,str(IS+1)]) in DICT.keys()) and DICT[":".join([chr,str(IS+1)])][1][0] == 3:
                        Nmajorpeak = len(DICT[":".join([chr,str(IS+1)])][1])
                        if Nminorpeak < Nmajorpeak * fraction or Nminorpeak < minimalMinorSupport:
                            drop_di[key] = value
                        else:
                            keep_di[key] = value
                    else:
                        keep_di[key] = value
            else:
                keep_di[key] = value 
    return drop_di, keep_di

def number_format_proviral(di):
    """
        caculate the number of split reads and average alignment score for each uniq locus
        { chr:pos: [Nsplit, AverScore], chr:pos: [], ...}
    """
    out_dict = {}
    for locus in di.keys():
        Nsplit = len(di[locus][1])
        Nscore = (sum(di[locus][0]))/Nsplit
        out_dict[locus] = [Nsplit, Nscore]
    return out_dict

def match_proviral(proviral_di, chr, IS, direction, threshold = 10000):
    out_pos_li = []
    out_number_li = []
    out_score_li = []
    if direction == 'left':
        for key, value in proviral_di.items():
            if key.split(":")[0] == chr and IS - 10000 < int(key.split(":")[1]) < IS:
                out_pos_li.append(int(key.split(":")[1]))
                out_number_li.append(int(value[0]))
                out_score_li.append(float(value[1]))
    elif direction == 'right':
        for key, value in proviral_di.items():
            if key.split(":")[0] == chr and IS < int(key.split(":")[1]) < IS + 10000:
                out_pos_li.append(int(key.split(":")[1]))
                out_number_li.append(int(value[0]))
                out_score_li.append(float(value[1]))
    else:
        sys.exit("undefined proviral split read direction!!")
    return [out_pos_li, out_number_li, out_score_li]
            
def ltr_proviral_merge(ltr_li, proviral_di, prime, plus=True):
    out_li = []
    for entry in ltr_li:
        chr = list(entry[0])[0].split(":")[0]
        IS = int(list(entry[0])[0].split(":")[1])
        if prime == '5LTR':
            if plus is True:
                proviral_li = match_proviral(proviral_di, chr, IS, 'left')
                out_li.append([entry, proviral_li])
            else:
                proviral_li = match_proviral(proviral_di, chr, IS, 'right')
                out_li.append([entry, proviral_li])
        elif prime == '3LTR':
            if plus is True:
                proviral_li = match_proviral(proviral_di, chr, IS, 'right')
                out_li.append([entry, proviral_li])
            else:
                proviral_li = match_proviral(proviral_di, chr, IS, 'left')
                out_li.append([entry, proviral_li])
    return out_li

def max_list(li):
    i = li[0]
    for j in li[1:]:
        if j > i:
            i = j
    return li.index(i)

def choose_most_support_IS(li):
    l = []
    for i in li:
        for value in i.values():
            l.append(int(value))
    return li[max_list(l)]

def set_IS(li):
    if len(li[0]) > 1:
        IS = list(choose_most_support_IS(li[0]))[0]
    elif len(li[0]) == 1:
        IS = list(list(li[0])[0])[0]
    return IS

def set_IS_number(li):
    if len(li[0]) > 1:
        for value in choose_most_support_IS(li[0]).values():
            num = value
    elif len(li[0]) == 1:
        for value in li[0][0].values():
            num = value
    return num

def possible_IS(li):
    out = []
    for i in li[0]:
        for key, value in i.items():
            IS = key.split(":")[1]
            Num = value
            out.append(":".join([IS, str(Num)]))
    return "/".join(out)

def plus_minus_merge(plus, minus):
    """
        final_list = [['chr:pos', pos/str, proviral_split, 'plus'], [], ...]
    """
    plus_minus_list = []
    for i in plus:
        IS = set_IS(i)
        IS_num = set_IS_number(i)
        possible_is = possible_IS(i)
        proviral_num = len(i[1][0])
        unit = [IS, IS_num, possible_is, proviral_num, 'plus']
        plus_minus_list.append(unit)
    for j in minus:
        IS = set_IS(j)
        IS_num = set_IS_number(j)
        possible_is = possible_IS(j)
        proviral_num = len(j[1][0])
        unit = [IS, IS_num, possible_is, proviral_num, 'minus']
        plus_minus_list.append(unit)
    return plus_minus_list

def plus_minus_dict(plus, minus):
    plus_minus_dict = {}
    for i in plus:
        value = i[1]
        IS = set_IS(i)
        proviral_num = len(value[0])
        value.append(proviral_num)
        plus_minus_dict[IS] = value
    for j in minus:
        value = j[1]
        IS = set_IS(j)
        proviral_num = len(value[0])
        value.append(proviral_num)
        plus_minus_dict[IS] = value
    return plus_minus_dict

def IS_SS_li(IS_SS_dict):
    out_li = []
    for IS, SS in IS_SS_dict.items():
        if len(SS[0]) == 0:
            out_li.append([IS, 'na', 0, 0, 0, 0])
        else:
            for i in range(len(SS[0])):
                dis=abs(SS[0][i] - int(IS.split(":")[1]))
                out_li.append([IS, SS[0][i], dis, SS[1][i], SS[2][i], SS[3]])
    return out_li

def print_help():
    if len(sys.argv) < 5:
        LT.fun_print_help("sample_name", "out.prefix", "split_folder", "out_folder", "prime")

def write_list(li, prefix):
    fo = LT.fun_open_file(prefix, 'w')
    for i in li:
        for j in range(len(i)-1):
            fo.write("{}\t".format(i[j]))
        fo.write("{}\n".format(i[-1]))
    fo.close()

def main():
    prime = sys.argv[5]
    print_help()
    if prime == "5LTR":
        # read 5LTR lib    
        left_ltr_5 = ltr_jc_clustering(number_format(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_5LTR_left_ltr_jc.tab'))))))
        right_ltr_5 = ltr_jc_clustering(number_format(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_5LTR_right_ltr_jc.tab'))))))
        drop_right, right_proviral_5 = proviral_jc_refine(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_5LTR_right_proviral_jc.tab')))), direction='right')
        right_proviral_5= number_format_proviral(right_proviral_5)
        drop_right = number_format_proviral(drop_right)
        drop_left, left_proviral_5 = proviral_jc_refine(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_5LTR_left_proviral_jc.tab')))), direction='left')
        left_proviral_5 = number_format_proviral(left_proviral_5)
        drop_left = number_format_proviral(drop_left)
        # 5LTR plus insertion
        plus_5 = ltr_proviral_merge(right_ltr_5, left_proviral_5, '5LTR', True)
        # 5LTR minus insertion
        minus_5 = ltr_proviral_merge(left_ltr_5, right_proviral_5, '5LTR', False)
        # merge plus add minus
        plus_minus_5 = plus_minus_merge(plus_5, minus_5)
        # write output
        write_list(plus_minus_5, os.path.join(sys.argv[4], sys.argv[1] + "_" + sys.argv[2] + '_5LTR.tab'))
        # write IS_SS
        write_list(IS_SS_li(plus_minus_dict(plus_5, minus_5)), os.path.join(sys.argv[4], sys.argv[1] + "_" + "IS_SS" + '_5LTR.tab'))
        # write dropped proviral pos
        write_list(drop_left, os.path.join(sys.argv[4], sys.argv[1] + "_" + "drop" + '_left_5LTR.tab'))
        write_list(drop_right, os.path.join(sys.argv[4], sys.argv[1] + "_" + "drop" + '_right_5LTR.tab'))
    else:
        # read 3LTR lib
        left_ltr_3 = ltr_jc_clustering(number_format(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_3LTR_left_ltr_jc.tab'))))))
        right_ltr_3 = ltr_jc_clustering(number_format(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_3LTR_right_ltr_jc.tab'))))))
        right_proviral_3 = number_format_proviral(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_3LTR_right_proviral_jc.tab')))))
        left_proviral_3 = number_format_proviral(de_multientries(filter_alignment_score(LT.fun_read(os.path.join(sys.argv[3], sys.argv[1] + '_3LTR_left_proviral_jc.tab')))))
        # 3LTR plus insertion
        plus_3 = ltr_proviral_merge(left_ltr_3, right_proviral_3, '3LTR', True)
        # 3LTR minus insertion
        minus_3 = ltr_proviral_merge(right_ltr_3, left_proviral_3, '3LTR', False)
        # plus add minus
        plus_minus_3 = plus_minus_merge(plus_3, minus_3)
        # write output
        write_list(plus_minus_3, os.path.join(sys.argv[4], sys.argv[1] + "_" + sys.argv[2] + '_3LTR.tab'))
        # write IS_SS
        write_list(IS_SS_li(plus_minus_dict(plus_3, minus_3)), os.path.join(sys.argv[4], sys.argv[1] + "_" + "IS_SS" + '_3LTR.tab'))

######### run ########
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        LT.fun_print_error("user interrupted, abort!!!")
        sys.exit(0)
