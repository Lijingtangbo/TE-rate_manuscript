#!/usr/bin/env python
# _*_ coding: utf-8 _*_

'''
This is a script analysing RRATE-seq reads.
Loading SAM file and generating four lists of junctions reads info.

Usage:
    python split_read_two_junctions.py in.sam out.prefix prime
'''

###### Load required modules ######

import LT_basic as LT
import os
import sys
from sam import parse_CIGAR_iter, parse_CIGAR
from Bio import Align

###### Define functions ######

def print_help():
    if len(sys.argv) < 3:
        LT.fun_print_help("in.sam", "out.prefix", "prime")

def fun_read_sam_to_list(sam):
    '''
        Read SAM file in a list and remove reads mapped to unplaced contigs
    '''
    fi = LT.fun_open_file(sam, "r")
    LIST_SAM_IN=[] 
    for l in fi:
        line=l.strip().split("\t")
        QNAME=line[0]
        FLAG=line[1]
        RNAME=line[2]
        POS=line[3]
        MAPQ=line[4]
        CIGAR=line[5]
        mateRNAME=line[6]
        matePOS=line[7]
        ISIZE=line[8]
        SEQ=line[9]
        QUAL=line[10]
        OTHER="/".join(line[11:])
        if RNAME.startswith('chr'): # Remove unplaced contigs
            LIST_SAM_IN.append([QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, mateRNAME, matePOS, ISIZE, SEQ, QUAL, OTHER])
    fi.close()
    return LIST_SAM_IN

def left_soft_clipped_list(LIST_SAM_IN):
    '''
        Recongize leftside of reads are softclipped and write the clipped part with moving 4 bases.
    '''
    left_soft_clipped_list=[]
    for read in LIST_SAM_IN:
        if read[5] != "*" and int(read[4]) != 255: # Cigar and MAPQ is not available.
            cigar_list = parse_CIGAR(read[5])
            if cigar_list[0][1] == "S":
                n = int(cigar_list[0][0])
                m = int(read[3])-1
                read = read[:] + [(m, n, read[9][0:n]), (m+1, n+1, read[9][0:n+1]), (m+2, n+2, read[9][0:n+2]), (m+3, n+3, read[9][0:n+3]), (m+4, n+4, read[9][0:n+4])]
                if n > 9 and int(read[4]) >= 40: # clipped length and MAPQ filter
                    left_soft_clipped_list.append(read)
    return left_soft_clipped_list

def right_consume_ref_length(cigar_str):
    '''
        Caculate the bases being consumed by mutations to adjust the position of rightclipped
    '''
    total_length = 0
    for counts, tag in parse_CIGAR(cigar_str):
        if tag in ('M','N','D', '=', 'X', 'EQ'):
            total_length += counts
    return total_length

def right_soft_clipped_list(LIST_SAM_IN):
    '''
        Recongize rightside of reads are softclipped and write the clipped part with moving 4 bases.
    '''
    right_soft_clipped_list=[]
    for read in LIST_SAM_IN:
        if read[5] != "*" and int(read[4]) != 255: # Cigar and MAPQ is not available.
            cigar_list=parse_CIGAR(read[5])
            if cigar_list[-1][1]=="S":
                m = int(right_consume_ref_length(read[5]) + int(read[3]))
                n = int(cigar_list[-1][0])
                read = read[:] + [(m, n, read[9][-n:]), (m-1, n+1, read[9][-(n+1):]), (m-2, n+2, read[9][-(n+2):]), (m-3, n+3, read[9][-(n+3):]), (m-4, n+4, read[9][-(n+4):])]
                if n > 9 and int(read[4]) >= 40: # clipped length and MAPQ filter
                    right_soft_clipped_list.append(read)
    return right_soft_clipped_list

Five_prime_PROVIRAL = 'TTCCAGACGAAGGCTCTGTTAACTTAGATGTCTGGGAAAA'
Five_prime_LTR = 'TGCGGGGAGCCGGTGAGGCATTCCACTCGTGACAAAGGTCATGAGGAAGGAGGCTCGGCATACGCAAAGGCGGGATCGAGCCTCAGGAGTCCCCCCGGATATTCTCGAGCATTTTCCCCCAAAAAAACCAGAGTCTGCCTACTTTATTGCTTTGTGCTCTCACCTCTGACTTTACTGGGGGCTGTC'
Three_prime_LTR = 'CCGCGTAAACCAAGCTACTCAGCTTCTTTTCTCCACTGAAATTTCCTACTGAGCTATCCTCATTCTATTGTTCTCTATATCCCTAATTAGCATATAAATAGTCGCCGACGCCGTCTCCCCTTCGAATACCCTGGATCAGCCGGGGCTGGTCCTCGGCA'
Three_prime_PROVIRAL = 'GATATGTGTAACACCCCTGAAAGTAAATGACACAGATTTTGAATGGGAAAAGATTAAAAACCATATTTCAGGTATTTGGAACAGCTCTGACATTAGCTTAGACTTAGGGAAACTTCACAATCAAATAGCAACCCTGGA'

def get_sense_alignment_ref(prime='', LR='', clipped_length=0):
    '''
        Return expercted sense alignment target reference.
    '''
    if prime == '5LTR':
        if LR == 'left':
            sense_ref = Five_prime_PROVIRAL[-clipped_length:]
        elif LR == 'right':
            sense_ref = Five_prime_LTR[:clipped_length]
        else:
            print('undefined clip direction')
            sys.exit()
    elif prime == '3LTR':
        if LR == 'left':
            sense_ref = Three_prime_LTR[-clipped_length:]
        elif LR == 'right':
            sense_ref = Three_prime_PROVIRAL[:clipped_length]
        else:
            print('undefined clip direction')
            sys.exit()
    else:
        print('undefined lib prime')
        sys.exit()
    return sense_ref

def get_antisense_alignment_ref(prime='', LR='', clipped_length=0):
    '''
        Return expercted antisense alignment target reference.
    '''
    if prime == '5LTR':
        if LR == 'left':
            antisense_ref = LT.fun_invert_dnaseq(Five_prime_LTR)[-clipped_length:]
        elif LR == 'right':
            antisense_ref = LT.fun_invert_dnaseq(Five_prime_PROVIRAL)[:clipped_length]
        else:
            print('undefined clip direction')
            sys.exit()
    elif prime == '3LTR':
        if LR == 'left':
            antisense_ref = LT.fun_invert_dnaseq(Three_prime_PROVIRAL)[-clipped_length:]
        elif LR == 'right':
            antisense_ref = LT.fun_invert_dnaseq(Three_prime_LTR)[:clipped_length]
        else:
            print('undefined clip direction')
            sys.exit()
    else:
        print('undefined lib prime')
        sys.exit()
    return antisense_ref
            
def pairwise_alignment(query='', ref='', mode = 'global'):
    '''
        Setting pairwise alignment parameters.
    '''
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.mismatch_score = -2
    aligner.target_end_gap_score = -2 
    aligner.query_end_gap_score = -2
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    alignment = aligner.align(ref, query)
    return alignment.score

def alignment_list_keep_seq(clipped_list, LR=''):
    '''
        Proceed pairwise alignment for softclipped part of reads.
    '''
    prime = sys.argv[3]
    li_out = []
    for read in clipped_list:
        sense_score_attach = []
        antisense_score_attach = []
        pos_attach = []
        for pos, length, seq in read[-5:]:
            sense_ref = get_sense_alignment_ref(prime, LR, length)
            antisense_ref = get_antisense_alignment_ref(prime, LR, length)
            sense_alignment_score = pairwise_alignment(seq, sense_ref)
            antisense_alignment_score = pairwise_alignment(seq, antisense_ref)
            sense_score_attach.append(float(sense_alignment_score/length))
            antisense_score_attach.append(float(antisense_alignment_score/length))
            pos_attach.append(int(pos))
        sense_high = max(sense_score_attach)
        antisense_high = max(antisense_score_attach)
        sense_pos = pos_attach[sense_score_attach.index(sense_high)]
        antisense_pos = pos_attach[antisense_score_attach.index(antisense_high)]
        read_out = read[0:10] 
        read_out.append(read[12+sense_score_attach.index(sense_high)]) 
        read_out.append(read[12+antisense_score_attach.index(antisense_high)])
        read_out.append(sense_score_attach.index(sense_high))
        read_out.append(sense_high)
        read_out.append(antisense_score_attach.index(antisense_high))
        read_out.append(antisense_high)
        li_out.append(read_out)
    return li_out

def get_LTR_JC_only(alignment_list, LR=''):
    '''
        Select the reads span host and LTR junctions.
    '''
    prime = sys.argv[3]
    LTR_JC_only = []
    if (prime == '5LTR' and LR == 'left') or (prime == '3LTR' and LR == 'right'):
        for read in alignment_list:
            if read[-1] > read[-3]:
                read_out = read[:-6] + read[-2:]
                read_out.append(read[-5])
                LTR_JC_only.append(read_out)
    else:
        for read in alignment_list:
            if read[-3] > read[-1]:
                read_out = read[:-6] + read[-4:-2]
                read_out.append(read[-6])
                LTR_JC_only.append(read_out)
    return LTR_JC_only

def LTR_JC_filter_without_seq(LTR_JC_only):
    LTR_JC_filtered = []
    if LTR_JC_only == []:
        return LTR_JC_filtered
    else:
        for read in LTR_JC_only:
            if read[-2] >0: # alignment_score filter
                read = [read[2], read[-1][0], read[-2], read[-3], read[-1][1]]
                LTR_JC_filtered.append(read)
    return LTR_JC_filtered

def get_PROVIRAL_JC_only(alignment_list, LR=''):
    '''
        Select reads span host-procirus junctions.
    '''
    prime = sys.argv[3]
    PROVIRAL_JC_only = []
    if (prime == '5LTR' and LR == 'left') or (prime == '3LTR' and LR == 'right'):
        for read in alignment_list:
            if read[-3] > read[-1]:
                read_out = read[:-6] + read[-4:-2]
                read_out.append(read[-6])
                PROVIRAL_JC_only.append(read_out)
    else:
        for read in alignment_list:
            if read[-1] > read[-3]:
                read_out = read[:-6] + read[-2:]
                read_out.append(read[-5])
                PROVIRAL_JC_only.append(read_out)
    return PROVIRAL_JC_only

def PROVIRAL_JC_refine_without_seq_5LTR(PROVIRAL_jc_only, LR = ''):
    PROVIRAL_JC_filtered = []
    if PROVIRAL_jc_only == []:
        return PROVIRAL_JC_filtered
    else:
        for read in PROVIRAL_jc_only:
            if LR == 'left' and str(read[-1][2]).endswith('GGGAAAA'):
                read = [read[2], read[-1][0], read[-2], read[-3], read[-1][1]]
                PROVIRAL_JC_filtered.append(read)
            elif LR == 'right' and str(read[-1][2]).startswith('TTTTCCC'):
                read = [read[2], read[-1][0], read[-2], read[-3], read[-1][1]]
                PROVIRAL_JC_filtered.append(read)
    return PROVIRAL_JC_filtered
               
def PROVIRAL_JC_refine_without_seq_3LTR(PROVIRAL_jc_only, LR = ''):
    PROVIRAL_JC_filtered = []
    if PROVIRAL_jc_only == []:
        return PROVIRAL_JC_filtered
    else:
        for read in PROVIRAL_jc_only:
            if LR == 'left' and str(read[-1][2]).endswith('ACATATC'):
                read = [read[2], read[-1][0], read[-2], read[-3], read[-1][1]]
                PROVIRAL_JC_filtered.append(read)
            elif LR == 'right' and str(read[-1][2]).startswith('GATATGT'):
                read = [read[2], read[-1][0], read[-2], read[-3], read[-1][1]]
                PROVIRAL_JC_filtered.append(read)
    return PROVIRAL_JC_filtered

def write_list(clipped_list, prefix):
    fo = LT.fun_open_file(prefix, 'w')
    for i in clipped_list:
        for j in range(len(i)-1):
            fo.write("{}\t".format(i[j]))
        fo.write("{}\n".format(i[-1]))
    fo.close()

def main():
    print_help()
    # read sam as a list
    LIST_SAM_IN = fun_read_sam_to_list(sys.argv[1])
    # write soft clipped reads with length more than 9 and MAPQ not less than 40 in itself and move 4 bases
    left_split = left_soft_clipped_list(LIST_SAM_IN)
    right_split = right_soft_clipped_list(LIST_SAM_IN)
    # align the soft_clipped part with corresponding ref in sense and antisense and write score and pos
    left_alignment = alignment_list_keep_seq(left_split, 'left')
    right_alignment = alignment_list_keep_seq(right_split, 'right')
    # write IS and SS separately and filter out reads alignment score less than 0
    LTR_JC_left = get_LTR_JC_only(left_alignment, 'left')
    LTR_JC_right = get_LTR_JC_only(right_alignment, 'right')
    LTR_JC_filtered_left = LTR_JC_filter_without_seq(LTR_JC_left)
    LTR_JC_filtered_right = LTR_JC_filter_without_seq(LTR_JC_right) 
    PROVIRAL_JC_left = get_PROVIRAL_JC_only(left_alignment, 'left')
    PROVIRAL_JC_right = get_PROVIRAL_JC_only(right_alignment, 'right')
    
    if sys.argv[3] == "5LTR":
        PROVIRAL_JC_filtered_left = PROVIRAL_JC_refine_without_seq_5LTR(PROVIRAL_JC_left, 'left')
        PROVIRAL_JC_filtered_right = PROVIRAL_JC_refine_without_seq_5LTR(PROVIRAL_JC_right, 'right')
    elif sys.argv[3] == "3LTR":
        PROVIRAL_JC_filtered_left = PROVIRAL_JC_refine_without_seq_3LTR(PROVIRAL_JC_left, 'left')
        PROVIRAL_JC_filtered_right = PROVIRAL_JC_refine_without_seq_3LTR(PROVIRAL_JC_right, 'right')
    
    write_list(LTR_JC_filtered_left, sys.argv[2] + '_' + sys.argv[3] + '_left_ltr_jc.tab')
    write_list(LTR_JC_filtered_right, sys.argv[2] + '_' + sys.argv[3] + '_right_ltr_jc.tab')
    write_list(PROVIRAL_JC_filtered_left, sys.argv[2] + '_' + sys.argv[3] + '_left_proviral_jc.tab')
    write_list(PROVIRAL_JC_filtered_right, sys.argv[2] + '_' + sys.argv[3] + '_right_proviral_jc.tab')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        LT.fun_print_error("user interrupted, abort!!!")
        sys.exit(0)