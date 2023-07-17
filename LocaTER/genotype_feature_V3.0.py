import pysam
import global_values as global_values
from multiprocessing import Pool
import os
import LT_basic as LT
import pybedtools
import sys


def unwrap_self_collect_clip_disc_features(arg, **kwarg):
    return XGenotyper.collect_features_one_site(*arg, **kwarg)

class XGenotyper():
    def __init__(self, sf_ref, working_folder, n_jobs, sf_rmsk):
        self.sf_reference = sf_ref
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.rmsk_annotation = sf_rmsk

    def call_genotype(self, sf_bam_list, sf_candidate_list, extnd, sf_out):
        l_bams = []
        with open(sf_bam_list) as fin_bam_list: #Reading list of BAM files
            for line in fin_bam_list:
                fields = line.split()
                l_bams.append(fields[0])
        # loop through Bams
        for sf_bam in l_bams:
            bam_id = os.path.basename(sf_bam)
            sf_out_features = sf_out + "{0}.txt".format(bam_id)
            l_chrm_records = []
            with open(sf_candidate_list) as fin_list: #Reading list of candidate loci
                # loop through candidates loci
                for line in fin_list:
                    fields = line.strip().split("\t")
                    pos = fields[0]
                    SUBFAM = fields[1]
                    fields = pos.split(":")
                    chrm = fields[0]
                    pos_l = int(fields[1].split("-")[0])
                    pos_r = int(fields[1].split("-")[1])
                    tmp_rcd = ((chrm, pos_l, pos_r, extnd, SUBFAM), sf_bam, self.working_folder)
                    l_chrm_records.append(tmp_rcd)
            pool = Pool(self.n_jobs)
            pool.map(unwrap_self_collect_clip_disc_features, list(zip([self]*len(l_chrm_records), l_chrm_records)), 1)
            pool.close()
            pool.join()

            with open(sf_out_features, 'w') as fout_all:
                # write one output feature file for each Bam
                for rcd in l_chrm_records:
                    chrm = rcd[0][0]
                    ins_pos = rcd[0][1]
                    s_pos_info = "{0}_{1}_{2}".format(chrm, ins_pos, bam_id)
                    sf_gntp_features = self.working_folder + s_pos_info + global_values.GNTP_FEATURE_SUFFIX

                    if os.path.isfile(sf_gntp_features) == False:
                        continue
                    with open(sf_gntp_features) as fin_site:
                        for line in fin_site:
                            fout_all.write(line)
                    os.remove(sf_gntp_features)


    def collect_features_one_site(self, record):
        # ((chrm, pos_l, pos_r, extnd), sf_bam, self.working_folder) ## record format

        chrm = record[0][0]
        ins_pos = record[0][1]
        ins_pos_r = record[0][2] 
        ins_m = (ins_pos + ins_pos_r)/2
        TSD_CUTOFF = (ins_pos_r - ins_pos)/2 + 5
        SUBFAM = record[0][4]
        extnd = record[0][3]
        start_pos = ins_pos - extnd
        if start_pos <= 0:
            start_pos = 1
        end_pos = ins_pos + extnd
        sf_bam = record[1]
        working_folder = record[2]
        bam_id = os.path.basename(sf_bam)

        chrm_in_bam = chrm
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename = self.sf_reference)
        m_chrm_id = self._get_chrm_id_name(samfile)
        l_disc_pairs = {}
        r_disc_pairs = {}
        n_l_disc_pairs = 0 #num of discordant read pairs at the left of breakpoint
        n_r_disc_pairs = 0 #num of discordant read pairs at the right of breakpoint
        l_disc_s = ""
        r_disc_s = ""
        n_concd_pairs = 0 #num of concordant read pairs bridge the breakpoint
        n_l_s_clip = 0 #num of left clipped reads within the region
        n_r_s_clip = 0 #num of right clipped reads within the region
        n_l_h_clip = 0 #num of left clipped reads for checking allel frequency
        n_r_h_clip = 0 #num of right clipped reads for cheking allel frequency
        n_full_map = 0
        n_l_cover_cp = 0 # pileup of reads in this given position (left of left-breakpoint)
        n_r_cover_cp = 0 # pileup of reads in this given position (right of right-breakpoint)
        set_cordant_pair_read_name = set() # read names of codant pairs to avoid counting same pair twice
        m_clip_qname = set() # clipped reads name
        s_pos_info = "{0}_{1}_{2}".format(chrm, ins_pos, bam_id)
        sf_gntp_features = working_folder + s_pos_info + global_values.GNTP_FEATURE_SUFFIX
        f_gntp_features = open(sf_gntp_features, "w")
        l_check_concord = []
        
        
        n_l_cover_cp = samfile.count(chrm_in_bam, (ins_pos -51), (ins_pos - 50) )
        n_r_cover_cp = samfile.count(chrm_in_bam, (ins_pos_r +50), (ins_pos_r + 51) )
        
        for algnmt in samfile.fetch(chrm_in_bam, start_pos, end_pos):
            if algnmt.is_duplicate == True: ## duplicate
                continue
            if algnmt.is_unmapped == True: # unmapped
                continue
            

            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            b_reverse = False
            if algnmt.is_reverse == True:
                b_reverse = True

            l_cigar = algnmt.cigar
            if len(l_cigar) < 1: # wrong alignment
                continue
            
            query_name = algnmt.query_name
            query_mapq = algnmt.mapping_quality
            map_pos = algnmt.reference_start
            map_end = algnmt.reference_end
            mate_chrm = '*'
            mate_pos = 0
            
            if (algnmt.next_reference_id in m_chrm_id) \
                and (algnmt.mate_is_unmapped == False) \
                and (algnmt.next_reference_id >= 0):
                mate_chrm = algnmt.next_reference_name
                mate_pos = algnmt.next_reference_start
####### reads fully mapped through insertion site
            b_fully_mapped = False
            i_map_start = map_pos
            i_map_end = -1
            if l_cigar[0][0] != 4 and l_cigar[-1][0] != 4 and l_cigar[0][0] != 5 and l_cigar[-1][0] != 5 and query_mapq != 0:
                b_fully_mapped = True
                i_map_end = map_end

            # check fully mapped reads cover the breakpoint
            if b_fully_mapped == True:
                #Read MUST cover 10 bp before and after IS
                if  i_map_start < (ins_pos-3) and (ins_pos_r+3) < i_map_end:
                    n_full_map += 1
                    #print(i_map_start, i_map_end)
                    
######## left soft clipped reads
            if l_cigar[0][0] == 4: 
                if algnmt.is_supplementary or algnmt.is_secondary:
                    continue
                if abs(map_pos - ins_m) <=  TSD_CUTOFF:
                    n_l_s_clip += 1
                    m_clip_qname.add(query_name)
######## right soft clipped
            if l_cigar[-1][0] == 4: 
                ## calculate the exact clip position
                ## only consider cigar types that consume ref bases
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:
                        continue
                    else:
                        map_pos += lenth
                if algnmt.is_supplementary or algnmt.is_secondary:
                    continue
                if abs(map_pos - ins_m) <= TSD_CUTOFF:
                    n_r_s_clip += 1
                    m_clip_qname.add(query_name)                       
######## left hard clipped reads
            if l_cigar[0][0] == 5: 
                if abs(map_pos - ins_m) <= TSD_CUTOFF:
                    n_l_h_clip += 1
                    m_clip_qname.add(query_name)
####### right hard clipped reads
            if l_cigar[-1][0] == 5: 
                ## calculate the exact clip position
                ## only consider cigar types that consume ref bases
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:
                        continue
                    else:
                        map_pos += lenth
                if abs(map_pos - ins_m) <= TSD_CUTOFF:
                    n_r_h_clip += 1
                    m_clip_qname.add(query_name)
####### Disconcordant pairs
            if mate_chrm == "*": ## mate unmapped reads are not interested!
                continue
            
            if query_mapq != 0 and self.is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos, global_values.DISC_THRESHOLD) == True:
                if self.is_in_LTR(mate_chrm, mate_pos):
                    subfam = self.is_in_LTR(mate_chrm, mate_pos)

                    if (b_first and b_reverse == False and (map_end <=  ins_pos)) or (b_first == False and b_reverse == False and (map_end <= ins_pos)):
                        try:
                            l_disc_pairs[subfam] += 1
                        except KeyError:
                            l_disc_pairs[subfam] = 1
                   
                    elif (b_first == False and b_reverse == True and (map_pos >= ins_pos_r)) or (b_first and b_reverse == True and (map_pos >= ins_pos_r)):
                        try:
                            r_disc_pairs[subfam] += 1
                        except KeyError:
                            r_disc_pairs[subfam] = 1

            elif algnmt.is_proper_pair:
                l_check_concord.append((query_name, chrm_in_bam, map_pos, mate_chrm, mate_pos, algnmt.inferred_length))
####### 
        for rcd_tmp in l_check_concord:
            if (rcd_tmp[0] not in m_clip_qname) and self.is_concrdant(rcd_tmp[1], rcd_tmp[2], rcd_tmp[3], rcd_tmp[4], ins_pos, ins_pos_r, rcd_tmp[5]) == True:
                set_cordant_pair_read_name.add(rcd_tmp[0])
                #print(rcd_tmp[0], rcd_tmp[2], rcd_tmp[5])
        n_concd_pairs = len(set_cordant_pair_read_name)

        if SUBFAM in l_disc_pairs:  
            n_l_disc_pairs = l_disc_pairs[SUBFAM]
           
        if SUBFAM in r_disc_pairs:
            n_r_disc_pairs = r_disc_pairs[SUBFAM]
        
        for key,value in l_disc_pairs.items():
            item = key + ":" + str(value) + "/"
            l_disc_s += item
        for key,value in r_disc_pairs.items():
            item = key + ":" + str(value) + "/"
            r_disc_s += item
        print(chrm, ins_pos)

        sinfo="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t".format(chrm, ins_pos, ins_pos_r, n_l_disc_pairs, n_r_disc_pairs, n_concd_pairs, n_l_s_clip, n_r_s_clip)
        sinfo1="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(n_l_h_clip, n_r_h_clip, n_full_map, n_l_cover_cp, n_r_cover_cp, l_disc_s, r_disc_s)

        f_gntp_features.write(sinfo+sinfo1)

        f_gntp_features.close()
        samfile.close()

    def _get_chrm_id_name(self, samfile):
        m_chrm = {}
        references = samfile.references
        for schrm in references:
            chrm_id = samfile.get_tid(schrm)
            m_chrm[chrm_id] = schrm
        m_chrm[-1] = "*"
        return m_chrm

    def is_discordant(self, chrm, map_pos, mate_chrm, mate_pos, is_threshold):
        if chrm != mate_chrm:
            return True
        else:
            if abs(mate_pos - map_pos) > is_threshold:
                return True
        return False

    # here require two reads are at the two sides of the insertion
    def is_concrdant(self, chrm, map_pos, mate_chrm, mate_pos, ins_pos, ins_pos_r, readLen): 
        if chrm == mate_chrm and map_pos < mate_pos:
            i_start = map_pos
            i_end = map_pos + readLen
            if i_end <= ins_pos and mate_pos >= ins_pos_r:
                return True
        return False

    def is_in_LTR(self, mate_chrm, mate_pos):
        '''
        This is the function to test if the reads mapping to a RepeatMasker annotated repeat.
        In this case only LTR elements, if yes from whichsubfamily.
        '''
        mate_end = mate_pos + global_values.READ_LENGTH
        rmsk = pybedtools.BedTool(self.rmsk_annotation)
        query_s = "{0} {1} {2}".format(mate_chrm, mate_pos, mate_end)
        query = pybedtools.BedTool(query_s, from_string=True)
        mate_interval = rmsk.intersect(query)
        if mate_interval:
            subfam,fam = mate_interval[0].name, mate_interval[0].score
            os.remove(query.fn)
            os.remove(mate_interval.fn)
            return subfam
            
        else:
            os.remove(query.fn)
            os.remove(mate_interval.fn)
            return False

def print_help():
    if len(sys.argv) < 3:
        print("May day! May day! May day!")
        LT.fun_print_help("bam.list", "targets.list", "working_folder")

def main():
    print_help()
    sf_ref = global_values.sf_ref
    extnd = global_values.extnd
    n_jobs = global_values.n_jobs
    sf_rmsk = global_values.sf_rmsk
    sf_bam_list = sys.argv[1]
    sf_candidate_list = sys.argv[2]
    sf_out = working_folder = sys.argv[3]


    ERVtype = XGenotyper(sf_ref, working_folder, n_jobs, sf_rmsk)
    ERVtype.call_genotype(sf_bam_list, sf_candidate_list, extnd, sf_out)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        LT.fun_print_error("user interrupted, abort!!!")
        sys.exit(0)
