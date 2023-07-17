args = commandArgs(TRUE)
inFile = args[1]
outFile = args[2]

if ( is.na(inFile) ) {
  stop( "Usage: Rscript inFile outFile")
}
stopifnot( file.exists(inFile) )

############### read genotype feature ##################
feature <- read.table(inFile, sep = "\t", header = F)
ID <- basename(inFile)
ID <- unlist(strsplit(ID, ".bam.txt"))
#feature <- read.table("File/82900120012_s.bam.txt", sep = "\t", header = F)
colnames(feature) <- c("chrm", "ins_pos_l", "ins_pos_r","n_l_disc_pairs", "n_r_disc_pairs", "n_concd_pairs", "n_l_s_clip", "n_r_s_clip","n_l_h_clip", "n_r_h_clip", "n_full_map", "n_l_cover_cp", "n_r_cover_cp", "l_disc_s", "r_disc_s")
## remove the chrX
feature <- feature[!feature$chrm == "chrX", ]

## remove the outliers
feature$n_disc_pairs <- feature$n_l_disc_pairs + feature$n_r_disc_pairs
feature$n_clipped_reads <- feature$n_l_s_clip + feature$n_r_s_clip + feature$n_l_h_clip + feature$n_r_h_clip
feature$IS <- paste(feature$chrm, feature$ins_pos_l, sep = ":")

high_limit_outlier_l = mean(feature$n_l_cover_cp) + sd(feature$n_l_cover_cp)
high_limit_outlier_r = mean(feature$n_r_cover_cp) + sd(feature$n_r_cover_cp)

OUTLIER <- feature[(feature$n_l_cover_cp == 0 & feature$n_r_cover_cp == 0 | feature$n_l_cover_cp > high_limit_outlier_l | feature$n_r_cover_cp > high_limit_outlier_r), ]

feature_v <- feature[!feature$IS %in% OUTLIER$IS, ]
high_limit_outlier_l_2 = mean(feature_v$n_l_cover_cp) + 5*sd(feature_v$n_l_cover_cp)
high_limit_outlier_r_2 = mean(feature_v$n_r_cover_cp) + 5*sd(feature_v$n_r_cover_cp)
low_limit_outlier_l_2 = mean(feature_v$n_l_cover_cp) - 3*sd(feature_v$n_l_cover_cp)
low_limit_outlier_r_2 = mean(feature_v$n_r_cover_cp) - 3*sd(feature_v$n_r_cover_cp)


outlier <- feature_v[(feature_v$n_l_cover_cp < low_limit_outlier_l_2 & feature_v$n_r_cover_cp < low_limit_outlier_r_2 | feature_v$n_l_cover_cp > high_limit_outlier_l_2 | feature_v$n_r_cover_cp > high_limit_outlier_r_2), ]

## WT
WT <- feature[feature$n_disc_pairs == 0 & feature$n_clipped_reads == 0, ]

## HOMO
HOMO <- feature[feature$n_concd_pairs == 0 & feature$n_full_map == 0 & feature$n_l_disc_pairs > 0 & feature$n_r_disc_pairs > 0,]

## define the lower bound for each feature
l_n_concd_pairs = unname(round(quantile(WT$n_concd_pairs,probs = c(0.05))))/2
l_n_full_map = unname(round(quantile(WT$n_full_map,probs = c(0.05))))/2
l_n_disc_pairs =  unname(round(quantile(HOMO$n_disc_pairs,probs = c(0.05))))/2
l_n_clipped_reads = unname(round(quantile(HOMO$n_clipped_reads,probs = c(0.05))))/2

HET <- feature[feature$n_concd_pairs >= l_n_concd_pairs & feature$n_full_map >= l_n_full_map & feature$n_disc_pairs >= l_n_disc_pairs & feature$n_clipped_reads >= l_n_clipped_reads,]

remain <- feature[!feature$IS %in% WT$IS, ]
remain <- remain[!remain$IS %in% HOMO$IS, ]
remain <- remain[!remain$IS %in% HET$IS, ]
remain <- remain[!remain$IS %in% outlier$IS, ]
remain <- remain[!remain$IS %in% OUTLIER$IS, ]

WT_r <- remain[remain$n_concd_pairs >= l_n_concd_pairs & remain$n_full_map >= l_n_full_map & remain$n_disc_pairs < l_n_disc_pairs & remain$n_clipped_reads < l_n_clipped_reads,]
HOMO_r <- remain[remain$n_concd_pairs == 0 & remain$n_full_map == 0 & remain$n_disc_pairs >= l_n_disc_pairs & remain$n_clipped_reads >= l_n_clipped_reads,]
HET_r <- remain[(remain$n_concd_pairs >= l_n_concd_pairs | remain$n_full_map >= l_n_full_map) & (remain$n_disc_pairs >= l_n_disc_pairs | remain$n_clipped_reads >= l_n_clipped_reads),]

## assign genotype
feature$geno <- ifelse(feature$IS %in% WT$IS, "WT",".")
feature$geno <- ifelse(feature$IS %in% WT_r$IS, "wt", feature$geno)
feature$geno <- ifelse(feature$IS %in% HOMO$IS, "HOMO", feature$geno)
feature$geno <- ifelse(feature$IS %in% HOMO_r$IS, "homo", feature$geno)
feature$geno <- ifelse(feature$IS %in% HET$IS, "HET", feature$geno)
feature$geno <- ifelse(feature$IS %in% HET_r$IS, "het", feature$geno)
feature$geno <- ifelse(feature$IS %in% outlier$IS, "outlier", feature$geno)
feature$geno <- ifelse(feature$IS %in% OUTLIER$IS, "OUTLIER", feature$geno)
feature <- subset(feature, select = -c(IS, n_disc_pairs, n_clipped_reads) )
feature$ID <- rep(ID, nrow(feature))
## write the genotypes
write.table(feature, file = outFile, sep = "\t", row.names = F, col.names = F, quote = F)


