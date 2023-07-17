
#### Running parameters; will be parsed 
sf_ref = "/home/u/s176815/Genomes_annotation/Fasta/ARS-UCD1.2.fa"
sf_rmsk = "/home/u/s176815/Genomes_annotation/Repeats/BosTau9_LTR_LIC_unplacedContigNames.bed"
extnd = 500 #window to extend for each side
n_jobs = 2 #integar; num of the cores 

#### Global parameters; default or set by running parameters
BPRE_FEATURE_SUFFIX = ".bpre_feature"
BPRE_SUFFIX = ".bpre"
GNTP_FEATURE_SUFFIX = ".gntp_feature"
BWA_HALF_READ_MIN_SCORE = 75
TSD_CUTOFF = 12
CLIP_EXACT_CLIP_SLACK = 10
DISC_THRESHOLD = 1500
DFT_IS = 300
LARGE_INDEL_IN_READ = 10
READ_LENGTH = 100
CK_POLYA_CLIP_WIN=25 #only check whether contain polyA for those clipped reads in a 25bp window
CK_POLYA_SEQ_MAX=20 #at most check 20 bases region for polyA
