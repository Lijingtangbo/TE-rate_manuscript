Scripts in this folder analysis adapting PCIP-seq reads to capture both polymorphic and *de novo* ERV insertions.

# Dependecies

- bwa (version 0.7.10)
- samtools (version 1.9)
- python (version 3.7)
	- Biopython

# Workflow

## 1. Alignment

The ensuing sequence reads were demultiplexed, quality assessed using [fastQC](https://github.com/s-andrews/FastQC), adapter sequences trimmed using [Cutadapt](https://github.com/marcelm/cutadapt), and trimmed reads mapped to the bovine reference genome using [BWA-MEM](https://github.com/lh3/bwa) and converted to BAM format using [SAMtools](https://github.com/samtools/). 
```
bwa mem -M \
${REFERENCE} \
${FASTQ}/${NAME}*_R1_001.fastq.gz ${FASTQ}/${NAME}*_R2_001.fastq.gz 2>${WD}/log/${NAME}.bwamem.log \
| samtools view -bS  > ${WD}/Bams/${NAME}.bam

samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 1800M -o ${WD}/Bams/${NAME}_sorted.bam ${WD}/Bams/${NAME}.bam

samtools index ${WD}/Bams/${NAME}_sorted.bam
```

## 2. Call split reads

Using a custom-made python script, we first identified clipped reads using CIGAR information.  We selected clipped reads with mapping quality ≥ 40 and a minimum of 10 clipped bases.  We then mapped the clipped reads to the segments of the ERVK[2-1_LTR] genome corresponding to the ERV-Tags in Fig. 1 (1 and 2 for the 5’LTR libraries and 3 and 4 for the 3’LTR libraries).  We demanded an alignment score ≥ 0.6 to declare a hit, and label the read as either an insertion site (IS) or a shearing site (SS) read.  When possible, we extended the alignment with ERVK[2-1-LTR] into non-clipped bases to refine the positions of the SS and IS. 
```
python split_read_two_junctions.py \
${SAM} \
${PRIME} # "5LTR" or "3LTR"
```

## 3. Cluster split reads

We then merged SS and IS site with same “breakpoint”, thereby identifying candidate SS and IS supported by multiple concordant reads.  We then paired IS with their cognate SS.  The pairing was based on orientation (f.i. a 5’ SS should be located upstream of the 5’ IS for an ERV element in “sense” orientation, and downstream of the 5’ IS for an ERV in “antisense” orientation; Fig. 1) and distance (the maximum distance between IS and SS was set at 5 Kb). 
```
python clustering_prediction_2.py \
${NAME} \
${SPLIT} \
${PRIME} # "5LTR" or "3LTR"
```

## 4. Summary

This module writes out an UMI table for polymorphic ERVs cross samples and a table for *de novo* ERVs filteration. 
```
python  summary_inherited_ERV.py \
meta.tab \
constitute_IS.list \
samplesheet \
known_IS \
output_path \
IS_depth \
bams_list\
WGS_cov
```

