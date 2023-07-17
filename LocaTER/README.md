# LocaTER

*LocaTER*(Localization of Transposable Retroviral elements) is designed to identify ERV insertions from pair-end next generation whole genome sequencing data. 

# Dependencies

- scala (version 2.12.0)
- jdk (version 1.8.0_66)
- htsjdk-2.0.1.jar
- python (version 3.7)
- R (version 3.5.1+)

# Input

- BAM files (Pair-end short reads were aligned to bovine reference genome ARS-UCD1.2 using [BAW-MEM](https://github.com/lh3/bwa))
- LTR repeats annotation (Download from [RepeatMasker](https://www.repeatmasker.org/) annotation)

# Run *LocaTER*

## Step 1. Call potential polymorphic ERVs from individual BAMs
```
JAVA_OPTS=-Xmx12g \
scala \
-cp htsjdk-2.0.1.jar:. findMobileElements \
${NAME}_${PREFIX} \
${BAMS} \
${PATH_TO_REPEATS} 
```

## Step 2. Insertions sites refinement by merging multiple calls
```
scala -J-Xmx8g \
-cp ./commons-io-2.4.jar \
ERVsPop_v4.scala \
> output.tab
```

## Step 3. Call features of refined insertions sites individually
```
python genotype_feature_V3.0.py \
${BAMS} \
${InsertionSites} \
${OUTPUT}
```

## Step 4. Call genotypes based on features
```
Rscripts from_features_to_genotypes_v2.R \
${Feature} \
${Genotypes}
```
