# PAH-noncoding-variant-modifiers-analysis
Analysis of regulatory variation in pulmonary arterial hypertension and lung-development pathway genes.

Code associated with the manuscript “Variable expressivity of TBX4 loss-of-function: cataloging noncoding variants as potential modifiers of lung development.”

The pipeline evaluates variants within topologically associated domains
(TADs) encompassing genes involved in lung development and pulmonary
vascular biology and incorporates functional genomic annotations,
population allele frequencies, variant quality metrics, and regulatory
features.

# TBX4 pathway TAD variant analysis: commands used

Commands used to extract, annotate and compare variants in the TADs of TBX4 and 12 lung-development pathway genes in the TBX4-PAH and TBX4-lethal lung cohorts, with 1000 Genomes EUR as a reference.

# Genome assembly: GRCh38

### 01. Merge per-sample VCFs
```bash
bcftools merge -m none --file-list sample_IDs.txt -Oz -o <cohort>.vcf.gz
tabix <cohort>.vcf.gz
```

### 02. Regulatory tracks and TAD regions
```bash
# H3K4me3 NHLF (promoters) and ATAC-seq IMR-90 (enhancers), ENCODE
wget https://www.encodeproject.org/files/ENCFF343DPH/@@download/ENCFF343DPH.bed.gz
wget https://www.encodeproject.org/files/ENCFF383XOM/@@download/ENCFF383XOM.bed.gz
zcat ENCFF343DPH.bed.gz | cut -f1-4 | sort -k1,1 -k2,2n | bgzip > H3K4Me3_NHLF.bed.gz
zcat ENCFF383XOM.bed.gz | cut -f1-4 | sort -k1,1 -k2,2n | bgzip > ATACseq_IMR90.bed.gz

sort -k1,1V -k2,2n pathway_genes_tad.bed | bgzip > pathway_genes_tad.bed.gz
tabix pathway_genes_tad.bed.gz
bedtools intersect -wb -a pathway_genes_tad.bed.gz -b H3K4Me3_NHLF.bed.gz > H3K4Me3_NHLF_tad.bed
bedtools intersect -wb -a pathway_genes_tad.bed.gz -b ATACseq_IMR90.bed.gz > ATACseq_IMR90_tad.bed
```

### 03. Extract TAD variants and normalize
```bash
bcftools view --samples-file sample_IDs.txt --regions-file pathway_genes_tad.bed.gz -Oz -o <cohort>_filtered.vcf.gz <cohort>.vcf.gz
bcftools +fill-tags <cohort>_filtered.vcf.gz -- -t AC,AN,AF | bcftools view --min-ac=1 -Oz -o <cohort>_filtered.vcf.gz

vt decompose -s <cohort>_filtered.vcf.gz \
  | vt normalize -r GRCh38.fa.gz - \
  | awk '$5!="*" && $5!="." && $10!="./."' \
  | bcftools annotate --set-id '%CHROM:%POS\_%REF>%FIRST_ALT' \
  | bcftools +fill-tags -- -t AC,AN,AF \
  | bgzip -c > <cohort>_filtered.vt.vcf.gz
tabix <cohort>_filtered.vt.vcf.gz
```

### 04. 1000 Genomes EUR reference
16 unrelated EUR individuals were randomly selected from the 1000 Genomes high-coverage call set (relatedness and QC flags from the gnomAD HGDP-1KG metadata).
```bash
shuf -n 16 TGP_EUR_unrelated.txt > selected_samples.txt
bcftools view -S selected_samples.txt -Oz -o chr.subset.vcf.gz <1KG_chr>.vcf.gz
bcftools concat -Oz *.subset.vcf.gz | bcftools sort -Oz \
  | bcftools +fill-tags -- -t AC,AN,AF | bcftools view --min-ac=1 -Oz -o selectedEUR_TGP.vcf.gz
```

### 05. VEP annotation
Ensembl VEP release 115.
```bash
vep --cache --offline --assembly GRCh38 --fasta GRCh38.fa.gz \
  --input_file <cohort>_filtered.vt.vcf.gz --output_file <cohort>_filtered.vt.tsv --tab \
  --canonical --xref_refseq --variant_class --symbol --numbers --mane --protein \
  --sift b --polyphen b --humdiv --biotype --show_ref_allele --total_length --hgvs --check_existing \
  --flag_pick_allele_gene --pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length \
  --af_gnomad --max_af \
  --custom H3K4Me3_NHLF.bed.gz,H3K4Me3_NHLF,bed \
  --custom ATACseq_IMR90.bed.gz,ATACseq_IMR90,bed \
  --plugin CADD,snv=whole_genome_SNVs.tsv.gz,indels=gnomad.genomes.r4.0.indel.tsv.gz \
  --dont_skip --individual all --exclude_predicted
```


### 06. Collapse variants and assign tiers
`05_collapse_duplicate_variants.R`: one row per variant, pathway-gene TAD assignment, and tiers (Tier 1: MAX_AF ≤ 0.1%, in a promoter or enhancer peak, CADD ≥ 15; Tier 2: same with CADD 10–15; Tier 3: all others).

### 07. Regulatory enrichment
`enrichment.R`: variant density in promoter/enhancer peaks vs TAD background (excluding peaks and GENCODE v46 gene bodies); one-sided Fisher's exact and Poisson rate-ratio tests, Benjamini–Hochberg across the 13 genes

### 09. Cohort comparison
- `batch_effects_plots.R`: per-sample QC metrics, read depth and PCA (Figures S1–S3).
- `pah_lldd_comparison.R`: rare-variant and regulatory-region burden per TAD, LLDD vs PAH Wilcoxon test with Benjamini–Hochberg correction (Figures S4–S5).
