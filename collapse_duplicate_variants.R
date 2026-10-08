
library(dplyr)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
vcf  <- args[2]
if (is.na(vcf)) stop("Usage: Rscript 05_collapse_duplicate_variants.R <name> <vt.vcf.gz>")

# VEP tsv output sits next to the VCF (04_annotate_vcf.sh): <x>.vt.vcf.gz -> <x>.vt.tsv
input <- sub("\\.vcf\\.gz$", ".tsv", vcf)
genes_input <- "data/genes.txt"
output_path <- "results/"
xlsx_file <- paste0(output_path, name, ".xlsx")
collapsed_output <- paste0(output_path, name, "_filtered_collapsed.tsv")

## ---- read main VEP/annotation TSV ----
lines <- readLines(input)
lines <- lines[!grepl("^##", lines)]
lines[1] <- sub("^#", "", lines[1])
df <- read.table(text = lines, header=TRUE, sep="\t")
genes <- read.table(genes_input, sep="\t", header=FALSE, col.names=c("chr","start","end","gene"))

## ---- parse CHR/POS ----
tmp <- t(sapply(strsplit(df$Location, "[:;-]"), `[`, 1:2))
df$CHR <- as.integer(sub("^chr", "", tmp[,1], ignore.case = TRUE))
df$POS <- as.integer(tmp[,2])
rm(tmp)

## ---- collapse duplicate variants ----
df <- df %>%
  group_by(Uploaded_variation) %>%
  summarise(
    across(everything(), ~ {
      u <- unique(.)
      if (length(u) == 1) u else paste(u, collapse = ";")
    }),
    .groups = "drop"
  )

## ---- add gene names to non-coding variants within the gene's TAD ----
df$GENE_TAD <- NA
genes$chr <- as.integer(sub("^chr", "", genes$chr))

for (i in seq_len(nrow(genes))) {
  idx <- df$CHR == genes$chr[i] &
    df$POS >= genes$start[i] &
    df$POS <= genes$end[i]
  df$GENE_TAD[idx] <- ifelse(
    is.na(df$GENE_TAD[idx]),
    genes$gene[i],
    paste(df$GENE_TAD[idx], genes$gene[i], sep = ";")
  )
}

df <- df[!is.na(df$Uploaded_variation), ]
## -------------------------------------------------------------------------
## TIER CLASSIFICATION
## -------------------------------------------------------------------------

## Tier 1: Rare (gnomAD MAF <0.001) variants 
##         fall within ATAC-seq peaks or H3K27ac/H3K4me3 marked regions
##         have high CADD (PHRED >= 15, top 3%)

## Tier 2: Rare (gnomAD MAF <0.001) variants
##         fall within ATAC-seq peaks or H3K27ac/H3K4me3 marked regions
##         moderate CADD (>=10, <15)

## Tier 3: Remaining

## --- Rarity ---
df$MAX_AF[df$MAX_AF == "-"] <- NA
df$MAX_AF <- as.numeric(df$MAX_AF)

## --- Regulatory region: ATAC-seq peak or H3K4me3 mark present ---
df$H3K4Me3_NHLF[df$H3K4Me3_NHLF == "-"] <- NA
df$ATACseq_IMR90[df$ATACseq_IMR90 == "-"] <- NA


## --- Tier assignment ---
## older CADD plugin versions name the column CADD_phred
if (!"CADD_PHRED" %in% names(df) & "CADD_phred" %in% names(df)) df <- df %>% rename(CADD_PHRED = CADD_phred)
df$MAX_AF <- as.numeric(df$MAX_AF)

if ("CADD_PHRED" %in% names(df)) {
  df$CADD_PHRED <- as.numeric(df$CADD_PHRED)
  df$Tier <- case_when(
    df$MAX_AF <= 0.001 & (!(is.na(df$H3K4Me3_NHLF)) | !(is.na(df$ATACseq_IMR90))) & df$CADD_PHRED >= 15 ~ 1,
    df$MAX_AF <= 0.001 & (!(is.na(df$H3K4Me3_NHLF)) | !(is.na(df$ATACseq_IMR90))) & df$CADD_PHRED >= 10 & df$CADD_PHRED < 15 ~ 2,
    TRUE ~ 3
  )
} else {
  ## no CADD annotation (e.g. 1KGP reference): tiers need CADD, so leave them empty
  print("No CADD_PHRED column: CADD_PHRED and Tier left empty")
  df$CADD_PHRED <- NA
  df$Tier <- NA
}

#remove columns where all values are -
names(df)[sapply(df, \(x) all(x == "-"))]
df <- df %>% select(-any_of(c("SOURCE", "CHECK_REF")))

df <- df %>% select(any_of(col_order), everything())

## Write collapsed df
write.table(df, collapsed_output,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            na="-")
print(paste("collapsed df written to", collapsed_output))

## Write each gene as an excel page
wb <- createWorkbook()
for(i in seq_len(nrow(genes))) {
  sub <- df[df$CHR == genes$chr[i] & df$POS >= genes$start[i] & df$POS <= genes$end[i], ]
  if (nrow(sub) > 0) {
    addWorksheet(wb, genes$gene[i])
    writeData(wb, sheet = genes$gene[i], sub)
  }
}
saveWorkbook(wb, xlsx_file, overwrite = T)

print(paste("Variants are written to ", xlsx_file))
