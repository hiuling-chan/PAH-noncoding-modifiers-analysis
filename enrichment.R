
pressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
})

name <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(name)) stop("Usage: Rscript 10_enrichment.R <name>  (e.g. isolated_pah, lldd, selectedEUR_TGP)")

results_dir <- "results"

variants_file <- file.path(results_dir, paste0(name, "_filtered_collapsed.tsv"))

today <- format(Sys.Date(), "%Y%m%d")
xlsx_file <- file.path(results_dir, paste0("TableS5_", name, "_burden_enrichment_HA", today, ".xlsx"))
output    <- file.path(results_dir, paste0("prom_enh_stats_", name, ".tsv"))

strip_chr <- function(x) sub("^chr", "", as.character(x))

# read and process input
tads  <- read.csv("data/tad_lengths.txt", sep = "\t", header = FALSE,
                     col.names = c("GENE","TAD_length"))
enh   <- read.csv("data/ATACseq_IMR90.bed", header = F, sep = "\t")
prom  <- read.csv("data/H3K4Me3_NHLF.bed", header = F, sep = "\t")
genes <- read.table("data/genes.txt", sep = "\t", header = FALSE,
                     col.names = c("chr","start","end","gene"))
gtf_file <- "data/gencode.v46.basic.annotation.gtf.gz"
df    <- read.csv(variants_file, header = T, sep = "\t")
wb    <- createWorkbook()

# remove chr from all coordinate sources
enh$V1  <- strip_chr(enh$V1)
prom$V1 <- strip_chr(prom$V1)
genes$chr <- strip_chr(genes$chr)

### A) Summarize
## 1. Count average number of variants per individual
inds <- strsplit(as.character(df$IND), ";")
vars <- df$Uploaded_variation
tmp <- do.call(rbind, Map(function(i, v) data.frame(IND = i, VAR = v), inds, vars))
df_count <- data.frame(
  IND = names(tapply(tmp$VAR, tmp$IND, length)),
  SUM_VAR = as.integer(tapply(tmp$VAR, tmp$IND, length)),
  row.names = NULL
)
addWorksheet(wb, "samples")
writeData(wb, sheet = "samples", df_count)

## 2. Count number of variants per pathway gene enhancer/promoter
df$promoter <- ifelse(!is.na(df$H3K4Me3_NHLF) & df$H3K4Me3_NHLF != "-", 1, NA)
df$enhancer <- ifelse(!is.na(df$ATACseq_IMR90) & df$ATACseq_IMR90 != "-", 1, NA)
df_count <- data.frame(
  gene = sort(unique(df$GENE_TAD)),
  promoter = as.integer(tapply(df$promoter, df$GENE_TAD, sum, na.rm = TRUE)),
  enhancer = as.integer(tapply(df$enhancer, df$GENE_TAD, sum, na.rm = TRUE)),
  row.names = NULL
)
addWorksheet(wb, "promoter_enhancer_counts")
writeData(wb, sheet = "promoter_enhancer_counts", df_count)

## 3. Count number of variants per gene
gene <- strsplit(as.character(df$GENE_TAD), ";")
tmp <- do.call(rbind, Map(function(i, v) data.frame(GENE_TAD = i, VAR = v), gene, vars))
df_count <- data.frame(
  GENE = names(tapply(tmp$VAR, tmp$GENE_TAD, length)),
  SUM_VAR = as.integer(tapply(tmp$VAR, tmp$GENE_TAD, length)),
  row.names = NULL
)
addWorksheet(wb, "variants_per_gene")
writeData(wb, sheet = "variants_per_gene", df_count)

### B) clustering
## 1. parse GTF -> GRanges of every gene body genome-wide (not just the 13
##    pathway genes)
gtf_raw <- readLines(gzfile(gtf_file, "rt"))
gtf_raw <- gtf_raw[!startsWith(gtf_raw, "#")]
gtf_split <- strsplit(gtf_raw, "\t")
gtf_gene_rows <- gtf_split[vapply(gtf_split, `[`, character(1), 3) == "gene"]
all_genes <- data.frame(
  chr   = strip_chr(vapply(gtf_gene_rows, `[`, character(1), 1)),
  start = as.integer(vapply(gtf_gene_rows, `[`, character(1), 4)),
  end   = as.integer(vapply(gtf_gene_rows, `[`, character(1), 5)),
  attr  = vapply(gtf_gene_rows, `[`, character(1), 9)
)
all_genes$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", all_genes$attr)
all_genes_gr <- GRanges(all_genes$chr, IRanges(all_genes$start, all_genes$end),
                         gene_name = all_genes$gene_name)

## 2. add number of peaks and peak length to each TAD, and count variants
##    against a background that excludes peaks AND every gene body (pathway
##    gene included)
enh$lengths  <- enh$V7 - enh$V6
prom$lengths <- prom$V7 - prom$V6

count_vars <- function(chr, start, end) sum(df$CHR == chr & df$POS >= start & df$POS <= end)
count_vars_gr <- function(gr) {
  if (length(gr) == 0) return(0)
  sum(mapply(count_vars, as.character(seqnames(gr)), start(gr), end(gr)))
}
enh$var  <- mapply(count_vars, enh$V1,  enh$V6 + 1, enh$V7)   # BED start is 0-based
prom$var <- mapply(count_vars, prom$V1, prom$V6 + 1, prom$V7)

for (col in c("number_of_enh_peaks","enh_length_sum","enh_var_sum",
              "number_of_prom_peaks","prom_length_sum","prom_var_sum",
              "pathway_gene_length","other_genic_length",
              "enh_background_length","prom_background_length",
              "enh_background_var","prom_background_var")) tads[[col]] <- 0

for (i in seq_len(nrow(tads))) {
  gene_id <- tads$GENE[i]
  tb <- genes[genes$gene == gene_id, ]   # genes.txt doubles as TAD coordinates
  if (nrow(tb) == 0) { warning("No TAD coords for ", gene_id, "; skipping"); next }
  tad_gr <- GRanges(tb$chr[1], IRanges(tb$start[1] + 1, tb$end[1]))

  e <- enh[enh$V4 == gene_id, ]
  p <- prom[prom$V4 == gene_id, ]
  tads$number_of_enh_peaks[i]  <- nrow(e)
  tads$enh_length_sum[i]       <- if (nrow(e)) sum(e$lengths) else 0
  tads$enh_var_sum[i]          <- if (nrow(e)) sum(e$var) else 0
  tads$number_of_prom_peaks[i] <- nrow(p)
  tads$prom_length_sum[i]      <- if (nrow(p)) sum(p$lengths) else 0
  tads$prom_var_sum[i]         <- if (nrow(p)) sum(p$var) else 0

  enh_gr  <- if (nrow(e)) reduce(GRanges(e$V1, IRanges(e$V6 + 1, e$V7))) else GRanges()
  prom_gr <- if (nrow(p)) reduce(GRanges(p$V1, IRanges(p$V6 + 1, p$V7))) else GRanges()

  genes_in_tad <- reduce(subsetByOverlaps(all_genes_gr, tad_gr))
  genes_in_tad <- GenomicRanges::intersect(genes_in_tad, tad_gr)   # clip to TAD
  pathway_gr <- GenomicRanges::intersect(genes_in_tad, all_genes_gr[all_genes_gr$gene_name == gene_id])
  tads$pathway_gene_length[i] <- sum(width(pathway_gr))
  tads$other_genic_length[i]  <- sum(width(genes_in_tad)) - tads$pathway_gene_length[i]

  enh_bg_gr  <- GenomicRanges::setdiff(tad_gr, reduce(c(enh_gr,  genes_in_tad)))
  prom_bg_gr <- GenomicRanges::setdiff(tad_gr, reduce(c(prom_gr, genes_in_tad)))
  tads$enh_background_length[i]  <- sum(width(enh_bg_gr))
  tads$prom_background_length[i] <- sum(width(prom_bg_gr))
  tads$enh_background_var[i]     <- count_vars_gr(enh_bg_gr)
  tads$prom_background_var[i]    <- count_vars_gr(prom_bg_gr)
}

## 3. Fisher exact test AND Poisson rate-ratio test: variant density in peaks vs background given expected genomic length
fisher_p <- function(vf, vb, lf, lb) mapply(function(a,b,c,d)
fisher.test(matrix(c(a,c,b,d), nrow = 2), alternative = "greater")$p.value, vf, vb, lf, lb)
rate_p <- function(vf, vb, lf, lb) mapply(function(a,b,c,d)
poisson.test(c(a, b), T = c(c, d), alternative = "greater")$p.value, vf, vb, lf, lb)

tads$enh_p_fisher   <- fisher_p(tads$enh_var_sum,  tads$enh_background_var,  tads$enh_length_sum,  tads$enh_background_length)
tads$prom_p_fisher  <- fisher_p(tads$prom_var_sum, tads$prom_background_var, tads$prom_length_sum, tads$prom_background_length)
tads$enh_p_poisson  <- rate_p(tads$enh_var_sum,  tads$enh_background_var,  tads$enh_length_sum,  tads$enh_background_length)
tads$prom_p_poisson <- rate_p(tads$prom_var_sum, tads$prom_background_var, tads$prom_length_sum, tads$prom_background_length)

tads$enh_p_fisher_adj   <- p.adjust(tads$enh_p_fisher,   method = "BH")
tads$prom_p_fisher_adj  <- p.adjust(tads$prom_p_fisher,  method = "BH")
tads$enh_p_poisson_adj  <- p.adjust(tads$enh_p_poisson,  method = "BH")
tads$prom_p_poisson_adj <- p.adjust(tads$prom_p_poisson, method = "BH")

tads[tads$prom_p_fisher_adj <= 0.05 | tads$prom_p_poisson_adj <= 0.05, ]$GENE   # genes with significant promoter enrichment
tads[tads$enh_p_fisher_adj <= 0.05 | tads$enh_p_poisson_adj <= 0.05, ]$GENE    # genes with significant enhancer enrichment

num_cols <- sapply(tads, is.numeric)
tads[num_cols] <- lapply(tads[num_cols], round, digits = 4)

# >>>>>> save output, named after `name`
write.table(tads, output, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote", output, "\n")
addWorksheet(wb, "features")
writeData(wb, sheet = "features", tads)

# >>>>>> legend sheet explaining every column across all sheets
legend <- data.frame(
  column = c(
    "GENE", "TAD_length", "number_of_enh_peaks", "enh_length_sum", "enh_var_sum",
    "number_of_prom_peaks", "prom_length_sum", "prom_var_sum",
    "pathway_gene_length", "other_genic_length",
    "enh_background_length", "prom_background_length",
    "enh_background_var", "prom_background_var",
    "enh_p_fisher", "prom_p_fisher", "enh_p_poisson", "prom_p_poisson",
    "enh_p_fisher_adj", "prom_p_fisher_adj", "enh_p_poisson_adj", "prom_p_poisson_adj"
  ),
  description = c(
    "Pathway gene",
    "Length (bp) of the gene's TAD",
    "Number of ATAC-seq enhancer peaks",
    "Total bp covered by enhancer peaks",
    "Total variants within the enhancer peaks",
    "Number of H3K4me3 promoter peaks",
    "Total bp covered by promoter peaks",
    "Total variants within the promoter peaks",
    "bp of the TAD covered by gene regions of the pathway gene",
    "bp of the TAD covered by gene regions of other genes, excluding the pathway gene",
    "bp of the TAD that is neither an enhancer peak nor any gene region",
    "bp of the TAD that is neither a promoter peak nor any gene region",
    "Variants falling in enh_background_length (outside enhancer peaks and outside all gene regions)",
    "Variants falling in prom_background_length (outside promoter peaks and outside all gene regions)",
    "Fisher's exact test p-value: enhancer variant/length vs background variant/length",
    "Fisher's exact test p-value: promoter variant/length vs background variant/length",
    "Poisson rate-ratio test p-value: enhancer variant density vs background variant density",
    "Poisson rate-ratio test p-value: promoter variant density vs background variant density",
    "BH-adjusted enh_p_fisher across the 13 genes",
    "BH-adjusted prom_p_fisher across the 13 genes",
    "BH-adjusted enh_p_poisson across the 13 genes",
    "BH-adjusted prom_p_poisson across the 13 genes"
  )
)
addWorksheet(wb, "legend")
writeData(wb, sheet = "legend", legend)
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("Wrote", xlsx_file, "\n")
