library(dplyr)
library(ggplot2)

lldd_file <- "results/lldd_filtered_collapsed.tsv"
pah_file  <- "results/isolated_pah_filtered_collapsed.tsv"
kg_file   <- "results/selectedEUR_TGP_filtered_collapsed.tsv"
output_path <- "figures/"

cohort_colors <- c("LLDD" = "#51C6E6", "PAH" = "#FF8833", "1KG_EUR" = "#DF536B")
exclude <- character(0)   # sample IDs to exclude (one LLDD sample with calls only in the TBX4 TAD; ID not shown)

# read and process input
lldd <- read.csv(lldd_file, header = T, sep = "\t", colClasses = "character")
pah  <- read.csv(pah_file,  header = T, sep = "\t", colClasses = "character")
kg   <- read.csv(kg_file,   header = T, sep = "\t", colClasses = "character")

lldd$cohort <- "LLDD"
pah$cohort  <- "PAH"
kg$cohort   <- "1KG_EUR"

cols <- c("Uploaded_variation", "IND", "GENE_TAD", "MAX_AF", "H3K4Me3_NHLF", "ATACseq_IMR90",
          "FILTER", "DifficultRegion", "cohort")
df <- rbind(lldd[, cols], pah[, cols], kg[, cols])

# high-quality: all FILTER values (one per duplicate VCF record) are PASS or ".", not a difficult region
filter_ok <- sapply(strsplit(df$FILTER, ";"), function(f) all(f %in% c("PASS", ".")))
df <- df[filter_ok & df$DifficultRegion != "YES", ]

df$MAX_AF[df$MAX_AF == "-"] <- NA
df$MAX_AF <- as.numeric(df$MAX_AF)
df$rare       <- is.na(df$MAX_AF) | df$MAX_AF < 0.01
df$regulatory <- df$H3K4Me3_NHLF != "-" | df$ATACseq_IMR90 != "-"

### A) One row per variant per individual per TAD
inds <- strsplit(df$IND, ";")
long <- data.frame(
  IND        = unlist(inds),
  cohort     = rep(df$cohort, lengths(inds)),
  VAR        = rep(df$Uploaded_variation, lengths(inds)),
  GENE_TAD   = rep(df$GENE_TAD, lengths(inds)),
  rare       = rep(df$rare, lengths(inds)),
  regulatory = rep(df$regulatory, lengths(inds))
)
tads <- strsplit(long$GENE_TAD, ";")
long <- long[rep(seq_len(nrow(long)), lengths(tads)), ]
long$GENE_TAD <- unlist(tads)
long <- long[!long$IND %in% exclude, ]

# every individual in every TAD, so individuals with 0 variants are counted
samples <- unique(long[, c("cohort", "IND")])
genes   <- sort(unique(long$GENE_TAD))
grid    <- merge(samples, data.frame(GENE_TAD = genes))

### B) Burden per TAD and LLDD vs PAH test
burden_per_tad <- function(class) {
  counts <- long[long[[class]], ] %>%
    distinct(cohort, IND, VAR, GENE_TAD) %>%
    count(cohort, IND, GENE_TAD, name = "n_var")
  counts <- left_join(grid, counts, by = c("cohort", "IND", "GENE_TAD"))
  counts$n_var[is.na(counts$n_var)] <- 0

  res <- counts %>%
    group_by(GENE_TAD) %>%
    summarise(
      mean_LLDD    = mean(n_var[cohort == "LLDD"]),
      mean_PAH     = mean(n_var[cohort == "PAH"]),
      mean_1KG_EUR = mean(n_var[cohort == "1KG_EUR"]),
      p = wilcox.test(n_var[cohort == "LLDD"], n_var[cohort == "PAH"], exact = FALSE)$p.value,
      .groups = "drop"
    )
  res$p[is.na(res$p)] <- 1           # identical counts in both cohorts
  res$q <- p.adjust(res$p, method = "BH")   # over the 13 pathway-gene TADs
  res
}

plot_burden <- function(res, xlab, title, file) {
  res <- res[order(res$mean_PAH), ]
  res$GENE_TAD <- factor(res$GENE_TAD, levels = res$GENE_TAD)
  res$label <- ifelse(res$q < 0.05, paste0("q=", signif(res$q, 2)), "")

  points <- data.frame(
    GENE_TAD = rep(res$GENE_TAD, 2),
    cohort   = rep(c("LLDD", "PAH"), each = nrow(res)),
    mean     = c(res$mean_LLDD, res$mean_PAH)
  )

  p <- ggplot(res, aes(y = GENE_TAD)) +
    geom_segment(aes(x = pmin(mean_LLDD, mean_PAH), xend = pmax(mean_LLDD, mean_PAH), yend = GENE_TAD),
                 color = "gray85", linewidth = 2) +
    geom_point(aes(x = mean_1KG_EUR, color = "1KG_EUR"), shape = "|", size = 6) +
    geom_point(data = points, aes(x = mean, color = cohort), size = 3) +
    geom_text(aes(x = pmax(mean_LLDD, mean_PAH), label = label), hjust = -0.3, size = 3, color = "gray30") +
    scale_color_manual(values = cohort_colors, breaks = c("LLDD", "PAH", "1KG_EUR"), name = NULL) +
    labs(x = xlab, y = NULL, title = title) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(paste0(output_path, file), p, width = 7, height = 5.5, dpi = 200)
}

### C) Plots
rare_res <- burden_per_tad("rare")
rare_res
plot_burden(rare_res, "Mean rare (<1%) variants per individual in TAD",
            "Rare-variant burden per TAD", "FigureS4_rare_burden_per_TAD.png")

reg_res <- burden_per_tad("regulatory")
reg_res
plot_burden(reg_res, "Mean variants in regulatory regions per individual in TAD",
            "Regulatory-region variant burden per TAD", "FigureS5_regulatory_burden_per_TAD.png")
