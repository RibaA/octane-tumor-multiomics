##------------------------------------------------------------
## 02. OCTANE EPIC v2 SAMPLE QC AND NOOB NORMALIZATION
##
## Project:
##   C:/octane-plasma-buffy-coat
##
## Input:
##   data/epic/proc/
##     OCTANE_EPICv2_raw_RGChannelSet.rds
##     OCTANE_EPICv2_targets_validated.csv
##
## Output:
##   data/epic/proc/
##   data/epic/proc/qc/
##
## Purpose:
##   1. Load validated raw EPIC v2 data
##   2. Calculate detection P-values
##   3. Summarize sample-level detection performance
##   4. Assess raw methylated/unmethylated signal intensity
##   5. Flag samples requiring review
##   6. Normalize all samples using Noob
##   7. Generate Beta-value and M-value matrices
##   8. Perform post-Noob density QC
##   9. Perform PCA using variable probes
##  10. Save processed objects and QC summaries
##
## Important:
##   - No sample is automatically excluded in this script.
##   - QC flags are descriptive only.
##   - Sample exclusion, if needed, will be decided after review.
##------------------------------------------------------------

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

user_lib <- "C:/Users/farno/R/win-library/4.5"

.libPaths(
  c(
    user_lib,
    .libPaths()
  )
)

############################################################
## 1. LIBRARIES
############################################################

suppressPackageStartupMessages({
  library(minfi)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(GenomicRanges)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(IlluminaHumanMethylationEPICv2manifest)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})

############################################################
## 2. DIRECTORIES
############################################################

project_dir <- getwd()

epic_dir <- file.path(
  project_dir,
  "data",
  "epic"
)

proc_dir <- file.path(
  epic_dir,
  "proc"
)

qc_dir <- file.path(
  proc_dir,
  "qc"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(epic_dir),
  dir.exists(proc_dir),
  dir.exists(qc_dir)
)

############################################################
## 3. INPUT FILES
############################################################

rgset_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_raw_RGChannelSet.rds"
)

targets_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_targets_validated.csv"
)

stopifnot(
  file.exists(rgset_file),
  file.exists(targets_file)
)

############################################################
## 4. LOAD DATA
############################################################

cat("\n============================================\n")
cat("LOADING VALIDATED EPIC v2 DATA\n")
cat("============================================\n")

rgSet_OCTANE_EPICv2 <- readRDS(
  rgset_file
)

targets <- read.csv(
  targets_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "RGChannelSet dimensions:",
  paste(
    dim(rgSet_OCTANE_EPICv2),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Targets rows:",
  nrow(targets),
  "\n"
)

print(
  annotation(
    rgSet_OCTANE_EPICv2
  )
)

stopifnot(
  inherits(
    rgSet_OCTANE_EPICv2,
    "RGChannelSet"
  ),
  ncol(
    rgSet_OCTANE_EPICv2
  ) == 23L,
  nrow(targets) == 23L,
  identical(
    sampleNames(
      rgSet_OCTANE_EPICv2
    ),
    targets$Sample_Name
  )
)

############################################################
## 5. RAW DETECTION P-VALUES
############################################################

cat("\n============================================\n")
cat("CALCULATING DETECTION P-VALUES\n")
cat("============================================\n")

detP <- minfi::detectionP(
  rgSet_OCTANE_EPICv2
)

cat(
  "Detection P-value matrix:",
  paste(
    dim(detP),
    collapse = " x "
  ),
  "\n"
)

stopifnot(
  identical(
    colnames(detP),
    sampleNames(
      rgSet_OCTANE_EPICv2
    )
  )
)

############################################################
## 6. SAMPLE-LEVEL DETECTION QC
############################################################

detection_cutoff <- 0.01

sample_qc <- data.frame(
  Sample_Name = colnames(detP),

  n_probes = nrow(detP),

  n_failed = colSums(
    detP > detection_cutoff,
    na.rm = TRUE
  ),

  prop_failed = colMeans(
    detP > detection_cutoff,
    na.rm = TRUE
  ),

  stringsAsFactors = FALSE
)

sample_qc$percent_failed <-
  100 * sample_qc$prop_failed

sample_qc <- sample_qc %>%
  arrange(
    desc(percent_failed)
  )

############################################################
## 7. ADD SAMPLE METADATA
############################################################

sample_qc <- sample_qc %>%
  left_join(
    targets %>%
      select(
        Sample_Name,
        Sample_Well,
        Sample_Plate,
        Sentrix_ID,
        Sentrix_Position
      ),
    by = "Sample_Name"
  )

############################################################
## 8. DETECTION QC FLAGS
############################################################

sample_qc <- sample_qc %>%
  mutate(

    detection_flag = case_when(

      prop_failed > 0.05 ~
        "Review: >5% failed",

      prop_failed > 0.01 ~
        "Review: 1-5% failed",

      TRUE ~
        "Pass"
    )
  )

cat("\nSample detection QC:\n")

print(
  sample_qc %>%
    select(
      Sample_Name,
      Sentrix_ID,
      Sentrix_Position,
      n_failed,
      percent_failed,
      detection_flag
    )
)

cat("\nDetection QC summary:\n")

print(
  table(
    sample_qc$detection_flag,
    useNA = "ifany"
  )
)

############################################################
## 9. DETECTION FAILURE FIGURE
############################################################

pdf(
  file.path(
    qc_dir,
    "OCTANE_EPICv2_detection_failure_by_sample.pdf"
  ),
  width = 12,
  height = 7
)

barplot(
  sample_qc$percent_failed,

  names.arg =
    sample_qc$Sample_Name,

  las = 2,

  ylab =
    "Failed probes (%)",

  main =
    "OCTANE EPIC v2: detection P-value failures",

  cex.names = 0.6
)

abline(
  h = 1,
  lty = 2
)

abline(
  h = 5,
  lty = 3
)

dev.off()

############################################################
## 10. RAW METHYLATED / UNMETHYLATED SIGNAL QC
############################################################

cat("\n============================================\n")
cat("ASSESSING RAW SIGNAL INTENSITY\n")
cat("============================================\n")

mSet_raw <- preprocessRaw(
  rgSet_OCTANE_EPICv2
)

qc_signal <- getQC(
  mSet_raw
)

print(
  qc_signal
)

qc_signal_df <- as.data.frame(
  qc_signal
) %>%
  rownames_to_column(
    "Sample_Name"
  )

sample_qc <- sample_qc %>%
  left_join(
    qc_signal_df,
    by = "Sample_Name"
  )

############################################################
## 11. RAW SIGNAL QC FIGURE
############################################################

pdf(
  file.path(
    qc_dir,
    "OCTANE_EPICv2_raw_signal_QC.pdf"
  ),
  width = 7,
  height = 6
)

plotQC(
  qc_signal
)

dev.off()

############################################################
## 12. SAVE RAW QC RESULTS
############################################################

write.csv(
  sample_qc,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_sample_QC_raw.csv"
  ),
  row.names = FALSE
)

saveRDS(
  detP,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_detectionP_raw.rds"
  )
)

saveRDS(
  qc_signal,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_raw_signal_QC.rds"
  )
)

############################################################
## 13. NOOB NORMALIZATION
############################################################

cat("\n============================================\n")
cat("RUNNING NOOB NORMALIZATION\n")
cat("============================================\n")

## At this stage all 23 samples are retained.
## If QC reveals a clear failure, we will revisit this step.

mSet_noob <- preprocessNoob(
  rgSet_OCTANE_EPICv2,
  dyeCorr = TRUE,
  verbose = TRUE
)

cat(
  "Noob MethylSet dimensions:",
  paste(
    dim(mSet_noob),
    collapse = " x "
  ),
  "\n"
)

stopifnot(
  identical(
    sampleNames(mSet_noob),
    sampleNames(
      rgSet_OCTANE_EPICv2
    )
  )
)

############################################################
## 14. GENERATE BETA AND M VALUES
############################################################

cat("\n============================================\n")
cat("GENERATING BETA AND M MATRICES\n")
cat("============================================\n")

beta_noob <- getBeta(
  mSet_noob
)

mval_noob <- getM(
  mSet_noob
)

cat(
  "Beta dimensions:",
  paste(
    dim(beta_noob),
    collapse = " x "
  ),
  "\n"
)

cat(
  "M-value dimensions:",
  paste(
    dim(mval_noob),
    collapse = " x "
  ),
  "\n"
)

stopifnot(
  identical(
    dim(beta_noob),
    dim(mval_noob)
  ),
  identical(
    rownames(beta_noob),
    rownames(mval_noob)
  ),
  identical(
    colnames(beta_noob),
    colnames(mval_noob)
  ),
  all(
    beta_noob >= 0 &
      beta_noob <= 1,
    na.rm = TRUE
  )
)

############################################################
## 15. POST-NOOB SAMPLE SIGNAL QC
############################################################

qc_signal_noob <- getQC(
  mSet_noob
)

qc_signal_noob_df <- as.data.frame(
  qc_signal_noob
) %>%
  rownames_to_column(
    "Sample_Name"
  )

names(
  qc_signal_noob_df
)[
  names(
    qc_signal_noob_df
  ) != "Sample_Name"
] <- paste0(
  names(
    qc_signal_noob_df
  )[
    names(
      qc_signal_noob_df
    ) != "Sample_Name"
  ],
  "_Noob"
)

sample_qc <- sample_qc %>%
  left_join(
    qc_signal_noob_df,
    by = "Sample_Name"
  )

pdf(
  file.path(
    qc_dir,
    "OCTANE_EPICv2_Noob_signal_QC.pdf"
  ),
  width = 7,
  height = 6
)

plotQC(
  qc_signal_noob
)

dev.off()

############################################################
## 16. POST-NOOB BETA DENSITY
############################################################

cat("\n============================================\n")
cat("CREATING POST-NOOB DENSITY QC\n")
cat("============================================\n")

## Use a fixed probe subset to keep plotting fast.
set.seed(20261002)

density_n_probes <- min(
  100000L,
  nrow(beta_noob)
)

density_probe_index <- sample(
  seq_len(
    nrow(beta_noob)
  ),
  density_n_probes
)

pdf(
  file.path(
    qc_dir,
    "OCTANE_EPICv2_Noob_beta_density.pdf"
  ),
  width = 8,
  height = 6
)

densityPlot(
  beta_noob[
    density_probe_index,
    ,
    drop = FALSE
  ],

  sampGroups = rep(
    "OCTANE",
    ncol(beta_noob)
  ),

  main =
    "OCTANE EPIC v2: Noob-normalized Beta distributions",

  xlab =
    "Beta value",

  legend = FALSE
)

dev.off()

############################################################
## 17. PCA FUNCTION
############################################################

run_pca <- function(
    m_values,
    metadata,
    n_probes = 20000L
) {

  stopifnot(
    identical(
      colnames(m_values),
      metadata$Sample_Name
    )
  )

  probe_var <- apply(
    m_values,
    1,
    var,
    na.rm = TRUE
  )

  probe_var[
    !is.finite(
      probe_var
    )
  ] <- NA_real_

  keep_var <- which(
    !is.na(
      probe_var
    )
  )

  if (length(keep_var) == 0L) {
    stop(
      "No probes with finite variance available for PCA."
    )
  }

  ordered_probes <- keep_var[
    order(
      probe_var[
        keep_var
      ],
      decreasing = TRUE
    )
  ]

  top_probes <- ordered_probes[
    seq_len(
      min(
        n_probes,
        length(
          ordered_probes
        )
      )
    )
  ]

  m_use <- m_values[
    top_probes,
    ,
    drop = FALSE
  ]

  ## Remove any probe still containing non-finite values
  valid_rows <- apply(
    m_use,
    1,
    function(x) {
      all(
        is.finite(x)
      )
    }
  )

  m_use <- m_use[
    valid_rows,
    ,
    drop = FALSE
  ]

  pca <- prcomp(
    t(m_use),
    center = TRUE,
    scale. = FALSE
  )

  variance_explained <-
    100 *
    pca$sdev^2 /
    sum(
      pca$sdev^2
    )

  pca_df <- data.frame(
    Sample_Name =
      rownames(
        pca$x
      ),

    PC1 =
      pca$x[, 1],

    PC2 =
      pca$x[, 2],

    stringsAsFactors = FALSE
  ) %>%
    left_join(
      metadata,
      by = "Sample_Name"
    )

  list(
    pca = pca,
    data = pca_df,
    variance_explained =
      variance_explained,
    n_probes_used =
      nrow(m_use)
  )
}

############################################################
## 18. POST-NOOB PCA
############################################################

cat("\n============================================\n")
cat("RUNNING POST-NOOB PCA\n")
cat("============================================\n")

targets_ordered <- targets[
  match(
    colnames(mval_noob),
    targets$Sample_Name
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    targets_ordered$Sample_Name,
    colnames(mval_noob)
  )
)

pca_result <- run_pca(
  m_values = mval_noob,
  metadata = targets_ordered,
  n_probes = 20000L
)

cat(
  "PCA probes used:",
  pca_result$n_probes_used,
  "\n"
)

cat(
  "PC1 variance explained:",
  round(
    pca_result$variance_explained[1],
    2
  ),
  "%\n"
)

cat(
  "PC2 variance explained:",
  round(
    pca_result$variance_explained[2],
    2
  ),
  "%\n"
)

############################################################
## 19. PCA FIGURES
############################################################

pca_plot_slide <- ggplot(
  pca_result$data,
  aes(
    x = PC1,
    y = PC2,
    shape = factor(
      Sentrix_ID
    )
  )
) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "OCTANE EPIC v2: PCA after Noob normalization",

    x = paste0(
      "PC1 (",
      round(
        pca_result$variance_explained[1],
        1
      ),
      "%)"
    ),

    y = paste0(
      "PC2 (",
      round(
        pca_result$variance_explained[2],
        1
      ),
      "%)"
    ),

    shape =
      "Sentrix ID"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    qc_dir,
    "OCTANE_EPICv2_Noob_PCA_by_SentrixID.pdf"
  ),
  plot = pca_plot_slide,
  width = 8,
  height = 6
)

############################################################
## 20. SAVE PCA TABLE
############################################################

write.csv(
  pca_result$data,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_Noob_PCA_coordinates.csv"
  ),
  row.names = FALSE
)

############################################################
## 21. SAVE FINAL SAMPLE QC TABLE
############################################################

write.csv(
  sample_qc,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_sample_QC.csv"
  ),
  row.names = FALSE
)

############################################################
## 22. SAVE NORMALIZED OBJECTS
############################################################

cat("\n============================================\n")
cat("SAVING NORMALIZED DATA\n")
cat("============================================\n")

saveRDS(
  mSet_noob,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_Noob_MethylSet.rds"
  )
)

saveRDS(
  beta_noob,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_Noob_beta.rds"
  )
)

saveRDS(
  mval_noob,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_Noob_M.rds"
  )
)

saveRDS(
  pca_result,
  file.path(
    proc_dir,
    "OCTANE_EPICv2_Noob_PCA.rds"
  )
)

############################################################
## 23. FINAL SUMMARY
############################################################

cat("\n============================================\n")
cat("SCRIPT 02 OCTANE EPIC v2 QC + NOOB COMPLETE\n")
cat("============================================\n")

cat(
  "Samples evaluated:",
  ncol(
    rgSet_OCTANE_EPICv2
  ),
  "\n"
)

cat(
  "Samples normalized:",
  ncol(
    beta_noob
  ),
  "\n"
)

cat(
  "Probes in normalized matrix:",
  nrow(
    beta_noob
  ),
  "\n"
)

cat(
  "Samples with >1% failed probes:",
  sum(
    sample_qc$prop_failed > 0.01
  ),
  "\n"
)

cat(
  "Samples with >5% failed probes:",
  sum(
    sample_qc$prop_failed > 0.05
  ),
  "\n"
)

cat(
  "PC1 variance explained:",
  round(
    pca_result$variance_explained[1],
    2
  ),
  "%\n"
)

cat(
  "PC2 variance explained:",
  round(
    pca_result$variance_explained[2],
    2
  ),
  "%\n"
)

cat(
  "\nIMPORTANT:\n",
  "No OCTANE sample was excluded automatically.\n",
  "Review detection P-values, signal QC, density plots, and PCA ",
  "before probe filtering.\n",
  sep = ""
)

cat(
  "\nProcessed data directory:\n",
  proc_dir,
  "\n"
)

cat(
  "\nQC figure directory:\n",
  qc_dir,
  "\n"
)
