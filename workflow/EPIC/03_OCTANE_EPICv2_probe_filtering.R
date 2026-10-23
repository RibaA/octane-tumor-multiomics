##------------------------------------------------------------
## 03. OCTANE EPIC v2 PROBE-LEVEL FILTERING
##
## Project:
##   C:/octane-plasma-buffy-coat
##
## Inputs:
##   data/epic/proc/
##     OCTANE_EPICv2_Noob_beta.rds
##     OCTANE_EPICv2_Noob_M.rds
##     OCTANE_EPICv2_detectionP_raw.rds
##     OCTANE_EPICv2_targets_validated.csv
##
## Outputs:
##   data/epic/proc/
##     OCTANE_EPICv2_filtered_beta.rds
##     OCTANE_EPICv2_filtered_M.rds
##     OCTANE_EPICv2_filtered_annotation.rds
##     OCTANE_EPICv2_probe_filtering_summary.csv
##     OCTANE_EPICv2_probe_QC_table.csv
##
## Purpose:
##   1. Apply probe-level detection P-value filtering
##   2. Remove probes with missing/non-finite normalized values
##   3. Attach EPIC v2 hg38 annotation
##   4. Remove probes without usable genomic coordinates
##   5. Flag sex-chromosome probes, but retain them for now
##   6. Save filtered Beta/M matrices for downstream projection
##
## Notes:
##   - OCTANE contains EPIC v2 only.
##   - No EPIC v1/v2 harmonization is required.
##   - Detection P-value threshold = 0.01.
##   - A probe must pass in >=90% of samples.
##   - With 23 samples, this corresponds to at least 21/23
##     samples passing detection P < 0.01.
##   - Sex-chromosome probes are retained at this stage.
##------------------------------------------------------------

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

############################################################
## 1. USER LIBRARY
############################################################

user_lib <- "C:/Users/farno/R/win-library/4.5"

.libPaths(
  c(
    user_lib,
    .libPaths()
  )
)

############################################################
## 2. LIBRARIES
############################################################

suppressPackageStartupMessages({
  library(minfi)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
  library(IlluminaHumanMethylationEPICv2manifest)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
})

############################################################
## 3. DIRECTORIES
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

stopifnot(
  dir.exists(epic_dir),
  dir.exists(proc_dir)
)

############################################################
## 4. INPUT FILES
############################################################

beta_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_Noob_beta.rds"
)

m_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_Noob_M.rds"
)

detp_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_detectionP_raw.rds"
)

targets_file <- file.path(
  proc_dir,
  "OCTANE_EPICv2_targets_validated.csv"
)

stopifnot(
  file.exists(beta_file),
  file.exists(m_file),
  file.exists(detp_file),
  file.exists(targets_file)
)

############################################################
## 5. LOAD DATA
############################################################

cat("\n============================================\n")
cat("LOADING OCTANE EPIC v2 DATA\n")
cat("============================================\n")

beta_noob <- readRDS(
  beta_file
)

mval_noob <- readRDS(
  m_file
)

detP <- readRDS(
  detp_file
)

targets <- read.csv(
  targets_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
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
  "M dimensions:",
  paste(
    dim(mval_noob),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Detection P dimensions:",
  paste(
    dim(detP),
    collapse = " x "
  ),
  "\n"
)

############################################################
## 6. BASIC MATRIX VALIDATION
############################################################

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
  identical(
    colnames(beta_noob),
    colnames(detP)
  ),
  identical(
    colnames(beta_noob),
    targets$Sample_Name
  )
)

cat(
  "Samples:",
  ncol(beta_noob),
  "\n"
)

cat(
  "Normalized EPIC v2 probes:",
  nrow(beta_noob),
  "\n"
)

############################################################
## 7. ALIGN DETECTION P TO NORMALIZED PROBES
############################################################

cat("\n============================================\n")
cat("ALIGNING DETECTION P-VALUES\n")
cat("============================================\n")

common_probes <- intersect(
  rownames(beta_noob),
  rownames(detP)
)

cat(
  "Normalized probes:",
  nrow(beta_noob),
  "\n"
)

cat(
  "Detection-P probes:",
  nrow(detP),
  "\n"
)

cat(
  "Common probe IDs:",
  length(common_probes),
  "\n"
)

if (length(common_probes) == 0L) {
  stop(
    "No common probe IDs between normalized matrices and detection P matrix."
  )
}

beta_work <- beta_noob[
  common_probes,
  ,
  drop = FALSE
]

mval_work <- mval_noob[
  common_probes,
  ,
  drop = FALSE
]

detP_work <- detP[
  common_probes,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(beta_work),
    rownames(detP_work)
  )
)

############################################################
## 8. DETECTION P-VALUE FILTER
############################################################

cat("\n============================================\n")
cat("PROBE-LEVEL DETECTION P FILTER\n")
cat("============================================\n")

detection_cutoff <- 0.01
minimum_pass_fraction <- 0.90

n_samples <- ncol(
  detP_work
)

minimum_pass_samples <- ceiling(
  minimum_pass_fraction *
    n_samples
)

maximum_failed_samples <-
  n_samples -
  minimum_pass_samples

cat(
  "Detection P threshold:",
  detection_cutoff,
  "\n"
)

cat(
  "Minimum passing fraction:",
  minimum_pass_fraction,
  "\n"
)

cat(
  "Number of samples:",
  n_samples,
  "\n"
)

cat(
  "Minimum samples required to pass:",
  minimum_pass_samples,
  "\n"
)

cat(
  "Maximum failed samples allowed:",
  maximum_failed_samples,
  "\n"
)

n_pass_detection <- rowSums(
  detP_work < detection_cutoff,
  na.rm = TRUE
)

n_fail_detection <- rowSums(
  detP_work >= detection_cutoff |
    is.na(detP_work)
)

keep_detection <-
  n_pass_detection >= minimum_pass_samples

cat(
  "Probes passing detection filter:",
  sum(keep_detection),
  "\n"
)

cat(
  "Probes removed by detection filter:",
  sum(!keep_detection),
  "\n"
)

############################################################
## 9. APPLY DETECTION FILTER
############################################################

beta_det <- beta_work[
  keep_detection,
  ,
  drop = FALSE
]

mval_det <- mval_work[
  keep_detection,
  ,
  drop = FALSE
]

cat(
  "Dimensions after detection filter:",
  paste(
    dim(beta_det),
    collapse = " x "
  ),
  "\n"
)

############################################################
## 10. MISSING / NON-FINITE VALUE FILTER
############################################################

cat("\n============================================\n")
cat("FILTERING MISSING / NON-FINITE VALUES\n")
cat("============================================\n")

beta_complete <- apply(
  beta_det,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

m_complete <- apply(
  mval_det,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

keep_complete <-
  beta_complete &
  m_complete

cat(
  "Probes with complete finite Beta and M values:",
  sum(keep_complete),
  "\n"
)

cat(
  "Probes removed for NA/non-finite values:",
  sum(!keep_complete),
  "\n"
)

beta_complete_mat <- beta_det[
  keep_complete,
  ,
  drop = FALSE
]

mval_complete_mat <- mval_det[
  keep_complete,
  ,
  drop = FALSE
]

############################################################
## 11. LOAD EPIC v2 hg38 ANNOTATION
############################################################

cat("\n============================================\n")
cat("LOADING EPIC v2 hg38 ANNOTATION\n")
cat("============================================\n")

annotation_epicv2 <- minfi::getAnnotation(
  IlluminaHumanMethylationEPICv2anno.20a1.hg38
)

cat(
  "Annotation rows:",
  nrow(annotation_epicv2),
  "\n"
)

cat(
  "Annotation columns:",
  ncol(annotation_epicv2),
  "\n"
)

cat(
  "\nAnnotation column names:\n"
)

print(
  colnames(annotation_epicv2)
)

cat(
  "\nFirst annotation row names:\n"
)

print(
  head(
    rownames(annotation_epicv2)
  )
)

############################################################
## 12. MATCH FILTERED PROBES TO ANNOTATION
############################################################

probe_ids <- rownames(
  beta_complete_mat
)

annotation_ids <- rownames(
  annotation_epicv2
)

cat("\n============================================\n")
cat("MATCHING FILTERED PROBES TO ANNOTATION\n")
cat("============================================\n")

cat(
  "Filtered probe IDs:",
  length(probe_ids),
  "\n"
)

cat(
  "Annotation probe IDs:",
  length(annotation_ids),
  "\n"
)

direct_match <- match(
  probe_ids,
  annotation_ids
)

cat(
  "Direct matches:",
  sum(!is.na(direct_match)),
  "\n"
)

cat(
  "Direct non-matches:",
  sum(is.na(direct_match)),
  "\n"
)

############################################################
## 13. CANONICAL EPIC v2 IDs
############################################################

canonical_probe_id <- function(x) {

  sub(
    "_.*$",
    "",
    x
  )
}

probe_ids_canonical <- canonical_probe_id(
  probe_ids
)

annotation_ids_canonical <- canonical_probe_id(
  annotation_ids
)

cat(
  "Unique canonical filtered CpGs:",
  length(
    unique(
      probe_ids_canonical
    )
  ),
  "\n"
)

cat(
  "Unique canonical annotation CpGs:",
  length(
    unique(
      annotation_ids_canonical
    )
  ),
  "\n"
)

############################################################
## 14. MATCH USING DIRECT IDs FIRST, THEN CANONICAL IDs
############################################################

## Start with direct matches where possible.

annotation_index <- direct_match

## For probes without a direct match, attempt matching using
## the canonical CpG ID.

missing_direct <- which(
  is.na(
    annotation_index
  )
)

if (length(missing_direct) > 0L) {

  ## Annotation may contain multiple EPIC v2 replicate IDs
  ## corresponding to the same canonical CpG.
  ## Use the first annotation entry for coordinate assignment.

  first_canonical_annotation <-
    !duplicated(
      annotation_ids_canonical
    )

  canonical_lookup_ids <-
    annotation_ids_canonical[
      first_canonical_annotation
    ]

  canonical_lookup_index <-
    which(
      first_canonical_annotation
    )

  canonical_match <- match(
    probe_ids_canonical[
      missing_direct
    ],
    canonical_lookup_ids
  )

  matched_canonical <-
    !is.na(
      canonical_match
    )

  annotation_index[
    missing_direct[
      matched_canonical
    ]
  ] <-
    canonical_lookup_index[
      canonical_match[
        matched_canonical
      ]
    ]
}

cat(
  "Total probes matched after direct/canonical matching:",
  sum(!is.na(annotation_index)),
  "\n"
)

cat(
  "Unmatched probes:",
  sum(is.na(annotation_index)),
  "\n"
)

############################################################
## 15. BUILD FILTERED ANNOTATION TABLE
############################################################

usable_annotation <- !is.na(
  annotation_index
)

annotation_filtered <- annotation_epicv2[
  annotation_index[
    usable_annotation
  ],
  ,
  drop = FALSE
]

annotation_df <- as.data.frame(
  annotation_filtered
)

annotation_df$Probe_ID <-
  probe_ids[
    usable_annotation
  ]

annotation_df$Canonical_CpG <-
  probe_ids_canonical[
    usable_annotation
  ]

## Check coordinate-column names before proceeding.
cat(
  "\nAvailable coordinate-related columns:\n"
)

print(
  intersect(
    c(
      "chr",
      "pos",
      "AddressA",
      "AddressB"
    ),
    colnames(
      annotation_df
    )
  )
)

############################################################
## 15B. STANDARDIZE hg38 COORDINATES
############################################################

if (!all(
  c(
    "chr",
    "pos"
  ) %in%
    colnames(
      annotation_df
    )
)) {

  stop(
    "Expected annotation columns 'chr' and 'pos' were not found. ",
    "Please inspect colnames(annotation_df)."
  )
}

annotation_df$chr <-
  as.character(
    annotation_df$chr
  )

annotation_df$position <-
  as.integer(
    annotation_df$pos
  )

usable_coordinate <-
  !is.na(
    annotation_df$chr
  ) &
  annotation_df$chr != "" &
  !is.na(
    annotation_df$position
  )

cat(
  "\nProbes successfully matched to annotation:",
  nrow(annotation_df),
  "\n"
)

cat(
  "Probes with usable hg38 coordinates:",
  sum(usable_coordinate),
  "\n"
)

cat(
  "Probes without usable hg38 coordinates:",
  sum(!usable_coordinate),
  "\n"
)

############################################################
## 15C. FILTER MATRICES TO ANNOTATED PROBES
############################################################

matched_probe_ids <-
  probe_ids[
    usable_annotation
  ]

beta_annotated <-
  beta_complete_mat[
    matched_probe_ids,
    ,
    drop = FALSE
  ]

mval_annotated <-
  mval_complete_mat[
    matched_probe_ids,
    ,
    drop = FALSE
  ]

## Apply usable-coordinate filter

beta_coord <-
  beta_annotated[
    usable_coordinate,
    ,
    drop = FALSE
  ]

mval_coord <-
  mval_annotated[
    usable_coordinate,
    ,
    drop = FALSE
  ]

annotation_coord <-
  annotation_df[
    usable_coordinate,
    ,
    drop = FALSE
  ]

stopifnot(
  nrow(beta_coord) ==
    nrow(mval_coord),

  nrow(beta_coord) ==
    nrow(annotation_coord),

  identical(
    rownames(beta_coord),
    annotation_coord$Probe_ID
  )
)

cat(
  "\nFinal annotated matrix dimensions:",
  paste(
    dim(beta_coord),
    collapse = " x "
  ),
  "\n"
)

############################################################
## 16. FLAG SEX-CHROMOSOME PROBES
############################################################

annotation_coord$sex_chromosome <-
  annotation_coord$chr %in%
  c(
    "chrX",
    "chrY",
    "X",
    "Y"
  )

cat(
  "\nSex-chromosome probes:",
  sum(
    annotation_coord$sex_chromosome
  ),
  "\n"
)

cat(
  "Autosomal probes:",
  sum(
    !annotation_coord$sex_chromosome
  ),
  "\n"
)

cat(
  "\nSex-chromosome probes are RETAINED in Script 03.\n"
)

############################################################
## 17. BUILD PROBE QC TABLE
############################################################

probe_qc <- data.frame(
  Probe_ID = rownames(beta_work),
  Canonical_CpG =
    canonical_probe_id(
      rownames(beta_work)
    ),
  n_pass_detection =
    n_pass_detection,
  n_fail_detection =
    n_fail_detection,
  pass_detection =
    keep_detection,
  stringsAsFactors = FALSE
)

complete_lookup <- data.frame(
  Probe_ID =
    rownames(beta_det),
  complete_finite =
    keep_complete,
  stringsAsFactors = FALSE
)

probe_qc <- probe_qc %>%
  left_join(
    complete_lookup,
    by = "Probe_ID"
  )

probe_qc$complete_finite[
  is.na(
    probe_qc$complete_finite
  )
] <- FALSE

coordinate_lookup <- data.frame(
  Probe_ID =
    annotation_df$Probe_ID,
  usable_coordinate =
    usable_coordinate,
  stringsAsFactors = FALSE
)

probe_qc <- probe_qc %>%
  left_join(
    coordinate_lookup,
    by = "Probe_ID"
  )

probe_qc$usable_coordinate[
  is.na(
    probe_qc$usable_coordinate
  )
] <- FALSE

probe_qc$final_keep <-
  probe_qc$pass_detection &
  probe_qc$complete_finite &
  probe_qc$usable_coordinate

############################################################
## 18. FILTERING SUMMARY
############################################################

filter_summary <- data.frame(
  step = c(
    "Noob-normalized probes",
    "Common with detection-P matrix",
    "Pass detection P filter",
    "Complete finite Beta/M values",
    "Usable hg38 coordinates",
    "Final retained probes"
  ),

  n_probes = c(
    nrow(beta_noob),
    nrow(beta_work),
    nrow(beta_det),
    nrow(beta_complete_mat),
    nrow(beta_coord),
    nrow(beta_coord)
  ),

  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("FILTERING SUMMARY\n")
cat("============================================\n")

print(
  filter_summary
)

############################################################
## 19. FINAL MATRIX VALIDATION
############################################################

stopifnot(
  identical(
    dim(beta_coord),
    dim(mval_coord)
  ),
  identical(
    rownames(beta_coord),
    rownames(mval_coord)
  ),
  identical(
    colnames(beta_coord),
    colnames(mval_coord)
  ),
  all(
    beta_coord >= 0 &
      beta_coord <= 1
  )
)

############################################################
## 20. SAVE FILTERED MATRICES
############################################################

beta_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_filtered_beta.rds"
)

m_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_filtered_M.rds"
)

annotation_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_filtered_annotation.rds"
)

summary_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_probe_filtering_summary.csv"
)

probe_qc_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_probe_QC_table.csv"
)

saveRDS(
  beta_coord,
  beta_output
)

saveRDS(
  mval_coord,
  m_output
)

saveRDS(
  annotation_coord,
  annotation_output
)

write.csv(
  filter_summary,
  summary_output,
  row.names = FALSE
)

write.csv(
  probe_qc,
  probe_qc_output,
  row.names = FALSE
)

############################################################
## 21. SAVE SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    proc_dir,
    "OCTANE_EPICv2_script03_sessionInfo.txt"
  )
)

############################################################
## 22. FINAL SUMMARY
############################################################

cat("\n============================================\n")
cat("SCRIPT 03 OCTANE EPIC v2 FILTERING COMPLETE\n")
cat("============================================\n")

cat(
  "Samples retained:",
  ncol(beta_coord),
  "\n"
)

cat(
  "Starting normalized probes:",
  nrow(beta_noob),
  "\n"
)

cat(
  "Final retained probes:",
  nrow(beta_coord),
  "\n"
)

cat(
  "Removed probes:",
  nrow(beta_noob) -
    nrow(beta_coord),
  "\n"
)

cat(
  "Retention rate:",
  round(
    100 *
      nrow(beta_coord) /
      nrow(beta_noob),
    2
  ),
  "%\n"
)

cat(
  "Sex-chromosome probes retained:",
  sum(
    annotation_coord$sex_chromosome
  ),
  "\n"
)

cat(
  "\nOutputs:\n"
)

cat(
  " ",
  beta_output,
  "\n"
)

cat(
  " ",
  m_output,
  "\n"
)

cat(
  " ",
  annotation_output,
  "\n"
)

cat(
  " ",
  summary_output,
  "\n"
)

cat(
  " ",
  probe_qc_output,
  "\n"
)

cat(
  "\nNext step:\n",
  "Review retained probe counts, then project filtered EPIC CpGs ",
  "to OCTANE plasma 300-bp regions.\n",
  sep = ""
)
