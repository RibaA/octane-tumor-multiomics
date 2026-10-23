##------------------------------------------------------------
## 01_OCTANE_EPIC_input_validation
##
## Project:
##   C:/octane-plasma-buffy-coat
##
## Input:
##   data/epic/
##     210246760049/
##     210246760071/
##     210246760072/
##     MethylationEPIC_Sample_Sheet_Wang_260902.csv
##
## Output:
##   data/epic/proc/
##
## Purpose:
##   1. Read the Illumina EPIC v2 sample sheet
##   2. Find all paired Green/Red IDAT files
##   3. Match sample-sheet entries to physical IDAT files
##   4. Validate Sentrix IDs and positions
##   5. Load the raw RGChannelSet
##   6. Confirm EPIC v2 annotation/platform
##   7. Save validated targets and raw RGChannelSet
##
## Important:
##   - No sample exclusion
##   - No normalization
##   - No probe filtering
##------------------------------------------------------------

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

############################################################
## 1. LIBRARIES
############################################################

suppressPackageStartupMessages({
  library(minfi)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
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

sample_sheet_file <- file.path(
  epic_dir,
  "MethylationEPIC_Sample_Sheet_Wang_260902.csv"
)

dir.create(
  proc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  dir.exists(project_dir),
  dir.exists(epic_dir),
  file.exists(sample_sheet_file),
  dir.exists(proc_dir)
)

cat("\nProject directory:\n", project_dir, "\n")
cat("\nEPIC raw directory:\n", epic_dir, "\n")
cat("\nEPIC processed directory:\n", proc_dir, "\n")
cat("\nSample sheet:\n", sample_sheet_file, "\n")

############################################################
## 3. READ ILLUMINA SAMPLE SHEET
############################################################

read_illumina_sample_sheet <- function(file) {

  x <- readLines(
    file,
    warn = FALSE
  )

  ## Remove possible UTF-8 BOM
  x <- sub(
    "^\ufeff",
    "",
    x
  )

  ## Illumina sample sheets may have:
  ## [Data],,,,,,
  ## rather than simply [Data]
  data_line <- grep(
    "^\\s*\\[Data\\]",
    x
  )

  if (length(data_line) != 1L) {
    stop(
      "Could not uniquely identify the [Data] section. ",
      "Number of matches: ",
      length(data_line)
    )
  }

  ## Skip the [Data] line itself.
  ## The following line becomes the column header.
  dat <- read.csv(
    file,
    skip = data_line,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  ## Remove completely empty columns, which commonly result
  ## from trailing commas in Illumina sample sheets
  keep_columns <- vapply(
    dat,
    function(z) {
      any(
        !is.na(z) &
          trimws(as.character(z)) != ""
      )
    },
    logical(1)
  )

  dat <- dat[
    ,
    keep_columns,
    drop = FALSE
  ]

  dat
}

############################################################
## 4. LOAD SAMPLE SHEET
############################################################

sample_sheet <- read_illumina_sample_sheet(
  sample_sheet_file
)

cat("\n============================================\n")
cat("SAMPLE SHEET\n")
cat("============================================\n")

cat("Rows:", nrow(sample_sheet), "\n")

cat("Columns:\n")
print(colnames(sample_sheet))

print(head(sample_sheet))

############################################################
## 5. IDENTIFY REQUIRED COLUMNS
############################################################

normalize_column_name <- function(x) {

  tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      x
    )
  )
}

normalized_names <- normalize_column_name(
  colnames(sample_sheet)
)

find_column <- function(
    candidates,
    required = TRUE
) {

  candidate_norm <- normalize_column_name(
    candidates
  )

  idx <- match(
    candidate_norm,
    normalized_names
  )

  idx <- idx[
    !is.na(idx)
  ]

  if (length(idx) == 0L) {

    if (required) {

      stop(
        "Could not identify column. Tried: ",
        paste(
          candidates,
          collapse = ", "
        )
      )
    }

    return(NULL)
  }

  colnames(sample_sheet)[
    idx[[1]]
  ]
}

sentrix_id_col <- find_column(
  c(
    "Sentrix_ID",
    "Sentrix ID",
    "SentrixID"
  )
)

sentrix_position_col <- find_column(
  c(
    "Sentrix_Position",
    "Sentrix Position",
    "SentrixPosition"
  )
)

sample_name_col <- find_column(
  c(
    "Sample_Name",
    "Sample Name",
    "Sample_ID",
    "Sample ID",
    "SampleID"
  ),
  required = FALSE
)

cat("\nDetected sample-sheet columns:\n")

cat(
  "Sentrix ID:",
  sentrix_id_col,
  "\n"
)

cat(
  "Sentrix position:",
  sentrix_position_col,
  "\n"
)

cat(
  "Sample name:",
  ifelse(
    is.null(sample_name_col),
    "Not detected",
    sample_name_col
  ),
  "\n"
)

############################################################
## 6. STANDARDIZE SAMPLE SHEET
############################################################

targets_sheet <- sample_sheet %>%
  mutate(

    Sentrix_ID =
      trimws(
        as.character(
          .data[[sentrix_id_col]]
        )
      ),

    Sentrix_Position =
      toupper(
        trimws(
          as.character(
            .data[[sentrix_position_col]]
          )
        )
      )
  )

if (!is.null(sample_name_col)) {

  targets_sheet$Sample_Name <- trimws(
    as.character(
      targets_sheet[[sample_name_col]]
    )
  )

} else {

  targets_sheet$Sample_Name <- paste(
    targets_sheet$Sentrix_ID,
    targets_sheet$Sentrix_Position,
    sep = "_"
  )
}

targets_sheet$Array_Key <- paste(
  targets_sheet$Sentrix_ID,
  targets_sheet$Sentrix_Position,
  sep = "_"
)

############################################################
## 7. FIND IDAT FILES
############################################################

idat_files <- list.files(
  epic_dir,
  pattern = "_(Grn|Red)\\.idat$",
  full.names = TRUE,
  recursive = TRUE
)

cat("\n============================================\n")
cat("IDAT FILE DISCOVERY\n")
cat("============================================\n")

cat(
  "Total IDAT files:",
  length(idat_files),
  "\n"
)

if (length(idat_files) == 0L) {

  stop(
    "No IDAT files were found under data/epic."
  )
}

############################################################
## 8. PARSE IDAT INFORMATION
############################################################

idat_table <- tibble(
  full_path = idat_files,
  filename = basename(
    idat_files
  )
) %>%
  mutate(

    channel = case_when(

      str_detect(
        filename,
        "_Grn\\.idat$"
      ) ~ "Grn",

      str_detect(
        filename,
        "_Red\\.idat$"
      ) ~ "Red",

      TRUE ~ NA_character_
    ),

    Basename = str_remove(
      full_path,
      "_(Grn|Red)\\.idat$"
    ),

    file_stem = basename(
      Basename
    ),

    Sentrix_ID = str_extract(
      file_stem,
      "^[0-9]+"
    ),

    Sentrix_Position = str_extract(
      file_stem,
      "R[0-9]{2}C[0-9]{2}$"
    ),

    Array_Key = paste(
      Sentrix_ID,
      Sentrix_Position,
      sep = "_"
    )
  )

print(
  head(idat_table)
)

############################################################
## 9. CHECK GREEN/RED PAIRS
############################################################

idat_pair_summary <- idat_table %>%
  count(
    Array_Key,
    Sentrix_ID,
    Sentrix_Position,
    channel
  ) %>%
  pivot_wider(
    names_from = channel,
    values_from = n,
    values_fill = 0
  )

cat("\nIDAT pair summary:\n")

print(
  idat_pair_summary
)

bad_pairs <- idat_pair_summary %>%
  filter(
    Grn != 1L |
      Red != 1L
  )

if (nrow(bad_pairs) > 0L) {

  print(
    bad_pairs
  )

  stop(
    "Some array positions do not have exactly one ",
    "Green and one Red IDAT file."
  )
}

cat(
  "\nComplete Green/Red pairs:",
  nrow(idat_pair_summary),
  "\n"
)

############################################################
## 10. CREATE ONE RECORD PER ARRAY POSITION
############################################################

idat_targets <- idat_table %>%
  distinct(
    Array_Key,
    Sentrix_ID,
    Sentrix_Position,
    Basename
  )

stopifnot(
  !anyDuplicated(
    idat_targets$Array_Key
  ),
  !anyDuplicated(
    idat_targets$Basename
  )
)

############################################################
## 11. MATCH SAMPLE SHEET TO IDAT FILES
############################################################

targets <- targets_sheet %>%
  left_join(
    idat_targets %>%
      select(
        Array_Key,
        Basename
      ),
    by = "Array_Key"
  )

############################################################
## 12. VALIDATE MATCHING
############################################################

missing_idat <- targets %>%
  filter(
    is.na(
      Basename
    )
  )

if (nrow(missing_idat) > 0L) {

  cat(
    "\nSamples in sample sheet without matching IDATs:\n"
  )

  print(
    missing_idat %>%
      select(
        Sample_Name,
        Sentrix_ID,
        Sentrix_Position
      )
  )

  stop(
    "Some sample-sheet entries could not be matched to IDAT files."
  )
}

idat_not_in_sheet <- setdiff(
  idat_targets$Array_Key,
  targets$Array_Key
)

if (length(idat_not_in_sheet) > 0L) {

  cat(
    "\nIDAT positions not represented in sample sheet:\n"
  )

  print(
    idat_not_in_sheet
  )

  stop(
    "Some IDAT positions are not represented in the sample sheet."
  )
}

############################################################
## 13. DUPLICATE CHECKS
############################################################

stopifnot(
  !anyDuplicated(
    targets$Array_Key
  ),
  !anyDuplicated(
    targets$Basename
  )
)

if (anyDuplicated(
  targets$Sample_Name
)) {

  warning(
    "Duplicated Sample_Name values were found."
  )
}

############################################################
## 14. EXPECTED SAMPLE COUNT
############################################################

expected_n_samples <- 23L

cat("\n============================================\n")
cat("MATCHING SUMMARY\n")
cat("============================================\n")

cat(
  "Sample-sheet rows:",
  nrow(targets_sheet),
  "\n"
)

cat(
  "Complete IDAT pairs:",
  nrow(idat_targets),
  "\n"
)

cat(
  "Successfully matched samples:",
  nrow(targets),
  "\n"
)

if (nrow(targets) != expected_n_samples) {

  warning(
    "Expected ",
    expected_n_samples,
    " samples, but found ",
    nrow(targets),
    "."
  )
}

############################################################
## 15. VERIFY PHYSICAL FILES
############################################################

stopifnot(
  all(
    file.exists(
      paste0(
        targets$Basename,
        "_Grn.idat"
      )
    )
  )
)

stopifnot(
  all(
    file.exists(
      paste0(
        targets$Basename,
        "_Red.idat"
      )
    )
  )
)

############################################################
## 16. LOAD RAW RGCHANNELSET
############################################################

cat("\n============================================\n")
cat("READING RAW EPIC v2 DATA\n")
cat("============================================\n")

rgSet_OCTANE_EPICv2 <- read.metharray.exp(
  targets = targets,
  extended = TRUE,
  verbose = TRUE
)

############################################################
## 17. SET SAMPLE NAMES
############################################################

sampleNames(
  rgSet_OCTANE_EPICv2
) <- targets$Sample_Name

stopifnot(
  identical(
    sampleNames(
      rgSet_OCTANE_EPICv2
    ),
    targets$Sample_Name
  )
)

############################################################
## 18. REVIEW ARRAY ANNOTATION
############################################################

cat("\n============================================\n")
cat("RGCHANNELSET SUMMARY\n")
cat("============================================\n")

cat(
  "RGChannelSet dimensions:",
  paste(
    dim(
      rgSet_OCTANE_EPICv2
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Number of samples:",
  ncol(
    rgSet_OCTANE_EPICv2
  ),
  "\n"
)

cat(
  "\nminfi annotation:\n"
)

print(
  annotation(
    rgSet_OCTANE_EPICv2
  )
)

############################################################
## 19. SAMPLES BY SENTRIX ID
############################################################

sentrix_summary <- targets %>%
  count(
    Sentrix_ID,
    name = "n_samples"
  ) %>%
  arrange(
    Sentrix_ID
  )

cat(
  "\nSamples by Sentrix ID:\n"
)

print(
  sentrix_summary
)

############################################################
## 20. SAVE CURATED INPUT OBJECTS
############################################################

targets_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_targets_validated.csv"
)

write.csv(
  targets,
  targets_output,
  row.names = FALSE
)

rgset_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_raw_RGChannelSet.rds"
)

saveRDS(
  rgSet_OCTANE_EPICv2,
  rgset_output
)

############################################################
## 21. SAVE SAMPLE-SHEET COPY
############################################################

sample_sheet_output <- file.path(
  proc_dir,
  "OCTANE_EPICv2_sample_sheet_curated.csv"
)

write.csv(
  targets_sheet,
  sample_sheet_output,
  row.names = FALSE
)

############################################################
## 22. FINAL SUMMARY
############################################################

cat("\n============================================\n")
cat("SCRIPT 01 OCTANE EPIC v2 INPUT VALIDATION COMPLETE\n")
cat("============================================\n")

cat(
  "Sample-sheet samples:",
  nrow(
    targets_sheet
  ),
  "\n"
)

cat(
  "Complete IDAT pairs:",
  nrow(
    idat_targets
  ),
  "\n"
)

cat(
  "Samples loaded into RGChannelSet:",
  ncol(
    rgSet_OCTANE_EPICv2
  ),
  "\n"
)

cat(
  "Array annotation:",
  paste(
    annotation(
      rgSet_OCTANE_EPICv2
    ),
    collapse = " / "
  ),
  "\n"
)

cat(
  "\nOutputs:\n"
)

cat(
  " ",
  targets_output,
  "\n"
)

cat(
  " ",
  rgset_output,
  "\n"
)

cat(
  " ",
  sample_sheet_output,
  "\n"
)

cat(
  "\nNo normalization or filtering was performed.\n"
)
