#-------------------------------------------------------------------------------
# Summary: helper to write one Source Data sheet (the exact data frame passed to ggplot) as CSV
#-------------------------------------------------------------------------------
# Used by every figure script. source_data_build/03_assemble_SourceData.R (Zenodo archive) collects the CSVs into
# SourceData_R1.xlsx (one sheet per CSV, sheet name = file name).
# Author: Kuldeep Kumar; 
# Date: 21 October, 2026
#-------------------------------------------------------------------------------

fn_write_SourceData_sheet <- function(df_input, sheet_name_C, out_dir_C) {
  # writes the exact data frame passed to ggplot as <sheet_name_C>.csv (full precision)
  if (nchar(sheet_name_C) > 31) stop("sheet name longer than 31 characters: ", sheet_name_C)
  dir.create(out_dir_C, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(as.data.frame(df_input), file.path(out_dir_C, paste0(sheet_name_C, ".csv")), row.names = FALSE, na = "NA")
  invisible(file.path(out_dir_C, paste0(sheet_name_C, ".csv")))
}
