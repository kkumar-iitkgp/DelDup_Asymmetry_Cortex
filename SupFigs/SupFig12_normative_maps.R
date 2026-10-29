#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 12 (20 normative cortical maps used for spatial decoding)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> One brain map per normative map (left hemisphere, Glasser-180, lateral above medial); colour = rank percentile across the 180 regions
#          -> Maps grouped by domain as Fig. 3a; within domain ordered by rho(cDEC, map); -T1w/T2w shown sign-reversed (= H)
#          -> Statistics are read from the Supplementary Tables (final numbering); the data frame passed to each
#             ggplot call is written to the Source Data sheet(s) listed in SHEET_NAMES.
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: 1
#
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (+ data/inputs/inputs_*.csv, see data/README.md)
#
#-----
# Output:
#       1) SupFigs/output/SupFig12.png, SupFig12.svg (+ checks/SupFig12_checks.txt)
#       2) data/source_data_sheets/<sheet>.csv
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Kuldeep's local system; Else: run from the repository root
flag_hires = 0  # 0: PNG at 300 dpi only (GitHub); 1: PNG at 600 dpi + SVG (publication files, Zenodo)
if( flag_local == 1){
  project_dir <- "/Users/kuldeepkumar/Documents/Projects/Temp_CNV_GWAS_oct2022/Release_NG_R1_29oct2026/GitHub"
} else {
  project_dir <- getwd()
}

#-----------------------------------------------------------------------------
# Packages, style and helpers
#-----------------------------------------------------------------------------

source(file.path(project_dir, "functions", "00_figure_style.R"))
source(file.path(project_dir, "functions", "fn_write_SourceData_sheet.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(data.table); library(readxl)
  library(ggseg); library(ggsegGlasser); library(cowplot)
})

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

FIG_ID <- "SupFig12"
SHEET_NAMES <- c(main = "SupFig12")

FIG_HEIGHT_IN <- 8.40
MAP_OPTION <- "viridis"

inputs_dir <- file.path(project_dir, "data", "inputs"); source_data_dir <- file.path(project_dir, "data", "source_data_sheets")
OUT_DIR <- file.path(project_dir, "SupFigs", "output"); CHECKS_DIR <- file.path(OUT_DIR, "checks"); dir.create(CHECKS_DIR, showWarnings = FALSE, recursive = TRUE)
XL   <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
MAPS <- file.path(inputs_dir, "inputs_main_neuromaps20_Glasser180.csv")                  # 20 parcellated maps, joined by region name
CMP  <- file.path(inputs_dir, "inputs_SupFig12_13_composite_maps_Glasser180.csv")         # Fig. 3 composites (row-order guard)
RHO_CDEC <- file.path(inputs_dir, "inputs_Fig3_reference_cDEC_vs_maps.csv")               # rho(cDEC, map): not held in ST9 (29 Oct 2026)
for (f in c(XL, MAPS, CMP, RHO_CDEC)) if (!file.exists(f)) stop("Missing required input file: ", f, call. = FALSE)
chk_log <- character()
chk <- function(ok, msg) { s <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", msg); message(s); chk_log <<- c(chk_log, s); if (!isTRUE(ok)) stop(msg, call. = FALSE) }
pct <- function(v) 100 * (rank(v) - 1) / (length(v) - 1)   # rank percentile 0-100 (as Fig. 3c,d)

# ------------------------------- DATA -----------------------------------------
st4 <- as.data.table(read_excel(XL, "ST4"))[gene_set_exp == "All_data_analysis_Fig2"]
roi <- st4[cnv_type == "DEL", gene_set_name]
stopifnot(length(roi) == 180)
cmp <- fread(CMP)[match(roi, region)]
chk(identical(cmp$region, roi), "composite file joined to ST4 regions by name (180/180)")
maps20 <- fread(MAPS)[match(roi, gene_set_name)]; chk(identical(maps20$gene_set_name, roi), "map file joined to ST4 regions by name (180/180)")
maps20[, c("gene_set_name", "glasser_region") := NULL]
chk(nrow(maps20) == 180 && ncol(maps20) == 20 && !anyNA(maps20), "map file: 180 regions x 20 maps, no NA")
# guards: composites and H recomputed from the map file equal the Fig. 3 composite file
rz <- function(v) { r <- rank(v); (r - mean(r)) / sd(r) }
chk(max(abs((rz(maps20$cmr02) + rz(maps20$cmruglu)) / 2 - cmp$energy_composite),
        abs((rz(maps20$Evolution.Expansion) + rz(maps20$devexp)) / 2 - cmp$expansion_composite)) < 1e-9,
    "row order guard: composites recomputed from the map file == composite file (|d| < 1e-9)")
chk(abs(cor(-maps20$T1T2ratio, cmp$H, method = "spearman") - 1) < 1e-12, "row order guard: H = -T1w/T2w ranks identical to the composite file's H")

ST9 <- as.data.table(read_excel(XL, "ST9"))[map_set == "neuromaps_20"]
ST9 <- merge(ST9, fread(RHO_CDEC)[, .(source_column, rho_cDEC, pspin_cDEC)], by = "source_column")
chk(nrow(ST9) == 20 && !anyNA(ST9$rho_cDEC), "ST9 has 20 neuromaps; rho(cDEC, map) joined by source_column")

# ------------------------------- ORDER (as Fig. 3a) ---------------------------
groups <- tribble(
  ~source_column,        ~label,                    ~map_group,
  "T1T2ratio",           "−T1w/T2w",                "Microstructure",
  "PC1.AHBA",            "Gene expression PC1",     "Transcriptomic",
  "Cortical.Thickness",  "Cortical thickness",      "Microstructure",
  "G1.fMRI",             "Functional gradient",     "Functional connectivity",
  "intersubjvar",        "Functional variability",  "Functional connectivity",
  "PC1.Neurosynth",      "NeuroSynth PC1",          "Functional connectivity",
  "cmr02",               "Oxygen metabolism",       "Energy metabolism",
  "cmruglu",             "Glucose metabolism",      "Energy metabolism",
  "Evolution.Expansion", "Evolutionary expansion",  "Cortical expansion",
  "devexp",              "Developmental expansion", "Cortical expansion",
  "CBF",                 "Cerebral blood flow",     "Hemodynamics",
  "cbv",                 "Cerebral blood volume",   "Hemodynamics",
  "scalingnih",          "Allometric scaling NIH",  "Allometric scaling",
  "scalingpnc",          "Allometric scaling PNC",  "Allometric scaling",
  "meg_delta",           "MEG delta",               "Electrophysiology",
  "meg_theta",           "MEG theta",               "Electrophysiology",
  "meg_alpha",           "MEG alpha",               "Electrophysiology",
  "meg_beta",            "MEG beta",                "Electrophysiology",
  "meg_gamma1",          "MEG low gamma",           "Electrophysiology",
  "meg_gamma2",          "MEG high gamma",          "Electrophysiology")
GROUP_LEVELS <- c("Transcriptomic", "Microstructure", "Functional connectivity", "Energy metabolism", "Cortical expansion", "Hemodynamics", "Allometric scaling", "Electrophysiology")
FLIP <- c(T1T2ratio = -1)   # -T1w/T2w (= H)
mapv <- function(k) maps20[[k]] * (if (k %in% names(FLIP)) FLIP[[k]] else 1)
chk(setequal(groups$source_column, ST9$source_column) && setequal(groups$source_column, names(maps20)), "20 group assignments == ST9 maps == map-file columns (no NT maps)")
ST9 <- copy(ST9)[source_column %in% names(FLIP), rho_cDEC := -rho_cDEC]   # rho(cDEC, -T1w/T2w) = +0.368
rows <- ST9[groups, on = "source_column"] %>% as_tibble() %>%
  mutate(map_group = factor(map_group, levels = GROUP_LEVELS)) %>%
  group_by(map_group) %>% arrange(desc(rho_cDEC), .by_group = TRUE) %>% ungroup()
chk(abs(rows$rho_cDEC[rows$source_column == "T1T2ratio"] - 0.368) < 0.001, sprintf("rho(cDEC, -T1w/T2w) = %+.3f (= rho(cDEC, H))", rows$rho_cDEC[rows$source_column == "T1T2ratio"]))

# ------------------------------- BRAIN TILES ----------------------------------
glasser_region <- sub("^L_", "", gsub(".", "-", roi, fixed = TRUE))
chk(all(glasser_region %in% names(glasser$palette)), "all 180 regions found in the ggseg Glasser atlas (name-based)")
pc <- lapply(setNames(rows$source_column, rows$source_column), function(s) pct(mapv(s)))   # -T1w/T2w: rank percentile of H
one_view <- function(v, view, legend = FALSE)
  ggseg(tibble(region = glasser_region, value = v), atlas = "glasser", hemisphere = "left", view = view,
        mapping = aes(fill = value), colour = "grey92", size = 0.02) +
  theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
  scale_fill_viridis_c(option = MAP_OPTION, limits = c(0, 100), breaks = c(0, 50, 100), name = "Rank percentile", na.value = "grey85") +
  guides(fill = guide_colourbar(barwidth = unit(90, "pt"), barheight = unit(6, "pt"), title.position = "top", title.hjust = 0.5)) +
  theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = if (legend) "bottom" else "none",
        legend.title = element_text(size = FONT_LEGEND_TITLE_PT), legend.text = element_text(size = FONT_LEGEND_TEXT_PT),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(0, 0, 0, 0))
wrap_name <- function(s) if (nchar(s) <= 16) s else sub("^(.*\\S) (\\S+)$", "\\1\n\\2", s)   # long names on two lines (words unchanged)
tile <- function(s, lab)
  plot_grid(ggdraw() + draw_label(wrap_name(lab), size = FONT_AXIS_TEXT_PT, lineheight = 0.9, colour = "black"),
            one_view(pc[[s]], "lateral"), one_view(pc[[s]], "medial"), ncol = 1, rel_heights = c(0.34, 1, 1))
group_title <- function(txt, fill = NA) {   # plain black; tint behind the name only (Energy / Expansion)
  tg <- grid::textGrob(txt, gp = grid::gpar(fontsize = FONT_AXIS_TITLE_PT, fontfamily = FONT_FAMILY, fontface = "plain", col = "black"))
  if (is.na(fill)) return(ggdraw() + draw_grob(tg))
  rg <- grid::rectGrob(width = grid::grobWidth(tg) + grid::unit(12, "pt"), height = grid::grobHeight(tg) + grid::unit(7, "pt"), gp = grid::gpar(fill = fill, col = NA))
  ggdraw() + draw_grob(grid::grobTree(rg, tg))
}
TINT <- c("Energy metabolism" = COL_BAND_ENERGY, "Cortical expansion" = COL_BAND_EXPANSION)
block <- function(g) {
  r <- filter(rows, map_group == g)
  plot_grid(group_title(g, if (g %in% names(TINT)) TINT[[g]] else NA),
            plot_grid(plotlist = Map(tile, r$source_column, r$label), nrow = 1), ncol = 1, rel_heights = c(0.13, 1))
}
GAP <- 0.12
row1 <- plot_grid(block("Transcriptomic"), NULL, block("Microstructure"), NULL, block("Functional connectivity"), nrow = 1, rel_widths = c(1, GAP, 2, GAP, 3))
row2 <- plot_grid(block("Energy metabolism"), NULL, block("Cortical expansion"), NULL, block("Hemodynamics"), nrow = 1, rel_widths = c(2, GAP, 2, GAP, 2))
cbar <- get_legend(one_view(pc[[1]], "lateral", legend = TRUE))
row3 <- plot_grid(block("Allometric scaling"), NULL, plot_grid(NULL, cbar, NULL, ncol = 1, rel_heights = c(0.45, 0.3, 0.25)), nrow = 1, rel_widths = c(2, GAP, 4))
row4 <- block("Electrophysiology")
fig <- plot_grid(row1, NULL, row2, NULL, row3, NULL, row4, ncol = 1, rel_heights = c(1, 0.05, 1, 0.05, 1, 0.05, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_figure(fig, FIG_ID, FIG_WIDTH_IN, FIG_HEIGHT_IN, out_dir = OUT_DIR)

# ------------------------------- SOURCE DATA / CHECKS -------------------------
sd <- rbindlist(lapply(seq_len(nrow(rows)), function(i) { s <- rows$source_column[i]
  data.table(order = i, map_group = as.character(rows$map_group[i]), map_name = rows$label[i], source_column = s,
             region = roi, value = mapv(s), sign_reversed = s %in% names(FLIP), rank_percentile = pc[[s]]) }))
fn_write_SourceData_sheet(sd, SHEET_NAMES[["main"]], source_data_dir)
chk(nrow(sd) == 20 * 180 && all(sd[, .(r = range(rank_percentile)), by = source_column][, r] %in% c(0, 100)), "source data: 20 maps x 180 regions; each map spans rank percentile 0..100")
txt <- c(rows$label, GROUP_LEVELS, "Rank percentile")
chk(!any(grepl("×|β|↓|↑|→|pending|ρ|p-spin", txt)), "text audit: only map names, group names and the colour-bar title (no x-sign, beta, arrows, statistics)")
chk(!file.exists(file.path(OUT_DIR, paste0(FIG_ID, ".pdf"))), "no PDF written")
chk_log <- c(chk_log, sprintf("[note] ggseg: region 10pp is in the data (%s); the known ggseg Glasser merge leaves it undrawn", "10pp" %in% glasser_region))
chk(isTRUE(all.equal(pc[["T1T2ratio"]], pct(-maps20$T1T2ratio))), "-T1w/T2w tile = rank percentile of H = -T1w/T2w (high = association cortex)")
writeLines(chk_log, file.path(CHECKS_DIR, paste0(FIG_ID, "_checks.txt")))
message("Saved ", FIG_ID, " (svg/png) and source data.")
