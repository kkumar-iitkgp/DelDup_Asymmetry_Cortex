#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 13 (energy metabolism and cortical expansion: unadjusted associations, single maps and replication)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a: unadjusted version of Fig. 3c,d: composite rank-percentile maps; Del and Dup effect size vs composite rank percentile
#          -> b: partial rho
#          -> H for the 19 normative maps other than T1w/T2w (Del / Dup)
#          -> c: replication of Del: Energy
#          -> H, Dup: Expansion
#          -> H and the interaction I_H in UK Biobank and 8 internal strata
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
#       1) SupFigs/output/SupFig13.png, SupFig13.svg (+ checks/SupFig13_checks.txt)
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
  library(cowplot); library(patchwork); library(ggseg); library(ggsegGlasser)
})

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

FIG_ID <- "SupFig13"
SHEET_NAMES <- c(a = "SupFig13a", a_stats = "SupFig13a_stats", b = "SupFig13b", c = "SupFig13c")

FIG_HEIGHT_IN <- 11.0
inputs_dir <- file.path(project_dir, "data", "inputs"); source_data_dir <- file.path(project_dir, "data", "source_data_sheets")
OUT_DIR <- file.path(project_dir, "SupFigs", "output"); CHECKS_DIR <- file.path(OUT_DIR, "checks"); dir.create(CHECKS_DIR, showWarnings = FALSE, recursive = TRUE)
XL <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
COL_I <- COL_DELTA_BAR
chk_log <- character()
chk <- function(ok, msg) { s <- sprintf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", msg); message(s); chk_log <<- c(chk_log, s); if (!isTRUE(ok)) stop(msg, call. = FALSE) }

# ------------------------------- DATA -----------------------------------------
st4 <- as.data.table(read_excel(XL, "ST4"))[gene_set_exp == "All_data_analysis_Fig2"]
roi <- st4[cnv_type == "DEL", gene_set_name]
es <- list(DEL = st4[cnv_type == "DEL"][match(roi, gene_set_name), Estimate], DUP = st4[cnv_type == "DUP"][match(roi, gene_set_name), Estimate])
maps20 <- fread(file.path(inputs_dir, "inputs_main_neuromaps20_Glasser180.csv"))[match(roi, gene_set_name)]; stopifnot(identical(maps20$gene_set_name, roi))
maps20[, c("gene_set_name", "glasser_region") := NULL]; stopifnot(nrow(maps20) == 180, !anyNA(maps20))
H <- -maps20$T1T2ratio
rz <- function(v) { r <- rank(v); (r - mean(r)) / sd(r) }
comp <- list(Energy = (rz(maps20$cmr02) + rz(maps20$cmruglu)) / 2, Expansion = (rz(maps20$Evolution.Expansion) + rz(maps20$devexp)) / 2)
cmpf <- fread(file.path(inputs_dir, "inputs_SupFig12_13_composite_maps_Glasser180.csv"))
chk(identical(cmpf$region, roi) && max(abs(cmpf$energy_composite - comp$Energy), abs(cmpf$expansion_composite - comp$Expansion)) < 1e-9, "composites recomputed == composite CSV; region order == ST4")
SD3 <- fread(file.path(inputs_dir, "inputs_SupFig13_Fig3_dimension_statistics.csv")); ADJ <- fread(file.path(inputs_dir, "inputs_SupFig13_partial_rho_adjH_neuromaps20.csv"))   # not held in the Supplementary Tables (29 Oct 2026)
mainD <- SD3[analysis == "main_composite" & sample == "All_data_analysis_Fig2"]
intI  <- SD3[analysis == "interaction_I" & sample == "All_data_analysis_Fig2"]

# ------------------------------- a: unadjusted (former main c/d) --------------
pct <- function(v) 100 * (rank(v) - 1) / (length(v) - 1)
pc_list <- lapply(comp, pct)
glasser_region <- sub("^L_", "", gsub(".", "-", roi, fixed = TRUE))
chk(all(glasser_region %in% names(glasser$palette)), "a: all 180 regions found in the ggseg Glasser atlas (name-based)")
brain_comp <- function(v) {
  one <- function(view, legend) ggseg(tibble(region = glasser_region, value = v), atlas = "glasser", hemisphere = "left", view = view,
        mapping = aes(fill = value), colour = "grey92", size = 0.02) +
    theme_brain(text.family = FONT_FAMILY, text.size = FONT_BRAIN_PT) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 100), breaks = c(0, 100), name = "Rank percentile", na.value = "grey85") +
    guides(fill = guide_colourbar(barwidth = unit(48, "pt"), barheight = unit(4, "pt"), title.position = "top", title.hjust = 0.5)) +
    theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = if (legend) "bottom" else "none",
          legend.title = element_text(size = FONT_LEGEND_TEXT_PT), legend.text = element_text(size = FONT_LEGEND_TEXT_PT),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(0, 0, 0, 0))
  plot_grid(plot_grid(one("lateral", FALSE), one("medial", FALSE), nrow = 1), get_legend(one("lateral", TRUE)), NULL, ncol = 1, rel_heights = c(1, 0.62, 0.12))
}
ylims <- lapply(es, function(v) { r <- range(v); d <- diff(r); c(r[1] - 0.12 * d, r[2] + 0.36 * d) })   # shared within CNV type; strip + headroom
uline <- function(r, p) c(sprintf("ρ = %.2f", r), paste0("p-spin = ", format_pspin(p)))
printed <- list(); scat <- list()
raw_scatter <- function(d, t) {
  m <- mainD[dimension == d]; r <- m[[paste0("rho_", tolower(t))]]; p <- m[[paste0("pspin_", tolower(t))]]
  l <- uline(r, p); printed[[paste(d, t)]] <<- l; lim <- ylims[[t]]; span <- diff(lim)
  g <- ggplot(tibble(x = pc_list[[d]], y = es[[t]]), aes(x, y)) +
    geom_tile(data = tibble(x = seq(0, 100, length.out = 101)), aes(x = x, y = lim[1] + 0.035 * span, fill = x), height = 0.045 * span, width = 1.01, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 100), guide = "none") +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_point(colour = COL_POINT, alpha = 0.7, size = 1.0) +
    annotate("text", x = 3, y = lim[2], label = paste(l, collapse = "\n"), hjust = 0, vjust = 1, lineheight = 0.95, size = pt_to_geom(FONT_ANNOTATION_PT)) +
    scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(add = 2)) +
    coord_cartesian(ylim = lim, expand = FALSE, clip = "off") +
    labs(x = NULL, y = paste("Effect size", if (t == "DEL") "Del" else "Dup")) + theme_ng() + theme(plot.margin = margin(2, 11, 2, 2))
  scat[[paste(d, t)]] <<- g; g
}
title_strip <- function(txt, fill) ggdraw() + draw_grob(grid::rectGrob(gp = grid::gpar(fill = fill, col = NA))) +
  draw_label(txt, colour = "black", size = FONT_AXIS_TITLE_PT, x = 0.07, hjust = 0)
dim_block <- function(d, txt, fill) {
  map <- plot_grid(NULL, brain_comp(pc_list[[d]]), NULL, nrow = 1, rel_widths = c(0.2, 0.6, 0.2))
  sc <- plot_grid(raw_scatter(d, "DEL"), raw_scatter(d, "DUP"), nrow = 1, align = "h", axis = "tb")
  plot_grid(title_strip(txt, fill), map, sc, ggdraw() + draw_label("Composite rank percentile", size = FONT_AXIS_TITLE_PT), ncol = 1, rel_heights = c(0.07, 0.50, 1, 0.07))
}
panel_a <- plot_grid(dim_block("Energy", "Energy metabolism", COL_BAND_ENERGY), dim_block("Expansion", "Cortical expansion", COL_BAND_EXPANSION),
                     nrow = 1, labels = "a", label_size = FONT_PANEL_TAG_PT, label_fontface = "bold", hjust = -0.1, vjust = 1.1)
SD_a <- data.table(region = roi, energy_composite_percentile = pc_list$Energy, expansion_composite_percentile = pc_list$Expansion,
                   effect_size_DEL = es$DEL, effect_size_DUP = es$DUP)
SD_a_stats <- rbindlist(lapply(c("Energy", "Expansion"), function(d) rbindlist(lapply(c("DEL", "DUP"), function(t) { m <- mainD[dimension == d]
  data.table(dimension = d, cnv_type = t, rho = m[[paste0("rho_", tolower(t))]], p_spin = m[[paste0("pspin_", tolower(t))]]) }))))
for (d in c("Energy", "Expansion")) for (t in c("DEL", "DUP")) { s <- SD_a_stats[dimension == d & cnv_type == t]
  chk(abs(cor(es[[t]], pc_list[[d]], method = "spearman") - s$rho) < 1e-12 && identical(printed[[paste(d, t)]], uline(s$rho, s$p_spin)),
      sprintf("a: %s-%s Spearman from plotted data == table; printed '%s' / '%s'", d, t, printed[[paste(d, t)]][1], printed[[paste(d, t)]][2])) }
lay <- unlist(lapply(scat, function(g) sapply(g$layers, function(l) class(l$geom)[1])))
chk(!any(lay == "GeomSmooth"), "a: no fit line")

# ------------------------------- b: 19 maps adjusted for H --------------------
groups <- tribble(~source_column, ~label, ~map_group,
  "PC1.AHBA", "Gene expression PC1", "Transcriptomic", "Cortical.Thickness", "Cortical thickness", "Microstructure",
  "G1.fMRI", "Functional gradient", "Functional connectivity", "intersubjvar", "Functional variability", "Functional connectivity",
  "PC1.Neurosynth", "NeuroSynth PC1", "Functional connectivity",
  "cmr02", "Oxygen metabolism", "Energy metabolism", "cmruglu", "Glucose metabolism", "Energy metabolism",
  "Evolution.Expansion", "Evolutionary expansion", "Cortical expansion", "devexp", "Developmental expansion", "Cortical expansion",
  "CBF", "Cerebral blood flow", "Hemodynamics", "cbv", "Cerebral blood volume", "Hemodynamics",
  "scalingnih", "Allometric scaling NIH", "Allometric scaling", "scalingpnc", "Allometric scaling PNC", "Allometric scaling",
  "meg_delta", "MEG delta", "Electrophysiology", "meg_theta", "MEG theta", "Electrophysiology", "meg_alpha", "MEG alpha", "Electrophysiology",
  "meg_beta", "MEG beta", "Electrophysiology", "meg_gamma1", "MEG low gamma", "Electrophysiology", "meg_gamma2", "MEG high gamma", "Electrophysiology")
GL <- c("Transcriptomic", "Microstructure", "Functional connectivity", "Energy metabolism", "Cortical expansion", "Hemodynamics", "Allometric scaling", "Electrophysiology")
GS <- c("Transcriptomic", "Microstructure", "Functional\nconnectivity", "Energy\nmetabolism", "Cortical\nexpansion", "Hemo-\ndynamics", "Allometric\nscaling", "Electro-\nphysiology")
ST9c <- fread(file.path(inputs_dir, "inputs_Fig3_reference_cDEC_vs_maps.csv"))[, .(source_column, rho_cDEC)]   # rho(cDEC, map) (Fig. 3a order)
adj19 <- ADJ[source_column != "T1T2ratio"]
chk(nrow(adj19) == 19 && setequal(adj19$source_column, groups$source_column) && !anyNA(adj19$rho_del_adjH), "b: 19 maps (T1w/T2w dropped), all with partial rho; no NT maps")
bd <- as_tibble(adj19[groups, on = "source_column"][ST9c, on = "source_column", nomatch = NULL]) %>%
  mutate(map_group = factor(map_group, levels = GL, labels = GS)) %>% group_by(map_group) %>% arrange(desc(rho_cDEC), .by_group = TRUE) %>% ungroup() %>%
  mutate(lab = factor(label, levels = rev(label)))   # Fig. 3a,b row order
lb <- bind_rows(bd %>% transmute(map_group, lab, CNV = "Del", rho = rho_del_adjH, p = pspin_del_adjH),
                bd %>% transmute(map_group, lab, CNV = "Dup", rho = rho_dup_adjH, p = pspin_dup_adjH))
band <- tibble(map_group = factor(GS[4:5], levels = GS), fill = c(COL_BAND_ENERGY, COL_BAND_EXPANSION))
panel_b <- ggplot() +
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = I(fill))) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_segment(data = bd, aes(x = rho_del_adjH, xend = rho_dup_adjH, y = lab, yend = lab), colour = COL_LINE, linewidth = 0.6) +
  geom_point(data = lb, aes(rho, lab, colour = CNV), size = 2) +
  geom_point(data = filter(lb, p < 0.05), aes(rho, lab), shape = 21, size = 2, stroke = 0.7, colour = "black") +
  scale_colour_manual(values = c(Del = COL_DEL, Dup = COL_DUP), name = NULL) +
  facet_grid(map_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) + coord_cartesian(xlim = c(-0.6, 0.6)) +
  labs(x = "Partial ρ | H (Del / Dup)", y = NULL, tag = "b") + theme_ng(grid_x = TRUE) +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, margin = margin(r = 8)),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.spacing.y = unit(3, "pt"),
        legend.position = "inside", legend.position.inside = c(0.9, 0.06), legend.key.size = unit(9, "pt"), legend.background = element_rect(fill = "white", colour = NA))
SD_b <- lb %>% transmute(map_group = GL[as.integer(map_group)], map_name = as.character(lab), cnv_type = CNV, rho_Hadj = rho, p_spin_Hadj = p)
chk(isTRUE(all.equal(sort(lb$rho), sort(c(adj19$rho_del_adjH, adj19$rho_dup_adjH)))), "b: plotted partial rho == partial rho input (38 values)")

# ------------------------------- c: replication -------------------------------
smp <- c(UKBB = "UKBB", Discovery = "Discovery", Replication = "Replication", Sex_Female = "Female", Sex_Male = "Male",
         Age_lessthan_60 = "Age < 60", Age_greaterthan_60 = "Age > 60", Ancestry_EUR = "EUR ancestry", UKBB_FI_center_visit = "UKBB FI centre visit")
sh <- unique(SD3[grepl("^replication_", analysis) & grepl("split", sample, ignore.case = TRUE), sample])   # split-half rows, when present
if (length(sh)) smp <- c(smp, setNames(sh, sh))
rp <- SD3[sample %in% names(smp) & grepl("^replication_", analysis)]
chk(all(setdiff(names(smp), sh) %in% rp$sample) && !("Other_5cohort" %in% rp$sample), "c: UKBB + 8 strata present; Other 5 cohorts not plotted")
Q <- c("Del: Energy | H", "Dup: Expansion | H", "Interaction")
dl <- rbind(rp[dimension == "Energy", .(sample, quantity = Q[1], value = rho_del_adjH, p = pspin_del_adjH)],
            rp[dimension == "Expansion", .(sample, quantity = Q[2], value = rho_dup_adjH, p = pspin_dup_adjH)],
            rp[grepl("interaction", analysis), .(sample, quantity = Q[3], value = delta_adjH, p = pspin_delta_adjH)])
full <- data.table(quantity = factor(Q, levels = Q), value = c(mainD[dimension == "Energy", rho_del_adjH], mainD[dimension == "Expansion", rho_dup_adjH], intI$delta_adjH))
dl[, `:=`(lab = factor(smp[sample], levels = rev(smp)), quantity = factor(quantity, levels = Q))]
panel_c <- ggplot(dl, aes(value, lab)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_vline(data = full, aes(xintercept = value), colour = "grey30", linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(colour = quantity), size = 2) +
  geom_point(data = dl[p < 0.05], shape = 21, size = 2, stroke = 0.7, colour = "black") +
  scale_colour_manual(values = setNames(c(COL_DEL, COL_DUP, COL_I), Q), guide = "none") +
  facet_grid(. ~ quantity, scales = "free_x") + scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  labs(x = "Partial ρ | H (Interaction: I_H)", y = NULL, tag = "c") + theme_ng(grid_x = TRUE) +
  theme(strip.background = element_blank(), panel.spacing.x = unit(14, "pt"))
SD_c <- dl[, .(sample, quantity = as.character(quantity), value, p_spin = p)]
SD_c <- rbind(SD_c, data.table(sample = "All_data_analysis_Fig2 (dashed line)", quantity = Q, value = full$value,
                               p_spin = c(mainD[dimension == "Energy", pspin_del_adjH], mainD[dimension == "Expansion", pspin_dup_adjH], intI$pspin_delta_adjH)))
chk(abs(dl[sample == "UKBB" & quantity == Q[3], value] - SD3[sample == "UKBB" & analysis == "replication_interaction_I", delta_adjH]) < 1e-12 &&
    abs(full$value[3] - intI$delta_adjH) < 1e-12, sprintf("c: UKBB I_H %.3f and full-sample I_H %.3f == ST_Fig3_dimensions", dl[sample == "UKBB" & quantity == Q[3], value], full$value[3]))

# ------------------------------- COMBINE / SAVE -------------------------------
fig <- plot_grid(panel_a, ggplotGrob(panel_b), ggplotGrob(panel_c), ncol = 1, rel_heights = c(1.15, 1.25, 0.82))
save_figure(fig, FIG_ID, FIG_WIDTH_IN, FIG_HEIGHT_IN, out_dir = OUT_DIR)
chk(!file.exists(file.path(OUT_DIR, paste0(FIG_ID, ".pdf"))), "no PDF written")
fn_write_SourceData_sheet(SD_a, SHEET_NAMES[["a"]], source_data_dir); fn_write_SourceData_sheet(SD_a_stats, SHEET_NAMES[["a_stats"]], source_data_dir)
fn_write_SourceData_sheet(SD_b, SHEET_NAMES[["b"]], source_data_dir); fn_write_SourceData_sheet(SD_c, SHEET_NAMES[["c"]], source_data_dir)
txt <- c("Energy metabolism", "Cortical expansion", unlist(printed), "Effect size Del", "Effect size Dup", "Composite rank percentile", "Rank percentile",
         levels(bd$lab), GS, Q, smp, "Partial ρ | H (Del / Dup)", "Partial ρ | H (Interaction: I_H)")
chk(!any(grepl("×|β|↓|pending|not applicable|Independent datasets|Internal stability|Other 5|T1w/T2w|Dashed line|H-adj|adj\\. for H|fit", txt)),
    "text audit: no x-sign, beta, arrow, 'pending', 'not applicable', dataset headings, Other 5 cohorts, T1w/T2w row, caption, fit notes")
chk_log <- c(chk_log, sprintf("[legend] Interaction I_H (full sample) = %.2f, p-spin = %s; I_H goes to supplementary source data c", intI$delta_adjH, format_pspin(intI$pspin_delta_adjH)))
writeLines(chk_log, file.path(CHECKS_DIR, paste0(FIG_ID, "_checks.txt")))
message("Saved ", FIG_ID)
