Fig3
================

## Fig. 3: Spatial decoding of deletion and duplication effects and their association with energy metabolism and cortical expansion
#### -- Figure legend -- ####

a) Spatial correlation between the cortical dosage-effect contrast (cDEC) and 20 normative cortical maps. The maps are grouped by domain: hierarchy and microstructure, energy metabolism, cortical expansion, hemodynamics, allometric scaling, and magnetoencephalography (MEG). Bars show Spearman's ρ. Asterisks mark two-sided spin-permutation p < 0.05. Positive ρ indicates that the map is high where duplication effects are relatively stronger than deletion effects. T1w/T2w is the inverse of the cortical hierarchy (H) used in Fig. 2. b) Spatial Spearman correlation of each map with the deletion (red) and duplication (blue) effect-size profiles. The grey line joining the two points shows their difference, the topographic contrast, T = ρ(Dup,map) - ρ(Del,map). Rows are in the same order as in a. Because effect sizes are negative, a negative ρ indicates larger effects where the map is high. A dark outline marks two-sided spin-permutation p < 0.05. Shaded rows mark the two dimensions analysed in c and d. c,d) Deletion (red) and duplication (blue) effects versus energy metabolism (c; oxygen and glucose metabolism) and cortical expansion (d; evolutionary and developmental expansion), adjusted for H. Partial Spearman ρ and two-sided p-spin (composite and H rotated jointly). Unadjusted results in Supplementary Fig. 13. e) Mean deletion versus duplication effect size across the 17 Yeo resting-state networks (± s.e.m.); colours match the inset network map. Spearman ρ and two-sided p-spin.

Spin-permutation tests used 10,000 null rotations and were two-sided. Source data are provided with this paper. Del, deletion; Dup, duplication; cDEC, cortical dosage-effect contrast; p-spin, spin-permutation p-value; ρ, Spearman correlation.
The code below is read from `Fig3.R` (same folder). Run it from the repository root, or knit this file. Statistics are read from the Supplementary Tables (`data/inputs/`); the plotted values are written to `data/source_data_sheets/` (Source Data).

#### Setup: paths, packages and style ####

``` r
#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Kuldeep's local system; Else: run from the repository root
if( flag_local == 1){
  project_dir <- "/Users/kuldeepkumar/Documents/Projects/Temp_CNV_GWAS_oct2022/Release_NG_R1_29oct2026/GitHub"
} else {
  project_dir <- getwd()
}

data_dir <- paste0(project_dir,"/data")
inputs_dir <- paste0(data_dir,"/inputs")
source_data_dir <- paste0(data_dir,"/source_data_sheets")
OUT <- paste0(project_dir,"/Figs/Fig3")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(data.table); library(readxl); library(tools)
  library(ggseg); library(ggsegGlasser); library(ggsegYeo2011); library(cowplot); library(patchwork); library(ggrepel); library(svglite)
})
source(file.path(project_dir, "functions/fn_write_SourceData_sheet.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters, style and input files
#-----------------------------------------------------------------------------

#----- style (from _style_R1.R, 19 Oct 2026; text sizes of Fig. 3 v3)
MM <- 1 / 25.4                                   # mm -> inch
FONT <- "Arial"
PT_TEXT <- 9; PT_SMALL <- 8; PT_TITLE <- 9.5; PT_TAG <- 11
LW_MIN <- 0.25 * 0.3528 / 0.75                   # 0.25 pt in ggplot linewidth units
pt2mm <- function(pt) pt / .pt                   # text size for geoms
COL_DEL <- "#B7222E"; COL_DUP <- "#2878B8"; COL_CDEC <- "#8E4C97"
COL_BAND_ENERGY <- "#E9F3E8"; COL_BAND_EXPANSION <- "#FDF0E3"; COL_POINT <- "grey60"; COL_LINE <- "grey70"
COL_TITLE_ENERGY <- "#CFE6CC"; COL_TITLE_EXPANSION <- "#FADDBE"
P_LABEL <- "p-spin"
SEED <- 20261019
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 1e-3) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
theme_R1 <- function(grid_x = FALSE) {
  theme_classic(base_size = PT_TEXT, base_family = FONT) +
    theme(axis.line = element_line(linewidth = LW_MIN * 1.6, colour = "black"), axis.ticks = element_line(linewidth = LW_MIN * 1.6, colour = "black"),
          axis.ticks.length = unit(1.2, "pt"), axis.text = element_text(size = PT_SMALL + 0.5, colour = "black"), axis.title = element_text(size = PT_TEXT),
          legend.text = element_text(size = PT_SMALL + 0.5), legend.title = element_text(size = PT_TEXT), legend.key.size = unit(6, "pt"),
          strip.text = element_text(size = PT_TEXT), strip.background = element_blank(), plot.margin = margin(2, 2, 2, 2),
          panel.grid.major.x = if (grid_x) element_line(colour = "grey88", linewidth = LW_MIN, linetype = "dashed") else element_blank())
}
pv2 <- function(obs, nul) (1 + sum(abs(nul) >= abs(obs))) / (length(nul) + 1)   # two-sided spin p
set.seed(SEED)
STEM <- "Fig3"

filename_ST_C <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
filename_maps_C <- file.path(inputs_dir, "inputs_main_neuromaps20_Glasser180.csv")
filename_rot_C <- file.path(inputs_dir, "inputs_main_spin_rotations_Glasser180.rds")
LOG <- character(); logf <- function(...) { s <- sprintf(...); message(s); LOG <<- c(LOG, s) }
chk <- function(ok, msg) { logf("[%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", msg); if (!isTRUE(ok)) stop(msg, call. = FALSE) }

#-----------------------------------------------------------------------------
# PART B: Load data
#-----------------------------------------------------------------------------
```

#### Data: effect sizes (Supplementary Table 4), 20 normative maps, spin rotations ####

``` r
st4 <- as.data.table(read_excel(filename_ST_C, "ST4"))[gene_set_exp == "All_data_analysis_Fig2"]
roi <- st4[cnv_type == "DEL", gene_set_name]
chk(identical(sub("^L_", "", gsub(".", "-", roi, fixed = TRUE)), names(glasser$palette)[1:180]), "ST4 region order = Glasser index order 1-180 (used by the fsaverage5 labels in panel e)")
DEL <- st4[cnv_type == "DEL"][match(roi, gene_set_name), Estimate]; DUP <- st4[cnv_type == "DUP"][match(roi, gene_set_name), Estimate]
maps <- fread(filename_maps_C); maps <- maps[match(roi, maps$gene_set_name)]; chk(!anyNA(maps$gene_set_name), "20 maps joined by region name")
H <- -maps$T1T2ratio
z <- function(v) as.numeric(scale(v))
chk(abs(cor(z(DEL) - z(DUP), H, method = "spearman") - 0.368) < 0.001, "guard: rho(cDEC, H) = +0.368 (inputs as Fig. 2)")
P <- readRDS(filename_rot_C); NR <- ncol(P)
chk(NR == 9999 && startsWith(unname(md5sum(filename_rot_C)), "a3ac277d"), "cleaned Glasser rotations: 10,000 generated, 9,999 retained (md5 a3ac277d...)")
ST8 <- as.data.table(read_excel(filename_ST_C, "ST9"))[map_set == "neuromaps_20"]   # Supplementary Table 9 (20 normative maps)

#-----------------------------------------------------------------------------
# PART C: Panels
#-----------------------------------------------------------------------------
```

#### Panels a and b ####

``` r
#----- Panels a, b: 20 normative maps (Fig. 3a,b)
groups <- tribble(
  ~source_column, ~label, ~map_group,
  "T1T2ratio", "−T1w/T2w", "Microstructure", "PC1.AHBA", "Gene expression PC1", "Transcriptomic",
  "Cortical.Thickness", "Cortical thickness", "Microstructure", "G1.fMRI", "Functional gradient", "Functional connectivity",
  "intersubjvar", "Functional variability", "Functional connectivity", "PC1.Neurosynth", "NeuroSynth PC1", "Functional connectivity",
  "cmr02", "Oxygen metabolism", "Energy metabolism", "cmruglu", "Glucose metabolism", "Energy metabolism",
  "Evolution.Expansion", "Evolutionary expansion", "Cortical expansion", "devexp", "Developmental expansion", "Cortical expansion",
  "CBF", "Cerebral blood flow", "Hemodynamics", "cbv", "Cerebral blood volume", "Hemodynamics",
  "scalingnih", "Allometric scaling NIH", "Allometric scaling", "scalingpnc", "Allometric scaling PNC", "Allometric scaling",
  "meg_delta", "MEG delta", "Electrophysiology", "meg_theta", "MEG theta", "Electrophysiology", "meg_alpha", "MEG alpha", "Electrophysiology",
  "meg_beta", "MEG beta", "Electrophysiology", "meg_gamma1", "MEG low gamma", "Electrophysiology", "meg_gamma2", "MEG high gamma", "Electrophysiology")
GROUP_LEVELS <- c("Transcriptomic", "Microstructure", "Functional connectivity", "Energy metabolism", "Cortical expansion", "Hemodynamics", "Allometric scaling", "Electrophysiology")
GROUP_STRIP <- c("Transcriptomic", "Microstructure", "Functional\nconnectivity", "Energy\nmetabolism", "Cortical\nexpansion", "Hemo-\ndynamics", "Allometric\nscaling", "Electro-\nphysiology")
chk(setequal(groups$source_column, ST8$source_column), "a/b: the 20 normative maps of ST9 (no neurotransmitter maps)")
# rho/p-spin(cDEC, map) (not stored in ST9): computed with the reference map rotated
rotR <- function(v) matrix(rank(v)[P], nrow = 180)
cdec <- z(DEL) - z(DUP)
FLIP <- c(T1T2ratio = -1)                                   # T1w/T2w shown as -T1w/T2w (= H)
mapv <- function(k) maps[[k]] * (if (k %in% names(FLIP)) FLIP[[k]] else 1)
cd <- rbindlist(lapply(groups$source_column, function(k) { m <- mapv(k); o <- cor(cdec, m, method = "spearman"); nul <- drop(cor(rank(cdec), rotR(m)))
  data.table(source_column = k, rho_cDEC = o, pspin_cDEC = pv2(o, nul)) }))
ref_cd <- fread(file.path(inputs_dir, "inputs_Fig3_reference_cDEC_vs_maps.csv"))
chk(isTRUE(all.equal(cd$rho_cDEC, ifelse(cd$source_column %in% names(FLIP), -1, 1) * ref_cd[match(cd$source_column, source_column), rho_cDEC])) &&
    isTRUE(all.equal(cd$pspin_cDEC, ref_cd[match(cd$source_column, source_column), pspin_cDEC])), "a: computed rho/p-spin(cDEC, map) == reference values of Oct 2026 (T1w/T2w sign-reversed)")
ST8f <- copy(ST8)[source_column %in% names(FLIP), `:=`(rho_del = -rho_del, rho_dup = -rho_dup, delta_r_dup_minus_del = -delta_r_dup_minus_del)]   # spin p two-sided: unchanged
rows <- ST8f[groups, on = "source_column"][cd, on = "source_column"] %>% as_tibble() %>%
  mutate(map_group = factor(map_group, levels = GROUP_LEVELS, labels = GROUP_STRIP)) %>% group_by(map_group) %>% arrange(desc(rho_cDEC), .by_group = TRUE) %>% ungroup() %>%
  mutate(lab = factor(label, levels = rev(label)))
s6 <- as.data.table(read_excel(filename_ST_C, "ST7"))[gene_set_exp == "All_data_analysis_Fig2"]; t1 <- rows[rows$source_column == "T1T2ratio", ]
chk(abs(t1$rho_cDEC - s6$rho_cDEC_H) < 1e-12 && abs(t1$rho_del - s6$rho_DEL_H) < 1e-12 && abs(t1$rho_dup - s6$rho_DUP_H) < 1e-12,
    sprintf("-T1w/T2w row == H in ST7: rho(cDEC,H) %+.3f (p-spin %s), rho(Del,H) %+.3f (p-spin %s), rho(Dup,H) %+.3f (p-spin %s)", t1$rho_cDEC, fmt_p(t1$pspin_cDEC),
            t1$rho_del, fmt_p(t1$pspin_del), t1$rho_dup, fmt_p(t1$pspin_dup)))
STARS <- c("Gene expression PC1", "Cortical thickness", "Functional variability", "Functional gradient", "−T1w/T2w", "Oxygen metabolism", "Evolutionary expansion", "MEG high gamma", "MEG delta")
chk(setequal(rows$label[rows$pspin_cDEC < .05], STARS), sprintf("a: %d asterisks: %s", sum(rows$pspin_cDEC < .05), paste(rows$label[rows$pspin_cDEC < .05], collapse = ", ")))

grp_plain <- setNames(GROUP_LEVELS, GROUP_STRIP)
df_source_ab <- rows %>% transmute(map = label, source_column, domain = grp_plain[as.character(map_group)],
                                   rho_cDEC, pspin_cDEC, star_a = pspin_cDEC < .05,
                                   rho_del, pspin_del, rho_dup, pspin_dup, outline_del = pspin_del < .05, outline_dup = pspin_dup < .05,
                                   delta_rho_dup_minus_del = delta_r_dup_minus_del, pspin_delta = pspin_delta_r, n_rot = NR)
fn_write_SourceData_sheet(df_source_ab, "Fig3a-b", source_data_dir)

band <- tibble(map_group = factor(GROUP_STRIP[4:5], levels = GROUP_STRIP), fill = c(COL_BAND_ENERGY, COL_BAND_EXPANSION))
th_rows <- theme_R1(grid_x = TRUE) + theme(panel.spacing.y = unit(2, "pt"), strip.placement = "outside", strip.background = element_blank(),
  strip.text.y.left = element_text(angle = 0, hjust = 1, size = PT_TEXT, face = "bold", margin = margin(r = 4)), axis.ticks.y = element_blank(), axis.line.y = element_blank(), plot.margin = margin(12, 2, 2, 24))
key_a <- tibble(map_group = factor(GROUP_STRIP[7], levels = GROUP_STRIP), x = 0.6, y = 1)   # key on the Allometric scaling PNC row
pa <- ggplot(rows, aes(y = lab)) +
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = I(fill)), inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = LW_MIN * 1.4) +
  geom_col(aes(x = rho_cDEC), fill = COL_CDEC, width = .68, orientation = "y") +
  geom_text(data = ~ filter(.x, pspin_cDEC < .05), aes(x = rho_cDEC + ifelse(rho_cDEC >= 0, .03, -.03), label = "*"), vjust = .75, size = pt2mm(PT_TITLE + 2), family = FONT) +
  geom_text(data = key_a, aes(x = x, y = y, label = paste0("* ", P_LABEL, " < 0.05")), inherit.aes = FALSE, hjust = 1, vjust = 0, size = pt2mm(PT_SMALL + .5), family = FONT) +
  facet_grid(map_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(-.4, 0, .4)) + coord_cartesian(xlim = c(-.55, .55), clip = "off") +
  labs(x = "ρ with cDEC (Del − Dup contrast)", y = NULL) + th_rows
long_b <- bind_rows(rows %>% transmute(map_group, lab, CNV = "Del", rho = rho_del, p = pspin_del), rows %>% transmute(map_group, lab, CNV = "Dup", rho = rho_dup, p = pspin_dup))
pb <- ggplot() +
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = I(fill))) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = LW_MIN * 1.4) +
  geom_segment(data = rows, aes(x = rho_del, xend = rho_dup, y = lab, yend = lab), colour = COL_LINE, linewidth = LW_MIN * 4) +
  geom_point(data = long_b, aes(x = rho, y = lab, colour = CNV), size = 2) +
  geom_point(data = filter(long_b, p < .05), aes(x = rho, y = lab, shape = "sig"), size = 2, stroke = .6, colour = "black", fill = NA) +
  scale_colour_manual(values = c(Del = COL_DEL, Dup = COL_DUP), name = NULL) +
  scale_shape_manual(values = c(sig = 21), labels = c(sig = paste0(P_LABEL, " < 0.05")), name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2)), shape = guide_legend(order = 2, override.aes = list(size = 2, stroke = .6))) +
  facet_grid(map_group ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = c(-.5, 0, .5)) + coord_cartesian(xlim = c(-.62, .62)) +
  labs(x = "ρ with Del / Dup effect profile", y = NULL) + th_rows +
  theme(strip.text = element_blank(), axis.text.y = element_blank(), legend.position = "inside", legend.position.inside = c(.84, .06),
        legend.background = element_rect(fill = "white", colour = NA), legend.spacing.y = unit(0, "pt"), legend.margin = margin(1, 1, 1, 1), legend.box = "vertical")
chk(all((long_b$p < .05) == (c(ST8f[match(rows$source_column, source_column), pspin_del], ST8f[match(rows$source_column, source_column), pspin_dup]) < .05)),
    sprintf("b: %d outlines == ST9 pspin_del/pspin_dup < 0.05", sum(long_b$p < .05)))
```

#### Panels c and d ####

``` r
#----- Panels c, d: energy metabolism and cortical expansion, partial rho | H (Fig. 3c,d)
rz <- function(v) { r <- rank(v); (r - mean(r)) / sd(r) }
dims <- list(Energy = (rz(maps$cmr02) + rz(maps$cmruglu)) / 2, Expansion = (rz(maps$Evolution.Expansion) + rz(maps$devexp)) / 2)
partial_obs <- function(y, x, h) { ry <- rank(y); rx <- rank(x); rh <- rank(h); cor(resid(lm(ry ~ rh)), resid(lm(rx ~ rh))) }
partial_null <- function(y, x, h) { RX <- rotR(x); RH <- rotR(h); ry <- rank(y); rxh <- cor(rank(x), rank(h)); ryx <- drop(cor(ry, RX)); ryh <- drop(cor(ry, RH))
  list(obs = partial_obs(y, x, h), null = (ryx - ryh * rxh) / sqrt((1 - ryh^2) * (1 - rxh^2))) }
raw_null <- function(y, x) list(obs = cor(y, x, method = "spearman"), null = drop(cor(rank(y), rotR(x))))
PR <- lapply(dims, function(x) list(DEL = list(adj = partial_null(DEL, x, H), raw = raw_null(DEL, x)), DUP = list(adj = partial_null(DUP, x, H), raw = raw_null(DUP, x))))
PS <- rbindlist(lapply(names(dims), function(d) rbindlist(lapply(c("DEL", "DUP"), function(t) { a <- PR[[d]][[t]]$adj; b <- PR[[d]][[t]]$raw
  data.table(dimension = d, cnv = t, partial_rho_H = a$obs, pspin_partial = pv2(a$obs, a$null), rho_raw = b$obs, pspin_raw = pv2(b$obs, b$null), n_rot = NR) }))))
IHf <- function(k, nul = FALSE) { g <- function(d, t) if (nul) PR[[d]][[t]][[k]]$null else PR[[d]][[t]][[k]]$obs; g("Energy", "DEL") - g("Energy", "DUP") - (g("Expansion", "DEL") - g("Expansion", "DUP")) }
IH <- data.table(statistic = c("I_H (adjusted for H)", "I (unadjusted)"), value = c(IHf("adj"), IHf("raw")), pspin = c(pv2(IHf("adj"), IHf("adj", TRUE)), pv2(IHf("raw"), IHf("raw", TRUE))), n_rot = NR)
ref <- fread(file.path(inputs_dir, "inputs_Fig3_reference_dimensions.csv"))
for (d in names(dims)) { r <- ref[analysis == "main_composite" & dimension == d]
  chk(abs(PS[dimension == d & cnv == "DEL", partial_rho_H] - r$rho_del_adjH) < 1e-12 && abs(PS[dimension == d & cnv == "DEL", pspin_partial] - r$pspin_del_adjH) < 1e-12 &&
      abs(PS[dimension == d & cnv == "DUP", partial_rho_H] - r$rho_dup_adjH) < 1e-12 && abs(PS[dimension == d & cnv == "DUP", pspin_partial] - r$pspin_dup_adjH) < 1e-12,
      sprintf("%s: partial rho | H Del %.3f (p %s), Dup %.3f (p %s) == reference", d, PS[dimension == d & cnv == "DEL", partial_rho_H], fmt_p(PS[dimension == d & cnv == "DEL", pspin_partial]),
              PS[dimension == d & cnv == "DUP", partial_rho_H], fmt_p(PS[dimension == d & cnv == "DUP", pspin_partial]))) }
ri <- ref[analysis == "interaction_I"]
chk(abs(IH$value[1] - ri$delta_adjH) < 1e-12 && abs(IH$pspin[1] - ri$pspin_delta_adjH) < 1e-12 && abs(IH$value[2] - ri$delta_raw) < 1e-12,
    sprintf("I_H = %.3f (p-spin %s); unadjusted %.3f (p-spin %s) == reference", IH$value[1], fmt_p(IH$pspin[1]), IH$value[2], fmt_p(IH$pspin[2])))
QA <- data.table(panel = c("c", "c", "d", "d"), test = c("Del ~ energy | H", "Dup ~ energy | H", "Del ~ expansion | H", "Dup ~ expansion | H"),
                 expected_rho = c(-.33, -.04, -.05, -.34), expected_p = c(.037, .644, .828, 2e-4))
QA[, `:=`(rho = PS[c(1, 2, 3, 4), partial_rho_H], pspin = PS[c(1, 2, 3, 4), pspin_partial])]
chk(all(abs(round(QA$rho, 2) - QA$expected_rho) < 0.006) && all(abs(QA$pspin - QA$expected_p) < 6e-4), "c/d reproduce the printed values")
# added-variable coordinates (rank percentile residualised on rank percentile of H)
hp <- 100 * rank(H) / 180; rres <- function(v) unname(resid(lm(100 * rank(v) / 180 ~ hp)))
X <- lapply(dims, rres); Y <- list(DEL = rres(DEL), DUP = rres(DUP))
for (d in names(dims)) for (t in c("DEL", "DUP")) chk(abs(cor(X[[d]], Y[[t]]) - PS[dimension == d & cnv == t, partial_rho_H]) < 1e-12, sprintf("%s-%s: cor(plotted x, plotted y) == partial rho | H", d, t))
XLIM <- c(-75, 75); YLIM <- c(-75, 110)
chk(all(abs(unlist(c(X, Y))) < 75), "all residuals within the shared axis limits")

pct <- function(v) 100 * (rank(v) - 1) / (length(v) - 1); PCT <- lapply(dims, pct)
df_source_cd <- data.frame(gene_set_name = roi, H = H, rank_H = rank(H),
                           energy_composite = dims$Energy, energy_composite_rank_percentile = PCT$Energy, energy_rank_resid_H = X$Energy,
                           expansion_composite = dims$Expansion, expansion_composite_rank_percentile = PCT$Expansion, expansion_rank_resid_H = X$Expansion,
                           del_rank_resid_H = Y$DEL, dup_rank_resid_H = Y$DUP)
fn_write_SourceData_sheet(df_source_cd, "Fig3c-d", source_data_dir)
df_stats_cd <- rbind(data.frame(panel = ifelse(PS$dimension == "Energy", "c", "d"), statistic = paste0("partial rho ", PS$cnv, " ~ ", PS$dimension, " | H"), value = PS$partial_rho_H, pspin = PS$pspin_partial, n_rot = NR),
                     data.frame(panel = "text", statistic = paste0("rho ", PS$cnv, " ~ ", PS$dimension, " (unadjusted)"), value = PS$rho_raw, pspin = PS$pspin_raw, n_rot = NR),
                     data.frame(panel = "text", statistic = IH$statistic, value = IH$value, pspin = IH$pspin, n_rot = NR))
fn_write_SourceData_sheet(df_stats_cd, "Fig3c-d_stats", source_data_dir)

scat <- function(d, t, ylab_on, xlab_on) {
  s <- PS[dimension == d & cnv == t]
  ggplot(tibble(x = X[[d]], y = Y[[t]]), aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = LW_MIN * 1.4) + annotate("segment", x = 0, xend = 0, y = YLIM[1], yend = 66, colour = "grey60", linewidth = LW_MIN * 1.4) +
    geom_point(colour = COL_POINT, alpha = .75, size = 1.1, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = if (t == "DEL") COL_DEL else COL_DUP, linewidth = LW_MIN * 4.5) +
    annotate("text", x = XLIM[1] + 3, y = YLIM[2], hjust = 0, vjust = 1, lineheight = .95, size = pt2mm(PT_SMALL + 0.5), family = FONT,
             label = sprintf("partial ρ | H = %.2f\n%s = %s", s$partial_rho_H, P_LABEL, fmt_p(s$pspin_partial))) +
    scale_x_continuous(breaks = c(-50, 0, 50)) + scale_y_continuous(breaks = c(-50, 0, 50)) +
    coord_cartesian(xlim = XLIM, ylim = YLIM, expand = FALSE, clip = "off") +
    labs(x = if (xlab_on) c(Energy = "Energy metabolism\n(rank, adjusted for H)", Expansion = "Cortical expansion\n(rank, adjusted for H)")[[d]] else NULL,
         y = if (ylab_on) paste0(if (t == "DEL") "Del" else "Dup", " effect\n(rank, adjusted for H)") else NULL) +
    theme_R1() + theme(axis.text.y = if (ylab_on) element_text(size = PT_SMALL + .5, colour = "black") else element_blank(),
                       axis.text.x = if (xlab_on) element_text(size = PT_SMALL + .5, colour = "black") else element_blank(), plot.margin = margin(2, 4, 2, 2))
}
SC <- list(cD = scat("Energy", "DEL", TRUE, FALSE), dD = scat("Expansion", "DEL", FALSE, FALSE), cU = scat("Energy", "DUP", TRUE, TRUE), dU = scat("Expansion", "DUP", FALSE, TRUE))

# composite brain maps (rank percentile, viridis), ggseg Glasser left hemisphere
glasser_region <- sub("^L_", "", gsub(".", "-", roi, fixed = TRUE)); chk(all(glasser_region %in% names(glasser$palette)), "180 regions in the ggseg Glasser atlas (name-based)")
one_view <- function(v, view, legend = FALSE) ggseg(tibble(region = glasser_region, value = v), atlas = "glasser", hemisphere = "left", view = view,
    mapping = aes(fill = value), colour = "grey92", size = .02) + theme_brain(text.family = FONT, text.size = PT_SMALL) +
  scale_fill_viridis_c(limits = c(0, 100), breaks = c(0, 50, 100), name = "Rank percentile", na.value = "grey85") +
  guides(fill = guide_colourbar(barwidth = unit(80, "pt"), barheight = unit(4, "pt"), title.position = "top", title.hjust = .5)) +
  theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = if (legend) "bottom" else "none",
        legend.title = element_text(size = PT_SMALL + .5, family = FONT), legend.text = element_text(size = PT_SMALL + .5, family = FONT), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0), plot.margin = margin(0, 0, 0, 0))
brain_comp <- function(v) plot_grid(one_view(v, "lateral"), one_view(v, "medial"), nrow = 1)
pal <- ggsegYeo2011::yeo17$palette
yeo_view <- function(v) ggseg(atlas = yeo17, hemisphere = "left", view = v, mapping = aes(fill = region), colour = "white", size = .03) +
  scale_fill_manual(values = pal, na.value = "grey90", guide = "none") + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
yeo_map <- plot_grid(yeo_view("lateral"), yeo_view("medial"), nrow = 1)
title_box <- function(txt, fill) { tg <- grid::textGrob(txt, gp = grid::gpar(fontsize = PT_TITLE, fontfamily = FONT, col = "black"))
  if (is.na(fill)) return(ggdraw() + draw_grob(tg))
  rg <- grid::rectGrob(width = grid::grobWidth(tg) + unit(8, "pt"), height = grid::grobHeight(tg) + unit(5, "pt"), gp = grid::gpar(fill = fill, col = NA))
  ggdraw() + draw_grob(grid::grobTree(rg, tg)) }
```

#### Panel e ####

``` r
#----- Panel e: Yeo-17 networks (Fig. 3e)
g <- read.table(file.path(inputs_dir, "inputs_Fig3_Glasser360_regions_mapped_to_fsaverage5.txt"))[, 1]; y <- read.table(file.path(inputs_dir, "inputs_Fig3_yeo17_networks_mapped_to_fsaverage5.txt"))[, 1]
g <- g[seq_len(length(g) / 2)]; y <- y[seq_len(length(y) / 2)]; keep <- g >= 1 & g <= 180 & y > 0
chk(length(g) == 10242, sprintf("fsaverage5 left hemisphere: 10,242 vertices; kept %d (dropped %d with Yeo = 0 [medial wall], %d with no Glasser region)", sum(keep), sum(y == 0), sum(g < 1 | g > 180)))
VT <- data.table(net = y[keep], reg = g[keep], DEL = DEL[g[keep]], DUP = DUP[g[keep]])   # Glasser index 1-180 = ST4 region order (checked above)
W <- as.matrix(table(factor(VT$net, levels = 1:17), factor(VT$reg, levels = 1:180))); Wn <- W / rowSums(W)
ye <- VT[, .(n_vertices = .N, n_glasser_regions = uniqueN(reg), mean_del = mean(DEL), sem_del_vertices = sd(DEL) / sqrt(.N), mean_dup = mean(DUP), sem_dup_vertices = sd(DUP) / sqrt(.N)), by = net][order(net)]
wsem <- function(xr, w) { w <- w / sum(w); m <- sum(w * xr); k <- sum(w > 0); if (k < 2) return(NA_real_); v <- sum(w * (xr - m)^2) / (1 - sum(w^2)); sqrt(v * sum(w^2)) }
ye[, `:=`(sem_del_regions = sapply(net, function(n) wsem(DEL, W[n, ])), sem_dup_regions = sapply(net, function(n) wsem(DUP, W[n, ])))]
chk(max(abs(drop(Wn %*% DEL) - ye$mean_del)) < 1e-12, "network mean = vertex-weighted mean of region values (W %*% profile)")
ye[, `:=`(network = sub("^Limibic", "Limbic", names(pal)[net]), colour = unname(pal)[net])]   # ggsegYeo2011 spells "Limibic"
# label colour: the network colour, darkened only where it is too pale to read on white (relative luminance > 0.6)
lum <- function(h) { r <- grDevices::col2rgb(h) / 255; r <- ifelse(r <= .03928, r / 12.92, ((r + .055) / 1.055)^2.4); drop(c(.2126, .7152, .0722) %*% r) }
dark <- function(h, f = .45) { r <- grDevices::col2rgb(h) * (1 - f); grDevices::rgb(r[1], r[2], r[3], maxColorValue = 255) }
ye[, label_colour := ifelse(sapply(colour, lum) > .6, sapply(colour, dark), colour)]
rho_y <- cor(ye$mean_del, ye$mean_dup, method = "spearman"); p_par <- suppressWarnings(cor.test(ye$mean_del, ye$mean_dup, method = "spearman", exact = FALSE)$p.value)
nx <- sapply(seq_len(NR), function(i) cor(drop(Wn %*% DEL[P[, i]]), ye$mean_dup, method = "spearman"))
ny <- sapply(seq_len(NR), function(i) cor(ye$mean_del, drop(Wn %*% DUP[P[, i]]), method = "spearman"))
p_spin_y <- mean(c(pv2(rho_y, nx), pv2(rho_y, ny)))
yref <- fread(file.path(inputs_dir, "inputs_Fig3_reference_Yeo17_spin.csv"))
chk(abs(rho_y - yref$rho_DEL_DUP_networks) < 1e-12 && abs(p_spin_y - yref$pspin_profile_pair_rule) < 1e-12 && abs(p_par - yref$p_parametric_old) < 1e-12,
    sprintf("Yeo-17: rho = %.3f, parametric P = %.3f (t approximation, n = 17), p-spin = %.3f (profile-pair rule, n_rot = %d) == reference", rho_y, p_par, p_spin_y, NR))

df_source_e <- ye[, .(network_index = net, network, colour_hex = colour, n_vertices, n_glasser_regions, mean_del, sem_del_vertices, sem_del_regions,
                      mean_dup, sem_dup_vertices, sem_dup_regions, rho_networks = rho_y, p_parametric = p_par, pspin_profile_pair = p_spin_y, n_rot = NR)]
fn_write_SourceData_sheet(df_source_e, "Fig3e", source_data_dir)

pe <- ggplot(ye, aes(mean_dup, mean_del, colour = factor(net))) +
  geom_linerange(aes(ymin = mean_del - sem_del_vertices, ymax = mean_del + sem_del_vertices), linewidth = LW_MIN * 4) +
  geom_linerange(aes(xmin = mean_dup - sem_dup_vertices, xmax = mean_dup + sem_dup_vertices), linewidth = LW_MIN * 4) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = network), colour = ye$label_colour, size = pt2mm(PT_SMALL), family = FONT, max.overlaps = Inf, min.segment.length = 0, segment.size = LW_MIN * 1.5,
                  box.padding = .18, force = 2, seed = SEED, show.legend = FALSE) +
  scale_colour_manual(values = unname(pal), guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, label = sprintf("ρ = %.2f", rho_y), size = pt2mm(PT_SMALL + .5), family = FONT) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) + labs(x = "Mean effect size, Dup", y = "Mean effect size, Del") + theme_R1()

#-----------------------------------------------------------------------------
# PART D: Assemble (183 mm x 225 mm) and save
#-----------------------------------------------------------------------------
```

#### Assemble and save ####

``` r
W_C <- 58; W_D <- 50; W_E <- 75
top <- (pa | pb) + plot_layout(widths = c(.44, .56))
top <- plot_grid(patchworkGrob(top)) + draw_label("a", x = .005, y = .995, hjust = 0, vjust = 1, size = PT_TAG, fontface = "bold", fontfamily = FONT) +
  draw_label("b", x = .585, y = .995, hjust = 0, vjust = 1, size = PT_TAG, fontface = "bold", fontfamily = FONT)
titles <- plot_grid(title_box("Energy metabolism", COL_TITLE_ENERGY), title_box("Cortical expansion", COL_TITLE_EXPANSION), title_box("Yeo-17 networks", NA), nrow = 1, rel_widths = c(W_C, W_D, W_E))
pad <- function(p, f = .12) plot_grid(NULL, p, NULL, nrow = 1, rel_widths = c(f, 1, f))
mapsrow <- plot_grid(plot_grid(NULL, pad(brain_comp(PCT$Energy)), nrow = 1, rel_widths = c(.12, 1)), pad(brain_comp(PCT$Expansion)), pad(yeo_map, .1), nrow = 1, rel_widths = c(W_C, W_D, W_E))
cbar <- plot_grid(plot_grid(NULL, get_legend(one_view(PCT$Energy, "lateral", TRUE)), nrow = 1, rel_widths = c(.12, 1)), NULL, nrow = 1, rel_widths = c(W_C + W_D, W_E))
sc4 <- plot_grid(SC$cD, SC$dD, SC$cU, SC$dU, ncol = 2, align = "hv", axis = "tblr", rel_widths = c(W_C, W_D))
lower <- plot_grid(sc4, plot_grid(pe, NULL, ncol = 1, rel_heights = c(1, .08)), nrow = 1, rel_widths = c(W_C + W_D, W_E))
bottom <- plot_grid(titles, mapsrow, cbar, lower, ncol = 1, rel_heights = c(7, 24, 12, 84))
bottom <- bottom + draw_label("c", x = .005, y = .995, hjust = 0, vjust = 1, size = PT_TAG, fontface = "bold", fontfamily = FONT) +
  draw_label("d", x = W_C / (W_C + W_D + W_E) + .005, y = .995, hjust = 0, vjust = 1, size = PT_TAG, fontface = "bold", fontfamily = FONT) +
  draw_label("e", x = (W_C + W_D) / (W_C + W_D + W_E) + .005, y = .995, hjust = 0, vjust = 1, size = PT_TAG, fontface = "bold", fontfamily = FONT)
fig <- plot_grid(top, bottom, ncol = 1, rel_heights = c(100, 125)) + theme(plot.background = element_rect(fill = "white", colour = NA))
f_out <- file.path(OUT, paste0(STEM, c(".svg", ".png")))   # PNG + SVG only
svglite::svglite(f_out[1], width = 183 * MM, height = 225 * MM, bg = "white"); print(fig); invisible(dev.off())
ragg::agg_png(f_out[2], width = 183 * MM, height = 225 * MM, units = "in", res = 600, background = "white"); print(fig); invisible(dev.off())
message("Saved Fig. 3 to: ", OUT)
```

#### Fig. 3 ####
![Fig. 3](Fig3.png)
