#-------------------------------------------------------------------------------
# Summary: R-script to plot Supplementary Fig. 1 a-c, g, h (summary of CNV carriers; unique genes per cortical gene set)
#-------------------------------------------------------------------------------

#-----
# # Summary:
#          -> a) unique genes (% of the 13,678 genes of the 180 AHBA gene sets), b) total genes, c) CNV carriers, Del and Dup
#             (Supplementary Table 3, overall burden, main analysis)
#          -> g) number of unique deleted / duplicated genes per region (Glasser-180 maps; Supplementary Table 4)
#          -> h) unique deleted vs duplicated genes per region (180 regions; dot-dash line = linear fit, as in the published panel)
#          -> Panels d-f (SST example: Z-score map, threshold, regions with relative high expression) are an illustration;
#             no code or source data. In the output, row 2 is left empty for these panels.
#
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 21 October, 2026
# Version: 1
#
#-----
# UPDATE log:
#       v1 (29 Oct 2026): release version of code/RscriptMac_mainFig1_count_nGenesPerROI_11May2023.R (panels a-c, g, h);
#                         counts read from Supplementary Tables 3 and 4 instead of the internal CC-GSBA result files;
#                         brain-map function copied from code/Rscript_Mac_BasePlotFunctions_PanelFigs_ggsegScatterBarBoxTile_95CI_14mar2022_v1.R
#
#-----
# Input files:
#       1) data/inputs/SupplementaryTables_NatGenRevision_29oct2026.xlsx (ST3, ST4)
#
#-----
# Output:
#       1) SupFigs/output/SupFig01.png, SupFig01.svg
#       2) data/source_data_sheets/SupFig1a-c.csv, SupFig1g-h.csv
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
inputs_dir <- file.path(project_dir, "data", "inputs")
source_data_dir <- file.path(project_dir, "data", "source_data_sheets")
output_dir <- file.path(project_dir, "SupFigs", "output")

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(ggseg); library(ggsegGlasser); library(dplyr); library(tibble); library(readxl)
  library(cowplot); library(ggpubr); library(svglite)
})
source(file.path(project_dir, "functions", "fn_write_SourceData_sheet.R"))

#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

filename_ST_C <- file.path(inputs_dir, "SupplementaryTables_NatGenRevision_29oct2026.xlsx")
in_fig_lab_size <- 20
in_txt_font_size <- 18
in_leg_count_font <- 16
in_atlas_name <- "glasser"
total_geneset_AHBA <- 13678    # genes in the union of the 180 regional gene sets (Z > 0.5)

#---- Brain map of counts (from Rscript_Mac_BasePlotFunctions_..._95CI_14mar2022_v1.R; Glasser branch)
fBrainPlot_ggseg_Cortex_array_count_single_plot <- function(array_count,cbar_min=0,cbar_max=10000,in_legtitle,in_atlas_name){
  pallete_ggseg = glasser$palette
  ggseg_region_names <- names(pallete_ggseg)[c(1:180)] # keep LH only
  # Max value limit: keep estimates within cbar limit
  array_count[array_count > cbar_max] <- cbar_max
  array_count[array_count < cbar_min] <- cbar_min
  temp_df <- as.data.frame(cbind(array_count))
  colnames(temp_df) <- c("count")
  someData <- tibble(region = ggseg_region_names, effect = temp_df[,"count"])
  p_left <-  ggseg(someData,
                   atlas = in_atlas_name,
                   position="dispersed",
                   hemisphere = "left",
                   mapping =aes(fill = effect,size=I(0.30)),
                   colour = "grey",
                   adapt_scales = TRUE) +
    theme_brain(text.family = "sans",text.size = 14) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c("yellow","orange","brown"))(n = 600)),0.8) +
    theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(title = in_legtitle, title.position = "top", barwidth = 7, barheight =0.7, title.hjust =0.5, title.vjust = 0)) +
    theme(legend.text = element_text(size = in_leg_count_font)) +
    theme(legend.title = element_text(size = in_leg_count_font))  +
    theme(legend.margin=margin(t = 0, unit='cm'))
  cowplot_legend <- cowplot::get_legend(p_left)
  p_left <- p_left + theme(legend.position = "none")
  p_single_plot <- p_left + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))
  return(list(p_single_plot,cowplot_legend))
}

#-----------------------------------------------------------------------------
# PART B: Counts from the Supplementary Tables
#-----------------------------------------------------------------------------

df_ST3 <- as.data.frame(read_excel(filename_ST_C, "ST3"))
df_ST3 <- df_ST3[df_ST3$gene_set_exp == "All_data_analysis_Fig2", ]
df_overall <- data.frame(n_unique_genes_DEL = df_ST3$n_unique_genes[df_ST3$cnv_type == "DEL"], n_unique_genes_DUP = df_ST3$n_unique_genes[df_ST3$cnv_type == "DUP"],
                         N_total_genes_DEL = df_ST3$N_total_genes[df_ST3$cnv_type == "DEL"], N_total_genes_DUP = df_ST3$N_total_genes[df_ST3$cnv_type == "DUP"],
                         total_carriers_DEL = df_ST3$total_carriers[df_ST3$cnv_type == "DEL"], total_carriers_DUP = df_ST3$total_carriers[df_ST3$cnv_type == "DUP"])
stopifnot(df_ST3$gene_list_n[1] == total_geneset_AHBA)

df_ST4 <- as.data.frame(read_excel(filename_ST_C, "ST4"))
df_ST4 <- df_ST4[df_ST4$gene_set_exp == "All_data_analysis_Fig2", ]
ggseg_region_names <- names(glasser$palette)[c(1:180)]
gene_list_name_order <- paste0("L_", gsub("-", ".", ggseg_region_names))     # region order of the ggseg atlas (joined by name)
df_reorder <- data.frame(gene_list_name = gene_list_name_order,
                         n_unique_genes_DEL = df_ST4$n_unique_genes[df_ST4$cnv_type == "DEL"][match(gene_list_name_order, df_ST4$gene_set_name[df_ST4$cnv_type == "DEL"])],
                         n_unique_genes_DUP = df_ST4$n_unique_genes[df_ST4$cnv_type == "DUP"][match(gene_list_name_order, df_ST4$gene_set_name[df_ST4$cnv_type == "DUP"])])
stopifnot(!anyNA(df_reorder))

#-----------------------------------------------------------------------------
# PART C: Panels a-c (bar plots)
#-----------------------------------------------------------------------------

fn_barplot_input <- function(count_del, count_dup, denom) {
  df_barplot_input <- data.frame(array_TYPE = factor(c("Del", "Dup"), levels = c("Del", "Dup")), array_count_summary = c(count_del, count_dup))
  df_barplot_input$array_count_summary_percent <- 100 * df_barplot_input$array_count_summary / denom
  df_barplot_input
}
in_legend_title <- " "

#----- Panel a: unique genes (% of the 13,678 genes) (Supplementary Fig. 1a)
df_a <- fn_barplot_input(df_overall$n_unique_genes_DEL, df_overall$n_unique_genes_DUP, total_geneset_AHBA)
p_barplot_unq <- ggplot(data=df_a,aes(x=array_TYPE,y=array_count_summary_percent,fill=array_TYPE))+
  geom_bar(stat = "identity") +
  geom_text(aes(label=array_count_summary), vjust=1.6, color="black", size=7)+
  theme_minimal()+ xlab(NULL)+ ylab("% genes")+ ylim(0,100)+ labs(fill=in_legend_title) +
  scale_fill_manual(values=c("#DC0000FF","#3C5488FF"),guide = 'none')+
  theme(text = element_text(size = in_txt_font_size))+ theme(axis.text.x = element_text(angle=45))

#----- Panel b: total genes (Supplementary Fig. 1b)
df_b <- fn_barplot_input(df_overall$N_total_genes_DEL, df_overall$N_total_genes_DUP, 15633)
p_barplot_Total <- ggplot(data=df_b,aes(x=array_TYPE,y=array_count_summary,fill=array_TYPE))+
  geom_bar(stat = "identity") +
  geom_text(aes(label=array_count_summary), vjust=1.6, color="black", size=7)+
  theme_minimal()+ xlab(NULL)+ ylab("# genes")+ ylim(0,80000)+ labs(fill=in_legend_title) +
  scale_fill_manual(values=c("#DC0000FF","#3C5488FF"),guide = 'none')+
  theme(text = element_text(size = in_txt_font_size))+ theme(axis.text.x = element_text(angle=45))

#----- Panel c: CNV carriers (Supplementary Fig. 1c)
df_c <- fn_barplot_input(df_overall$total_carriers_DEL, df_overall$total_carriers_DUP, 15633)
p_barplot_CNV_Carriers <- ggplot(data=df_c,aes(x=array_TYPE,y=array_count_summary,fill=array_TYPE))+
  geom_bar(stat = "identity") +
  geom_text(aes(label=array_count_summary), vjust=1.6, color="black", size=7)+
  theme_minimal()+ xlab(NULL)+ ylab("# individuals")+ ylim(0,30000)+ labs(fill=in_legend_title) +
  scale_fill_manual(values=c("#DC0000FF","#3C5488FF"),guide = 'none')+
  theme(text = element_text(size = in_txt_font_size))+ theme(axis.text.x = element_text(angle=45))

df_SD_ac <- rbind(data.frame(panel = "a", quantity = "unique genes", cnv_type = df_a$array_TYPE, count = df_a$array_count_summary, percent_of_13678_genes = df_a$array_count_summary_percent),
                  data.frame(panel = "b", quantity = "total genes", cnv_type = df_b$array_TYPE, count = df_b$array_count_summary, percent_of_13678_genes = NA),
                  data.frame(panel = "c", quantity = "CNV carriers", cnv_type = df_c$array_TYPE, count = df_c$array_count_summary, percent_of_13678_genes = NA))
fn_write_SourceData_sheet(df_SD_ac, "SupFig1a-c", source_data_dir)

#-----------------------------------------------------------------------------
# PART D: Panels g-h (unique genes per region)
#-----------------------------------------------------------------------------

#----- Panel g: unique deleted / duplicated genes per region (Supplementary Fig. 1g)
fn_count_map <- function(array_count) {
  in_cbar_min = floor(min(array_count)/100)*100
  in_cbar_max = ceiling(max(array_count)/100)*100
  list_plots <- fBrainPlot_ggseg_Cortex_array_count_single_plot(array_count,in_cbar_min,in_cbar_max,'# unique genes',in_atlas_name)
  ggarrange(list_plots[[1]],NULL,list_plots[[2]], ncol = 1,heights = c(0.3,-0.055,0.15))
}
p_map_Count_DEL <- fn_count_map(df_reorder[,"n_unique_genes_DEL"])
p_map_Count_DUP <- fn_count_map(df_reorder[,"n_unique_genes_DUP"])
p_map_DEL_DUP_unq <- ggarrange(p_map_Count_DEL,p_map_Count_DUP,nrow = 1,labels = c("Del", "Dup"),heights = c(0.8,0.8), label.x=0.43, label.y = c(0.8,0.8))

#----- Panel h: unique deleted vs duplicated genes per region (Supplementary Fig. 1h)
df_scatter_input <- data.frame(array_count_DEL = df_reorder[,"n_unique_genes_DEL"], array_count_DUP = df_reorder[,"n_unique_genes_DUP"])
fn_write_SourceData_sheet(data.frame(region = df_reorder$gene_list_name, n_unique_genes_DEL = df_reorder$n_unique_genes_DEL, n_unique_genes_DUP = df_reorder$n_unique_genes_DUP),
                          "SupFig1g-h", source_data_dir)
p_scatter_unq = ggplot(df_scatter_input, aes(x=array_count_DUP, y=array_count_DEL)) +
  theme_minimal()+
  geom_point(color = "grey80",size=3,alpha = 1) +
  geom_smooth(method = "lm",se=FALSE,linetype="dotdash",fullrange=TRUE,color = "black")+
  theme(legend.position = "none")+
  theme(text = element_text(size = in_txt_font_size))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("# genes Dup")+
  ylab("# genes Del")

#-----------------------------------------------------------------------------
# PART E: Assemble (row 2 = panels d-f, illustration, left empty); save PNG + SVG
#-----------------------------------------------------------------------------

p_row1 <- ggarrange(p_barplot_unq + labs(tag = "Unique genes") + theme(plot.tag.position = "top") ,
                    NULL,p_barplot_Total + labs(tag = "Total genes") + theme(plot.tag.position = "top") ,
                    NULL,p_barplot_CNV_Carriers + labs(tag = "CNV carriers") + theme(plot.tag.position = "top") ,
                    nrow = 1,labels = c("a","","b","","c"),font.label = list(size = in_fig_lab_size),
                    widths = c(0.8,0.02,0.8,0.02,0.8)) + theme(text = element_text(size = 32))
p_row3 <- ggarrange(p_map_DEL_DUP_unq,NULL,p_scatter_unq,nrow = 1,
                    labels = c("g","","h"),font.label = list(size = in_fig_lab_size),
                    widths = c(0.7,0.02,0.3)) + theme(text = element_text(size = 32))
p_row3_arrow <- p_row3 + geom_segment(aes(x = 0.68, y = 0.5, xend =0.71, yend = 0.5),size = 2,arrow = arrow(length = unit(0.5, "cm")))
p_fig1_ext <- ggarrange(NULL,p_row1, NULL,NULL, NULL,p_row3_arrow, NULL,ncol=1,heights = c(0.01,0.8,0.02,0.8,0.05,0.8,0.005)) +
  theme(text = element_text(size = 32), plot.background = element_rect(fill = "white", colour = NA))

in_height <- 1.5*1600; in_width <- 1.5*2400   # pixels at 300 dpi, as the published panel
png(filename = file.path(output_dir, "SupFig01.png"), width = in_width, height = in_height, res = 300)
print(p_fig1_ext)
invisible(dev.off())
if (flag_hires == 1) {   # SVG only for the publication files (the PNG is 300 dpi, as the published panel)
  svglite::svglite(file.path(output_dir, "SupFig01.svg"), width = in_width / 300, height = in_height / 300, bg = "white")
  print(p_fig1_ext)
  invisible(dev.off())
}
