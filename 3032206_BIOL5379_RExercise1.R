## Introduction ##########
# Loading all the necessary libraries 
#install.packages("rlang")
#BiocManager::install("org.Hs.eg.db", force=TRUE)

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(eulerr)
library(reshape2)
library(amap)
library(BiocManager)
library(clusterProfiler) 
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(STRINGdb) 
library(devEMF)
library(cowplot)
library(rlang)

# loading the functions
source("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/assignments_functions.R") 
output_path = "C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/plots/"


# Loading the data
# ss - sample_sheet; EM - ; de1, de2, de3 -
em = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/EM.csv", header = TRUE, row.names = 1, sep = "\t")

de_S_vs_P = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/DE_Senes_vs_Prolif.csv", header = TRUE, row.names = 1, sep = "\t")
de_Md_vs_S = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/DE_Senes_MtD_vs_Senes.csv", header = TRUE, row.names = 1, sep = "\t")
de_Md_vs_P = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/DE_Senes_MtD_vs_Prolif.csv", header = TRUE, row.names = 1, sep = "\t")

ss = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/sample_sheet.csv", header = TRUE, row.names = 1, sep = "\t")

annotations = read.table("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/Assignment/Human_Background_GRCh38.p13.csv", header = TRUE, row.names = 1, sep = "\t")

# Part 1: master table for each comparison ####
# merging the data
# firstly, I want to have a separate table for each 1vs1 analysis. Then I want to have an ultramaster table with all groups
master_SP = merge_master(de_S_vs_P, annotations, em)
master_MdS = merge_master(de_Md_vs_S, annotations, em)
master_MdP = merge_master(de_Md_vs_P, annotations, em)

# I'm not using all the columns, and I do not want to keep those
master_SP = master_SP[,-c(13:15)]
master_MdS = master_MdS[,-c(13:15)]
master_MdP = master_MdP[,-c(13:15)]

# I also do not want any empty values
master_SP = na.omit(master_SP)
master_MdS = na.omit(master_MdS)
master_MdP = na.omit(master_MdP)

# changing the row names into symbols
master_SP = row_rename(master_SP, 11, "gene_symbol") 
master_MdS = row_rename(master_MdS, 11, "gene_symbol")
master_MdP = row_rename(master_MdP, 11, "gene_symbol")

names(master_SP)[1] = "ENSEMBL"
names(master_MdS)[1] = "ENSEMBL"
names(master_MdP)[1] = "ENSEMBL"

# making a new column in your master table for mean expression
master_means = rowMeans(master_SP[2:10])
master_SP$expression_mean = master_means

master_means = rowMeans(master_MdS[2:10])
master_MdS$expression_mean = master_means

master_means = rowMeans(master_MdP[2:10])
master_MdP$expression_mean = master_means

# Task 4. Add a column for -log10p to the master table.
log_values = -log10(master_SP$p)         # calc the p values as vector
master_SP$mlog10p = log_values           # adding the p values to the table

log_values = -log10(master_MdS$p) 
master_MdS$mlog10p = log_values

log_values = -log10(master_MdP$p) 
master_MdP$mlog10p = log_values

# some have a value of Inf bc of the p = 0

# creating the description columns based on the significance
master_SP$sig = as.factor(master_SP$p.adj < 0.001 & abs(master_SP$log2fold) > 2.0)
master_MdS$sig = as.factor(master_MdS$p.adj < 0.001 & abs(master_MdS$log2fold) > 2.0)
master_MdP$sig = as.factor(master_MdP$p.adj < 0.001 & abs(master_MdP$log2fold) > 2.0)

# Part 2: creating the ultramaster with all comparisons ####
# merging the tables
ultramaster_de = merge(de_Md_vs_S, de_Md_vs_P, by.x=0, by.y=0, suffixes=c(".MvS",""))
ultramaster_de = merge(ultramaster_de, de_S_vs_P, by.x=1, by.y=0, suffix=c(".MvP", ".SvP"))
row.names(ultramaster_de) = ultramaster_de[,1]
ultramaster_de = ultramaster_de[,-1]
ultramaster = merge(em,annotations, by.x=0, by.y=0)
ultramaster = merge(ultramaster, ultramaster_de, by.x=1, by.y=0)

# renaming the row indexes and deleting the spare columns
row.names(ultramaster) = ultramaster[,1]
ultramaster = ultramaster[,-c(13:15)]

# changing the row names of the ultramaster
ultramaster = row_rename(ultramaster, 11, "gene_symbol")
names(ultramaster)[1] = "ENSEMBL"

# creating the description columns based on the significance
ultramaster$sig_SvP = as.factor(ultramaster$p.adj.SvP < 0.001 & abs(ultramaster$log2fold.SvP) > 2.0)
ultramaster$sig_MvS = as.factor(ultramaster$p.adj.MvS < 0.001 & abs(ultramaster$log2fold.MvS) > 2.0)
ultramaster$sig_MvP = as.factor(ultramaster$p.adj.MvP < 0.001 & abs(ultramaster$log2fold.MvP) > 2.0)

# scaled versions of the tables
em_SP_scaled = scale_omit(master_SP, 2:10)
em_MdS_scaled = scale_omit(master_MdS, 2:10)
em_MdP_scaled = scale_omit(master_MdP, 2:10)

row.names(em_SP_scaled) = row.names(master_SP)
row.names(em_MdS_scaled) = row.names(master_MdS)
row.names(em_MdP_scaled) = row.names(master_MdP)

# Narrowing the data tables to only the significant genes
# and sorting by the p value
# and grouping by the direction of the change
# grouping the full master table for the sake of the volcano plot

master_SP_sign = subset(master_SP, master_SP$sig == TRUE)
master_SP_sign = sort_p(master_SP_sign)
master_SP_sign = group_sign_only(master_SP_sign)
master_SP = group_sign_full(master_SP)

master_MdS_sign = subset(master_MdS, master_MdS$sig == TRUE)
master_MdS_sign = sort_p(master_MdS_sign)
master_MdS_sign = group_sign_only(master_MdS_sign)
master_MdS = group_sign_full(master_MdS)

master_MdP_sign = subset(master_MdP, master_MdP$sig == TRUE)
master_MdP_sign = sort_p(master_MdP_sign)
master_MdP_sign = group_sign_only(master_MdP_sign)
master_MdP = group_sign_full(master_MdP)

# lists of the top sign genes for each comparison

top10up_SP = top10up(master_SP_sign, 0.001, 2)
top10up_MdS = top10up(master_MdS_sign, 0.001, 2)
top10up_MdP = top10up(master_MdP_sign, 0.001, 2)

top10down_SP = top10down(master_SP_sign, 0.001, 2)
top10down_MdS = top10down(master_MdS_sign, 0.001, 2)
top10down_MdP = top10down(master_MdP_sign, 0.001, 2)

# for ultramaster: any of the comparisons significant

ultramaster_AnySign = subset(ultramaster, ultramaster$sig_SvP == TRUE | ultramaster$sig_MvS == TRUE | ultramaster$sig_MvP == TRUE)

# columns for the direction of changes
  # for the SvP
ultramaster_AnySign_up = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.SvP > 1.0)
ultramaster_AnySign_down = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.SvP < -1.0)
ultramaster_AnySign_up$direction.SvP = "a_up"
ultramaster_AnySign_down$direction.SvP = "b_down"
ultramaster_AnySign =  rbind(ultramaster_AnySign_up, ultramaster_AnySign_down)

  # for the MvS
ultramaster_AnySign_up = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.MvS > 1.0)
ultramaster_AnySign_down = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.MvS < -1.0)
ultramaster_AnySign_up$direction.MvS = "a_up"
ultramaster_AnySign_down$direction.MvS = "b_down"
ultramaster_AnySign =  rbind(ultramaster_AnySign_up, ultramaster_AnySign_down)

  # for the MvP
ultramaster_AnySign_up = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.MvP > 1.0)
ultramaster_AnySign_down = subset(ultramaster_AnySign, ultramaster_AnySign$log2fold.MvP < -1.0)
ultramaster_AnySign_up$direction.MvP = "a_up"
ultramaster_AnySign_down$direction.MvP = "b_down"
ultramaster_AnySign =  rbind(ultramaster_AnySign_up, ultramaster_AnySign_down)

# grouping the ultramaster_AnySign based on the significant comparisons
# as function because I do not need those tables and they are creating chaos in my environment
ultramaster_AnySign =  group_direction(ultramaster)

# where the change was only true when senescent were without mitochondria - so when mitochondria made the difference
ultramaster_SPvsMdP_true = subset(ultramaster, ultramaster$sig_SvP == FALSE & ultramaster$sig_MvS == TRUE)
ultramaster_SPvsMdP_true_list = row.names(ultramaster_SPvsMdP_true[,- c(2:10)])

# reverse: lack of mitochondria made senescent cells resemble proliferating cells for those genes
ultramaster_SPvsMdP_false = subset(ultramaster, ultramaster$sig_SvP == TRUE & ultramaster$sig_MvS == FALSE)
ultramaster_SPvsMdP_false_list = row.names(ultramaster_SPvsMdP_false[,- c(2:10)])

# for ultramaster: all of the comparisons significant
ultramaster_AllSign = subset(ultramaster, ultramaster$sig_SvP == TRUE & ultramaster$sig_MvS == TRUE & ultramaster$sig_MvP == TRUE)

# top 10 significant by each comparison
sort_order = order(ultramaster_AllSign[,"p.adj.SvP"], decreasing=FALSE)
ultramaster_AllSign_by_SvP = ultramaster_AllSign[sort_order,]
AllSign_top10_SvP = row.names(ultramaster_AllSign_by_SvP[1:10,])

sort_order = order(ultramaster_AllSign[,"p.adj.MvP"], decreasing=FALSE)
ultramaster_AllSign_by_MvP = ultramaster_AllSign[sort_order,]
AllSign_top10_MvP = row.names(ultramaster_AllSign_by_MvP[1:10,])

sort_order = order(ultramaster_AllSign[,"p.adj.MvS"], decreasing=FALSE)
ultramaster_AllSign_by_MvS = ultramaster_AllSign[sort_order,]
AllSign_top10_MvS = row.names(ultramaster_AllSign_by_MvS[1:10,])

## Making the plots ########
# Plot 1: Volcano plots ####

SP_volcano2 = plot_volcano2(master_SP, 0.001, 2, 
                            my_title="Proliferating vs senescent") + theme(legend.position = "none")
MdP_volcano2 = plot_volcano2(master_MdP, 0.001, 2, 
                             my_title="Proliferating vs senescent without mitochondria") + theme(legend.position = "none")
MdS_volcano2 = plot_volcano2(master_MdS, 0.001, 2, 
                             my_title="Senescent with vs without mitochondria") + theme(legend.position = "right")
SP_volcano2
MdP_volcano2
MdS_volcano2
#save_plot_png(SP_volcano2, plot_name="VolcanoPlot_SP", output_path, height = 700, width = 700) 
#save_plot_png(MdP_volcano2, plot_name="VolcanoPlot_MdP", output_path, height = 700, width = 700) 
#save_plot_png(MdS_volcano2, plot_name="VolcanoPlot_MdS", output_path, height = 700, width = 700) 

# using cowplot to arrange them into multipanel figure
# not facet since the up/down-regulated genes all come from different data tables and it's easier that way

combined_volcano <- plot_grid(SP_volcano2, MdP_volcano2, MdS_volcano2, ncol = 3, labels = c("A.", "B.", "C."), rel_widths = c(0.75,0.75,1))
combined_volcano

save_plot_png(combined_volcano, plot_name="VolcanoPlot", output_path, height = 500, width = 1300) 

# Plot 2: MA plots ####
SP_MA = plot_MA2(master_SP, 0.001, 2, my_title="Proliferating vs senescent") + 
  theme(legend.position = "none")

MdP_MA = plot_MA2(master_MdP, 0.001, 2, my_title="Proliferating vs senescent without mitochondria") + 
  theme(legend.position = "none")

MdS_MA = plot_MA2(master_MdS, 0.001, 2, my_title="Senescent with vs without mitochondria") + 
  theme(legend.position = "right")

SP_MA
MdP_MA
MdS_MA

#save_plot_png(SP_MA, plot_name="MAPlot_SP", output_path, height = 700, width = 700) 
#save_plot_png(MdP_MA, plot_name="MAPlot_MdP", output_path, height = 700, width = 700) 
#save_plot_png(MdS_MA, plot_name="MAPlot_MdS", output_path, height = 700, width = 700) 

#combined_MA <- plot_grid(SP_MA, MdP_MA, MdS_MA, ncol = 3, labels = c("A.", "B.", "C."), rel_widths = c(0.75,0.75,1))
#combined_MA

#save_plot_emf(combined_MA, plot_name="MAPlot", output_path, height = 700, width = 1300) 

# Plot 3: Venn diagram. Overlapping genes ####

# Genes that signify only senescent without mitochondria, not with
# = significant for MdS comparison
SeneNoMito_genes = row.names(master_MdS_sign) #SENEscent NO MITOchondria

# Genes that differentiate senescent from proliferating
Sene_genes = row.names(master_SP_sign)

# Genes that differentiate proliferating and senescent without mitochondria, need the list
Prol_SeneNoMito_genes = row.names(master_MdP_sign)

# with the Venn diagram I will illustrate:
# - genes that are unique for senescent cells and independent of mitochondria
# - genes that are unique for senescent cells and linked to mitochondria (it's absence or presence)
#     - genes that are unique for senescent cells but only in the absence of mito (compensation?)
#     - genes that are unique for senescent cells but only in the presence of mito (mitochondria mediated senescence?)
# - Due to the absence of mitochondria: genes that are expressed like in proliferating but not like senescent (rescue effect?)
# - Due to the absence of mitochondria: genes that are expressed like in senescent but not like proliferating (mitochondria mediated senescence?)

venn_data = list("Prolif vs senes" = Sene_genes, "Prolif vs senes no mito" = Prol_SeneNoMito_genes)
ggp = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, fill=c("azure", "cornflowerblue")) 
ggp
save_plot_png(ggp, "Venn1", output_path)

venn_data = list("Senes no mito vs senes" = SeneNoMito_genes, "Senes no mito vs prolif" = Prol_SeneNoMito_genes)
ggp = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, fill=c("cornsilk", "cornflowerblue"))
ggp
save_plot_png(ggp, "Venn2", output_path)

# 3D Venn diagram
venn_data = list("Prolif vs senes" = Sene_genes, "Prolif vs senes no mito" = Prol_SeneNoMito_genes, "Senes no mito vs senes" = SeneNoMito_genes)
Venn3 = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, fill=c("azure", "cornflowerblue", "cornsilk")) 
Venn3
save_plot_png(Venn3, "Venn_complete", output_path)

# Exploring the data: how many genes in each; later supported by Venn
count_sig_SvP = sum(ultramaster$sig_SvP == "TRUE", na.rm=TRUE)
count_sig_MvS = sum(ultramaster$sig_MvS == "TRUE", na.rm=TRUE)
count_sig_MvP = sum(ultramaster$sig_MvP == "TRUE", na.rm=TRUE)

count_sig_all = sum(ultramaster$sig_MvP == "TRUE" & ultramaster$sig_MvS == "TRUE" & ultramaster$sig_SvP == "TRUE", na.rm=TRUE)
count_sig_MvPS = sum((ultramaster$sig_MvP == "TRUE" | ultramaster$sig_MvS == "TRUE") & ultramaster$sig_SvP != "TRUE", na.rm=TRUE)
count_sig_MvPS2 = sum((ultramaster$sig_MvP == "TRUE" & ultramaster$sig_MvS == "TRUE") & ultramaster$sig_SvP != "TRUE", na.rm=TRUE)

# Hypergeometric test for p value ####
#. . . 1: for the sene regardless of mito and such ####
hypometric_SvP_MdP = hypometric(group1 = Venn3[["data"]][["fitted.values"]][["Prolif vs senes no mito"]],
                                group2 = Venn3[["data"]][["fitted.values"]][["Prolif vs senes"]],
                                both_groups = Venn3[["data"]][["fitted.values"]][["Prolif vs senes&Prolif vs senes no mito"]],
                                total_table = ultramaster)

#. . . 2: for the sene regardless of mito ####
hypometric_MdS_MdP = hypometric(group1 = Venn3[["data"]][["fitted.values"]][["Prolif vs senes no mito"]],
                                group2 = Venn3[["data"]][["fitted.values"]][["Senes no mito vs senes"]],
                                both_groups = Venn3[["data"]][["fitted.values"]][["Prolif vs senes no mito&Senes no mito vs senes"]],
                                total_table = ultramaster)

# Plot 4: Scatter plot of the log2fold change for both groups ####
#. . . 1: for the sene regardless of mito and such ####

scatter1a = ggplot(ultramaster, aes(x= log2fold.SvP, y= log2fold.MvP)) + geom_point(aes(fill = p.SvP)) +
  my_theme_basic() +
  scale_color_gradient(low = "blue", high = "red") + 
  labs(title="Scatter plot", x="Log2fold change: Prolif vs Sene", y="Log2fold change: Prolif vs Sene no mito")
scatter1a

scatter1 = ggplot(ultramaster, aes(x= log2fold.SvP, y= log2fold.MvP)) + geom_point(aes(fill = p.SvP, size = p.MvP, alpha = p.MvP), shape=21) +
  theme_bw() +
  scale_fill_gradient(low = "skyblue4", high = "white") + 
  scale_size_continuous(range = c(3, 1)) +
  scale_alpha_continuous(range = c(0.5, 0.2)) +
  labs(title="Comparison of expression change matrix by groups", x="Log2fold change: Prolif vs Sene", y="Log2fold change: Prolif vs Sene no mito")
scatter1

#save_plot_png(scatter1, "scatter_SvP_MvP", output_path)


#. . . 2: for the sene regardless of mito ####
scatter2 = ggplot(ultramaster, aes(x= log2fold.MvS, y= log2fold.MvP)) + geom_point(aes(fill = p.MvS, size = p.MvP, alpha = p.MvP), shape=21) +
  theme_bw() +
  scale_fill_gradient(low = "skyblue4", high = "white") + 
  scale_size_continuous(range = c(3, 1)) +
  scale_alpha_continuous(range = c(0.5, 0.2)) +
  labs(title="Comparison of expression change matrix by groups", x="Log2fold change: Sene vs Sene no mito", y="Log2fold change: Prolif vs Sene no mito")
scatter2

#save_plot_png(scatter2, "scatter_MvS_MvP", output_path)

#making a multipanel plot
scatter1 = scatter1 + theme(legend.position = "none")
scatter2 = scatter2 + theme(legend.position = "right")

combined_scatter <- plot_grid(scatter1, scatter2, ncol=2, labels = c("B.", "C."), rel_widths = c(0.8, 1))
combined_scatter
combined_VannScatter <- plot_grid(Venn3, combined_scatter, ncol = 1, labels = c("A.", ""), rel_heights = c(1,1))
combined_VannScatter

save_plot_png(combined_VannScatter, "CombinedVannScatter", output_path, 1000, 1000)

# .  correlations ####
correlation_SvP_MvP = cor.test(ultramaster$log2fold.SvP, ultramaster$log2fold.MvP)
correlation_MvS_MvP = cor.test(ultramaster$log2fold.MvS, ultramaster$log2fold.MvP)
correlation_SvP_MvP_cor = correlation_SvP_MvP[["estimate"]][["cor"]]
correlation_MvS_MvP_cor = correlation_MvS_MvP[["estimate"]][["cor"]]

# Plot 5: PCA plot ####
# they are the same 
ggp = plot_pca1(em_SP_scaled)
ggp
save_plot_png(ggp, "PCA_all", output_path)

# Plot 6: Expression density plot ####
ggp = plot_expressions(em, 1:9)
ggp
save_plot_png(ggp, "ExpressionDensities", output_path, 200, 1000)

# Plot 7: Box and violin plots
em_symbols = ultramaster[,2:10]
top10down_MdP
#  [1] "TPX2"   "BIRC5"  "KIF4A"  "SFRP1"  "COL1A1" "KIF20A" "DHCR24" "DLGAP5" "LDLR"   "CCNB1" 
# SFRP1, COL1A1, DHCR24, LDLR interesting

top10up_MdP
#[1] "SLC2A3"  "PGK1"    "NDRG1"   "ENO2"    "PFKFB4"  "SLC2A1"  "PPP1R3C" "MXI1"    "LDHA"    "IGFBP3" 
# SLC2A1, PPP1R3C, MXI1

AllSign_top10_SvP
#  [1] "ARRDC4"  "HAS2"    "RIPOR3"  "IGFBP3"  "A2M"     "C3"      "FOSB"    "PLXDC2"  "COLEC10" "LMOD1"  

AllSign_top10_MvP
#  [1] "IGFBP3"  "COL3A1"  "MT-ND5"  "MT-ND1"  "MT-RNR2" "APLN"    "AK4"     "MT-ND6"  "MT-RNR1" "FABP3"  

AllSign_top10_MvS
#  [1] "RIPOR3"  "IGFBP3"  "AK4"     "COL3A1"  "MT-ND6"  "MT-ND5"  "MT-ND1"  "MT-RNR2" "MT-RNR1" "APLN"   

gene = "SLC7A5"

gene_data = em_symbols[gene,] 
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP

# changing the order  of coloums
gene_data = gene_data[c(2,1)]   
names(gene_data) = c("sample_group", "expression")

# ordering the data in rows by sample group
levels(gene_data$sample_group) #current data order is NULL
gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
my_title = paste(gene,"Expression", sep = " ")

ggp = ggplot(gene_data, aes(x=sample_group, y=expression)) + 
  geom_violin(colour="skyblue4", fill="white") +
  my_theme_darkblue() +
  labs(title=my_title, x="Sample group", y="Expression")
ggp

save_plot_png(ggp, my_title, output_path, 500, 500)

# Plot 7: top10 boxplots ####
# MvS 15, MvP 18 SvP 21
em_scaled = em_SP_scaled
top10_MvS = top10_boxplot_scaled(ultramaster, 15, my_title="Expression of the top 10 most significant genes: MvS")
top10_MvS
save_plot_png(top10_MvS, "top10_MvS", output_path, 400, 1000)

top10_MvP = top10_boxplot_scaled(ultramaster, 18, my_title="Expression of the top 10 most significant genes: MvP")
top10_MvP
save_plot_png(top10_MvP, "top10_MvP", output_path, 400, 1000)

top10_SvP = top10_boxplot_scaled(ultramaster, 21, my_title="Expression of the top 10 most significant genes: SvP")
top10_SvP
save_plot_png(top10_SvP, "top10_SvP", output_path, 400, 1000)

top10_SvP = top10_SvP + theme(legend.position="right")
top10_MvP = top10_MvP + theme(legend.position="right")
top10_MvS = top10_MvS + theme(legend.position="right")

combined_top10 <- plot_grid(top10_SvP, top10_MvP, top10_MvS, ncol=1, labels = c("A.", "B.", "C."))
combined_top10

save_plot_png(combined_top10, "combined_top10", output_path, 700, 1000)

# Plot 8: Heatmaps ####
heatmap_SP = plot_heatmap(master_SP, 2:10)
heatmap_SP
save_plot_png(heatmap_SP, "heatmap_SvP", output_path, 1000, 600)

heatmap_MdS = plot_heatmap(master_MdS, 2:10)
heatmap_MdS
save_plot_png(heatmap_MdS, "heatmap_MdS", output_path, 1000, 600)

heatmap_MdP = plot_heatmap(master_MdP, 2:10)
heatmap_MdP
save_plot_png(heatmap_MdP, "heatmap_MdP", output_path, 1000, 600)

heatmap_SP = heatmap_SP + theme(legend.position = "bottom")
heatmap_MdP = heatmap_MdP  + theme(legend.position = "bottom")
heatmap_MdS = heatmap_MdS + theme(legend.position = "bottom")

combined_heatmaps <- plot_grid(heatmap_SP, heatmap_MdS, heatmap_MdP, ncol=3, labels = c("A.", "B.", "C."))
combined_heatmaps

save_plot_png(combined_heatmaps, "combined_heatmaps", output_path, 1200, 1000)

# Plot 9: ORA
sig_ENSEMBL_SP = master_SP_sign[,1]
sig_ENSEMBL_MdP = master_MdS_sign[,1]
sig_ENSEMBL_MdS = master_MdS_sign[,1]

sig_ENTREZ_SP = bitr(sig_ENSEMBL_SP, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_ENTREZ_MdP = bitr(sig_ENSEMBL_MdP, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_ENTREZ_MdS = bitr(sig_ENSEMBL_MdS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ora_results_SP = enrichGO(gene = sig_ENTREZ_SP$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdP = enrichGO(gene = sig_ENTREZ_MdP$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdS = enrichGO(gene = sig_ENTREZ_MdS$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)

ggp = dotplot(ora_results_SP, showCategory=10)+ my_theme_basic() + labs(x="Gene ratio")
ggp
save_plot_png(ggp, "ora_SP", output_path, 600, 600)

ggp = dotplot(ora_results_MdP, showCategory=10)+ my_theme_basic() + labs(x="Gene ratio")
ggp
save_plot_png(ggp, "ora_MdP", output_path, 600, 600)

ggp = dotplot(ora_results_MdS, showCategory=10)+ my_theme_basic() + labs(x="Gene ratio")
ggp
save_plot_png(ggp, "ora_MdS", output_path, 600, 600)

# separately for upregulated and downregulated genes
sig_ENSEMBL_SP_up = subset(master_SP_sign, p.adj<0.001 & log2fold>1)[,1]
sig_ENSEMBL_SP_down = subset(master_SP_sign, p.adj<0.001 & log2fold< -1)[,1]

sig_ENSEMBL_MdP_up = subset(master_MdP_sign, p.adj<0.001 & log2fold>1)[,1]
sig_ENSEMBL_MdP_down = subset(master_MdP_sign, p.adj<0.001 & log2fold< -1)[,1]

sig_ENSEMBL_MdS_up = subset(master_MdS_sign, p.adj<0.001 & log2fold>1)[,1]
sig_ENSEMBL_MdS_down = subset(master_MdS_sign, p.adj<0.001 & log2fold< -1)[,1]
###
sig_ENTREZ_SP_up = bitr(sig_ENSEMBL_SP_up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_ENTREZ_SP_down = bitr(sig_ENSEMBL_SP_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

sig_ENTREZ_MdP_up = bitr(sig_ENSEMBL_MdP_up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_ENTREZ_MdP_down = bitr(sig_ENSEMBL_MdP_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

sig_ENTREZ_MdS_up = bitr(sig_ENSEMBL_MdS_up, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
sig_ENTREZ_MdS_down = bitr(sig_ENSEMBL_MdS_down, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
###
ora_results_SP_up = enrichGO(gene = sig_ENTREZ_SP_up$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdP_up = enrichGO(gene = sig_ENTREZ_MdP_up$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdS_up = enrichGO(gene = sig_ENTREZ_MdS_up$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)

ora_results_SP_down = enrichGO(gene = sig_ENTREZ_SP_down$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdP_down = enrichGO(gene = sig_ENTREZ_MdP_down$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
ora_results_MdS_down = enrichGO(gene = sig_ENTREZ_MdS_down$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)

ora_SP_up = dotplot(ora_results_SP_up, showCategory=10)+ theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes upregulated in Sene vs Prolif", x="Gene ratio")
ora_SP_up
save_plot_png(ora_SP_up, "ora_SP_up", output_path, 600, 600)

ora_MdP_up = dotplot(ora_results_MdP_up, showCategory=10)+ theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes upregulated in Sene no mito vs Prolif", x="Gene ratio")
ora_MdP_up
save_plot_png(ora_MdP_up, "ora_MdP_up", output_path, 600, 600)

ora_MdS_up = dotplot(ora_results_MdS_up, showCategory=10)+ theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes upregulated in Sene no mito vs Sene", x="Gene ratio")
ora_MdS_up
save_plot_png(ora_MdS_up, "ora_MdS_up", output_path, 600, 600)

ora_SP_down = dotplot(ora_results_SP_down, showCategory=10)+ theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes downregulated in Sene vs Prolif", x="Gene ratio")
ora_SP_down
save_plot_png(ora_SP_down, "ora_SP_down", output_path, 600, 600)

ora_MdP_down = dotplot(ora_results_MdP_down, showCategory=10)+theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes downregulated in Sene no mito vs Prolif", x="Gene ratio")
ora_MdP_down
save_plot_png(ora_MdP_down, "ora_MdP_down", output_path, 600, 600)

ora_MdS_down = dotplot(ora_results_MdS_down, showCategory=10)+ theme_bw() + theme(legend.position = "right") + labs(title="ORA of genes downregulated in Sene no mito vs Sene", x="Gene ratio")
ora_MdS_down
save_plot_png(ora_MdS_down, "ora_MdS_down", output_path, 600, 600)

ora_up_combined <- plot_grid(ora_SP_up, ora_MdP_up, ora_MdS_up, ncol=3, labels = c("A.", "B.", "C."))
ora_up_combined

ora_down_combined <- plot_grid(ora_SP_down, ora_MdP_down, ora_MdS_down, ncol=3, labels = c("D.", "E.", "F."))
ora_down_combined

ora_all_combined <- plot_grid(ora_up_combined, ora_down_combined, ncol=1)
ora_all_combined

save_plot_png(ora_all_combined, "combined_ora_updown", output_path, 900, 1800)

# Plot 11: GSEA plot ####

#gsea_SP = plot_gsea(master_SP)
#gsea_SP
#save_plot_png(gsea_SP, "gsea_SP", output_path, 1000, 1000)
#
#gsea_MdS = plot_gsea(master_MdS)
#gsea_MdS
#save_plot_png(gsea_MdS, "gsea_MdS", output_path, 1000, 1000)
#
#gsea_MdP = plot_gsea(master_MdP)
#gsea_MdP
#save_plot_png(gsea_MdP, "gsea_MdP", output_path, 1000, 1000)

# Plot 10: enriched genes box plot

gene_sets_SP = ora_results_SP$geneID
description_SP = ora_results_SP$Description
p.adj_SP = ora_results_SP$p.adjust

ora_results_table_SP = data.frame(cbind(description_SP, gene_sets_SP, p.adj_SP))
row.names(ora_results_table_SP) = ora_results_table_SP[,1]
ora_results_table_SP = ora_results_table_SP[,2:3]

# getting the most enriched genes as the list
enriched_gene_set_SP = as.character(ora_results_table_SP[1,1]) 
candidate_genes_SP = unlist(strsplit(enriched_gene_set_SP, "/"))
candidate_genes_SP

candidate_genes_table_SP = data.frame(candidate_genes_SP) 
names(candidate_genes_table_SP) = "gene" 

# load the database, 9606 is human, 10090 is mouse. Google the function for more information. 
string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="") 

# map our genes to the database 
# "gene" means gene symbols. We could put "ENSEMBL" here if needed. 
string_mapped = string_db$map(candidate_genes_table_SP, "gene", removeUnmappedRows = TRUE )

# now we plot 
new_ggp = string_db$plot_network(string_mapped) 
