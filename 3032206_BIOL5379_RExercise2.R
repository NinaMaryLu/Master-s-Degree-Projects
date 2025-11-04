library(ggplot2)
library(ggrepel)

# to call them run
#source("C:/Users/User/Documents/UoG 2024 2025/Data Vis in R/omics_functions.R") 

# ORGANISING TABLES ####

# creating a master table out of: the de table, annotation table, em table

merge_master = function(de, annotations, em)
{
  # Join all of the tables together to create a “Master” table
  master_temp = merge(em,annotations, by.x=0, by.y=0)
  master = merge(master_temp, de, by.x=1, by.y=0)
  return(master)
} 

# rename and delete duplicate columns for ENSEMBL
row_rename = function(de_table, column_index, name)
{
  # assigning the row names
  row.names(de_table) = de_table[,column_index]
  #renaming the column "str"
  names(de_table)[column_index] = name
  return(de_table)
}

# making the mlog10p column
append_mlog2p = function(de_table)
{
  log_values = -log10(de_table$p)         # calc the p values as vector
  de_table$mlog10p = log_values           # adding the p values to the table
  return(de_table)
}

# making the expression_mean column
#   de_table - table name, range - range of columns with expression data
append_expression_mean = function(de_table, range)
{
  de_table_means = rowMeans(de_table[range])
  de_table$expression_mean = de_table_means
  return(de_table)
}

issignificant = function(de_table)
  # evaluates significance based on columns "p.adj" and "log2fold"
  # returns significance classification as vector
{
  significance = as.factor(de_table$p.adj < 0.001 & abs(de_table$log2fold) > 2.0) # as.factor changes the variable to categorical
  return(significance)
}


# sorting by p value
#   requirement: p column has to be names "p.adj"
sort_p = function(de_table)
{
  sort_order = order(de_table[,"p.adj"], decreasing=FALSE)
  de_table = de_table[sort_order,]
}

# scaling and omitting the NA values
#   de_table - name of the table; start, end - range for the scaling
scale_omit = function(de_table, range)
{
  table_scaled = data.frame(t(scale(data.frame(t(de_table[,range])))))
  table_scaled = na.omit(table_scaled)
  return(table_scaled)
}

# returns the expression table of only significant genes
#   requirements: p column has to be names "p.adj" and log2 fold "log2fold"
table_sign = function(de_table, range)
{
  #Make a list of significant genes using subset
  table_sig = subset(de_table, p.adj < 0.001 & abs(log2fold) >2)
  # subset the table to have expression data only
  expression = table_sig[,range]
  # without the NA values
  expression = na.omit(expression)
  return(expression)
}

# returns the scaled expression table of only significant genes
table_sign_scaled = function(de_table, range)
{
  #Make a list of significant genes using subset
  table_sig = subset(de_table, p.adj < 0.001 & abs(log2fold) >2)
  # subset the table to have expression data only
  table_sig = table_sig[,range]
  expression = na.omit(data.frame(t(scale(data.frame(t(table_sig))))))
  return(expression)
}

top10up = function(de_table, p_threshold=0.001, log2fold_threshold = 2)
{
  de_table_sig_up = subset(de_table, de_table$p.adj<p_threshold &  de_table$log2fold > log2fold_threshold)
  de_table_sig_up = sort_p(de_table_sig_up)
  top10up = row.names(de_table_sig_up[1:10,])
  return(top10up)
}

top10down = function(de_table, p_threshold=0.001, log2fold_threshold = 2)
{
  de_table_sig_down = subset(de_table, de_table$p.adj<p_threshold &  de_table$log2fold < log2fold_threshold)
  de_table_sig_down = sort_p(de_table_sig_down)
  top10down = row.names(de_table_sig_down[1:10,])
  return(top10down)
}

top10 = function(de_table, p_threshold=0.001, log2fold_AbsThreshold = 2)
{
  de_table_sig = subset(de_table, de_table$p.adj<p_threshold &  abs(de_table$log2fold) > log2fold_threshold)
  de_table_sig = sort_p(de_table_sig)
  top10 = row.names(de_table_sig[1:10,])
  return(top10)
}

group_direction = function(de_table)
{
  ultramaster_AnySign_onlySP = subset(de_table, de_table$sig_SvP == TRUE & de_table$sig_MvS == FALSE & de_table$sig_MvP == FALSE)
  ultramaster_AnySign_onlyMvS = subset(de_table, de_table$sig_SvP == FALSE & de_table$sig_MvS == TRUE & de_table$sig_MvP == FALSE)
  ultramaster_AnySign_onlyMvP = subset(de_table, de_table$sig_SvP == FALSE & de_table$sig_MvS == FALSE & de_table$sig_MvP == TRUE)
  ultramaster_AnySign_SP_MvS = subset(de_table, de_table$sig_SvP == TRUE & de_table$sig_MvS == TRUE & de_table$sig_MvP == FALSE)
  ultramaster_AnySign_SP_MvP = subset(de_table, de_table$sig_SvP == TRUE & de_table$sig_MvS == FALSE & de_table$sig_MvP == TRUE)
  ultramaster_AnySign_MvS_MvSP = subset(de_table, de_table$sig_SvP == FALSE & de_table$sig_MvS == TRUE & de_table$sig_MvP == TRUE)
  
  ultramaster_AnySign_onlySP$group = "only_SP"
  ultramaster_AnySign_onlyMvS$group = "only_MvS"
  ultramaster_AnySign_onlyMvP$group = "only_MvP"
  ultramaster_AnySign_SP_MvS$group = "SP_MvS"
  ultramaster_AnySign_SP_MvP$group = "SP_MvP"
  ultramaster_AnySign_MvS_MvSP$group = "MvS_MvP"
  
  ultramaster_AnySign =  rbind(ultramaster_AnySign_onlySP, ultramaster_AnySign_onlyMvS, ultramaster_AnySign_onlyMvP, ultramaster_AnySign_SP_MvS, ultramaster_AnySign_SP_MvP, ultramaster_AnySign_MvS_MvSP)
  return(ultramaster_AnySign)
}

# creates a table with a new column of "direction" with groups : 
#                     "a" for"non-significant" / "b" for "downregulated" / "c" for "upregulated"
#   requirements: table column named "sig" with boolean values TRUE / FALSE for significance
#                 column names: "p.adj" and "log2fold"
group_sign_full = function(de_table)
{
  #ordering the data
  de_table$sig = factor(de_table$sig, levels = c("TRUE", "FALSE"))  # changing the order
  
  # creating 3 groups: non significant, down-reg and up-reg
  de_table_non_sig = subset(de_table, sig == FALSE)
  de_table_sig_up = subset(de_table, de_table$p.adj<0.001 &  de_table$log2fold > 2.0)
  de_table_sig_down = subset(de_table, de_table$p.adj<0.001 &  de_table$log2fold < -2.0)
  
  de_table_non_sig$direction = "a_nonsig"
  de_table_sig_up$direction = "b_up"
  de_table_sig_down$direction = "c_down"
  
  de_table =  rbind(de_table_non_sig, de_table_sig_up, de_table_sig_down)
  return(de_table)
}

group_sign_only = function(de_table)
{
  # creating 2 groups: down-reg and up-reg
  de_table_sig_up = subset(de_table, de_table$log2fold > 2.0)
  de_table_sig_down = subset(de_table, de_table$log2fold < -2.0)
  
  de_table_sig_up$direction = "a_up"
  de_table_sig_down$direction = "b_down"
  
  de_table =  rbind(de_table_sig_up, de_table_sig_down)
  return(de_table)
}


# themes
my_theme_basic = function()
{
  my_theme_basic = theme( 
  plot.title = element_text(size=15, face = "bold", margin = margin(t = 10, b = 10)), 
  axis.text.x = element_text(size=9, colour="slategrey"), 
  axis.text.y = element_text(size=9, colour="slategrey"), 
  axis.title.x = element_text(size=12, colour="slategrey", face="bold"), 
  axis.title.y = element_text(size=12, colour="slategrey", face="bold"),
  legend.title = element_text(size=10, colour="slategrey", face="bold"),
  legend.text = element_text(size=13, colour="grey1"),
  legend.background = element_rect(fill="grey92", colour="grey88"),
  legend.position = "right",
  #legend.justification = c(0, 0),
  legend.spacing.x = unit(0, 'cm'),
  panel.grid.major = element_line(linetype="solid"),
  panel.grid.minor = element_line(linetype="solid"))
  return(my_theme_basic)
}

my_theme_darkblue = function()
{
  my_theme_darkblue = theme( 
  plot.title = element_text(size=18, face = "bold", colour = "white", margin = margin(t = 10, b = 10)), 
  axis.text.x = element_text(size=9, colour="white"), 
  axis.text.y = element_text(size=9, colour="white"), 
  axis.title.x = element_text(size=12, colour="white", face="bold"), 
  axis.title.y = element_text(size=12, colour="white", face="bold"),
  legend.axis.line = element_line(colour = "azure2"),
  legend.title = element_text(size=10, colour="grey30", face="bold"),
  legend.text = element_text(size=8, colour="grey30"),
  legend.background = element_rect(fill="azure1", colour="azure3"),
  legend.position = "right",
  legend.justification = c(0, 1),
  panel.grid.major = element_line(linetype="solid", colour="azure2"),
  panel.grid.minor = element_line(linetype="solid", colour="azure2"),
  plot.background = element_rect(fill = "skyblue4", colour = "slategrey"), # on the edge
  panel.background = element_rect(fill = "azure1", colour = "slategrey"))# panel with the plot
  return(my_theme_darkblue)
}

my_theme_beige = function()
{
  my_theme_beige = theme( 
    plot.title = element_text(size=18, face = "bold", colour = "white", margin = margin(t = 10, b = 10)), 
    axis.text.x = element_text(size=9, colour="white"), 
    axis.text.y = element_text(size=9, colour="white"), 
    axis.title.x = element_text(size=12, colour="white", face="bold"), 
    axis.title.y = element_text(size=12, colour="white", face="bold"),
    legend.axis.line = element_line(colour = "burlywood4"),
    legend.title = element_text(size=10, colour="grey30", face="bold"),
    legend.text = element_text(size=8, colour="grey30"),
    legend.background = element_rect(fill="bisque3", colour="burlywood4"),
    legend.position = "right",
    legend.justification = c(0, 1),
    panel.grid.major = element_line(linetype="solid", colour="bisque3"),
    panel.grid.minor = element_line(linetype="solid", colour="bisque3"),
    plot.background = element_rect(fill = "burlywood4", colour = "chocolate4"), # on the edge
    panel.background = element_rect(fill = "beige", colour = "chocolate4"))# panel with the plot
  return(my_theme_beige)
}


#geom_point(labs(colour = "Cylinders"))) adds a specific title on teh legend

# PLOTTING GRAPHS ####

# save current plot

save_plot_png = function(plot=ggp, plot_name, output_path, height = 1000, width = 1000)
{
  file_path = paste(output_path, plot_name, ".png", sep="")
  png(file_path, height = height, width = width) 
  print(plot)
  dev.off()
}

save_plot_emf = function(plot=ggp, plot_name, output_path, height = 1000, width = 1000)
{
  file_path = paste(output_path, plot_name, ".emf", sep="")
  emf(file_path, height = height, width = width) 
  print(plot)
  dev.off()
}


# with data points as directions
# row index has to be in the first column
plot_volcano1 = function(de_table, p_threshold=0.001, fold_threshold=2)
{ 
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) >fold_threshold)  
  de_table_sig_up = subset(de_table_sig, log2fold>0) #log2fold>0 means it's upregulated
  de_table_sig_up = sort_p(de_table_sig_up)
  de_table_sig_down = subset(de_table_sig, log2fold<0) #log2fold<0 means it's downregulated
  de_table_sig_down = sort_p(de_table_sig_down)
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 =  de_table_sig_down[1:10,]
  
  label_nonsig = "non-significant"
  fill_nonsig = "lightgrey"
  colour_nonsig = "grey40"
  
  label_sigup = "upregulated"
  fill_sigup = "lightgreen"
  colour_sigup = "darkolivegreen3"
  colour_text_sigup = "chartreuse4"
  
  label_sigdown = "downregulated"
  fill_sigdown = "seashell"
  colour_sigdown = "tomato"
  colour_text_sigdown = "darkred"
  
  colour_threshold = "grey70"
  
  shape_nonsig = 21
  shape_sigup = 24
  shape_sigdown = 25

  
  ggp = ggplot(de_table, aes(x=log2fold, y=mlog10p, colour=direction)) + 
    geom_point(size = 1, fill = fill_nonsig, shape = shape_nonsig) + 
    geom_point(aes(colour="b_up"), fill=fill_sigup, data=de_table_sig_up, shape=shape_sigup) + 
    geom_point(aes(colour="c_down"), fill=fill_sigdown, data=de_table_sig_down, shape=shape_sigdown) +
    labs(title="Volcano plot", x="Log2fold change", y="-Log10 p") +
    my_theme_basic() + 
    geom_vline(xintercept = -1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_vline(xintercept = 1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour = colour_threshold, size=0.5) +
    xlim(c(-10,10)) +
    ylim(c(0,400)) +
    geom_text_repel(data=de_table_sig_up_top10, aes(label=row.names(de_table_sig_up_top10)), size=3, colour=colour_text_sigup) +
    geom_text_repel(data=de_table_sig_down_top10, aes(label=row.names(de_table_sig_down_top10)), size=3, colour=colour_text_sigdown) +
    scale_colour_manual(values = c(colour_nonsig, colour_sigup, colour_sigdown), labels=c(label_nonsig, label_sigup, label_sigdown)) # legend
  ggp 
  return(ggp) 
}

# unified shape of the datapoints
plot_volcano2 = function(de_table, p_threshold=0.001, fold_threshold=2, my_title="Volcano plot")
{ 
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) >fold_threshold)  
  de_table_sig_up = subset(de_table_sig, log2fold>0) #log2fold>0 means it's upregulated
  de_table_sig_up = sort_p(de_table_sig_up)
  de_table_sig_down = subset(de_table_sig, log2fold<0) #log2fold<0 means it's downregulated
  de_table_sig_down = sort_p(de_table_sig_down)
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 =  de_table_sig_down[1:10,]

  label_nonsig = "non-significant"
  fill_nonsig = "lightgrey"
  colour_nonsig = "grey40"
  
  label_sigup = "upregulated"
  fill_sigup = "lightgreen"
  colour_sigup = "darkolivegreen3"
  colout_text_sigup = "chartreuse4"
  
  label_sigdown = "downregulated"
  fill_sigdown = "seashell"
  colour_sigdown = "tomato"
  colout_text_sigdown = "darkred"
  
  colour_threshold = "grey70"
  
  shape_nonsig = 21
  shape_sigup = 21
  shape_sigdown = 21
  
  ggp = ggplot(de_table, aes(x=log2fold, y=mlog10p, colour=direction)) + 
    geom_point(size = 1, fill = fill_nonsig, shape = shape_nonsig) + 
    geom_point(aes(colour="b_up"), fill=fill_sigup, data=de_table_sig_up, shape=shape_sigup) + 
    geom_point(aes(colour="c_down"), fill=fill_sigdown, data=de_table_sig_down, shape=shape_sigdown) +
    labs(title=my_title, x="Log2fold change", y="-Log10 p") +
    my_theme_basic() +
    geom_vline(xintercept = -1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_vline(xintercept = 1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour = colour_threshold, size=0.5) +
    xlim(c(-10,10)) +
    ylim(c(0,400)) +
    geom_text_repel(data=de_table_sig_up_top10, aes(label=row.names(de_table_sig_up_top10)), size=3, colour=colout_text_sigup) +
    geom_text_repel(data=de_table_sig_down_top10, aes(label=row.names(de_table_sig_down_top10)), size=3, colour=colout_text_sigdown) +
    scale_colour_manual(values = c(colour_nonsig, colour_sigup, colour_sigdown), labels=c("non-significant", "upregulated", "downregulated")) # legend
  ggp 
  return(ggp) 
}

# with datapoint shape as direction
plot_MA1 = function(de_table, p_threshold=0.001, fold_threshold=2)
{
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) >fold_threshold)  
  de_table_sig_up = subset(de_table_sig, log2fold>0) #log2fold>0 means it's upregulated
  de_table_sig_up = sort_p(de_table_sig_up)
  de_table_sig_down = subset(de_table_sig, log2fold<0) #log2fold<0 means it's downregulated
  de_table_sig_down = sort_p(de_table_sig_down)
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 =  de_table_sig_down[1:10,]
  de_table = group_significance(de_table)
  
  label_nosig = "non-significant"
  fill_nonsig = "lightgrey"
  colour_nonsig = "grey40"
  
  label_sigup = "upregulated"
  fill_sigup = "lightgreen"
  colour_sigup = "darkolivegreen3"
  colour_text_sigup = "chartreuse4"
  
  label_sigdown = "downregulated"
  fill_sigdown = "seashell"
  colour_sigdown = "tomato"
  colour_text_sigdown = "darkred"
  
  colour_threshold = "grey2"
  
  shape_nonsig = 21
  shape_sigup = 24
  shape_sigdown = 25
  
  ggp = ggplot(de_table, aes(x=log10(expression_mean), y=log2fold, colour=direction)) + 
    geom_point(size = 1, fill = fill_nonsig, shape = shape_nonsig) + 
    geom_point(aes(colour="b"), fill=fill_sigup, data=de_table_sig_up, shape=shape_sigup) + 
    geom_point(aes(colour="c"), fill=fill_sigdown, data=de_table_sig_down, shape=shape_sigdown) +
    labs(title="MA plot", x="Log10 Mean expression", y="Log2fold change") +
    my_theme_basic() + 
    geom_hline(yintercept = -1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_hline(yintercept = 1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_text_repel(data=de_table_sig_up_top10, aes(label=row.names(de_table_sig_up_top10)), size=3, colour=colour_text_sigup) +
    geom_text_repel(data=de_table_sig_down_top10, aes(label=row.names(de_table_sig_down_top10)), size=3, colour=colour_text_sigdown) +
    scale_colour_manual(values = c(colour_nonsig, colour_sigup, colour_sigdown), labels=c(label_nosig, label_sigup, label_sigdown)) # legend
  ggp
  return(ggp)
}

# with unified datapoint shape
plot_MA2 = function(de_table, p_threshold=0.001, fold_threshold=2, my_title="MA plot")
{
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) >fold_threshold)  
  de_table_sig_up = subset(de_table_sig, log2fold>0) #log2fold>0 means it's upregulated
  de_table_sig_up = sort_p(de_table_sig_up)
  de_table_sig_down = subset(de_table_sig, log2fold<0) #log2fold<0 means it's downregulated
  de_table_sig_down = sort_p(de_table_sig_down)
  de_table_sig_up_top10 = de_table_sig_up[1:10,]
  de_table_sig_down_top10 =  de_table_sig_down[1:10,]

  label_nosig = "non-significant"
  fill_nonsig = "lightgrey"
  colour_nonsig = "grey40"
  
  label_sigup = "upregulated"
  fill_sigup = "lightgreen"
  colour_sigup = "darkolivegreen3"
  colout_text_sigup = "chartreuse4"
  
  label_sigdown = "downregulated"
  fill_sigdown = "seashell"
  colour_sigdown = "tomato"
  colout_text_sigdown = "darkred"
  
  colour_threshold = "grey2"
  
  shape_nonsig = 21
  shape_sigup = 21
  shape_sigdown = 21
  
  ggp = ggplot(de_table, aes(x=log10(expression_mean), y=log2fold, colour=direction)) + 
    geom_point(size = 1, fill = fill_nonsig, shape = shape_nonsig) + 
    geom_point(aes(colour="b_up"), fill=fill_sigup, data=de_table_sig_up, shape=shape_sigup) + 
    geom_point(aes(colour="c_down"), fill=fill_sigdown, data=de_table_sig_down, shape=shape_sigdown) +
    labs(title=my_title, x="Log10 Mean expression", y="Log2fold change") +
    my_theme_basic() + 
    geom_hline(yintercept = -1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_hline(yintercept = 1, linetype="dashed", colour = colour_threshold, size=0.5) +
    geom_text_repel(data=de_table_sig_up_top10, aes(label=row.names(de_table_sig_up_top10)), size=3, colour=colout_text_sigup) +
    geom_text_repel(data=de_table_sig_down_top10, aes(label=row.names(de_table_sig_down_top10)), size=3, colour=colout_text_sigdown) +
    scale_colour_manual(values = c(colour_nonsig, colour_sigup, colour_sigdown), labels=c(label_nosig, label_sigup, label_sigdown)) # legend
  ggp
  return(ggp)
}

# unable to add your own names to the conditions with this function
plot_venn = function(de_table, condition1, condition2)
{
  c1_count = row.names(subset(de_table, condition1)) 
  c2_count = row.names(subset(de_table, condition2))
  
  # creating the Venn diagram
  venn_data = list("Condition 1" = c1_count, "Condition 2" = c2_count)
  plot(euler(venn_data, shape = "ellipse"), quantities = TRUE) 
}

plot_pca1 = function(de_table, my_title="Sample group")
{
  pca_matrix = as.matrix(sapply(de_table, as.numeric))
  pca = prcomp(t(pca_matrix))
  pca_coordinates = data.frame(pca$x)
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour = row.names(pca_coordinates))) + 
    geom_point(size=3) + 
    labs(title = "PCA graph") + 
    geom_text_repel(aes(label=row.names(pca_coordinates))) +
    my_theme_basic() + 
    guides(colour=guide_legend(title=my_title)) +
    scale_color_manual(values=c("tan", "tomato3","sienna4", "thistle", "plum3", "purple1", "skyblue", "steelblue", "navyblue" )) 
  ggp
  return(ggp)
}

plot_expression = function(de_table, index)
{
  group = de_table[,index]
  ggp = ggplot(de_table, aes(x=log10(group))) + 
    geom_density(colour = "magenta", size = 1, alpha = 0.5)  + 
    labs(title = "Expression density plot")
  ggp
}

plot_expressions = function(de_table, range)
{
  data = de_table[,range]
  ncol = ncol(data)
  data$gene_name = row.names(data)
  data.m = melt(data) # message: using gene_names as id variables
  
  ggp = ggplot(data.m, aes(x=log10(value))) + # value stands for value in the table
    geom_density(colour = "seagreen", size = 1, alpha = 0.5)  + 
    labs(title = "Expression density plot") + 
    my_theme_darkblue() +
    facet_wrap(~variable, ncol=ncol) # variable stands for the column name
  ggp
}

plot_boxplot = function(de_table, gene, sample_table)
{
  gene_data = de_table[gene,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = sample_table$SAMPLE_GROUP
  
  # changing the order  of coloums
  gene_data = gene_data[c(2,1)]   
  names(gene_data) = c("sample_group", "expression")
  
  # ordering the data in rows by sample group
  levels(gene_data$sample_group) #current data orderis NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("gut","duct","node")) 
  my_title = paste(gene,"Expression", sep = " ")
  
  ggp = ggplot(gene_data, aes(x=sample_group, y=expression)) + 
    geom_boxplot(colour="skyblue4", fill="white") +
    my_theme_darkblue() +
    labs(title=my_title, x="Sample group", y="Expression")
  ggp
}

plot_violin = function(de_table, gene, sample_table)
{
  gene_data = de_table[gene,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = sample_table$SAMPLE_GROUP
  
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
  return(ggp)
}

plot_jitter = function(de_table, gene, sample_table)
{
  gene_data = de_table[gene,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = sample_table$SAMPLE_GROUP
  
  # changing the order  of coloums
  gene_data = gene_data[c(2,1)]   
  names(gene_data) = c("sample_group", "expression")
  
  # ordering the data in rows by sample group
  levels(gene_data$sample_group) #current data orderis NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
  my_title = paste(gene,"Expression", sep = " ")
  
  ggp = ggplot(gene_data, aes(x=sample_group, y=expression)) + 
    geom_jitter(colour="skyblue4", fill="white", size=2) +
    my_theme_darkblue() +
    labs(title=my_title, x="Sample group", y="Expression")
  ggp
}

plot_violin_jitter = function(de_table, gene, sample_table)
{
  gene_data = de_table[gene,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = sample_table$SAMPLE_GROUP
  
  # changing the order  of coloums
  gene_data = gene_data[c(2,1)]   
  names(gene_data) = c("sample_group", "expression")
  
  # ordering the data in rows by sample group
  levels(gene_data$sample_group) #current data orderis NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
  my_title = paste(gene,"Expression", sep = " ")
  
  ggp = ggplot(gene_data, aes(x=sample_group, y=expression)) + 
    geom_violin(colour="azure4", fill="white") +
    geom_jitter(colour="skyblue4", fill="white", size=2) +
    my_theme_darkblue() +
    labs(title=my_title, x="Sample group", y="Expression")
  ggp
}
# box plot in the argument?
plot_boxplot_jitter = function(de_table, gene, sample_table)
{
  gene_data = de_table[gene,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  
  # changing the order  of coloums
  gene_data = gene_data[c(2,1)]   
  names(gene_data) = c("sample_group", "expression")
  
  # ordering the data in rows by sample group
  levels(gene_data$sample_group) #current data orderis NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
  my_title = paste(gene,"Expression", sep = " ")
  
  ggp = ggplot(gene_data, aes(x=sample_group, y=expression)) + 
    geom_boxplot(colour="skyblue4", fill="white") +
    geom_jitter(colour="skyblue4", fill="white", size=2) +
    my_theme_darkblue() +
    labs(title=my_title, x="Sample group", y="Expression")
  ggp
}

plot_heatmap = function(de_table, range, my_title="Heatmap")
{
  # create the same but with scaled expression
  em_sig_scaled = sort_p(de_table)
  em_sig_scaled = em_sig_scaled[1:100,]
  em_sig_scaled = table_sign_scaled(em_sig_scaled, range)
  
  # data need to be a matrix before clustering
  em_sig_scaled_matrix = as.matrix(em_sig_scaled)
  
  # correlation here:
  y.dist = Dist(em_sig_scaled_matrix, method="spearman") 
  # building the dendrogram
  y.cluster = hclust(y.dist, method="average") 
  # extracting the dendorgram
  y.dd = as.dendrogram(y.cluster) 
  # Step4: we untangle the dendrogram (re-order) 
  y.dd.reorder = reorder(y.dd,0,FUN="average") 
  # Step5: we get the new gene order (as a vector of row indexes) 
  y.order = order.dendrogram(y.dd.reorder) 
  # Step6: we re-order the heatmap matrix to have the new, clustered gene order 
  hm.matrix_clustered = em_sig_scaled_matrix[y.order,] 
  # Step7: now we proceed with the melt and plot as before 
  
  em_sig_scaled_matrix = melt(hm.matrix_clustered)
  
  
  ggp = ggplot(em_sig_scaled_matrix, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() +
    my_theme_basic() +
    labs(x="Sample group", y="Gene", Title = my_title)
  ggp
  return(ggp)
}

# Over Representation Analysis (ORA)
#category: # bp = biological processes # cc = cellular components # mf = molecular functions # all = all of the above
# sig: all sig genes, up - sig upregulated, down - sig downregulated

plot_ora_sig = function(de_table, category="all")
{
  sig = subset(de_table, p.adj<0.001 & abs(log2fold) > 2)
  sig_ENSEMBL = sig[,1]
  sig_ENTREZ = bitr(sig_ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ora_results = enrichGO(gene = sig_ENTREZ$ENTREZID, OrgDb = org.Mm.eg.db, readable = T,  ont = category,  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
  
  ggp = barplot(ora_results, showCategory=10) + my_theme_basic() + labs(title="significant genes by category")
  ggp
  return(ggp)
}

plot_ora_up = function(de_table, category="all")
{
  sig_up = subset(de_table, p.adj<0.001 & log2fold >2)
  sig_up_ENSEMBL = sig_up[,1]
  sig_up_ENTREZ = bitr(sig_up_ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ora_results = enrichGO(gene = sig_up_ENTREZ$ENTREZID, OrgDb = org.Mm.eg.db, readable = T,  ont = category,  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
  ggp = barplot(ora_results, showCategory=10) + my_theme_basic() + labs(title="significantly upregulated genes by category")
  return(ggp)
}

plot_ora_down = function(de_table, category="all")
{
  sig_down = subset(de_table, p.adj<0.001 & log2fold < -2)
  sig_down_ENSEMBL = sig_down[,1]
  sig_down_ENTREZ = bitr(sig_down_ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ora_results = enrichGO(gene = sig_down_ENTREZ$ENTREZID, OrgDb = org.Mm.eg.db, readable = T,  ont = category,  pvalueCutoff = 0.001,  qvalueCutoff = 0.10)
  ggp = barplot(ora_results, showCategory=10) + my_theme_basic() + labs(title="significantly downregulated genes by category")
  return(ggp)
}

plot_gsea = function(de_table)
{
  gsea_input = de_table$log2fold
  names(gsea_input) = row.names(de_table)
  gsea_input = na.omit(gsea_input)
  gsea_input = sort(gsea_input, decreasing = TRUE) 
  gse_results = gseGO(geneList=gsea_input,  
                      ont ="all",
                      keyType = "SYMBOL", 
                      nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.001,
                      verbose = TRUE,
                      OrgDb = org.Hs.eg.db,  
                      pAdjustMethod = "none")
  ggp = ridgeplot(gse_results) + my_theme_basic() + guides(colour=guide_legend(title="adjusted p value"))
  return(ggp)
}

# COMPOUND PLOTS ####


top10_boxplot = function(de_table, p_col_indx)
  # change the column name for the p value: 17 here
  # and the names of the sample groups
{
  sort_order = order(de_table[,p_col_indx], decreasing=FALSE)
  de_table = de_table[sort_order,]
  
  candidate_genes = de_table[1:10,]
  candidate_genes = row.names(candidate_genes)
  
  # build the table
  gene_data = em_symbols[candidate_genes,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  
  # ordering the data 
  levels(gene_data$sample_group) #current data order is NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
  
  # melting the table
  gene_data.m = melt(gene_data, id.vars="sample_group")
  
  ggp = ggplot(gene_data.m, aes(x=variable, y=value, colour=sample_group, fill=sample_group)) + 
    geom_boxplot(alpha=0.8) +
    my_theme_darkblue() +
    scale_fill_manual(values=c("skyblue","lightgreen", "pink")) +
    scale_color_manual(values=c("darkblue", "darkgreen", "pink4")) +
    labs(title="Expression of the top 10 most significant genes", x="Gene", y="Expression")
  ggp
}

top10_boxplot_scaled = function(de_table, p_col_indx, my_title="Expression of the top 10 most significant genes: MvS")
  # change the column name for the p value: 17 here
  # and the names of the sample groups
{
  sort_order = order(de_table[,p_col_indx], decreasing=FALSE)
  de_table = de_table[sort_order,]
  
  candidate_genes = de_table[1:10,]
  candidate_genes = row.names(candidate_genes)
  
  # build the table
  gene_data = em_scaled[candidate_genes,] 
  gene_data = data.frame(t(gene_data))
  gene_data$sample_group = ss$SAMPLE_GROUP
  
  # ordering the data 
  levels(gene_data$sample_group) #current data order is NULL
  gene_data$sample_group = factor(gene_data$sample_group, levels=c("Prolif","Senes","Senes_MtD")) 
  
  # melting the table
  gene_data.m = melt(gene_data, id.vars="sample_group")
  
  ggp = ggplot(gene_data.m, aes(x=sample_group, y=value, facet= variable,colour=sample_group, fill=sample_group)) + 
    geom_boxplot(alpha=0.8) +
    facet_wrap(~variable, ncol=10) +
    my_theme_darkblue() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # must be placed after any other theme 
    theme(legend.position = "bottom") +
    scale_fill_manual(values=c("skyblue","lightgreen", "pink")) +
    scale_color_manual(values=c("darkblue", "darkgreen", "pink4")) +
    labs(title=my_title, x="Sample group", y="Expression")
  ggp
}

# STATISTICS ####
hypometric = function(group1_only, group2_only, both_groups, total_table)
{
  group1 = group1_only + both_groups
  group2 = group2_only + both_groups
  total = nrow(total_table)
  overlap = both_groups
  result = phyper(overlap-1, group2, total-group2, group1,lower.tail= FALSE)
  return(result)
}

