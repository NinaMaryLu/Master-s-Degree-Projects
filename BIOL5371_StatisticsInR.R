

# install.packages("dplyr")
# install.packages("FSA")
# install.packages("openxlsx")
#install.packages("caret")
install.packages("rstatix")
install.packages("tidyverse")
install.packages("ggpubr")
library(ggplot2)
library(dplyr)
library(FSA)
library(openxlsx)
library(caret)
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggpubr)

annotations <- read.table("C:/Users/User/Documents/UoG 2024 2025/Statistics in R/Assessment/Data for Report Assessment-20241026/Annotations.csv",
                        header=TRUE, row.names=1, sep = "\t")
DE_GOUT_vs_HC <- read.table("C:/Users/User/Documents/UoG 2024 2025/Statistics in R/Assessment/Data for Report Assessment-20241026/DE_GOUT_vs_HC.csv",
                            header=TRUE, row.names=1, sep = "\t")
DE_SA_vs_HC <- read.table("C:/Users/User/Documents/UoG 2024 2025/Statistics in R/Assessment/Data for Report Assessment-20241026/DE_SA_vs_HC.csv",
                          header=TRUE, row.names=1, sep = "\t")
Expression_Table <- read.table("C:/Users/User/Documents/UoG 2024 2025/Statistics in R/Assessment/Data for Report Assessment-20241026/Expression_Table.csv",
                               header=TRUE, row.names=1, sep = "\t")
Sample_Info <- read.table("C:/Users/User/Documents/UoG 2024 2025/Statistics in R/Assessment/Data for Report Assessment-20241026/Sample_Information.csv",
                                 header=TRUE, row.names=1, sep = "\t")
# lacking = subset(DE_SA_vs_HC, p.adj=="")
# lacking = subset(DE_SA_vs_HC, p.adj == "")


# Sample info analysis
summary(Sample_Info)
str(Sample_Info)
Sample_Info$SEX = as.factor(Sample_Info$SEX)
Sample_Info$SAMPLE_GROUP = as.factor(Sample_Info$SAMPLE_GROUP)
str(Sample_Info)
summary(Sample_Info)

samples_Females = subset(Sample_Info,SEX == "F")
samples_Males = subset(Sample_Info,SEX == "M")
summary(samples_Females)
summary(samples_Males)

    # Figure: sex distribution in groups
ggp = ggplot(Sample_Info, mapping = aes(x = SAMPLE_GROUP)) +
  geom_bar(aes(fill = SEX), position="stack") + xlab("Group") + scale_y_discrete(limits=seq(0,10,1))
ggp

    # Figure: neutrophil count in groups and by sex
ggp = ggplot(Sample_Info, aes(x=SAMPLE_GROUP,y = NEUTROPHILS)) +
  geom_boxplot(size=0.5, alpha=0.3, fill="grey") + xlab("Group") +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.4)
ggp

ggp = ggplot(Sample_Info, aes(x=SAMPLE_GROUP,y = NEUTROPHILS)) +
  geom_violin(trim=FALSE, fill="brown", alpha=0.3) + xlab("Group")
ggp

ggp = ggplot(Sample_Info, aes(x=SAMPLE_GROUP,y = NEUTROPHILS, fill=SEX)) +
  geom_boxplot(size=0.5, alpha=0.3) + xlab("Group")
ggp

    # Nice graph for neutrophil distribution; Neutrophil distribution within samples
ggp = ggplot(Sample_Info, aes(x=NEUTROPHILS)) + geom_histogram(colour="brown",  fill="white", size=0.5, alpha=0.8)
ggp
ggp = ggplot(Sample_Info, aes(x=NEUTROPHILS, fill=SAMPLE_GROUP)) + geom_density(size=0.5, alpha=0.8)
ggp
ggp = ggplot(Sample_Info, aes(x=NEUTROPHILS, fill=SEX)) + geom_density(size=0.5, alpha=0.8)
ggp

    # Neutrophil count per gender
samples_HC = subset(Sample_Info,SAMPLE_GROUP == "HC")
samples_GOUT = subset(Sample_Info,SAMPLE_GROUP == "GOUT")
samples_SEPSIS = subset(Sample_Info,SAMPLE_GROUP == "SEPSIS")

summary(samples_HC) # mean 5.04
sd(samples_HC$NEUTROPHILS)  # 2.58

summary(samples_GOUT)   # mean 7.48
sd(samples_GOUT$NEUTROPHILS)  # 3.45

summary(samples_SEPSIS)   # mean 12.78
sd(samples_SEPSIS$NEUTROPHILS)  # 2.69

ggp = ggplot(samples_HC, aes(x=NEUTROPHILS, fill=SEX)) + geom_density(size=0.5, alpha=0.8)      # looks like femals might have higher Neutrophil count
ggp
ggp = ggplot(samples_GOUT, aes(x=NEUTROPHILS, fill=SEX)) + geom_density(size=0.5, alpha=0.8)    # oh what is that. Looks intriquing. Low n?
ggp
ggp = ggplot(samples_SEPSIS, aes(x=NEUTROPHILS, fill=SEX)) + geom_density(size=0.5, alpha=0.8)  # ohhhhhh males have high neutrophil but not females
ggp

    # Females neutrophil by group
ggp = ggplot(samples_Females, aes(x=NEUTROPHILS, fill=SAMPLE_GROUP)) + geom_density(size=0.5, alpha=0.8)  # ohhhhhh males have high neutrophil but not females
ggp
# looks to me that: healthy low, sepsis high, gout varied with low / medium / high. "Halves" but for 3 groups would be cool but my sample is too low
ggp = ggplot(samples_Males, aes(x=NEUTROPHILS, fill=SAMPLE_GROUP)) + geom_density(size=0.5, alpha=0.8)  # ohhhhhh males have high neutrophil but not females
ggp
# males have a MUCH higher response to sepsis in a neutrophil count.

# Multiple T tests for gender differences in groups - no significant differences between sex neutrophil scores within groups
HC_fem = subset(samples_HC,SEX=="F")
HC_male = subset(samples_HC,SEX=="M")
wilcox.test(samples_HC$NEUTROPHILS ~ samples_HC$SEX, alternative = "two.sided")   # NOT SIGNIFICANT
HC_p = t.test(HC_fem$NEUTROPHILS, HC_male$NEUTROPHILS)     # NOT SIGNIFICANT, HC almost significant 
HC_p

GOUT_fem = subset(samples_GOUT,SEX=="F")
GOUT_male = subset(samples_GOUT,SEX=="M")
wilcox.test(samples_GOUT$NEUTROPHILS ~ samples_GOUT$SEX)   # cannot compute
GOUT_p = t.test(GOUT_fem$NEUTROPHILS, GOUT_male$NEUTROPHILS)     # NOT SIGNIFICANT
GOUT_p

SEPSIS_fem = subset(samples_SEPSIS,SEX=="F")
SEPSIS_male = subset(samples_SEPSIS,SEX=="M")
wilcox.test(samples_SEPSIS$NEUTROPHILS ~ samples_SEPSIS$SEX, alternative = "two.sided")   # NOT SIGNIFICANT
SEPSIS_p = t.test(SEPSIS_fem$NEUTROPHILS, SEPSIS_male$NEUTROPHILS)     # NOT SIGNIFICANT
SEPSIS_p

# neutrophil scores between groups
model1 = lm(Sample_Info$NEUTROPHILS ~ Sample_Info$SAMPLE_GROUP)   # SIGNIFICANT
anova(model1)

model2 = lm(Sample_Info$NEUTROPHILS ~ Sample_Info$SEX)          # NOT SIGNIFICANT
anova(model2)

# model3 = lm(samples_Females$NEUTROPHILS ~ Sample_Info$SAMPLE_GROUP)          # groups not even
# anova(model3)
# 
# model4 = lm(samples_Males$NEUTROPHILS ~ Sample_Info$SAMPLE_GROUP)          # groups not even
# anova(model4)

model3 = lm(Sample_Info$NEUTROPHILS ~ (Sample_Info$SEX + Sample_Info$SAMPLE_GROUP))
anova(model3)    # significant, but data variance is not equal bc of the sepsis sex groups        
summary(model3)

kruskal.test(Sample_Info$NEUTROPHILS ~ SEX, data = Sample_Info)           # NOT SIGNIFICANT
kruskal.test(Sample_Info$NEUTROPHILS ~ SAMPLE_GROUP, data = Sample_Info)  # SIGNIFICANT
dunnTest(NEUTROPHILS ~ SAMPLE_GROUP,
         data=Sample_Info,
         method="bonferroni")

# genes significant in GOUT 
sign_GOUT_vs_HC = subset(DE_GOUT_vs_HC, p.adj<0.05)   # 69
GOUT_up_all = subset(sign_GOUT_vs_HC, log2Fold>=0)    # 32
GOUT_down_all = subset(sign_GOUT_vs_HC, log2Fold<0)    # 37
sign_up_GOUT_vs_HC = subset(sign_GOUT_vs_HC, log2Fold>1)  # 10
sign_down_GOUT_vs_HC = subset(sign_GOUT_vs_HC, log2Fold<(-1)) # 5

# genes significant in SEPSIS - mnóstwo
sign_SA_vs_HC = subset(DE_SA_vs_HC, p.adj<0.05)   # mnóstwo, 13 046
SA_up_all = subset(sign_SA_vs_HC, log2Fold>=0)    # 5 879 
SA_down_all = subset(sign_SA_vs_HC, log2Fold<0)    # 7 167
sign_up_SA_vs_HC = subset(sign_SA_vs_HC, log2Fold>1)   # 2 158
sign_down_SA_vs_HC = subset(sign_SA_vs_HC, log2Fold<(-1))   # 4 112

# further groups
sign_GOUT_SA = merge(sign_GOUT_vs_HC, sign_SA_vs_HC, by=0)    # 47 genes significant for both
GOUT_SA_up_all = merge(GOUT_up_all, SA_up_all, by=0)          # 15
GOUT_SA_down_all = merge(GOUT_down_all, SA_down_all, by=0)          # 26
# sign_GOUT_SA_up = merge(sign_up_GOUT_vs_HC, sign_up_SA_vs_HC, by=0) # 4 genes upregulated for both
# sign_GOUT_SA_down = merge(sign_down_GOUT_vs_HC, sign_down_SA_vs_HC, by=0) # 2 genes downregulated for both

sign_GOUT_up_SA_down = merge(GOUT_up_all, SA_down_all, by=0)    # 4 genes 
sign_GOUT_down_SA_up = merge(GOUT_down_all, SA_up_all, by=0)    # 2 genes

GOUT_SA_vs_HC = merge(DE_GOUT_vs_HC, DE_SA_vs_HC, by=0)
row.names(GOUT_SA_vs_HC) = GOUT_SA_vs_HC[,1]
names(GOUT_SA_vs_HC)[1] = "GENE_ID"
GOUT_only = subset(GOUT_SA_vs_HC, p.adj.x<0.05 & p.adj.y>=0.05) # 22 genes sign changed for gout but not SA
SEPSIS_only = subset(GOUT_SA_vs_HC, p.adj.y>=0.05 & p.adj.y<0.05)  # no genes changed only for SA
GOUT_only_up = subset(GOUT_only, log2Fold.x>1)  # 4 genes upregulated for GOUT with no change in SA
GOUT_only_down = subset(GOUT_only, log2Fold.x<(-1)) # 1 gene downregulated for GOUT with no change in SA

# All genes unique for GOUT
GOUT_only = merge(GOUT_only, annotations, by=0)
row.names(GOUT_only) = GOUT_only[,1]
names(GOUT_only)[1] = "GENE_ID"
GOUT_only = GOUT_only[,-1]
GOUT_only = GOUT_only[,-5:-7]
new_names = c("log2Fold", "p", "p.adj.") 
names(GOUT_only)[2:4] = new_names
GOUT_only <- GOUT_only[order(GOUT_only$log2Fold, decreasing = TRUE),]
write.xlsx(GOUT_only, 'GOUT_only.xlsx')

# Genes unique for Gout and upregulated
DE_GOUT_only_up = merge(GOUT_only_up, annotations, by=0)
row.names(DE_GOUT_only_up) = DE_GOUT_only_up[,1]
names(DE_GOUT_only_up)[1] = "GENE_ID"
DE_GOUT_only_up = DE_GOUT_only_up[,-1]
DE_GOUT_only_up = DE_GOUT_only_up[,-5:-7]
new_names = c("log2Fold", "p", "p.adj.") 
names(DE_GOUT_only_up)[2:4] = new_names
DE_GOUT_only_up <- DE_GOUT_only_up[order(DE_GOUT_only_up$log2Fold, decreasing = TRUE),]
write.xlsx(DE_GOUT_only_up, 'DE_GOUT_only_up.xlsx')

# Genes unique for Gout and downregulated
DE_GOUT_only_down = merge(GOUT_only_down, annotations, by=0)
row.names(DE_GOUT_only_down) = DE_GOUT_only_down[,1]
names(DE_GOUT_only_down)[1] = "GENE_ID"
DE_GOUT_only_down = DE_GOUT_only_down[,-1]
DE_GOUT_only_down = DE_GOUT_only_down[,-5:-7]
new_names = c("log2Fold", "p", "p.adj.") 
names(DE_GOUT_only_down)[2:4] = new_names
write.xlsx(DE_GOUT_only_down, 'DE_GOUT_only_down.xlsx')

# Expression tables
Expression_Table = t(Expression_Table)
Expression_Table = data.frame(Expression_Table)

DE_Expression = merge(Sample_Info, Expression_Table, by=0)
row.names(DE_Expression) = DE_Expression[,1]
names(DE_Expression)[1] = "Sample_ID"

ggp = ggplot(DE_Expression, aes(x=NEUTROPHILS, y=ENSG00000240204, colour=SEX)) + geom_point() #x and y are two genes
ggp

cor(DE_Expression$NEUTROPHILS, DE_Expression$ENSG00000130540, method = c("spearman"))

expression_females = subset(DE_Expression, DE_Expression$SEX == "F")
expression_males = subset(DE_Expression, DE_Expression$SEX == "M")

ggp = ggplot(expression_females, aes(x=SAMPLE_GROUP,y = ENSG00000091513, fill=SEX)) +
  geom_boxplot(size=0.5, alpha=0.3) + xlab("Group")
ggp

ggp = ggplot(DE_Expression, aes(x=SAMPLE_GROUP,y = ENSG00000198759, fill=SEX)) +
  geom_boxplot(size=0.5, alpha=0.3) + xlab("Group")
ggp

ggp = ggplot(DE_Expression, aes(x=ENSG00000240204, fill=SAMPLE_GROUP)) +
  geom_density(size=0.5, alpha=0.8)
ggp

# GOUT Separating data from the expression table by the top 10 significant gout genes
Expression_Table = t(Expression_Table)
Expression_Table = data.frame(Expression_Table)

top10GOUT <- sign_GOUT_vs_HC[order(sign_GOUT_vs_HC$p.adj), ]
top10GOUT <- top10GOUT[1:10, ]

top10GOUT_expression = merge(top10GOUT,Expression_Table, by=0)
row.names(top10GOUT_expression) = top10GOUT_expression[,1]
GOUT_means = top10GOUT_expression[,-c(1,2,3,4)]
GOUT_means = t(GOUT_means)
GOUT_means = data.frame(GOUT_means)
DE_GOUT_means = merge(Sample_Info,GOUT_means,by=0)
row.names(DE_GOUT_means) = DE_GOUT_means[,1]
DE_GOUT_means = DE_GOUT_means[,-1]
summary(DE_GOUT_means)

# GOUT Normalising the expression of those genes to plot together

for (column in 4:ncol(DE_GOUT_means)){
  DE_GOUT_means[,column] = as.numeric(DE_GOUT_means[, column])
  process <- preProcess(as.data.frame(DE_GOUT_means[, column]), method=c("range"))
  norm_scale <- predict(process, as.data.frame(DE_GOUT_means[,column]))
  DE_GOUT_means[,column] = norm_scale
}
DE_GOUT_means = DE_GOUT_means[,-c(2,3)]

#multiple t tests for gender for GOUT only for GOUT vs SEPSIS
mydata <- subset(DE_GOUT_means, SAMPLE_GROUP=="GOUT" | SAMPLE_GROUP=="SEPSIS")
mydata.long <- mydata %>%
  pivot_longer(-SAMPLE_GROUP, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(18)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ SAMPLE_GROUP) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
write.xlsx(stat.test, 'stat.test_top10_GOUT.xlsx')

# for the visualisations
DE_GOUT_means <- melt(DE_GOUT_means ,  id.vars = 'SAMPLE_GROUP', variable.name = 'Gene')

ggp = ggplot(DE_GOUT_means) + geom_boxplot(aes(x = Gene, y = value, fill = SAMPLE_GROUP)) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Genes") +
  ylab("Normalized expression") +
  ggtitle("Top 10 Significant DEG in GOUT") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_grid(cols = vars(Gene), scales=scales = "free_x", space = "free")
ggp

ggp = ggplot(DE_GOUT_means) + geom_boxplot(aes(x = Gene, y = value, fill = SAMPLE_GROUP)) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Genes") +
  ylab("Normalized expression") +
  ggtitle("Top 10 Significant DEG in GOUT by sample group") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(vars(Gene), nrow=1, scales = "free_x")
ggp

ggp = ggplot(DE_GOUT_means) + geom_density(aes(x = value, fill = SAMPLE_GROUP), alpha=0.5) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("Density") +
  ggtitle("Distribution of the Top 10 Significant DEG in GOUT by group") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(vars(Gene), nrow=1, scales="free_x")
ggp

# SA Separating data from the expression table by the top 10 significant SA genes
top10SA <- sign_SA_vs_HC[order(sign_SA_vs_HC$p.adj), ]
top10SA <- top10SA[1:10, ]

top10SA_expression = merge(top10SA,Expression_Table, by=0)
row.names(top10SA_expression) = top10SA_expression[,1]
SA_means = top10SA_expression[,-c(1,2,3,4)]
SA_means = t(SA_means)
SA_means = data.frame(SA_means)
DE_SA_means = merge(Sample_Info,SA_means,by=0)
row.names(DE_SA_means) = DE_SA_means[,1]
DE_SA_means = DE_SA_means[,-1]
summary(DE_SA_means)

# SA Normalising the expression of those genes to plot together

for (column in 4:ncol(DE_SA_means)){
  DE_SA_means[,column] = as.numeric(DE_SA_means[,column])
  process <- preProcess(as.data.frame(DE_SA_means[,column]), method=c("range"))
  norm_scale <- predict(process, as.data.frame(DE_SA_means[,column]))
  DE_SA_means[,column] = norm_scale
}
DE_SA_means = DE_SA_means[,-c(2,3)]

#multiple t tests for gender for GOUT only for GOUT vs SEPSIS
mydata <- subset(DE_SA_means, SAMPLE_GROUP=="GOUT" | SAMPLE_GROUP=="SEPSIS")
mydata.long <- mydata %>%
  pivot_longer(-SAMPLE_GROUP, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(18)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ SAMPLE_GROUP) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
write.xlsx(stat.test, 'stat.test_top10_SA.xlsx')

# for the visualisation
DE_SA_means <- melt(DE_SA_means ,  id.vars = 'SEX', variable.name = 'Gene')

ggp = ggplot(DE_SA_means) + geom_boxplot(aes(x = Gene, y = value, fill = SEX)) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Genes") +
  ylab("Normalized expression") +
  ggtitle("Top 10 Significant DEG in SA by sex") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_grid(cols = vars(Gene), scales = "free_x", space = "free")
ggp

ggp = ggplot(DE_SA_means) + geom_density(aes(x = value, fill = SEX), alpha=0.5) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("Density") +
  ggtitle("Distribution of the Top 10 Significant DEG in SA by sex") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(vars(Gene), nrow=1, scales="free_x")
ggp

#multiple t tests for gender
mydata.long <- DE_SA_means %>%
  pivot_longer(-SEX, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(27)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ SEX) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

#multiple t tests for gender ON THE WHOLE EXPRESSION TABLE - too many genes
#multiple t tests for gender for GOUT 

mydata.long <- DE_GOUT_means %>%
  pivot_longer(-SEX, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(27)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ SEX) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

# ggp = ggplot(DE_SA_means) + geom_boxplot(aes(x = Gene, y = value, fill = SAMPLE_GROUP)) +
#   scale_fill_brewer(palette = "Set3") +
#   xlab("Genes") +
#   ylab("Normalized expression") +
#   ggtitle("Top 10 Significant DEG in SA by group") +
#   theme_bw() +
#   theme(axis.text.x = element_blank()) +
#   facet_wrap(vars(Gene), nrow=2, scales="free")
# ggp


top10SA_overlap_GOUT = merge(top10SA, top10GOUT, by=0) # all different

top10SA_symbols = merge(annotations, top10SA, by=0)
top10SA_symbols = top10SA_symbols[,-c(3)]
write.xlsx(top10SA_symbols, 'top10SA_symbols.xlsx')

top10GOUT_symbols = merge(annotations, top10GOUT, by=0)
top10GOUT_symbols = top10GOUT_symbols[,-c(3)]
write.xlsx(top10GOUT_symbols, 'top10GOUT_symbols.xlsx')

# ONLY GOUT Separating data from the expression table by the top 10 significant gout genes
onlyGOUT_expression = merge(GOUT_only,Expression_Table, by=0)
row.names(onlyGOUT_expression) = onlyGOUT_expression[,1]
onlyGOUT_means = onlyGOUT_expression[,-c(1,2,3,4)]
onlyGOUT_means = t(onlyGOUT_means)
onlyGOUT_means = data.frame(onlyGOUT_means)
DE_onlyGOUT_means = merge(Sample_Info,onlyGOUT_means,by=0)
row.names(DE_onlyGOUT_means) = DE_onlyGOUT_means[,1]
DE_onlyGOUT_means = DE_onlyGOUT_means[,-1]
summary(DE_onlyGOUT_means)


# onlyGOUT Normalising the expression of those genes to plot together

for (column in 4:ncol(DE_onlyGOUT_means)){
  DE_onlyGOUT_means[,column] = as.numeric(DE_onlyGOUT_means[,column])
  process <- preProcess(as.data.frame(DE_onlyGOUT_means[,column]), method=c("range"))
  norm_scale <- predict(process, as.data.frame(DE_onlyGOUT_means[,column]))
  DE_onlyGOUT_means[,column] = norm_scale
}
DE_onlyGOUT_means = DE_onlyGOUT_means[,-c(1,3)]
DE_onlyGOUT_means <- melt(DE_onlyGOUT_means ,  id.vars = 'SEX', variable.name = 'Gene')

#multiple t tests for gender for GOUT only for GOUT vs SEPSIS
mydata <- subset(DE_onlyGOUT_means, SAMPLE_GROUP=="GOUT" | SAMPLE_GROUP=="SEPSIS")
mydata.long <- mydata %>%
  pivot_longer(-SAMPLE_GROUP, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(18)

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ SAMPLE_GROUP) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
write.xlsx(stat.test, 'stat.test_onlyGOUT.xlsx')
# ONLY GOUT VISUALISATION
ggp = ggplot(DE_onlyGOUT_means) + geom_boxplot(aes(x = Gene, y = value, fill = SEX)) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Genes") +
  ylab("Normalized expression") +
  ggtitle("DEG unique for GOUT") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_grid(cols = vars(Gene), scales = "free_x", space = "free")
ggp

ggp = ggplot(DE_onlyGOUT_means) + geom_boxplot(aes(x = Gene, y = value, fill = SEX)) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Genes") +
  ylab("Normalized expression") +
  ggtitle("DEG unique for GOUT") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(vars(Gene), ncol=9, scales="free_x")
ggp

ggp = ggplot(DE_SA_means) + geom_density(aes(x = value, fill = SAMPLE_GROUP), alpha=0.5) +
  scale_fill_brewer(palette = "Set3") +
  xlab("") +
  ylab("Density") +
  ggtitle("Distribution of the Top 10 Significant DEG in SA") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  facet_wrap(vars(Gene), nrow=2, scales="free")
ggp

onlyGOUT_symbols = merge(annotations, GOUT_only, by=0)
onlyGOUT_symbols = onlyGOUT_symbols[,-c(3)]
write.xlsx(onlyGOUT_symbols, 'top10GOUT_symbols.xlsx')

# ONLY GOUT Separating data from the expression table significant gout genes
onlyGOUT_expression = merge(GOUT_only,Expression_Table, by=0)
row.names(onlyGOUT_expression) = onlyGOUT_expression[,1]
onlyGOUT_means = onlyGOUT_expression[,-c(1,2,3,4)]
onlyGOUT_means = t(onlyGOUT_means)
onlyGOUT_means = data.frame(onlyGOUT_means)
DE_onlyGOUT_means = merge(Sample_Info,onlyGOUT_means,by=0)
row.names(DE_onlyGOUT_means) = DE_onlyGOUT_means[,1]
DE_onlyGOUT_means = DE_onlyGOUT_means[,-1]
summary(DE_onlyGOUT_means)

# onlyGOUT Normalising the expression of those genes to plot together

for (column in 4:ncol(DE_onlyGOUT_means)){
  DE_onlyGOUT_means[,column] = as.numeric(DE_onlyGOUT_means[,column])
  process <- preProcess(as.data.frame(DE_onlyGOUT_means[,column]), method=c("range"))
  norm_scale <- predict(process, as.data.frame(DE_onlyGOUT_means[,column]))
  DE_onlyGOUT_means[,column] = norm_scale
}

# ggp = ggplot(DE_onlyGOUT_means) + geom_point(aes(x = NEUTROPHILS, y = ENSG00000130701, colour=SEX)) +
#   xlab("Neutrophil count") +
#   ylab("Normalized gene expression") +
#   ggtitle("ENSG0000013070 by gender") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1))
# ggp

DE_onlyGOUT_means = DE_onlyGOUT_means[,-c(1,2)]
DE_onlyGOUT_means <- melt(DE_onlyGOUT_means ,  id.vars = 'NEUTROPHILS', variable.name = 'Gene')

ggp = ggplot(DE_onlyGOUT_means) + geom_point(aes(x = NEUTROPHILS, y = value)) +
  xlab("Neutrophil count") +
  ylab("Normalized gene expression") +
  ggtitle("DEG unique for GOUT: correlation plots vs. neutrophil counts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  facet_wrap(vars(Gene), nrow=4, scales = "free_x")
ggp

DE_onlyGOUT_means_top10 = merge(GOUT_only, top10GOUT)
