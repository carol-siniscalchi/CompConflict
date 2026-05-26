library(phytools)
library(ggplot2)
library(dplyr)
library(apaTables)

# read data file 
triplets <- read.csv("~/LikelihoodMapping/alignment_info_with_triplets.csv")
triplets$FinalT <- as.factor(triplets$FinalT)

# check counts
triplets %>% 
  count(FinalT) 

# order factor levels
triplets$FinalT <- ordered(triplets$FinalT, levels = c("T1", "T2", "T3", "T12", "T13", "T23", "T123"))


## Plot barplots and boxplots ##

pdf("bar_topology_3.pdf", height = 5, width = 5)
ggplot(triplets, aes(x = FinalT, fill = FinalT)) +
  geom_bar(color = "black") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.3) +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "Number of alignments") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none")
dev.off()

pdf("violin_alignlength.pdf", height = 4, width = 6)
ggplot(triplets, aes(x = FinalT, y = Alignment_length, fill = FinalT)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill = "grey75", alpha = 0.5) +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "Alignments length (bp)") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none")
dev.off()

pdf("violin_alignPI.pdf", height = 4, width = 6)
ggplot(triplets, aes(x = FinalT, y = Parsimony_informative_sites, fill = FinalT)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill = "grey75", alpha = 0.5) +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "PI sites (bp)") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none")
dev.off()

pdf("violin_aligntaxa.pdf", height = 4, width = 6)
ggplot(triplets, aes(x = FinalT, y = No_of_taxa, fill = FinalT)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill = "grey75", alpha = 0.5) +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "Number of taxa in alignment") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none")
dev.off()

pdf("bp_alignPI.pdf", height = 3, width = 4)
ggplot(triplets, aes(x = FinalT, y = Parsimony_informative_sites, fill = FinalT)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "PI sites (bp)") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none") +
  coord_flip()
dev.off()

pdf("bp_aligntaxa.pdf", height = 3, width = 4)
ggplot(triplets, aes(x = FinalT, y = No_of_taxa, fill = FinalT)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#dd5e5e", "#e3afaf", "#c4875e", "#dddddd", "#dddddd", "#dddddd", "#ababab")) +
  theme_classic() +
  labs(x = "Chosen topology", y = "Number of taxa in alignment") +
  theme(axis.title=element_text(size=16)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position = "none") +
  coord_flip()
dev.off()

###### Means, ANOVA and Tukey HSD #######

## Number of taxa in alignment
triplets %>% 
  group_by(FinalT) %>% 
  summarise(mean = mean(No_of_taxa))

#FinalT  mean
#1 T1      39.3
#2 T2      37.9
#3 T3      39.9
#4 T12     39.7
#5 T13     40.4
#6 T23     39.7
#7 T123    35.6

taxa.aov <- aov(No_of_taxa ~ FinalT, data = triplets)
summary(taxa.aov)
plot(taxa.aov)

taxa.hsd <- TukeyHSD(taxa.aov,conf.level = 0.95)
mat <- as.matrix(taxa.hsd)
df_res <- as.data.frame(mat[1])
plot(taxa.hsd, col= ifelse(df_res[,4]<0.05, 'red', 'black'))

# significant comparisons:
# T123-T1  -3.669103e+00  -6.210947 -1.1272593 0.0004785
# T123-T3  -4.305556e+00  -8.378513 -0.2325985 0.0304930

## Alignment length
triplets %>% 
  group_by(FinalT) %>% 
  summarise(mean = mean(Alignment_length))

#FinalT  mean
#1 T1     2434.
#2 T2     1940.
#3 T3     2335.
#4 T12    2034.
#5 T13    2468.
#6 T23    2038.
#7 T123   1554.

len.aov <- aov(Alignment_length ~ FinalT, data = triplets)
summary(len.aov)
plot(len.aov)

len.hsd <- TukeyHSD(len.aov,conf.level = 0.95)
mat <- as.matrix(len.hsd)
df_res <- as.data.frame(mat[1])
plot(len.hsd, col= ifelse(df_res[,4]<0.05, 'red', 'black'))

# Significant comparison:
# T123-T1  -879.807505 -1553.6498 -205.9652 0.0024507

## PI Sites

triplets %>% 
  group_by(FinalT) %>% 
  summarise(mean = mean(Parsimony_informative_sites))

#FinalT  mean
#1 T1      356.
#2 T2      270.
#3 T3      332.
#4 T12     312.
#5 T13     309.
#6 T23     291.
#7 T123    237.

pis.aov <- aov(Parsimony_informative_sites ~ FinalT, data = triplets)
summary(pis.aov)
plot(pis.aov)

pis.hsd <- TukeyHSD(pis.aov,conf.level = 0.95)
mat <- as.matrix(pis.hsd)
df_res <- as.data.frame(mat[1])
plot(pis.hsd, col= ifelse(df_res[,4]<0.032, 'red', 'black'))

#T2-T1     -85.775381 -167.56453  -3.986228 0.0328541
#T123-T1  -119.204191 -188.43609 -49.972297 0.0000112

### ANOVA tables ###

table.len <- apa.aov.table(len.aov, filename = "len.doc")
table.taxa <- apa.aov.table(taxa.aov, filename = "taxa.doc")
table.pis <- apa.aov.table(pis.aov, filename = "pis.doc")

### Tukey HSD plots ####

pdf("tukeyHSD.pdf", width = 15, height = 7)
par(mfrow = c(1,3))
plot(len.hsd, col= ifelse(df_res[,4]<0.05, 'red', 'black'))
plot(pis.hsd, col= ifelse(df_res[,4]<0.032, 'red', 'black'))
plot(taxa.hsd, col= ifelse(df_res[,4]<0.05, 'red', 'black'))
dev.off()

###### Plot simulated topologies trees ########

topo1 <- read.tree("~/LikelihoodMapping/test_topo1.tre")
topo2 <- read.tree("~/LikelihoodMapping/test_topo2.tre")
topo3 <- read.tree("~/LikelihoodMapping/test_topo3.tre")

pdf("topologies.pdf", width = 15, height = 7)
par(mfrow = c(1,3))
plot.phylo(topo1)
plot.phylo(topo2)
plot.phylo(topo3)
dev.off()