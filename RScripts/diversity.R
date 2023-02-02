#####   Name: diversity.R
#####   Version: 1.0
#####   Description: Diversity plots by sample
#####   Date of Creation: Jan-2023
#####   Author: Natália Faraj Murad
#####   E-MAIL: nataliafmurad@gmail.com
#####   PROJECT: https://github.com/natmurad/amplicon16S


############################
###--- Load libraries ---###
############################
library(tidyverse) # data manipulation
library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(reshape2) # data manipulation
library(ggtree) # for visualizing phylogenetic trees
library(ape) # for manipulating phylogenetic trees
library(ggplot2)
library(ggthemes) # colors
library(ggpubr) # add stats on plots
library(qiime2R) # read qiime2 outputs
library(metagenomeSeq) # for differential abundance analysis - zero inflated model
library(microbiome)
library(Matrix) 
library(vegan)


# I have followed these 2 tutorials: https://rpubs.com/lconteville/713954 & https://rpubs.com/lconteville/714853

###--- Graph config ---###
theme_set(theme_bw())

###--- Set directory ---###
setwd("/Users/natmurad/Documents/")
outDir <- "/Users/natmurad/Documents/results/"

# Creating a Phyloseq Object
physeq<-qza_to_phyloseq(
  features="./DADA2_denoising_output/table.qza",
  tree="rooted_tree.qza",
  "taxonomy.qza",
  metadata = "q2metadata.tsv"
)

###--- Organizing treatments in an order ---###
sample_data(physeq)$Subgroup1 <- factor((sample_data(physeq)$Subgroup1), levels=c("Control","AD","NDGA","AD_NDGA"))

###--- Alpha diversity ---###
richness <- estimate_richness(physeq)
head(richness)

richnessg <- plot_richness(physeq, color = "Subgroup1") + 
  scale_x_discrete(labels = sample_data(physeq)$customer_label) +
  scale_colour_tableau("Classic 10") +
  ggtitle("Diversity Measures") + labs(color = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("richnessMeasures.pdf", plot = richnessg,
       device = "pdf", path = outDir, width = 20, height = 4,
       units = "in")

###--- Significance ---###
a_my_comparisons <- list( c("AD", "Control"), c("NDGA", "Control"), c("AD_NDGA", "AD"), c("AD_NDGA", "Control"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

alphadiv <- plot_richness(physeq, x="Subgroup1", measures="Shannon", color = "Subgroup1")+
  geom_boxplot(aes(fill=Subgroup1), alpha=0.2) + 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) +
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) +
  xlab("Condition") +
  scale_colour_tableau("Classic 10") +
  scale_fill_tableau("Classic 10") +
  ggtitle("Alpha Diversity") + labs(color = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("alphaDiversity.pdf", plot = alphadiv,
       device = "pdf", path = outDir, width = 2.8, height = 5.5,
       units = "in")
# A histogram can be created to check if the Shannon values estimated for the metagenomes are normally distributed: 
hist(richness$Shannon, main="Shannon index", xlab="")

# Since it is, we can run anova tests and check if the variable Group impacts the Shannon diversity:
anova.sh = aov(richness$Shannon ~ sample_data(physeq)$Subgroup1)
summary(anova.sh)

# Based on the anova tests results, it’s possible to compute the Tukey Honest Significant Differences:
TukeyHSD(anova.sh)

# For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test:
x <- kruskal.test(richness$Shannon ~ sample_data(physeq)$Subgroup1)

# We can also get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness$Shannon, sample_data(physeq)$Subgroup, p.adj = "BH")

###--- Beta diversity ---###
# transform the data to relative abundance and create a new phyloseq object:

relab_genera = transform_sample_counts(physeq, function(x) x / sum(x) * 100)
head(otu_table(relab_genera)[,1:6])

# Calculate Bray-Curtis distance among samples and convert the result to a matrix:
abrel_bray <- phyloseq::distance(relab_genera, method = "bray")
abrel_bray <- as.matrix(abrel_bray)
head(abrel_bray)[,1:6]
write.csv(abrel_bray, paste0(outDir, "brayscurtis.csv"))

# The “abrel_bray” matrix presents the distance between all samples, but since we wanna generate a boxplot with distances considering the metagenomes in each group separately, we need to filter this matrix. With the code below we end with a dataframe “df.bray” storing Bray Curtis Distances between metagenomes from the same groups.
sub_dist <- list()
groups_all <- sample_data(relab_genera)$Subgroup1

for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- sample_names(relab_genera)[row_group]
  sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups<- melt(sub_dist)
df.bray <- braygroups[complete.cases(braygroups), ]
df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))

head(df.bray)

###--- boxplot ---###
brayc <- ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(aes(fill=L1), alpha=0.2) + 
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) +
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12)) +
  scale_colour_tableau("Classic 10") +
  scale_fill_tableau("Classic 10") +
  ggtitle("Bray-Curtis Distance") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("braycurtis.pdf", plot = brayc,
       device = "pdf", path = outDir, width = 2.8, height = 5.5,
       units = "in")

###--- PcoA plot ---###

# To generate the PcoA plot, first get the ordination result with the “ordinate” command from phyloseq. It requires a phyloseq object, and it accepts diverse methods and distances. Next, generate the plot using the “plot_ordination()” function.

ord = ordinate(relab_genera, method="PCoA", distance = "bray")


pcoa <- plot_ordination(relab_genera, ord, color = "Subgroup1", shape="condition") + 
  geom_point(size=4) + 
  stat_ellipse(aes(group=Subgroup1)) +
  scale_colour_tableau("Classic 10") +
  scale_fill_tableau("Classic 10") +
  ggtitle("PcoA") +
  theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Condition", shape = "")

ggsave("pcoa.pdf", plot = pcoa,
       device = "pdf", path = outDir, width = 6, height = 5,
       units = "in")

###--- Significance ---###
# To test if the Control and the Treatment Groups are significantly different from each other we can run a Permanova test using the adonis function from the vegan package.
samples <- data.frame(sample_data(relab_genera))
adonis(abrel_bray ~ condition, data = samples)

# species network
#plot_net(physeq, color="Species", type="taxa")
# beta diversity network
net <- plot_net(physeq, color="Subgroup1", distance="bray") +
  scale_colour_tableau("Classic 10") +
  ggtitle("Beta Diversity Network - Bray-Curtis Distance") +
  theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Condition")

ggsave("betadiversityNet.pdf", plot = net,
       device = "pdf", path = outDir, width = 7, height = 5,
       units = "in")
