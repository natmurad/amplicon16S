#####   Name: composition.R
#####   Version: 1.0
#####   Description: Composition plots by sample
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
library(ggrepel)
library(ggtree) # phylogenetic trees visualization
library(ape) # phylogenetic trees manipulation
library(ggplot2)
library(ggthemes) # colors
library(RColorBrewer)
library(qiime2R) # read qiime2 outputs
library(phyloseq)
############################


###--- Graph config ---###
theme_set(theme_bw())

###--- Set directory ---###
setwd("/Users/natmurad/Documents/")
outDir <- "/Users/natmurad/Documents/results/"

###--- Read files ---###
# Create a Phyloseq Object
physeq<-qza_to_phyloseq(
  features="./DADA2_denoising_output/table.qza",
  tree="rooted_tree.qza",
  "taxonomy.qza",
  metadata = "q2metadata.tsv"
)

###--- Sample Composition ---###

###--- Filtering ---###
# Remove OTUs that do not show appear more than 5 times in more than half the samples
wh0 = genefilter_sample(physeq, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq))
physeq_filt = prune_taxa(wh0, physeq)

# Transform to even sampling depth.
physeq_filt = transform_sample_counts(physeq_filt, function(x) 1E6 * x/sum(x))

# Keep only the most abundant five phyla.
phylum.sum = tapply(taxa_sums(physeq_filt), tax_table(physeq_filt)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
physeq_filt = prune_taxa((tax_table(physeq_filt)[, "Phylum"] %in% top5phyla), physeq_filt)

### dendrogram
tree <- plot_tree(physeq_filt, color = "Subgroup1", label.tips = "Genus", size = "abundance",
                  plot.margin = 0.2, ladderize = TRUE, sizebase = 2) +
  scale_colour_tableau("Classic 10") +
  ggtitle("Phylogenetic Tree - Genus") +
  theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Condition")

ggsave("treeGenus.pdf", plot = tree,
       device = "pdf", path = outDir, width = 15, height = 15,
       units = "in")

###--- barplot composition ---###
# scaling
# set seed
set.seed(1)  

# scaling the mouse data to the smallest samples. Note: rngseed is similar to set.seed
physeq_scaled <- rarefy_even_depth(physeq_filt,sample.size=19000, replace=FALSE, rngseed = 1) 

###--- Organizing treatments in an order ---###
sample_data(physeq_scaled)$Subgroup1 <- factor((sample_data(physeq_scaled)$Subgroup1), levels=c("Control","AD","NDGA","AD_NDGA"))

# labels for the samples
rownames <- sample_data(physeq_scaled)[,"customer_label"]
dim(otu_table(physeq_scaled))
colourCount <- length(unique(tax_table(physeq_scaled)[,"Genus"])) # number of genus for colors

getPalette = colorRampPalette(brewer.pal(10, "Paired")) # create color palette

###--- bar plot ---###
composition <- plot_bar(physeq_scaled, x = "customer_label", fill = "Genus") +
  facet_grid(~Subgroup1, scales="free_x") +
  scale_fill_manual(values=getPalette(colourCount))+
  ggtitle("Sample Composition - Genus") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave("sampleComposition.pdf", plot = composition,
       device = "pdf", path = outDir, width = 12, height = 8.5,
       units = "in")


#p <- plot_composition(physeq, 
#                      "Family", 
#                      numberOfTaxa = 30, fill = Subgroup1)
## plot facetting
#p <- p + facet_wrap(~EnvType, 
#scales = "free_x"
#, nrow = 1)
#plot(p)


###--- Heatmap ---###

# from here: https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

metadata<-read_q2metadata("q2metadata.tsv")
countsTable<-read_qza("./DADA2_denoising_output/table.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data

countsTable<-apply(countsTable, 2, function(x) x/sum(x)*100) # convert to percent

countsToPlot<-  
  data.frame(MeanAbundance=rowMeans(countsTable)) %>% # find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(40, MeanAbundance) %>%
  pull(Feature.ID) # extract only the names from the table

heatmap <- countsTable %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% countsToPlot,  Feature.ID, "Remainder")) %>% # flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=customer_label, y=Taxon, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~Subgroup1, scales="free_x") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")

ggsave("heatmapSpecies.pdf", plot = heatmap,
       device = "pdf", path = outDir, width = 14, height = 8.5,
       units = "in")

