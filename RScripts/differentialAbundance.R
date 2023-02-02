#####   Name: differentalAbundance.R
#####   Version: 1.0
#####   Description: Differential abundance analysis across groups/samples
#####   Date of Creation: Jan-2023
#####   Author: Natália Faraj Murad
#####   E-MAIL: nataliafmurad@gmail.com
#####   PROJECT: https://github.com/natmurad/amplicon16S

###--- Load libraries ---###
library(tidyverse) # data manipulation
library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(EnhancedVolcano) # plot volcano
library(DT) # visualize tables
library(htmltools) # export html
library(ggrepel) # for offset labels
library(ggplot2)
library(ggthemes) # colors
library(ggpubr) # add stats on plots
library(qiime2R) # read qiime2 outputs
library(metagenomeSeq) # for differential abundance analysis - zero inflated model

###--- Graph config ---###
theme_set(theme_bw())

###--- Set directory ---###
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

###--- Conversions needed ---###
##--- metagenomeseq ---##
mtseq <- phyloseq_to_metagenomeSeq(physeq)

####################################
###--- Differential Abundance ---###
####################################

###--- Filtering ---###
#filterData(obj, present = 10, depth = 1000)

###--- Calculating normalization factors ---###
p = cumNormStatFast(mtseq)
mtseq = cumNorm(mtseq, p = p)

###--- Statistical testing ---###
###--- Zero-inflated Gaussian Mixture model ---###


normFactor = normFactors(mtseq)
normFactor = log2(normFactor/median(normFactor) + 1)

# maxit=1 is for demonstration purposes
settings = zigControl(maxit = 1, verbose = FALSE)
mod = model.matrix(~ 0 + Subgroup1, pData(mtseq))
#colnames(mod) = levels(as.factor(treat))
# fitting the ZIG model
res = fitZig(obj = mtseq, mod = mod, control = settings)

# The output of fitZig contains a list of various useful
# items. hint: names(res). Probably the most useful is the
# limma 'MLArrayLM' object called fit.
zigFit = slot(res, "fit")
finalMod = slot(res, "fit")$design

# Ad x Control; Ad_Drug x Ad; Drug x Control

contrasts <- combn(colnames(finalMod)[1:length(colnames(finalMod))-1], 2)
updn_cols <- c(RColorBrewer::brewer.pal(6, 'Blues')[3], RColorBrewer::brewer.pal(6, 'Reds')[3])

for(i in 1:ncol(contrasts)){
  x <- paste0(contrasts[2,i], "-", contrasts[1,i])
  contrast.matrix = makeContrasts(x, levels  = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2 = eBayes(fit2)
  table <- topTable(fit2, sort="none",n=Inf)
  table$seq <- rownames(table)
  table <- merge(table, fData(mtseq)[,c("Genus", "Species")], by = "row.names")
  write.csv(table, paste0(outDir, gsub("Subgroup1", "", paste0(contrasts[2,i], "vs", contrasts[1,i])), '.csv'))
  html <-  table %>%
    datatable(caption = gsub("Subgroup1", "", paste0(contrasts[2,i], " vs ", contrasts[1,i]))) %>%
    DT::formatStyle('logFC',
                    valueColumns = 'logFC',
                    backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
    DT::formatSignif(1:6, digits = 4)
  save_html(html, file = paste0(outDir, gsub("Subgroup1", "", paste0(contrasts[2,i], "vs", contrasts[1,i])), '.html'), background = "white", libdir = "lib", lang = "en") 
  
  ###--- Volcano ---###
  volcano <- EnhancedVolcano(table,
                             lab = paste0(table[,"Genus"], " ", table[,"Species"]),
                             x = 'logFC',
                             y = 'adj.P.Val',
                             labSize = 3,
                             xlim = c(-10, 10),
                             title = gsub("Subgroup1", "", paste0(contrasts[2,i], " vs ", contrasts[1,i])),
                             subtitle = "LogFC vs -Log10Pvalue",
                             pCutoff = 0.05,
                             col=c('#030000', '#15B01A', '#0343DF', '#E50000')
  )
  file <- gsub("Subgroup1", "", paste0("volcano", contrasts[2,i], "vs", contrasts[1,i], ".pdf"))
  ggsave(file, plot = volcano,
         device = "pdf", path = outDir, width = 7, height = 7,
         units = "in")
}

for(i in 1:ncol(contrasts)){
  x <- paste0(contrasts[1,i], "-", contrasts[2,i])
  contrast.matrix = makeContrasts(x, levels  = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2 = eBayes(fit2)
  table <- topTable(fit2, sort="none",n=Inf)
  table$seq <- rownames(table)
  table <- merge(table, fData(mtseq)[,c("Genus", "Species")], by = "row.names")
  write.csv(table, paste0(outDir, gsub("Subgroup1", "", paste0(contrasts[1,i], "vs", contrasts[2,i])), '.csv'))
  html <-  table %>%
    datatable(caption = gsub("Subgroup1", "", paste0(contrasts[1,i], " vs ", contrasts[2,i]))) %>%
    DT::formatStyle('logFC',
                    valueColumns = 'logFC',
                    backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
    DT::formatSignif(1:6, digits = 4)
  save_html(html, file = paste0(outDir, gsub("Subgroup1", "", paste0(contrasts[1,i], "vs", contrasts[2,i])), '.html'), background = "white", libdir = "lib", lang = "en") 
  
  ###--- Volcano ---###
  volcano <- EnhancedVolcano(table,
                             lab = paste0(table[,"Genus"], " ", table[,"Species"]),
                             x = 'logFC',
                             y = 'adj.P.Val',
                             labSize = 3,
                             xlim = c(-10, 10),
                             title = gsub("Subgroup1", "", paste0(contrasts[1,i], " vs ", contrasts[2,i])),
                             subtitle = "LogFC vs -Log10Pvalue",
                             pCutoff = 0.05,
                             col=c('#030000', '#15B01A', '#0343DF', '#E50000')
  )
  file <- gsub("Subgroup1", "", paste0("volcano", contrasts[1,i], "vs", contrasts[2,i], ".pdf"))
  ggsave(file, plot = volcano,
         device = "pdf", path = outDir, width = 7, height = 7,
         units = "in")
}