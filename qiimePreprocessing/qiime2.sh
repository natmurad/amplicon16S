#! usr/bin/bash

# download env config
curl -sL \
  "https://data.qiime2.org/distro/core/qiime2-2020.8-py36-osx-conda.yml" > \
  "qiime2.yml"

# create conda env
conda env create -n qiime2 --file qiime2.yml

# remove env file
rm qiime2.yml

# activate conda env
conda activate qiime2

# import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manifest1.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux.qza

# summariza data and create a visualization for the quality of the reads
qiime demux summarize --i-data demux.qza \
--o-visualization demux.qzv

# trimm primers off
qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GACTACHVGGGTATCTAATCC \
--p-discard-untrimmed \
--o-trimmed-sequences demux-trimmed.qza \
--verbose

# denoise with DADA2
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-trimmed.qza \
--p-trunc-len-f 280 \
--p-trunc-len-r 170 \
--p-n-threads 4 \
--output-dir DADA2_denoising_output \
--verbose

# summarize outputs
# table
qiime feature-table summarize \
--i-table /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/table.qza \
--o-visualization table \
--m-sample-metadata-file metadata.tsv

# rep_seqs
qiime feature-table tabulate-seqs \
--i-data /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/representative_sequences.qza \
--o-visualization rep_seqs

# denoising stats
qiime metadata tabulate \
--m-input-file /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/denoising_stats.qza \
--o-visualization denoising_stats

# build phylogenetic tree
# align with mafft
qiime alignment mafft  \
--i-sequences /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/representative_sequences.qza \
--o-alignment aligned_rep_seqs.qza

# mask alignment
qiime alignment mask --i-alignment aligned_rep_seqs.qza \
--o-masked-alignment masked_aligned_rep_seqs.qza

# build an unrooted tree
qiime phylogeny fasttree \
--i-alignment masked_aligned_rep_seqs.qza --o-tree unrooted_tree

# root the tree
qiime phylogeny midpoint-root \
--i-tree unrooted_tree.qza --o-rooted-tree rooted_tree

# diversity analyses
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted_tree.qza \
--i-table /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/table.qza \
--p-sampling-depth 19000 --m-metadata-file metadata.tsv \
--output-dir core_metrics_results

# create alpha div visualizations
qiime diversity alpha-group-significance \
--i-alpha-diversity core_metrics_results/faith_pd_vector.qza \
--m-metadata-file metadata.tsv \
--o-visualization core_metrics_results/faith_pd_group_significance

qiime diversity alpha-group-significance \
--i-alpha-diversity core_metrics_results/evenness_vector.qza \
--m-metadata-file metadata.tsv \
 --o-visualization core_metrics_results/evenness_group_significance

qiime diversity alpha-group-significance \
--i-alpha-diversity core_metrics_results/shannon_vector.qza --m-metadata-file metadata.tsv \
--o-visualization core_metrics_results/shannon_group_significance

# alpha rarefaction plots
qiime diversity alpha-rarefaction \
--i-table core_metrics_results/rarefied_table.qza \
--p-max-depth 30000 --m-metadata-file metadata.tsv \
--p-steps 25 --o-visualization alpha_rarefaction.qzv

# train classifier
# download silva files https://www.arb-silva.de/download/archive/qiime/

# import silva files
qiime tools import --type 'FeatureData[Sequence]' \
--input-path /Users/natmurad/Documents/Buck/Dipa16S/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path silva132_99

qiime tools import --type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path /Users/natmurad/Documents/Buck/Dipa16S/SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt \
--output-path silva132_99_ref_taxonomy

# extract ref reads
qiime feature-classifier extract-reads --i-sequences silva132_99.qza \
--p-f-primer GTGCCAGCMGCCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT \
--p-trunc-len 300 --p-n-jobs 5 --o-reads ref_seqs

# train classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref_seqs.qza \
--i-reference-taxonomy silva132_99_ref_taxonomy.qza \
--o-classifier classifier.qza

# assign taxonomy
qiime feature-classifier classify-sklearn \
--i-classifier /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/silva-138-99-nb-classifier.qza \
--i-reads /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/representative_sequences.qza \
--o-classification taxonomy.qza

# create visualization
qiime metadata tabulate --m-input-file taxonomy.qza \
--o-visualization taxonomy

# create barplots
qiime taxa barplot --i-table core_metrics_results/rarefied_table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file metadata.tsv --o-visualization taxa_bar_plots.qzv

# export to tsv
qiime tools export \
--input-path /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/DADA2_denoising_output/table.qza \
--output-path export/table

biom convert -i /Users/natmurad/Documents/Buck/Dipa16S/emp-paired-end-sequences/export/table/feature-table.biom \
-o table.tsv --to-ts