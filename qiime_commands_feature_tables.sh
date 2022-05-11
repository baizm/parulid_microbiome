module load python/3.6.8
module load anaconda3
source activate /storage/home/mdb5946/.conda/envs/qiime2-2020.8

#######-----2017-2019 samples (Batch 1)------######

#filtered out files for test samples and recaptures and re-import (raw_reads directory)
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path /gpfs/group/dut374/default/marcella/16s2_ny/raw_reads/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /gpfs/group/dut374/default/marcella/16s2_ny/demux-seqs.qza
#took ~15 mins

qiime dada2 denoise-paired \
  --p-n-threads 0 \
  --i-demultiplexed-seqs /gpfs/group/dut374/default/marcella/16s2_ny/demux-seqs.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --output-dir /gpfs/group/dut374/default/marcella/16s2_ny/denoised/ \
  --o-table feature_table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#use the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva_classifier.qza \
  --i-reads ../rep-seqs.qza \
  --p-n-jobs 4 \
  --o-classification silva_taxonomy.qza

qiime metadata tabulate \
  --m-input-file silva_taxonomy.qza \
  --o-visualization silva_taxonomy.qzv

#visualize classified sequences
qiime taxa barplot \
  --i-table ../feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --m-metadata-file ../metadata2_16s_ny.tsv \
  --o-visualization taxa-bar-plots.qzv

#filter out mithochondrial, chloroplast and unassigned ASVs, AND eukaryota
qiime taxa filter-table \
  --i-table ../feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --p-exclude mitochondria,chloroplast,unassigned,eukaryota \
  --o-filtered-table feature_table_no_MtCpUnEuk.qza

#summarize new feature tabulate
qiime feature-table summarize \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --o-visualization feature_table_no_MtCpUnEuk.qzv

#filter seqs to remove mt,cp,un for phylogeny
qiime feature-table filter-seqs \
  --i-data ../rep-seqs.qza \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --o-filtered-data rep-seqs_no_MtCpUnEuk.qza

#visualize
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_no_MtCpUnEuk.qza \
  --o-visualization rep-seqs_no_MtCpUnEuk.qzv

#generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_no_MtCpUnEuk.qza \
  --o-alignment aligned-rep-seqs_no_MtCpUnEuk.qza \
  --o-masked-alignment masked-aligned-rep-seqs_no_MtCpUnEuk.qza \
  --o-tree unrooted-tree_silva.qza \
  --o-rooted-tree rooted-tree_silva.qza

#produce alpha-rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --i-phylogeny rooted-tree_silva.qza \
  --p-max-depth 21000 \
  --m-metadata-file ../metadata2_16s_ny.tsv \
  --o-visualization alpha-rarefaction_silva.qzv

#for decontam:
#export feature table to otu table (puts it in a folder called out-table.biom)
qiime tools export \
  --input-path feature_table_no_MtCpUnEuk.qza \
  --output-path otu_table_silva.biom

#convert to tsv file for decontam
biom convert -i  ./otu_table_silva.biom/feature-table.biom -o ./otu_table_silva.biom/feature-table.tsv --to-tsv

#export tree for decontam
qiime tools export \
  --input-path rooted-tree_silva.qza \
  --output-path exported-tree_silva

####--before continuing to next step, work in decontam.R to ID contaminant ASVs to remove (listed in asv_to_remove.txt)---####

#filter out contaminant ASVs from decontam
qiime feature-table filter-features \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --p-exclude-ids \
  --m-metadata-file asv_to_remove.txt \
  --o-filtered-table feature_table_no_contam.qza

#summarize new feature tabulate
qiime feature-table summarize \
  --i-table feature_table_no_contam.qza \
  --o-visualization feature_table_no_contam.qzv

#filter seqs to remove contaminants
qiime feature-table filter-seqs \
  --i-data ../rep-seqs.qza \
  --i-table feature_table_no_contam.qza \
  --o-filtered-data rep-seqs_no_contam.qza

#visualize
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_no_contam.qza \
  --o-visualization rep-seqs_no_contam.qzv
  
 ####------2020 samples (Batch 2)-------#######
 
 qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path /gpfs/group/dut374/default/metabarcoding_data/03-12-2021/16s/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /gpfs/group/dut374/default/marcella/16s_2020samples/demux-seqs.qza
#took ~15 mins

qiime dada2 denoise-paired \
  --p-n-threads 0 \
  --i-demultiplexed-seqs /gpfs/group/dut374/default/marcella/16s_2020samples/demux-seqs.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --output-dir /gpfs/group/dut374/default/marcella/16s_2020samples/denoised/ \
  --o-table feature_table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#use the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier /gpfs/group/dut374/default/marcella/16s2_ny/silva/silva_classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 4 \
  --o-classification silva_taxonomy.qza

qiime metadata tabulate \
  --m-input-file silva_taxonomy.qza \
  --o-visualization silva_taxonomy.qzv

#visualize classified sequences
qiime taxa barplot \
  --i-table feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --m-metadata-file metadata_16s_2020samples \
  --o-visualization taxa-bar-plots.qzv

#for decontam:
#export feature table to otu table for decontam (puts it in a folder called otu-table)
qiime tools export \
  --input-path feature_table.qza \
  --output-path otu_table_forDecontam

#convert to tsv file for decontam
biom convert -i  ./otu_table_forDecontam/feature-table.biom -o ./otu_table_forDecontam/feature-table.tsv --to-tsv

#generate a tree for decontam
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree_forDecontam.qza \
  --o-rooted-tree rooted-tree_forDecontam.qza

#export tree for decontam
qiime tools export \
  --input-path rooted-tree_forDecontam.qza \
  --output-path exported-tree_forDecontam

####--before continuing to next step, work in decontam.R to ID contaminant ASVs to remove (listed in asv_to_remove.txt)---####

#pick up here after using R decontam to ID contaminant ASVs to remove
#filter out mithochondrial, chloroplast and unassigned ASVs, AND eukaryota
qiime taxa filter-table \
  --i-table feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --p-exclude mitochondria,chloroplast,unassigned,eukaryota \
  --o-filtered-table feature_table_no_MtCpUnEuk.qza

#filter out contaminant ASVs from decontam
qiime feature-table filter-features \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --p-exclude-ids \
  --m-metadata-file asv_to_remove.txt \
  --o-filtered-table feature_table_no_contam.qza

#summarize new feature tabulate
qiime feature-table summarize \
  --i-table feature_table_no_contam.qza \
  --o-visualization feature_table_no_contam.qzv

#filter seqs to remove mt,cp,un,euk,contam for phylogeny
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table feature_table_no_contam.qza \
  --o-filtered-data rep-seqs_no_contam.qza

#visualize
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_no_contam.qza \
  --o-visualization rep-seqs_no_contam.qzv

####----merge feature tables from batch 1 and batch 2, produce final feature table and other files----####

#merge filtererd data from 2017-2019 batch with 2020 batch
qiime feature-table merge \
  --i-tables /gpfs/group/dut374/default/marcella/16s_2020samples/feature_table_no_contam.qza \
  --i-tables /gpfs/group/dut374/default/marcella/16s2_ny/silva/feature_table_no_contam.qza \
  --o-merged-table feature_table_merged.qza

qiime feature-table merge-seqs \
  --i-data /gpfs/group/dut374/default/marcella/16s_2020samples/rep-seqs_no_contam.qza \
  --i-data /gpfs/group/dut374/default/marcella/16s2_ny/silva/rep-seqs_no_contam.qza \
  --o-merged-data rep-seqs_merged.qza

qiime feature-table merge-taxa \
  --i-data /gpfs/group/dut374/default/marcella/16s_2020samples/silva_taxonomy.qza \
  --i-data /gpfs/group/dut374/default/marcella/16s2_ny/silva/silva_taxonomy.qza \
  --o-merged-data silva_taxonomy_merged.qza

#summarize new feature tabulate
qiime feature-table summarize \
  --i-table feature_table_merged.qza \
  --o-visualization feature_table_merged.qzv

#visualize, should have  sequences (subtract cp,mt,un,euk,contam) - !
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_merged.qza \
  --o-visualization rep-seqs_merged.qzv

qiime metadata tabulate \
  --m-input-file silva_taxonomy_merged.qza \
  --o-visualization silva_taxonomy_merged.qzv
#don't forget to download tsv of silva_taxonomy.qzv from qiime2 view

#generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_merged.qza \
  --o-alignment aligned-rep-seqs_merged.qza \
  --o-masked-alignment masked-aligned-rep-seqs_merged.qza \
  --o-tree unrooted-tree_merged.qza \
  --o-rooted-tree rooted-tree_merged.qza

#export feature table to otu table (puts it in a folder called out-table.biom)
qiime tools export \
  --input-path feature_table_merged.qza \
  --output-path otu_table_merged

#convert to tsv file for phyloseq
biom convert -i  ./otu_table_merged/feature-table.biom -o ./otu_table_merged/feature-table.tsv --to-tsv

#export tree for phyloseq
qiime tools export \
  --input-path rooted-tree_merged.qza \
  --output-path exported-tree_merged
