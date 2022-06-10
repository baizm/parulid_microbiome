source activate qiime2-2020.8

##group by mean ceiling (average OTU count by species--rounds up so you don't have fractions of OTUs)
qiime feature-table group \
--i-table ./core-metrics-results_filtered/rarefied_table.qza \
--o-grouped-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
--p-mode mean-ceiling \
--m-metadata-file ../../../../metadata2_16s_ny.tsv \
--m-metadata-column Species \
--p-axis sample

qiime feature-table summarize \
  --i-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
  --o-visualization ./core-metrics-results_filtered/rarefied_table_grouped.qzv

#calculte diversity metrics on grouped table
qiime diversity beta-phylogenetic \
  --i-phylogeny rooted-tree_silva.qza \
  --i-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
  --p-metric unweighted_unifrac \
  --o-distance-matrix ./core-metrics-results_filtered/grouped_beta/grouped_unifrac.qza

qiime tools export \
  --input-path ./core-metrics-results_filtered/grouped_beta/grouped_unifrac.qza \
  --output-path ./core-metrics-results_filtered/grouped_beta/exported_unifrac

qiime diversity beta-phylogenetic \
  --i-phylogeny rooted-tree_silva.qza \
  --i-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
  --p-metric weighted_unifrac \
  --o-distance-matrix ./core-metrics-results_filtered/grouped_beta/grouped_weightedunifrac.qza

qiime tools export \
  --input-path ./core-metrics-results_filtered/grouped_beta/grouped_weightedunifrac.qza \
  --output-path ./core-metrics-results_filtered/grouped_beta/exported_weightedunifrac

qiime diversity beta \
  --i-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
  --p-metric braycurtis \
  --o-distance-matrix ./core-metrics-results_filtered/grouped_beta/grouped_braycurtis.qza

qiime tools export \
  --input-path ./core-metrics-results_filtered/grouped_beta/grouped_braycurtis.qza \
  --output-path ./core-metrics-results_filtered/grouped_beta/exported_braycurtis

qiime diversity beta \
  --i-table ./core-metrics-results_filtered/rarefied_table_grouped.qza \
  --p-metric jaccard \
  --o-distance-matrix ./core-metrics-results_filtered/grouped_beta/grouped_jaccard.qza

qiime tools export \
  --input-path ./core-metrics-results_filtered/grouped_beta/grouped_jaccard.qza \
  --output-path ./core-metrics-results_filtered/grouped_beta/exported_jaccard
