
# to run correlation analysis with pre-calculated null maps

module load anaconda3
source activate imt_env

imagingtranscriptomics -i .../imaging_transcriptomics/data/maps/scz_vs_hc_ut.nii.gz \
--nulls /imaging_transcriptomics/data/brainsmash_nulls_schaeferxiao/scz_vs_hc_122.npy \
-o outdir \
--geneset .../imaging_transcriptomics/data/jorstad_2023_top_celltype_markers_diffstable.gmt \
corr
