import os
import numpy as np
import nibabel as nib
from neuromaps import datasets, stats
from statsmodels.stats import multitest
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import csv

# Define paths and constants
path = os.path.expanduser("~") + "/local_data/github/cbfreceptors/"
parcellation_atlas = os.path.join(path, 'data/atlas/Schaefer2018_100Parcels_7Networks_Xiao_2019_SubCorSeg_resampled_asl.nii')

#supplementary analysis
#parcellation_atlas = path+'data/atlas/Schaefer2018_100Parcels_7Networks_2mm' #cortical only
#parcellation_atlas = path+'data/atlas/DKT_space-MNI152NLin6_res-2x2x2_no_brainstem.nii' #DK

scale = '122'
asl_maps = ['ssd_vs_hc', 'chr_p_vs_hc', 'chr_t_vs_chr_nt'] #main analysis
receptor_names = np.load(os.path.join(path, 'data/receptor_names_pet.npy'))

# Load data from CSV and numpy files
map_data = {map_name: [] for map_name in asl_maps}
with open(os.path.join(path, 'data/asl_parcellated_'+scale+'.csv'), newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        for i, map_name in enumerate(asl_maps):
            map_data[map_name].append(float(row[i]))

receptor_data = np.genfromtxt(os.path.join(path, 'data/receptors_parcellated_'+scale+'.csv'), delimiter=',')
parcellated_receptors = {name: data for name, data in zip(receptor_names, receptor_data.T)}

#colormap
sexyRdBu = matplotlib.colors.LinearSegmentedColormap.from_list("", [
    "#00B7DF",#blue
    "#0FD4FF",
    "#31DAFF",
    "#57E1FF",
    "#7BE7FF",
    "#B1F1FF",
    "#ffffff",
    "#FFBFB1",
    "#FF937B",
    "#FF7557",
    "#FF5631",
    "#FF3A0F",
    "#DF2800"#red
])

# Load nulls for correlations
nulls = {i: np.load(os.path.join(path, 'data/saved_nulls/', asl_maps[i]+'_'+scale+'.npy')) for i in range(len(asl_maps))}

# Run correlations
results = []
for i, map_name in enumerate(asl_maps):
    result_dict = {}
    for receptor_name, parcellated_map in parcellated_receptors.items():
        result = stats.compare_images(map_data[map_name], parcellated_map, metric='spearmanr', nulls=nulls[i])
        result_dict[receptor_name] = result
    results.append((map_name, result_dict))

corr_df_indices = [v[0] for v in results]
corr_df_columns = results[0][1].keys()
corr_df_data = [[r[1][c][0] for c in corr_df_columns] for r in results]
corr_table = pd.DataFrame(index=corr_df_indices, columns=corr_df_columns, data=corr_df_data)

# Extract the p-values
p_values_data = [[r[1][c][1] for c in corr_df_columns] for r in results]

# Perform FDR correction on the p-values
fdr_corrected = multitest.fdrcorrection(np.concatenate(p_values_data), alpha=0.05)
fdr_corrected_p_values = fdr_corrected[1]

# Reshape the corrected p-values back into the original DataFrame shape
fdr_corrected_p_values_matrix = np.array(fdr_corrected_p_values).reshape(len(corr_df_indices), len(corr_df_columns))

# Create the final DataFrame with corrected p-values
fdr_p_values = pd.DataFrame(index=corr_df_indices, columns=corr_df_columns, data=fdr_corrected_p_values_matrix)

# Determine significant indices based on FDR-corrected p-values
significant_indices = []
alpha = 0.05  # Adjust alpha as needed
for i, row in enumerate(fdr_p_values.values):
    significant_indices.append([j for j, p_value in enumerate(row) if p_value < alpha])
    
corr_table_subset = corr_table.iloc[:3, :]

# Create heatmap of correlations
sns.heatmap(corr_table_subset, cmap=sexyRdBu, annot=False, center=0, vmin=-1, vmax=1, linewidth=.5)

# Mark significance level on the heatmap
for i, indices in enumerate(significant_indices):
    for j in indices:
        plt.text(j + 0.5, i + 0.5, "*", ha="center", va="center", color="white", fontsize=14)

print(corr_table_subset)

# Save the heatmap as an image file
plt.savefig(path+'figures/scale'+scale+'/spearmans_correlation_heatmap.svg')
plt.title("Spearman's Correlation Heatmap")
plt.show()