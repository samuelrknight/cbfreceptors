import os
import numpy as np
import nibabel as nib
from neuromaps import datasets, stats
from statsmodels.stats import multitest
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
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

# Custom colormap
sexyRdBu = matplotlib.colors.LinearSegmentedColormap.from_list("", [
    "#00B7DF", "#0FD4FF", "#31DAFF", "#57E1FF", "#7BE7FF", "#B1F1FF",
    "#ffffff",
    "#FFBFB1", "#FF937B", "#FF7557", "#FF5631", "#FF3A0F", "#DF2800"
])

# Load nulls
nulls = {i: np.load(os.path.join(path, 'data/saved_nulls/', f'{asl_maps[i]}_{scale}.npy')) for i in range(len(asl_maps))}

# Run correlations
results = []
for i, map_name in enumerate(asl_maps):
    result_dict = {}
    for receptor_name, parcellated_map in parcellated_receptors.items():
        result = stats.compare_images(map_data[map_name], parcellated_map, metric='spearmanr', nulls=nulls[i])
        result_dict[receptor_name] = result
    results.append((map_name, result_dict))

# Convert to DataFrames
corr_df_indices = [v[0] for v in results]
corr_df_columns = results[0][1].keys()
corr_df_data = [[r[1][c][0] for c in corr_df_columns] for r in results]
corr_table = pd.DataFrame(index=corr_df_indices, columns=corr_df_columns, data=corr_df_data)

# FDR correction
p_values_data = [[r[1][c][1] for c in corr_df_columns] for r in results]
fdr_corrected = multitest.fdrcorrection(np.concatenate(p_values_data), alpha=0.05)
fdr_corrected_p_values = fdr_corrected[1]
fdr_corrected_p_values_matrix = np.array(fdr_corrected_p_values).reshape(len(corr_df_indices), len(corr_df_columns))
fdr_p_values = pd.DataFrame(index=corr_df_indices, columns=corr_df_columns, data=fdr_corrected_p_values_matrix)

# Identify significant values
significant_indices = []
alpha = 0.05
for i, row in enumerate(fdr_p_values.values):
    significant_indices.append([j for j, p in enumerate(row) if p < alpha])

# Subset of correlation table for plotting
corr_table_subset = corr_table.iloc[:3, :]
data = corr_table_subset.values.astype(float)

# Plot using matplotlib (no edges)
fig, ax = plt.subplots(figsize=(12, 6))
im = ax.imshow(data, cmap=sexyRdBu, vmin=-1, vmax=1, aspect='auto', interpolation='none')

# Set axis ticks and labels
ax.set_xticks(np.arange(len(corr_df_columns)))
ax.set_yticks(np.arange(len(corr_table_subset.index)))
ax.set_xticklabels(corr_df_columns, rotation=45, ha="right")
ax.set_yticklabels(corr_table_subset.index)

# Remove all spines (border lines around the plot)
for edge in ["top", "bottom", "left", "right"]:
    ax.spines[edge].set_visible(False)

# Remove tick marks
ax.tick_params(top=False, bottom=False, left=False, right=False)

# Add significance markers
for i, row in enumerate(data):
    for j, val in enumerate(row):
        color = "white" if abs(val) > 0.5 else "black"
        if j in significant_indices[i]:
            ax.text(j, i, "*", ha="center", va="center", color=color, fontsize=14)

# Add colorbar
cbar = ax.figure.colorbar(im, ax=ax)
cbar.outline.set_visible(False)

cbar.ax.set_ylabel("Spearman's ρ", rotation=-90, va="bottom")

plt.title("Spearman's Correlation Heatmap")
plt.tight_layout()

# Save the heatmap (no black border, white background)
output_path = os.path.join(path, f'figures/scale{scale}/spearmans_correlation_heatmap.svg')
os.makedirs(os.path.dirname(output_path), exist_ok=True)
plt.savefig(
    output_path,
    bbox_inches='tight',
    pad_inches=0.0,
    facecolor='white',
    edgecolor='none',
    dpi=600
)

plt.show()

# Print the correlation values
print(corr_table_subset)
