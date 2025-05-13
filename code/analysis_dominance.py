# -*- coding: utf-8 -*-
"""
ASL Patient vs Control difference map x PET maps Dominance analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netneurotools import datasets, stats, utils
from scipy.stats import zscore
from matplotlib.colors import ListedColormap
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
import matplotlib

def get_reg_r_sq(X, y):
    lin_reg = LinearRegression()
    lin_reg.fit(X, y)
    yhat = lin_reg.predict(X)
    SS_Residual = sum((y - yhat) ** 2)
    SS_Total = sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (float(SS_Residual)) / SS_Total
    adjusted_r_squared = 1 - (1 - r_squared) * \
        (len(y) - 1) / (len(y) - X.shape[1] - 1)
    return adjusted_r_squared

def get_perm_p(emp, null):
    return (1 + sum(abs(null - np.mean(null))
                    > abs(emp - np.mean(null)))) / (len(null) + 1)

def get_reg_r_pval(X, y, smash, nsmash):
    emp = get_reg_r_sq(X, y)
    null = np.zeros((nsmash, ))
    for s in range(nsmash):
        null[s] = get_reg_r_sq(X, smash[:, s])
    return (1 + sum(null > emp))/(nsmash + 1)

#colormap
sexyRdBu = matplotlib.colors.LinearSegmentedColormap.from_list("", [
    #"#00B7DF",#blue
    #"#0FD4FF",
    #"#31DAFF",
    #"#57E1FF",
    #"#7BE7FF",
    #"#B1F1FF",
    "#ffffff",
    "#FFBFB1",
    "#FF937B",
    "#FF7557",
    "#FF5631",
    "#FF3A0F",
    "#DF2800"#red
])

"""
set-up
"""

scale = '122'
path = os.path.expanduser("~") + "/local_data/cbfreceptors-main/"

#brainsmash implementation for cort+subcort maps. Ensure nulls match parcellated_asl_scale.
nsmash = 5000
smash = dict([])
smash[0] = np.load(path+'data/saved_nulls/ssd_vs_hc_122.npy')
smash[1] = np.load(path+'data/saved_nulls/chr_p_vs_hc_122.npy')
smash[2] = np.load(path+'data/saved_nulls/chr_t_vs_chr_nt_122.npy')

#supplementary analysis
#smash[0] = np.load(path+'data/saved_nulls/ch_t_vs_hc_122.npy')
#smash[1] = np.load(path+'data/saved_nulls/chr_nt_vs_hc_122.npy')
#smash[2] = np.load(path+'data/saved_nulls/ssd_vs_hc_age_sex_122.npy')
#smash[3] = np.load(path+'data/saved_nulls/ssd_vs_hc_age_sex_med_122.npy')
#smash[4] = np.load(path+'data/saved_nulls/chr_p_vs_hc_age_sex_122.npy')
#smash[5] = np.load(path+'data/saved_nulls/chr_p_vs_hc_age_sex_med_122.npy')


# load the receptor data
receptor_data = np.genfromtxt(path+'data/receptors_parcellated_'+scale+'.csv', delimiter=',')
receptor_names = np.load(path+'data/receptor_names_pet.npy')

# load the asl maps
#asl = np.genfromtxt(path+'data/ut/asl_'+scale+'.csv', delimiter=',')
asl = np.genfromtxt(path+'data/asl_parcellated_'+scale+'.csv', delimiter=',')
asl = np.delete(asl, (3,4,5,6,7,8), axis=1)
# 0 = 'SSD vs HC'
# 1 = 'CHR-P vs HC'
# 2 = 'CHR-T vs CHR-NT'

#supplementary analysis
# 3 = 'CHR-T vs HC'
# 4 = 'CHR-NT vs HC'
# 5 = 'SSD vs HC age, sex'
# 6 = 'SSD vs HC age, sex, med'
# 7 = 'CHR vs HC age, sex'
# 8 = 'CHR vs HC age, sex, med'

# list of disorders included in analysis. 
disorders = ['SSD vs HC',
             'CHR-P vs HC',
             'CHR-T vs CHR-NT',
             
             #supplementary analysis
             #'CHR-T vs HC',
             #'CHR-NT vs HC',
             #'SSD vs HC age, sex',
             #'SSD vs HC age, sex, med',
             #'CHR-P vs HC age, sex',
             #'CHR-P vs HC age, sex, med'
            ]

"""
Dominance analysis
"""

model_metrics = dict([])
model_pval = np.zeros((len(disorders), ))

for i in range(len(disorders)):
    print(i)
    m, _ = stats.get_dominance_stats(zscore(receptor_data),
                                     zscore(asl[:, i]))
    model_metrics[disorders[i]] = m

    # get p-value of model
    model_pval[i] = get_reg_r_pval(zscore(receptor_data),
                                   zscore(asl[:, i]), zscore(smash[i]), nsmash)
                                  # spins, nspins)

model_pval = multipletests(model_pval, method='fdr_bh')[1]

dominance = np.zeros((len(disorders), len(receptor_names)))

for i in range(len(model_metrics)):
    tmp = model_metrics[disorders[i]]
    dominance[i, :] = tmp["total_dominance"]
np.save(path+'results/dominance/scale_'+scale+'/dominance.npy', dominance)

# Bar plot of model importance
plt.ion()
plt.figure(figsize=(8, 4))
plt.barh(np.arange(len(disorders)), np.sum(dominance, axis=1), tick_label=disorders, color="#FF7557")
plt.xticks(rotation='horizontal')
plt.tight_layout()
plt.savefig(path + f'figures/scale_{scale}/bar_dominance.svg')

# Heatmap of dominance (normalized)
dominance[np.where(model_pval >= 0.05)[0], :] = 0
dominance_norm = dominance / np.sum(dominance, axis=1)[:, None]

fig, ax = plt.subplots(figsize=(12, 6))
im = ax.imshow(dominance_norm, cmap=sexyRdBu, vmin=0, vmax=dominance_norm.max(), aspect='auto', interpolation='none')

# Axis ticks
ax.set_xticks(np.arange(len(receptor_names)))
ax.set_yticks(np.arange(len(disorders)))
ax.set_xticklabels(receptor_names, rotation=45, ha="right")
ax.set_yticklabels(disorders)

# Remove spines and ticks
for edge in ["top", "bottom", "left", "right"]:
    ax.spines[edge].set_visible(False)
ax.tick_params(top=False, bottom=False, left=False, right=False)

# Colorbar
cbar = ax.figure.colorbar(im, ax=ax)
cbar.outline.set_visible(False)
cbar.ax.set_ylabel("Normalized Dominance", rotation=-90, va="bottom")
cbar.ax.tick_params(size=0)

plt.title("Dominance Analysis Heatmap")
plt.tight_layout()

# Save figure cleanly
output_path = path + f'figures/scale{scale}/heatmap_dominance.svg'
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
