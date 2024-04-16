import os
import nibabel as nib
import numpy as np
import neuromaps
from neuromaps import datasets
from neuromaps import stats
from neuromaps.parcellate import Parcellater
import matplotlib
import pandas as pd

"""
set-up
"""

#setup
# Define paths
path = os.path.expanduser("~") + "/local_data/github/cbfreceptors/"

# Choose parcellation atlas
parcellation_atlas = path+'data/atlas/Schaefer2018_100Parcels_7Networks_Xiao_2019_SubCorSeg_resampled_asl.nii'

#supplementary analysis
#parcellation_atlas = path+'data/atlas/Schaefer2018_100Parcels_7Networks_2mm' #cortical only
#parcellation_atlas = path+'data/atlas/DKT_space-MNI152NLin6_res-2x2x2_no_brainstem.nii' #DK

# Choose scale
scale = '122' #number of ROI in atlas
space = 'MNI152'
resolution = '2mm'

parcellater = Parcellater(parcellation_atlas, space)

# Define contrast maps 
asl_maps = [
    ('ssd_vs_hc', 'cbfbirn_asl_SSD_vs_HC.nii'),
    ('chr_p_vs_hc', 'neutop_prod_asl_CHR_P_vs_HC.nii'),
    ('chr_t_vs_chr_nt', 'neutop_prod_asl_CHR_T_vs_CHR_NT.nii'),
    ('chr_t_vs_hc', 'neutop_prod_asl_CHR_T_vs_HC.nii'),
    ('chr_nt_vs_hc', 'neutop_prod_asl_CHR_NT_vs_HC.nii'),
    ('ssd_vs_hc_age_sex', 'cbfbirn_asl_SSD_vs_HC_sex_age.nii'),
    ('ssd_vs_hc_age_sex_med', 'cbfbirn_asl_SSD_vs_HC_sex_age_med.nii'),
    ('chr_p_vs_hc_age_sex', 'neutop_prod_asl_CHR_P_vs_HC_age_sex.nii'),
    ('chr_p_vs_hc_age_sex_med', 'neutop_prod_asl_CHR_P_vs_HC_age_sex_med.nii'),
    # Add more maps here
]

# Define receptor maps 
receptor_maps = [
    ('serotonin1a', '5HT1a_way_hc36_savli.nii'),
    ('serotonin1b', '5HT1b_p943_hc65_gallezot.nii'),
    ('serotonin2a', '5HT2a_cimbi_hc29_beliveau.nii'),
    ('serotonin4', '5HT4_sb20_hc59_beliveau.nii'),
    ('serotonin6', '5HT6_gsk_hc30_radhakrishnan.nii'),
    ('serotoninT', '5HTT_dasb_hc100_beliveau.nii'),
    ('A4B2', 'A4B2_flubatine_hc30_hillmer.nii'),
    ('CB1', 'CB1_omar_hc77_normandin.nii'),
    ('D1', 'D1_SCH23390_hc13_kaller.nii'),
    ('D2', 'D2_flb457_hc55_sandiego.nii'),
    ('DAT', 'DAT_fpcit_hc174_dukart_spect.nii'),
    ('GABAa', 'GABAa-bz_flumazenil_hc16_norgaard.nii'),
    ('H3', 'H3_cban_hc8_gallezot.nii'),
    ('M1', 'M1_lsn_hc24_naganawa.nii'),
    ('mGluR5', 'mGluR5_abp_hc73_smart.nii'),
    ('MU', 'MU_carfentanil_hc204_kantonen.nii'),
    ('NET', 'NET_MRB_hc77_ding.nii'),
    ('NMDA', 'NMDA_ge179_hc29_galovic.nii'),
    ('VAChT', 'VAChT_feobv_hc18_aghourian_sum.nii')
    # Add more maps here
]

receptor_names = np.array(["5HT1a",
                           "5HT1b",
                           "5HT2a",
                           "5HT4",
                           "5HT6",
                           "5HTT",
                           "A4B2",
                           "CB1",
                           "D1",
                           "D2",
                           "DAT",
                           "GABAa",
                           "H3",
                           "M1",
                           "mGluR5",
                           "MOR",
                           "NET",
                           "NMDA",
                           "VAChT"])
np.save(path+'data/receptor_names_pet.npy', receptor_names)

"""
parcellate!
"""

#Parcellate asl maps 
parcellated_asl = {}

for map_name, map_path in asl_maps:
    full_map_path = os.path.join(path+'data/asl', map_path)
    parcellated_map = parcellater.fit_transform(full_map_path, space, ignore_background_data= True)
    parcellated_asl[map_name] = parcellated_map


# Create an empty DataFrame with columns as the asl names
asldf = pd.DataFrame(columns=parcellated_asl.keys())

# Populate the DataFrame with reshaped parcellated maps
for asl_name, parcellated_map in parcellated_asl.items():
    reshaped_map = parcellated_map.flatten()
    asldf[asl_name] = reshaped_map

# Save the DataFrame to a CSV file
csv_filename = path+'data/asl_parcellated_'+scale+'.csv'
asldf.to_csv(csv_filename, index=False, header = None)

#Parcellate receptor maps 
parcellated_receptors = {}

for map_name, map_path in receptor_maps:
    full_map_path = os.path.join(path+'data/pet_atlas', map_path)
    parcellated_map = parcellater.fit_transform(full_map_path, space, ignore_background_data= True)
    parcellated_receptors[map_name] = parcellated_map


# Create an empty DataFrame with columns as the receptor names
rdf = pd.DataFrame(columns=parcellated_receptors.keys())

# Populate the DataFrame with reshaped parcellated maps
for receptor_name, parcellated_map in parcellated_receptors.items():
    reshaped_map = parcellated_map.flatten()
    rdf[receptor_name] = reshaped_map

# Save the DataFrame to a CSV file
csv_filename = path+'data/receptors_parcellated_'+scale+'.csv'
rdf.to_csv(csv_filename, index=False, header = None)

print(f"Saved parcellated_receptors to {csv_filename}")    
    
print("Finished parcellating asl and receptor maps")

"""
generate nulls
"""

print("Generating nulls. Please wait!")

# Generate nulls with brainsmash burt2020
for i, (map_name, _) in enumerate(asl_maps):
    nulls = neuromaps.nulls.burt2020(globals()[asldf[{map_name}], atlas=space, density=resolution,
                                     n_perm=5000, seed=1212, parcellation=parcellation_atlas)
    np.save(path+'data/saved_nulls/'+f'{map_name}_{scale}.npy', nulls)
    print(f"Nulls for {map_name} generated and saved.")

print("Finished generating nulls")
