import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import zscore

sexyRdBu = matplotlib.colors.LinearSegmentedColormap.from_list("", [
    "#ffffff",
    "#FFBFB1",
    "#FF937B",
    "#FF7557",
    "#FF5631",
    "#FF3A0F",
    "#DF2800"
])

path = '/Users/k1801291/local_data/github/cbfreceptors/'
scale = '122'

X = pd.read_csv(path+'data/receptors_parcellated_'+scale+'.csv', header=None)

data = pd.read_csv(path+'data/asl_parcellated_'+scale+'.csv', header=None)
data = zscore(data)

# Load the sample names from file "atlas.csv"
info = pd.read_csv(path+'data/atlas/schaefer2018_xiao_2019.csv')

# Initialize a DataFrame to store Cook's distance for each dependent variable
cooks_distance_df = pd.DataFrame()

for col in range(0, 9):
    y = data[col]
    
    # Filter rows where scale is "scale122" and select the columns for label, and hemisphere
    sample_names = info.query('scale == "scale122"')[['label', 'hemisphere']]

    # Fit the linear regression model
    model = sm.OLS(y, X).fit()

    # Calculate Cook's distance
    influence = model.get_influence()
    cooks_distance = influence.cooks_distance[0]

    # Save Cook's distance as a column in the DataFrame
    cooks_distance_df[path+f'results/Cooks_Distance_{col + 1}'] = cooks_distance

    # Sort the samples by Cook's distance and extract the sample indices in order
    sorted_indices = np.argsort(cooks_distance)[::-1]

    # Create labels for each sample using its name, yeo_7, and hemisphere
    sample_labels = [f"{row.label}, {row.hemisphere})" for _, row in sample_names.iterrows()]

    # Print the top 20 most important samples
    print(f"Top 20 most important samples (column {col}):")
    for i in sorted_indices[:19]:
        print(f"{sample_labels[i]}: {cooks_distance[i]:.4f}")

    # Save the top 20 most important samples to a CSV file
    top_samples = [sample_labels[i] for i in sorted_indices[:19]]
    top_cooks = [cooks_distance[i] for i in sorted_indices[:19]]
    top_df = pd.DataFrame({'Sample': top_samples, "Cook's Distance": top_cooks})
    top_df.to_csv(path+f'results/Cooks_Distance/top_samples_col{col}.csv', index=False)

# Save the DataFrame to a CSV file
cooks_distance_df.to_csv(path+'results/cooks_distance/cooks_distance_'+scale+'.csv', index=False, header=False)
