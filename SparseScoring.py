# This code takes input data directory and output file as arguments
# It assumes input data directory has a set of files, once for each type of cell
# that contains ChiDist for fusion candidates, so these files have one line per fusion candidate
# These files are expected to be tab separated, containing at least two columns named FusionName and Zscore
# After reading these files, for each fusion candidate, zscores are combined for each fusion candidate
# Finally, it writes the results which is combined zscores for each fusion candidate into output file
import sys
import os
import subprocess

import pandas as pd
import numpy as np
from scipy import stats

# take a list of scores, numpy array
# calculate sparsity score, L1
def scoring_by_abs_mean(scores):
    if len(scores) == 0:
        return 0
    return np.mean(np.abs(scores))

def scoring_L2(scores):
    if len(scores) == 0:
        return 0
    return np.sqrt(np.sum(np.square(scores)))/len(scores)


# data_dir = "/Users/sadiaatique/PycharmProjects/4761Project/results/PValueData"
# output_file = "/Users/sadiaatique/PycharmProjects/4761Project/results/PValueData/sparsity_result.txt"

data_dir = sys.argv[1]
output_file = sys.argv[2]

all_files = [f for f in os.listdir(data_dir) if f.endswith('.txt') or f.endswith('.csv')]

columns_to_extract = ['FusionName', 'Zscore']

all_dfs = []
for file in all_files:
    filepath = os.path.join(data_dir, file)
    try:
        df = pd.read_csv(filepath, sep='\t')  # change sep if needed
        selected = df[columns_to_extract]
        all_dfs.append(selected)
    except Exception as e:
        print(f"Skipping {file}: {e}")

combined_df = pd.concat(all_dfs, ignore_index=True)

grouped_scores_by_fusion = grouped_df = combined_df.groupby('FusionName', as_index=False).agg({'Zscore': list})

print(grouped_scores_by_fusion)
fusion_candidate_count = grouped_scores_by_fusion.shape[0]

abs_mean_scores = np.zeros(fusion_candidate_count)
l2_scores = np.zeros(fusion_candidate_count)
l2_over_l1 = np.zeros(fusion_candidate_count)
for index, row in grouped_scores_by_fusion.iterrows():
    fusion_name = row['FusionName']
    z = np.array(row['Zscore'])
    print(fusion_name)
    print(z)
    abs_mean_scores[index] = scoring_by_abs_mean(z)
    l2_scores[index] = scoring_L2(z)
    l2_over_l1[index] = l2_scores[index]/abs_mean_scores[index]


grouped_scores_by_fusion['L1_score'] = abs_mean_scores
grouped_scores_by_fusion['L2_score'] = l2_scores
grouped_scores_by_fusion['L2_over_L1'] = l2_over_l1
print(grouped_scores_by_fusion)
grouped_scores_by_fusion.to_csv(output_file, sep='\t', index=False)

