import pandas as pd
import numpy as np
import random

def perform_sim_matching(hr_emd_mat, hrv_emd_mat, demographics, noise_factor=0.15):
    """
    Implements the SIM (Smartwatch-Informed Matching) procedure.
    
    Parameters:
    - hr_emd_mat: Symmetric matrix of EMD distances for heart rate distributions.
    - hrv_emd_mat: Symmetric matrix of EMD distances for HRV distributions.
    - demographics: DataFrame containing 'participant_num', 'Male', and 'age_group'.
    - noise_factor: Maximum scale of perturbation (default 15% of empirical std).
    """
    
    # 1. Combine HR and HRV distances using Euclidean distance
    # Dissimilarity = sqrt( d_HR^2 + d_HRV^2 )
    common_ids = hr_emd_mat.index.intersection(hrv_emd_mat.index)
    hr_sub = hr_emd_mat.loc[common_ids, common_ids]
    hrv_sub = hrv_emd_mat.loc[common_ids, common_ids]
    
    dissimilarity_mat = np.sqrt(hr_sub**2 + hrv_sub**2)
    
    # 2. Stochastic Component: Perturb the matrix prior to matching
    # Calculate the empirical standard deviation of the dissimilarity scores
    # We take only the upper triangle to avoid the diagonal and double-counting
    empirical_std = dissimilarity_mat.values[np.triu_indices_from(dissimilarity_mat, k=1)].std()
    
    # Draw noise from a distribution scaled to 15% of the empirical standard deviation
    noise_scale = noise_factor * empirical_std
    noise = np.random.normal(0, noise_scale, size=dissimilarity_mat.shape)
    
    # Ensure the noise matrix is symmetric to maintain consistent pairwise distance
    noise = (noise + noise.T) / 2
    perturbed_mat = dissimilarity_mat + noise
    
    # Set diagonal to infinity to prevent matching a participant with themselves
    np.fill_diagonal(perturbed_mat.values, np.inf)
    
    # 3. Stratification and Greedy Matching
    # Groups participants by sex (Male) and age_group (e.g., 18-55, 56-85)
    matched_results = []
    
    # Prepare demographics for iteration
    demo_relevant = demographics[demographics['participant_num'].isin(common_ids)]
    strata = demo_relevant.groupby(['gender', 'age_group'])
    
    for (is_male, age_grp), group in strata:
        ids = group['participant_num'].tolist()
        if len(ids) < 2:
            continue # Cannot match a single participant
            
        # Create a local copy for greedy matching within this subgroup
        current_pool = set(ids)
        sub_perturbed = perturbed_mat.loc[ids, ids].copy()
        
        while len(current_pool) > 1:
            # Nearest-neighbor greedy matching: Find the minimum perturbed distance
            temp_mat = sub_perturbed.loc[list(current_pool), list(current_pool)]
            p1, p2 = temp_mat.stack().idxmin()
            
            # Constrained Randomization: Randomly assign to Intervention or Control
            if random.random() > 0.5:
                pair = {'intervention': p1, 'control': p2}
            else:
                pair = {'intervention': p2, 'control': p1}
            
            pair.update({'strata_male': is_male, 'strata_age': age_grp})
            matched_results.append(pair)
            
            # Remove matched individuals from the pool
            current_pool.remove(p1)
            current_pool.remove(p2)
            
    return pd.DataFrame(matched_results)

# --- Usage Example ---
hr_matrix = pd.read_pickle("data/hr_dist_mat.pkl")
hrv_matrix = pd.read_pickle("data/hrv_dist_mat.pkl")
demographics_df = pd.read_pickle("data/demographics_df.pkl")
matched_df = perform_sim_matching(hr_matrix, hrv_matrix, demographics_df)
print(matched_df)