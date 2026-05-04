import pandas as pd
import numpy as np
import scipy.stats as st
from collections import namedtuple
from tqdm import tqdm
from statsmodels.stats.power import tt_solve_power

# --- 1. CONFIGURATION & DATA STRUCTURES ---
# Define a container for clinical trial parameters
# pi: symptomatic infection probability
# ve: vaccine efficacy
# r: observed correlation between matched pairs
Params = namedtuple('Params', ['pi', 've', 'r'])

# --- 2. CORE SIMULATION LOGIC ---

def run_simulation(n, pi, r, ve, b=0):
    """
    Generates simulated clinical trial outcomes for matched pairs.
    Following Appendix D, this incorporates dependencies based on observed 
    correlations in the PerMed dataset.
    """
    # Create n/2 matched pairs
    half_n = int(n / 2)
    
    # Sample latent probabilities
    unbiased = np.random.uniform(0, 1, half_n)
    biased = np.round(np.random.uniform(0, 1, half_n) + b, 7)
    
    # Split participants into control and intervention (vaccine) arms
    ub_1, ub_2 = np.split(unbiased, 2)
    b_1, b_2 = np.split(biased, 2)
    
    control_grp = np.concatenate([ub_1, b_1])
    case_grp = np.concatenate([ub_2, b_2])
    
    # Draw initial symptomatic infection status
    control_symptoms = (control_grp <= pi).astype(int)
    case_symptoms = (case_grp <= pi).astype(int)
    
    # Apply correlation (r) to reflect physiological similarity between pairs
    # This simulates the 'matching' effect where pairs share similar risk profiles
    correlated_indices = np.random.choice(half_n, int(r * half_n), replace=False)
    case_symptoms[correlated_indices] = control_symptoms[correlated_indices]
    
    # Apply vaccine efficacy (VE) to the intervention (case) group
    # Only symptomatic individuals in the case group are eligible for reduction
    symptomatic_mask = (case_symptoms == 1)
    vaccination_impact = np.random.binomial(1, p=ve, size=np.sum(symptomatic_mask))
    case_symptoms[symptomatic_mask] -= vaccination_impact
    
    return pd.DataFrame({
        "control_symptoms": control_symptoms,
        "case_symptoms": case_symptoms
    })

# --- 3. POWER ANALYSIS HELPERS ---

def get_ci(lst):
    """Calculates 95% Confidence Interval for estimated power."""
    if len(lst) < 2:
        return 0, 1
    return st.t.interval(confidence=0.95, df=len(lst)-1,
                         loc=np.mean(lst),
                         scale=st.sem(lst))

def get_power(j, simulation_results):
    """
    Computes statistical power for sample size 'j' using a paired t-test approach 
    derived from simulation outcomes.
    """
    j = int(j)
    # Global standard deviation of differences
    stds = (simulation_results.control_symptoms - simulation_results.case_symptoms).std()
    
    # Draw a sample of size j and compute effect size (Cohen's d)
    cur_sample = simulation_results.sample(n=j)
    sample_diff = cur_sample.control_symptoms - cur_sample.case_symptoms
    cohen_d = sample_diff.mean() / stds
    
    # Solve for power given alpha=0.05
    return tt_solve_power(effect_size=cohen_d, alpha=0.05, nobs=j, alternative='two-sided')


def optimize_sample_size(target_power, sim_results, ci_tol=0.03, start_n=1000):
    """
    Iteratively finds the minimum sample size required to maintain target power.
    """
    iters = 0
    current_x = start_n
    step_size = max(int(current_x * 0.1), 1)
    
    memo = {} # Cache results for n
    
    while iters < 250:
        iters += 1
        
        if current_x not in memo:
            memo[current_x] = []
            
        # Perform iterations to build a confidence interval for power at current n
        while len(memo[current_x]) < 25:
            memo[current_x].append(get_power(current_x, sim_results))
            
        est_power = np.mean(memo[current_x])
        lower_ci, upper_ci = get_ci(memo[current_x])
        ci_size = upper_ci - lower_ci
        
        # Check for convergence: Target power is within the CI and CI is narrow enough
        if target_power >= lower_ci and target_power <= upper_ci and ci_size <= ci_tol:
            return int(current_x), ci_size
            
        # Adjust search based on target
        if target_power > upper_ci:
            # Power too low, increase n
            step_size = int(step_size * 1.5)
            current_x += step_size
        elif target_power < lower_ci:
            # Power too high, decrease n
            step_size = int(step_size * 1.5)
            current_x = max(step_size // 2, current_x - step_size)
        else:
            # CI is too wide or we are close; refine search
            step_size = max(0.25 * step_size, 1)
            if est_power > target_power:
                current_x -= step_size
            else:
                current_x += step_size
                
    return int(current_x), ci_size

# --- 5. MAIN EXECUTION BLOCK ---

if __name__ == "__main__":
    n_iterations = 100
    pi_values = [0.003, 0.030, 0.100]
    ve_values = [0.40, 0.50, 0.60, 0.90]
    target_powers = [0.80, 0.95]
    correlations = {"Standard": 0.012, "SIM_Symptom": 0.176, "SIM_Physio": 0.245}
    
    raw_results = []

    print(f"Starting simulation: {n_iterations} iterations per cell...")
    
    # Outer loop for iterations
    for run_idx in range(n_iterations):
        print(f"Iteration {run_idx + 1}/{n_iterations}")
        
        for pi in pi_values:
            for ve in ve_values:
                for label, r in correlations.items():
                    # Generate fresh population data for this run
                    pop_data = run_simulation(n=1_000_000, pi=pi, r=r, ve=ve)
                    
                    for pwr in target_powers:
                        n_req = optimize_sample_size(pwr, pop_data, ci_tol=0.03, start_n=1000)
                        raw_results.append({
                            "Iteration": run_idx,
                            "pi": pi,
                            "VE": ve,
                            "Matching": label,
                            "Target_Power": pwr,
                            "N_Required": n_req
                        })

    # --- 4. DATA AGGREGATION & SAVING ---
    all_runs_df = pd.DataFrame(raw_results)
    
    # Save the raw data for all 100 runs in case you need to calculate SD later
    all_runs_df.to_csv("simulated_sample_sizes_RAW.csv", index=False)
    
    # Calculate the mean value for each cell
    mean_results_df = all_runs_df.groupby(['pi', 'VE', 'Matching', 'Target_Power'])['N_Required'].mean().reset_index()
    
    # Save the final mean values
    mean_results_df.to_csv("simulated_sample_sizes_MEAN.csv", index=False)
    
    print("Simulation complete.")
    print("Raw results saved to simulated_sample_sizes_RAW.csv")
    print("Mean values saved to simulated_sample_sizes_MEAN.csv")