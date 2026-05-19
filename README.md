# Enhancing Randomized Controlled Trials Through Smartwatch-Guided Participant Matching
## Overview
Randomized controlled trials (RCTs) are the gold standard for assessing medical interventions, but they often face challenges regarding high costs and recruitment burdens. Traditional stratification usually relies on limited demographic variables like age and sex.  We developed Smartwatch-Informed Matching (SIM), a pre-randomization framework that leverages continuous physiological measurements from consumer smartwatches specifically heart rate (HR) and heart rate variability (HRV) to group similar participants and optimize trial efficiency. Using a prospective cohort of 4,795 individuals, we demonstrated that SIM can reduce required sample sizes by 9–18% while maintaining statistical power.  

## Repository Structure
The code is organized into two primary components to allow for reproduction of the paper's results:
1. **Matching Algorithm** **(matching.py)* - This script implements the core SIM framework.
2. **Power Analysis Simulation** *(simulation.py)* - This script reproduces the results found in Table 3 of the paper.

## Citation
If you use this code or the SIM framework in your research, please cite:

Shahmoon, E., Yechezkel, M., Snir, S. et al. Enhancing randomized controlled trials through smartwatch-guided participant matching for infectious disease outcomes. Sci Rep (2026). https://doi.org/10.1038/s41598-026-52579-4
