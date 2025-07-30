#!/usr/bin/env python3
"""
Simple PMF Complete Workflow Runner
Author: PRISM TEAM
"""

import pmf
import os

def main():
    print("PMF Complete Workflow Runner")
    print("=" * 40)
    
    # Configuration
    config = {
        'reference_group': 'Protein',
        'moving_group': 'LIG',
        'smd': {
            'pull_rate': 0.005,
            'pull_k': 1000.0
        },
        'distance': {
            'start': 0.3,
            'end': 2.0
        },
        'umbrella': {
            'production_time_ps': 22000,
            'sampling_interval_ps': 10
        },
        'analysis': {
            'begin_time_ps': 2000,
            'bootstrap_iterations': 50
        }
    }
    
    # Initialize PMF system
    pmf.setup(config, "./GMX_PROLIG_PMF", "./gaff_model")
    
    # Step 1: SMD Preparation
    print("\nStep 1: SMD Preparation")
    results1 = pmf.pmf_step1_smd_preparation()
    print(f"SMD setup completed: {results1['smd_dir']}")
    print(f"Run command: cd {results1['smd_dir']} && bash run_smd.sh")
    
    input("Press Enter after SMD simulation completes...")
    
    # Step 2: Umbrella Preparation
    print("\nStep 2: Umbrella Sampling Preparation")
    results2 = pmf.pmf_step2_umbrella_preparation()
    print(f"Generated {results2['n_windows']} umbrella windows")
    print(f"Run command: cd {results2['umbrella_dir']} && bash run_all_umbrella.sh parallel")
    
    input("Press Enter after umbrella sampling completes...")
    
    # Step 3: WHAM Analysis
    print("\nStep 3: WHAM Analysis")
    results3 = pmf.pmf_step3_wham_analysis()
    
    # Generate detailed analysis
    detailed = pmf.pmf_analyze_pmf()
    
    # Final Results
    print("\nResults:")
    print(f"Binding Energy: {results3['binding_energy']['value']:.2f} {results3['binding_energy']['unit']}")
    print(f"PMF Plot: {results3['plots']['pmf_plot']}")
    print(f"Report: {results3['report']}")
    print(f"Detailed Analysis: {detailed['plots_dir']}")
    print("\nWorkflow completed successfully!")

if __name__ == "__main__":
    main()