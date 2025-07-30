# PMF Calculation Script - Complete Usage Guide

## Environment Setup

### 1. Software Dependencies
```bash
# Install GROMACS (Required)
# Ubuntu/Debian:
sudo apt-get install gromacs

# CentOS/RHEL:
sudo yum install gromacs

# Or compile from source
```

### 2. Python Environment Configuration
```bash
# Create virtual environment
conda create -n pmf_env python=3.8
conda activate pmf_env

# Install Python dependencies
pip install numpy matplotlib pathlib pyyaml
```

## File Structure Preparation

### Required Input File Structure
```
your_project/
├── gaff_model/                    # Modeling task directory
│   ├── GMX_PROLIG_MD/
│   │   ├── prod/
│   │   │   └── md.gro             # Required: Final MD structure
│   │   ├── topol.top              # Required: Topology file
│   │   └── posre.itp              # Required: Position restraint file
│   └── LIG.amb2gmx/               # Required: Ligand force field directory
│       ├── LIG.itp
│       ├── LIG.gro
│       └── other_force_field_files
└── pmf.py                         # PMF script file
```

## Usage Steps

### Step 1: Script Setup and Configuration

```python
# 1. Import PMF module
import pmf

# 2. Get example configuration and customize
config = pmf.get_example_config()

# 3. Adjust configuration for your system
config.update({
    'reference_group': 'Protein',      # Reference group (usually protein)
    'moving_group': 'LIG',             # Moving group (ligand name)
    'smd': {
        'pull_rate': 0.005,            # Pull rate (nm/ps)
        'pull_k': 1000.0               # Pull force constant (kJ/mol/nm²)
    },
    'distance': {
        'start': 0.3,                  # Starting distance (nm)
        'end': 2.0                     # End distance (nm)
    }
})

# 4. Initialize PMF system
pmf.setup(config, "./GMX_PROLIG_PMF", "./gaff_model")

# 5. View workflow overview
pmf.print_workflow_summary()
```

### Step 2: Step 1 - SMD Preparation

```python
# Execute Step 1: SMD preparation
results1 = pmf.pmf_step1_smd_preparation()
print(f"SMD setup completed. Files in: {results1['smd_dir']}")
print(f"Run script: {results1['run_script']}")

# View generated files
import os
os.listdir("./GMX_PROLIG_PMF/smd/")
# Output: ['md.gro', 'topol.top', 'index.ndx', 'smd.mdp', 'run_smd.sh', 'analyze_smd.sh', 'LIG.amb2gmx/']
```

### Step 3: Manual SMD Simulation

```bash
# Change to SMD directory
cd GMX_PROLIG_PMF/smd

# Run SMD simulation script
bash run_smd.sh

# Wait for simulation completion (may take several hours)
# Upon completion, generates:
# - results/smd.gro
# - results/smd.xtc  
# - results/smd_pullf.xvg
# - results/smd_pullx.xvg
```

### Step 4: Analyze SMD Results

```python
# Analyze SMD results
smd_analysis = pmf.pmf_analyze_smd()
print(f"SMD analysis plots: {smd_analysis['plots']}")

# View force and distance curves
from matplotlib import pyplot as plt
import numpy as np

# Read and display pull force curve
pullf_data = np.loadtxt("./GMX_PROLIG_PMF/smd/results/smd_pullf.xvg", comments=['#', '@'])
plt.figure(figsize=(10, 6))
plt.plot(pullf_data[:, 0], pullf_data[:, 1])
plt.xlabel('Time (ps)')
plt.ylabel('Force (kJ/mol/nm)')
plt.title('SMD Pull Force')
plt.show()
```

### Step 5: Step 2 - Umbrella Sampling Preparation

```python
# Execute Step 2: Umbrella sampling preparation
results2 = pmf.pmf_step2_umbrella_preparation()
print(f"Generated {results2['n_windows']} umbrella windows")
print(f"Run script: {results2['run_script']}")

# View generated windows
for window in results2['windows'][:5]:  # Show first 5 windows
    print(f"Window {window['id']}: {window['distance']:.3f} nm")
```

### Step 6: Manual Umbrella Sampling

```bash
# Change to umbrella sampling directory
cd GMX_PROLIG_PMF/umbrella

# Choose execution method:

# Method 1: Parallel execution (recommended)
bash run_all_umbrella.sh parallel

# Method 2: Sequential execution (for testing)
bash run_all_umbrella.sh sequential

# Method 3: Generate cluster job scripts
bash run_all_umbrella.sh cluster
# Then submit jobs: sbatch job_window_*.sh

# Or manually run individual windows
cd window_000
bash run_window.sh
```

### Step 7: Step 3 - WHAM Analysis

```python
# Check umbrella sampling completion status
status = pmf.get_workflow_status()
print(f"Completed windows: {status['step2_details']['n_windows_completed']}")

# Execute Step 3: WHAM analysis
results3 = pmf.pmf_step3_wham_analysis()
print(f"Binding energy: {results3['binding_energy']['value']:.2f} {results3['binding_energy']['unit']}")
print(f"Analysis report: {results3['report']}")
```

### Step 8: Detailed Analysis

```python
# Generate detailed PMF analysis
detailed_analysis = pmf.pmf_analyze_pmf()
print(f"Detailed plots: {detailed_analysis['plots']}")

# View PMF curve
import matplotlib.pyplot as plt
import numpy as np

# Read PMF data
pmf_data = np.loadtxt("./GMX_PROLIG_PMF/analysis/wham/pmf.xvg", comments=['#', '@'])
distances = pmf_data[:, 0]
pmf_values = pmf_data[:, 1]

# Plot PMF curve
plt.figure(figsize=(10, 6))
plt.plot(distances, pmf_values, 'b-', linewidth=2)
plt.xlabel('Distance (nm)')
plt.ylabel('PMF (kcal/mol)')
plt.title('Potential of Mean Force')
plt.grid(True)
plt.show()

# Calculate binding energy
binding_energy = np.max(pmf_values) - np.min(pmf_values)
print(f"Binding energy: {binding_energy:.2f} kcal/mol")
```

## Monitoring and Troubleshooting

### Monitor Workflow Status
```python
# Check overall status
status = pmf.get_workflow_status()
print("Workflow Status:")
for step, completed in status.items():
    if isinstance(completed, bool):
        print(f"  {step}: {'✓' if completed else '✗'}")

# Get runtime estimates
estimates = pmf.get_expected_runtime_estimates()
print(f"Expected SMD runtime: {estimates['smd_simulation']['typical_runtime']}")
print(f"Expected umbrella runtime: {estimates['umbrella_per_window']['typical_runtime']}")
```

### Common Issues and Solutions

#### 1. Group Name Errors
```python
# If group name errors occur, check index.ndx file
with open("./GMX_PROLIG_PMF/smd/index.ndx", 'r') as f:
    content = f.read()
    print("Available groups:")
    for line in content.split('\n'):
        if line.startswith('[') and line.endswith(']'):
            print(f"  {line}")

# Then update configuration
config['reference_group'] = 'correct_protein_group_name'
config['moving_group'] = 'correct_ligand_group_name'
```

#### 2. Clean and Restart
```python
# Clean specific steps
pmf.clean_workflow(['step1'])  # Only clean Step 1
pmf.clean_workflow(['step2', 'step3'])  # Clean Steps 2 and 3

# Clean all (requires confirmation)
pmf.clean_workflow()
```

#### 3. File Verification
```python
# Verify required files exist
import os
required_files = [
    "./gaff_model/GMX_PROLIG_MD/prod/md.gro",
    "./gaff_model/GMX_PROLIG_MD/topol.top", 
    "./gaff_model/GMX_PROLIG_MD/posre.itp"
]

for file_path in required_files:
    if os.path.exists(file_path):
        print(f"Found: {file_path}")
    else:
        print(f"Missing: {file_path}")
```

## Results Interpretation

### PMF Curve Interpretation
- **Binding energy < 0**: Favorable binding, more negative = stronger binding
- **Binding energy > 0**: Unfavorable binding, possible issues
- **Minimum position**: Most stable binding distance
- **Curve smoothness**: Reflects sampling quality

### Output File Description
```
GMX_PROLIG_PMF/
├── smd/results/                   # SMD simulation results
│   ├── smd_pullf.xvg             # Pull force vs time
│   └── smd_pullx.xvg             # Distance vs time
├── umbrella/window_xxx/results/   # Each umbrella window results
│   └── umbrella_pullf.xvg        # Pull force for each window
├── analysis/
│   ├── wham/
│   │   ├── pmf.xvg               # Main PMF data
│   │   ├── pmferror.xvg          # Error estimates
│   │   └── bsprofile.xvg         # Bootstrap convergence analysis
│   ├── pmf_curve.png             # PMF curve plot
│   ├── pmf_analysis_report.txt   # Detailed analysis report
│   └── detailed_plots/           # Additional analysis plots
```

## Performance Optimization Tips

### Computational Resource Planning
```python
# Get computation time estimates
estimates = pmf.get_expected_runtime_estimates()

# Adjust based on system size
if your_system_atoms > 100000:
    print("Large system detected:")
    print("- Use GPU for SMD simulation")
    print("- Use HPC cluster for umbrella sampling")
    print("- Consider reducing production time for testing")
```

### Parallelization Strategy
```bash
# Parallel umbrella sampling (on cluster)
# 1. Generate job scripts
cd GMX_PROLIG_PMF/umbrella
bash run_all_umbrella.sh cluster

# 2. Submit jobs in batch
for script in job_window_*.sh; do
    sbatch $script
done

# 3. Monitor job status
squeue -u $USER
```

## Complete Example Script

```python
#!/usr/bin/env python3
"""
Complete PMF calculation workflow example
"""

import pmf
import os

def main():
    print("=== PMF Calculation Workflow ===")
    
    # 1. Configuration setup
    config = pmf.get_example_config()
    config.update({
        'reference_group': 'Protein',
        'moving_group': 'LIG',
        'distance': {'start': 0.3, 'end': 2.0},
        'smd': {'pull_rate': 0.005, 'pull_k': 1000.0}
    })
    
    # 2. Initialization
    pmf.setup(config, "./GMX_PROLIG_PMF", "./gaff_model")
    pmf.print_workflow_summary()
    
    # 3. Step 1: SMD preparation
    print("\n--- Step 1: SMD Preparation ---")
    try:
        results1 = pmf.pmf_step1_smd_preparation()
        print(f"SMD setup completed: {results1['run_script']}")
        print("Please run SMD simulation manually:")
        print(f"   cd {results1['smd_dir']} && bash run_smd.sh")
        
        # Wait for user confirmation of SMD completion
        input("\nPress Enter after SMD simulation completes...")
        
        # SMD analysis
        smd_analysis = pmf.pmf_analyze_smd()
        print(f"SMD analysis completed: {smd_analysis['plots']}")
        
    except Exception as e:
        print(f"Step 1 failed: {e}")
        return
    
    # 4. Step 2: Umbrella sampling preparation
    print("\n--- Step 2: Umbrella Preparation ---")
    try:
        results2 = pmf.pmf_step2_umbrella_preparation()
        print(f"Generated {results2['n_windows']} umbrella windows")
        print("Please run umbrella sampling:")
        print(f"   cd {results2['umbrella_dir']} && bash run_all_umbrella.sh parallel")
        
        # Wait for user confirmation of umbrella sampling completion
        input("\nPress Enter after umbrella sampling completes...")
        
    except Exception as e:
        print(f"Step 2 failed: {e}")
        return
    
    # 5. Step 3: WHAM analysis
    print("\n--- Step 3: WHAM Analysis ---")
    try:
        results3 = pmf.pmf_step3_wham_analysis()
        print(f"WHAM analysis completed")
        print(f"Binding energy: {results3['binding_energy']['value']:.2f} {results3['binding_energy']['unit']}")
        print(f"Report: {results3['report']}")
        
        # Detailed analysis
        detailed = pmf.pmf_analyze_pmf()
        print(f"Detailed analysis: {detailed['plots_dir']}")
        
    except Exception as e:
        print(f"Step 3 failed: {e}")
        return
    
    print("\nPMF calculation workflow completed successfully!")
    print(f"Final binding energy: {results3['binding_energy']['value']:.2f} {results3['binding_energy']['unit']}")

if __name__ == "__main__":
    main()
```

## Advanced Configuration Options

### Detailed Configuration Parameters
```python
# Complete configuration example
config = {
    'reference_group': 'Protein',           # Reference group name
    'moving_group': 'LIG',                  # Moving group name
    'smd': {
        'pull_rate': 0.005,                 # Pull rate (nm/ps)
        'pull_k': 1000.0                    # Pull force constant (kJ/mol/nm²)
    },
    'distance': {
        'start': 0.3,                       # Starting distance (nm)
        'end': 2.0                          # End distance (nm)
    },
    'simulation': {
        'dt': 0.002,                        # Time step (ps)
        'temperature': 310.0,               # Temperature (K)
        'pressure': 1.0                     # Pressure (bar)
    },
    'umbrella': {
        'production_time_ps': 22000,        # Production time per window (ps)
        'sampling_interval_ps': 10,         # Sampling interval (ps)
        'sample_interval_near': 0.1,        # Sampling interval near binding (nm)
        'sample_interval_far': 0.2,         # Sampling interval far from binding (nm)
        'cutoff_distance': 1.5              # Distance cutoff for adaptive sampling (nm)
    },
    'analysis': {
        'begin_time_ps': 2000,              # Analysis start time (ps)
        'bootstrap_iterations': 50,         # Bootstrap iterations for error estimation
        'energy_unit': 'kCal'               # Energy unit (kCal or kJ)
    }
}
```

### System-Specific Adjustments
```python
# For large systems (>100k atoms)
config['umbrella']['production_time_ps'] = 10000  # Shorter production time
config['analysis']['bootstrap_iterations'] = 20   # Fewer bootstrap iterations

# For small systems (<10k atoms)
config['umbrella']['production_time_ps'] = 50000  # Longer production time
config['analysis']['bootstrap_iterations'] = 100  # More bootstrap iterations

# For weak binding systems
config['distance']['end'] = 3.0                   # Extend to larger distances
config['umbrella']['sample_interval_far'] = 0.3   # Coarser sampling at large distances
```

## Best Practices

### 1. Preparation Phase
- Validate all input files before starting
- Test with shorter simulations first
- Ensure adequate disk space (10-100 GB)
- Check group names in index file

### 2. Execution Phase
- Monitor SMD for reasonable force profiles
- Use parallel execution for umbrella sampling
- Check individual window convergence
- Save intermediate results regularly

### 3. Analysis Phase
- Examine PMF curve for smoothness
- Check bootstrap error estimates
- Validate results against experimental data
- Consider multiple independent runs

### 4. Troubleshooting
- Check log files for errors
- Verify GROMACS installation
- Monitor system resources
- Use workflow status functions

## Key Points Summary

### Critical Preparation

1. **Environment Requirements**:
   - GROMACS installed and in PATH
   - Python 3.x + numpy, matplotlib
   - Sufficient disk space (GB-level)

2. **Required Files**:
   ```
   gaff_model/
   ├── GMX_PROLIG_MD/prod/md.gro     # Final MD structure
   ├── GMX_PROLIG_MD/topol.top       # Topology file  
   ├── GMX_PROLIG_MD/posre.itp       # Position restraints
   └── LIG.amb2gmx/                  # Ligand force field
   ```

### Three-Step Workflow

1. **Step 1 - SMD Preparation**:
   ```python
   pmf.pmf_step1_smd_preparation()
   # Then manually run: bash run_smd.sh
   ```

2. **Step 2 - Umbrella Sampling**:
   ```python
   pmf.pmf_step2_umbrella_preparation() 
   # Then manually run: bash run_all_umbrella.sh parallel
   ```

3. **Step 3 - WHAM Analysis**:
   ```python
   pmf.pmf_step3_wham_analysis()
   # Automatically generates PMF curve and binding energy
   ```

### Important Notes

1. **Group Configuration**: Ensure `reference_group` and `moving_group` exist in your system
2. **Computation Time**: SMD takes 1-6 hours, umbrella sampling 2-12 hours per window
3. **Parallelization**: Strongly recommended for umbrella sampling windows
4. **Disk Space**: Ensure sufficient space for trajectory files

### Quick Start

```python
import pmf

# 1. Configuration
config = pmf.get_example_config()
config['reference_group'] = 'Protein'  # Adjust for your system
config['moving_group'] = 'LIG'         # Adjust for your ligand

# 2. Initialize
pmf.setup(config, "./GMX_PROLIG_PMF", "./gaff_model")

# 3. View overview
pmf.print_workflow_summary()

# 4. Start first step
results = pmf.pmf_step1_smd_preparation()
```

This script offers **high automation**, **comprehensive error handling**, and **detailed result analysis**, making it particularly suitable for research projects requiring batch PMF calculations. Following this guide, you should be able to successfully calculate protein-ligand binding free energies!