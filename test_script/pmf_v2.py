#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Optimized PMF (Potential of Mean Force) module for PRISM

This module provides streamlined PMF calculation capabilities including:
- Steered Molecular Dynamics (SMD) preparation and execution
- Umbrella Sampling window generation and execution
- WHAM analysis workflow
- Automated result visualization

Author: PRISM Team (Optimized)
Version: 2.0
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import yaml
import subprocess
import shutil
from typing import Dict, List, Tuple, Optional, Union
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PMFCalculator:
    """
    Streamlined PMF (Potential of Mean Force) calculator for protein-ligand systems
    
    Combines SMD simulation, umbrella sampling, and WHAM analysis in a unified workflow
    """
    
    def __init__(self, config: Dict, output_dir: str, md_results_dir: str):
        """
        Initialize PMF calculator
        
        Parameters:
        -----------
        config : Dict
            PMF configuration dictionary
        output_dir : str
            Output directory for PMF results
        md_results_dir : str
            Directory containing MD simulation results
        """
        self.config = config
        self.output_dir = Path(output_dir).resolve()
        self.md_results_dir = Path(md_results_dir).resolve()
        
        # Extract all configuration in one go
        self._extract_all_config()
        
        # Setup directory structure
        self._setup_directories()
        
        logger.info("PMF Calculator initialized successfully")
        self._print_summary()
    
    def _extract_all_config(self):
        """Extract and validate all configuration parameters at once"""
        # Basic settings
        self.method = self.config.get('method', 'umbrella')
        self.reference_group = self.config.get('reference_group', 'Protein')
        self.moving_group = self.config.get('moving_group', 'LIG')
        
        # SMD parameters
        smd = self.config.get('smd', {})
        self.smd_pull_rate = smd.get('pull_rate', 0.005)
        self.smd_pull_k = smd.get('pull_k', 1000.0)
        self.smd_pull_direction = smd.get('pull_direction', [-1, 0, 0])
        
        # Distance parameters
        dist = self.config.get('distance', {})
        self.distance_start = dist.get('start', 0.3)
        self.distance_end = dist.get('end', 2.0)
        self.force_constant = dist.get('force_constant', 1000.0)
        
        # Simulation parameters
        sim = self.config.get('simulation', {})
        self.dt = sim.get('dt', 0.002)
        self.temperature = sim.get('temperature', 310.0)
        self.pressure = sim.get('pressure', 1.0)
        
        # Umbrella sampling
        umbrella = self.config.get('umbrella', {})
        self.n_windows = umbrella.get('n_windows', 20)
        self.production_time = umbrella.get('production_time_ps', 22000)  # 22ns
        self.sampling_interval = umbrella.get('sampling_interval_ps', 10)
        
        # Analysis parameters
        analysis = self.config.get('analysis', {})
        self.analysis_method = analysis.get('method', 'wham')
        self.wham_begin_time = analysis.get('begin_time_ps', 2000)
        self.bootstrap_iterations = analysis.get('bootstrap_iterations', 50)
    
    def _setup_directories(self):
        """Create necessary directory structure"""
        self.pmf_dir = self.output_dir / "pmf"
        self.smd_dir = self.pmf_dir / "smd"
        self.umbrella_dir = self.pmf_dir / "umbrella"
        self.analysis_dir = self.pmf_dir / "analysis"
        self.mdp_dir = self.output_dir / "mdps"
        
        for directory in [self.pmf_dir, self.smd_dir, self.umbrella_dir, self.analysis_dir, self.mdp_dir]:
            directory.mkdir(parents=True, exist_ok=True)
    
    def _print_summary(self):
        """Print configuration summary"""
        logger.info("=== PMF Configuration Summary ===")
        logger.info(f"Method: {self.method}")
        logger.info(f"Groups: {self.reference_group} → {self.moving_group}")
        logger.info(f"Distance range: {self.distance_start} - {self.distance_end} nm")
        logger.info(f"SMD pull rate: {self.smd_pull_rate} nm/ps")
        logger.info(f"Windows: {self.n_windows}")
        logger.info(f"Production time: {self.production_time/1000:.1f} ns")
    
    def generate_mdp_files(self):
        """Generate all required MDP files"""
        logger.info("Generating MDP files...")
        
        # SMD MDP
        smd_mdp = self._create_smd_mdp()
        with open(self.mdp_dir / "smd.mdp", 'w') as f:
            f.write(smd_mdp)
        
        # Umbrella MDP template
        umbrella_mdp = self._create_umbrella_mdp()
        with open(self.mdp_dir / "umbrella.mdp", 'w') as f:
            f.write(umbrella_mdp)
        
        logger.info("✓ MDP files generated")
        return str(self.mdp_dir / "smd.mdp"), str(self.mdp_dir / "umbrella.mdp")
    
    def _create_smd_mdp(self):
        """Create SMD MDP content"""
        distance_range = self.distance_end - self.distance_start
        smd_time_ps = distance_range / self.smd_pull_rate
        nsteps = int(smd_time_ps / self.dt)
        output_interval = int(5.0 / self.dt)  # Every 5 ps
        
        return f"""; SMD (Steered Molecular Dynamics) Parameters
title               = SMD simulation ({distance_range:.2f} nm in {smd_time_ps:.1f} ps)
integrator          = md
nsteps              = {nsteps}
dt                  = {self.dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}

; Bonds and constraints
continuation        = yes
constraint_algorithm = lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4

; Neighbor searching and short-range nonbonded interactions
cutoff-scheme       = Verlet
nstlist             = 10
rlist               = 1.0
rcoulomb            = 1.0
rvdw                = 1.0

; Electrostatics
coulombtype         = PME
pme_order           = 4
fourierspacing      = 0.16

; Temperature coupling
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1     0.1
ref_t               = {self.temperature}     {self.temperature}

; Pressure coupling
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 1.0
ref_p               = {self.pressure}
compressibility     = 4.5e-5

; Periodic boundary conditions
pbc                 = xyz
DispCorr            = EnerPres

; Velocity generation
gen_vel             = no

; COM motion removal
refcoord_scaling    = com
comm-mode           = Linear
comm-grps           = Protein Non-Protein

; Pull code for SMD
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {self.reference_group}
pull_group2_name    = {self.moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y N N
pull-coord1-vec     = {' '.join(map(str, self.smd_pull_direction))}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = {self.smd_pull_rate}
pull_coord1_k       = {self.smd_pull_k}
pull-pbc-ref-prev-step-com = yes
"""
    
    def _create_umbrella_mdp(self):
        """Create umbrella sampling MDP template"""
        nsteps = int(self.production_time / self.dt)
        output_interval = int(self.sampling_interval / self.dt)
        
        return f"""; Umbrella Sampling Parameters
title               = Umbrella sampling (22.0 ns)
integrator          = md
nsteps              = {nsteps}
dt                  = {self.dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}

; Bonds and constraints
continuation        = yes
constraint_algorithm = lincs
constraints         = h-bonds
lincs_iter          = 1
lincs_order         = 4

; Neighbor searching and short-range nonbonded interactions
cutoff-scheme       = Verlet
nstlist             = 10
rlist               = 1.0
rcoulomb            = 1.0
rvdw                = 1.0

; Electrostatics
coulombtype         = PME
pme_order           = 4
fourierspacing      = 0.16

; Temperature coupling
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1     0.1
ref_t               = {self.temperature}     {self.temperature}

; Pressure coupling
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 1.0
ref_p               = {self.pressure}
compressibility     = 4.5e-5

; Periodic boundary conditions
pbc                 = xyz
DispCorr            = EnerPres

; Velocity generation
gen_vel             = no

; COM motion removal
refcoord_scaling    = com
comm-mode           = Linear
comm-grps           = Protein Non-Protein

; Pull code for umbrella sampling (static constraint)
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {self.reference_group}
pull_group2_name    = {self.moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y N N
pull-coord1-vec     = {' '.join(map(str, self.smd_pull_direction))}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = 0.0
pull_coord1_k       = {self.smd_pull_k}
pull_coord1_r0      = DISTANCE_PLACEHOLDER
pull-pbc-ref-prev-step-com = yes
"""
    
    def generate_complete_workflow(self):
        """Generate complete PMF workflow script"""
        logger.info("Generating complete PMF workflow...")
        
        distance_range = self.distance_end - self.distance_start
        smd_time_ps = distance_range / self.smd_pull_rate
        n_frames = 300
        frame_interval = int(smd_time_ps / n_frames)
        
        script_content = f"""#!/bin/bash
######################################################
# COMPLETE PMF WORKFLOW SCRIPT
######################################################

set -e  # Exit on any error

echo "=== Starting Complete PMF Workflow ==="

# Configuration
REFERENCE_GROUP="{self.reference_group}"
MOVING_GROUP="{self.moving_group}"
PULL_RATE="{self.smd_pull_rate}"
PULL_K="{self.smd_pull_k}"
N_FRAMES="{n_frames}"
FRAME_INTERVAL="{frame_interval}"

# Step 1: Setup and validation
echo "Step 1: Setting up PMF environment..."

# Create necessary directories
mkdir -p smd trajectory_frames distance_analysis umbrella analysis

# Create index file
if [ ! -f index.ndx ]; then
    echo "Creating index file..."
    echo "r $MOVING_GROUP\\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx
    echo "✓ Index file created"
fi

# Step 2: SMD Simulation
echo "Step 2: Running SMD simulation..."
if [ ! -f smd/smd.gro ]; then
    echo "Starting SMD simulation..."
    echo "  Pull rate: $PULL_RATE nm/ps"
    echo "  Force constant: $PULL_K kJ/mol/nm²" 
    echo "  Simulation time: {smd_time_ps:.1f} ps"
    
    gmx grompp -f ../mdps/smd.mdp -c solv_ions.gro -n index.ndx -p topol.top -o smd/smd.tpr -maxwarn 10
    gmx mdrun -s smd/smd.tpr -deffnm smd/smd -ntmpi 1 -ntomp 10 -gpu_id 0 -v
    echo "✓ SMD simulation completed"
else
    echo "✓ SMD already completed"
fi

# Step 3: Analyze SMD trajectory
echo "Step 3: Analyzing SMD trajectory..."

# Extract frames
if [ ! -f trajectory_frames/frame_0.gro ]; then
    echo "Extracting $N_FRAMES trajectory frames..."
    echo "0" | gmx trjconv -f smd/smd.xtc -s smd/smd.tpr -o trajectory_frames/frame.gro -sep -dt $FRAME_INTERVAL
    echo "✓ Trajectory frames extracted"
fi

# Calculate distances
if [ ! -f distance_analysis/summary_distances.dat ]; then
    echo "Calculating distances for each frame..."
    echo "# Frame Distance(nm)" > distance_analysis/summary_distances.dat
    
    for ((i=0; i<$N_FRAMES; i++)); do
        if [ -f trajectory_frames/frame_${{i}}.gro ]; then
            echo "Processing frame ${{i}}..."
            gmx distance -s smd/smd.tpr -f trajectory_frames/frame_${{i}}.gro -n index.ndx \\
                -select "com of group $REFERENCE_GROUP plus com of group $MOVING_GROUP" \\
                -oall distance_analysis/dist_${{i}}.xvg 2>/dev/null
            
            d=$(tail -n 1 distance_analysis/dist_${{i}}.xvg | awk '{{print $2}}')
            echo "${{i}} ${{d}}" >> distance_analysis/summary_distances.dat
            rm distance_analysis/dist_${{i}}.xvg
        fi
    done
    echo "✓ Distance analysis completed"
fi

# Step 4: Generate umbrella windows
echo "Step 4: Generating umbrella sampling windows..."
python3 << 'EOF'
import numpy as np
import os

# Read distances
data = np.loadtxt('distance_analysis/summary_distances.dat')
frames = data[:, 0].astype(int)
distances = data[:, 1]

# Adaptive sampling: dense near protein, sparse far away
sampled_indices = []
current_idx = 0
sampled_indices.append(current_idx)

while current_idx < len(distances):
    current_dist = distances[current_idx]
    
    # Adaptive interval based on distance
    if current_dist < 1.5:
        target_interval = 0.1  # 0.1 nm for close distances
    else:
        target_interval = 0.2  # 0.2 nm for far distances
    
    target_dist = current_dist + target_interval
    
    # Find closest frame to target distance
    remaining_dists = distances[current_idx:]
    diffs = np.abs(remaining_dists - target_dist)
    next_relative_idx = np.argmin(diffs)
    next_idx = current_idx + next_relative_idx
    
    if next_idx == current_idx:
        break
    
    sampled_indices.append(next_idx)
    current_idx = next_idx

# Save umbrella windows
with open('umbrella/windows_list.dat', 'w') as f:
    f.write("# Window Frame Distance(nm)\\n")
    for i, idx in enumerate(sampled_indices):
        f.write(f"{{i:03d}} {{frames[idx]:d}} {{distances[idx]:.3f}}\\n")

print(f"Generated {{len(sampled_indices)}} umbrella windows")
EOF

# Step 5: Setup umbrella sampling windows
echo "Step 5: Setting up umbrella sampling windows..."
if [ -f umbrella/windows_list.dat ]; then
    while read -r window frame distance; do
        if [[ $window == \#* ]]; then continue; fi
        
        window_dir="umbrella/window_$window"
        mkdir -p "$window_dir"
        
        # Create window-specific MDP file
        sed "s/DISTANCE_PLACEHOLDER/$distance/g" ../mdps/umbrella.mdp > "$window_dir/umbrella.mdp"
        sed -i "s/Umbrella sampling (22.0 ns)/Umbrella sampling at $distance nm (frame $frame)/g" "$window_dir/umbrella.mdp"
        
        # Copy starting structure from SMD frame
        cp "trajectory_frames/frame_$frame.gro" "$window_dir/start.gro"
        
        # Create run script for this window
        cat > "$window_dir/run_umbrella.sh" << WINDOW_EOF
#!/bin/bash
set -e
echo "Running umbrella window $window at distance $distance nm..."

gmx grompp -f umbrella.mdp -c start.gro -n ../../index.ndx -p ../../topol.top -o umbrella.tpr -maxwarn 10
gmx mdrun -s umbrella.tpr -deffnm umbrella -ntmpi 1 -ntomp 10 -gpu_id 0 -v

echo "✓ Window $window completed"
WINDOW_EOF
        chmod +x "$window_dir/run_umbrella.sh"
        
        echo "✓ Setup window $window (frame $frame, distance $distance nm)"
    done < umbrella/windows_list.dat
fi

# Step 6: Run all umbrella sampling (optional - can be run separately)
echo "Step 6: Umbrella sampling windows are ready"
echo "To run all windows:"
echo "  cd umbrella && find . -name 'run_umbrella.sh' -execdir bash run_umbrella.sh \\;"

# Step 7: Generate WHAM analysis script
echo "Step 7: Generating WHAM analysis script..."
cat > analysis/run_wham.sh << 'WHAM_EOF'
#!/bin/bash
set -e

echo "=== Starting WHAM Analysis ==="

cd {self.analysis_dir}
mkdir -p wham
cd wham

# Generate file lists for WHAM
echo "Collecting umbrella sampling results..."

# Find all completed umbrella windows
find ../../umbrella -name "umbrella.tpr" | sort > tpr-files.dat
find ../../umbrella -name "umbrella_pullf.xvg" | sort > pullf-files.dat

n_windows=$(wc -l < tpr-files.dat)
echo "Found $n_windows completed umbrella windows"

if [ $n_windows -lt 5 ]; then
    echo "Error: Need at least 5 windows for WHAM analysis"
    exit 1
fi

# Run WHAM analysis
echo "Running gmx wham..."
gmx wham \\
    -it tpr-files.dat \\
    -if pullf-files.dat \\
    -b {self.wham_begin_time} \\
    -o pmf.xvg \\
    -unit kCal \\
    -bsres pmferror.xvg \\
    -bsprof zerrorprofile.xvg \\
    -nBootstrap {self.bootstrap_iterations} \\
    -bs-method b-hist

echo "✓ WHAM analysis completed!"
echo "Results:"
echo "  - PMF curve: pmf.xvg"  
echo "  - Bootstrap error: pmferror.xvg"
echo "  - Error profile: zerrorprofile.xvg"
WHAM_EOF

chmod +x analysis/run_wham.sh

echo "=== PMF Workflow Setup Complete ==="
echo ""
echo "Generated files:"
echo "  - SMD setup: smd/"
echo "  - Trajectory frames: trajectory_frames/"
echo "  - Distance analysis: distance_analysis/"
echo "  - Umbrella windows: umbrella/"
echo "  - WHAM analysis: analysis/"
echo ""
echo "Next steps:"
echo "1. Run SMD: bash this script runs SMD automatically"
echo "2. Run umbrella sampling: cd umbrella && find . -name 'run_umbrella.sh' -execdir bash run_umbrella.sh \\;"
echo "3. Run WHAM analysis: cd analysis && bash run_wham.sh"
echo ""
echo "✓ Complete PMF workflow ready!"
"""
        
        script_path = self.pmf_dir / "complete_pmf_workflow.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_path, 0o755)
        
        logger.info(f"✓ Complete workflow generated: {script_path}")
        return str(script_path)
    
    def run_analysis(self, umbrella_results_dir=None):
        """Run WHAM analysis on completed umbrella sampling"""
        if umbrella_results_dir is None:
            umbrella_results_dir = self.umbrella_dir
        
        logger.info("Running WHAM analysis...")
        
        # Collect umbrella results
        window_dirs = list(Path(umbrella_results_dir).glob("window_*"))
        if not window_dirs:
            raise RuntimeError("No umbrella sampling results found")
        
        # Setup WHAM analysis
        wham_dir = self.analysis_dir / "wham"
        wham_dir.mkdir(exist_ok=True)
        
        # Generate file lists
        tpr_files = []
        pullf_files = []
        
        for window_dir in sorted(window_dirs):
            tpr_file = window_dir / "umbrella.tpr"
            pullf_file = window_dir / "umbrella_pullf.xvg"
            
            if tpr_file.exists() and pullf_file.exists():
                tpr_files.append(str(tpr_file.relative_to(wham_dir)))
                pullf_files.append(str(pullf_file.relative_to(wham_dir)))
        
        if len(tpr_files) < 5:
            raise RuntimeError(f"Insufficient umbrella windows: {len(tpr_files)} (need ≥5)")
        
        # Write file lists
        with open(wham_dir / "tpr-files.dat", 'w') as f:
            f.write('\n'.join(tpr_files))
        
        with open(wham_dir / "pullf-files.dat", 'w') as f:
            f.write('\n'.join(pullf_files))
        
        # Run WHAM
        cmd = [
            "gmx", "wham",
            "-it", "tpr-files.dat",
            "-if", "pullf-files.dat", 
            "-b", str(self.wham_begin_time),
            "-o", "pmf.xvg",
            "-unit", "kCal",
            "-bsres", "pmferror.xvg",
            "-bsprof", "zerrorprofile.xvg",
            "-nBootstrap", str(self.bootstrap_iterations),
            "-bs-method", "b-hist"
        ]
        
        try:
            result = subprocess.run(cmd, cwd=wham_dir, capture_output=True, text=True, check=True)
            logger.info("✓ WHAM analysis completed successfully")
            
            # Generate visualization
            self._plot_pmf_results(wham_dir)
            
            return str(wham_dir)
        
        except subprocess.CalledProcessError as e:
            logger.error(f"WHAM analysis failed: {e.stderr}")
            raise
    
    def _plot_pmf_results(self, wham_dir):
        """Generate PMF plots from WHAM results"""
        pmf_file = Path(wham_dir) / "pmf.xvg"
        error_file = Path(wham_dir) / "pmferror.xvg"
        
        if not pmf_file.exists():
            logger.warning("PMF file not found, skipping plot generation")
            return
        
        # Read PMF data
        pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
        distances = pmf_data[:, 0]
        pmf_values = pmf_data[:, 1]
        
        # Read error data if available
        errors = None
        if error_file.exists():
            error_data = np.loadtxt(error_file, comments=['#', '@'])
            errors = error_data[:, 1]
        
        # Create PMF plot
        plt.figure(figsize=(12, 8))
        
        if errors is not None:
            plt.errorbar(distances, pmf_values, yerr=errors, fmt='o-', capsize=3, linewidth=2)
        else:
            plt.plot(distances, pmf_values, 'o-', linewidth=2)
        
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.ylabel('PMF (kcal/mol)', fontsize=14)
        plt.title('Potential of Mean Force', fontsize=16)
        plt.grid(True, alpha=0.3)
        
        # Add binding energy annotation
        min_pmf = np.min(pmf_values)
        max_pmf = np.max(pmf_values)
        binding_energy = max_pmf - min_pmf
        
        plt.text(0.05, 0.95, f'Binding Energy: {binding_energy:.2f} kcal/mol', 
                transform=plt.gca().transAxes, fontsize=12, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = Path(wham_dir) / "pmf_curve.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"✓ PMF plot saved: {plot_file}")
        logger.info(f"✓ Estimated binding energy: {binding_energy:.2f} kcal/mol")
    
    def run_complete_workflow(self):
        """Run the complete PMF calculation workflow"""
        logger.info("=== Starting Complete PMF Workflow ===")
        
        try:
            # Step 1: Generate MDP files
            smd_mdp, umbrella_mdp = self.generate_mdp_files()
            
            # Step 2: Generate complete workflow script
            workflow_script = self.generate_complete_workflow()
            
            logger.info("=== PMF Workflow Setup Complete ===")
            logger.info(f"Workflow script: {workflow_script}")
            logger.info("Execute the workflow script to run complete PMF calculation")
            
            return {
                'workflow_script': workflow_script,
                'smd_mdp': smd_mdp,
                'umbrella_mdp': umbrella_mdp,
                'pmf_dir': str(self.pmf_dir)
            }
            
        except Exception as e:
            logger.error(f"PMF workflow setup failed: {e}")
            raise


def run_pmf_analysis(config: Dict, output_dir: str, md_results_dir: str) -> str:
    """
    Main entry point for PMF analysis
    
    Parameters:
    -----------
    config : Dict
        PMF configuration dictionary
    output_dir : str
        Output directory for PMF results
    md_results_dir : str
        Directory containing MD simulation results
    
    Returns:
    --------
    str : Path to the PMF results directory
    """
    if not config.get('enabled', False):
        logger.info("PMF analysis is disabled in configuration")
        return None
    
    try:
        calculator = PMFCalculator(config, output_dir, md_results_dir)
        results = calculator.run_complete_workflow()
        return results['pmf_dir']
    except Exception as e:
        logger.error(f"PMF analysis failed: {e}")
        raise


# Usage example
if __name__ == "__main__":
    # Example configuration
    config = {
        'enabled': True,
        'method': 'umbrella',
        'reference_group': 'Protein',
        'moving_group': 'LIG',
        'smd': {
            'pull_rate': 0.005,
            'pull_k': 1000.0,
            'pull_direction': [-1, 0, 0]
        },
        'distance': {
            'start': 0.3,
            'end': 2.0,
            'force_constant': 1000.0
        },
        'simulation': {
            'dt': 0.002,
            'temperature': 310.0,
            'pressure': 1.0
        },
        'umbrella': {
            'n_windows': 20,
            'production_time_ps': 22000,
            'sampling_interval_ps': 10
        },
        'analysis': {
            'method': 'wham',
            'begin_time_ps': 2000,
            'bootstrap_iterations': 50
        }
    }
    
    # Run PMF analysis
    pmf_results = run_pmf_analysis(
        config=config,
        output_dir="./pmf_output",
        md_results_dir="./md_results"
    )
    
    print(f"PMF analysis results: {pmf_results}")