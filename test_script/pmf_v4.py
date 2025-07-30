#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple API PMF (Potential of Mean Force) calculation system

Usage:
    import pmf
    
    # Setup global configuration
    pmf.setup(config, output_dir, md_results_dir)
    
    # Run individual modules
    pmf.smd_preparation()
    pmf.smd()
    pmf.umbrella_sampling()
    pmf.pmf_analysis()

Author: PRISM Team (Simple API)
Version: 4.0
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

# Global configuration storage
_CONFIG = None
_OUTPUT_DIR = None
_MD_RESULTS_DIR = None
_PMF_DIR = None
_MDP_DIR = None


def setup(config: Dict, output_dir: str, md_results_dir: str = None):
    """
    Setup global PMF configuration
    
    Parameters:
    -----------
    config : Dict
        PMF configuration dictionary
    output_dir : str
        Output directory for PMF results
    md_results_dir : str, optional
        Directory containing MD simulation results (required for smd_preparation)
    
    Example:
    --------
    >>> import pmf
    >>> config = {
    ...     'reference_group': 'Protein',
    ...     'moving_group': 'LIG',
    ...     'smd': {'pull_rate': 0.005, 'pull_k': 1000.0}
    ... }
    >>> pmf.setup(config, "./pmf_output", "./md_results")
    """
    global _CONFIG, _OUTPUT_DIR, _MD_RESULTS_DIR, _PMF_DIR, _MDP_DIR
    
    _CONFIG = config
    _OUTPUT_DIR = Path(output_dir).resolve()
    _MD_RESULTS_DIR = Path(md_results_dir).resolve() if md_results_dir else None
    
    # Setup directory structure
    _PMF_DIR = _OUTPUT_DIR / "pmf"
    _MDP_DIR = _OUTPUT_DIR / "mdps"
    
    for directory in [_PMF_DIR, _MDP_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    
    logger.info("PMF system configured successfully")
    logger.info(f"Output directory: {_OUTPUT_DIR}")
    logger.info(f"PMF directory: {_PMF_DIR}")


def _check_setup():
    """Check if PMF system has been setup"""
    if _CONFIG is None:
        raise RuntimeError("PMF system not configured. Please call pmf.setup() first.")


def _get_config(section: str = None, key: str = None, default=None):
    """Get configuration value"""
    _check_setup()
    if section is None:
        return _CONFIG
    if key is None:
        return _CONFIG.get(section, {})
    return _CONFIG.get(section, {}).get(key, default)


def _create_common_mdp_sections():
    """Create common MDP file sections"""
    temp = _get_config('simulation', 'temperature', 310.0)
    pressure = _get_config('simulation', 'pressure', 1.0)
    
    return f"""
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
ref_t               = {temp}     {temp}

; Pressure coupling
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 1.0
ref_p               = {pressure}
compressibility     = 4.5e-5

; Periodic boundary conditions
pbc                 = xyz
DispCorr            = EnerPres

; Velocity generation
gen_vel             = no

; COM motion removal
refcoord_scaling    = com
comm-mode           = Linear
comm-grps           = Protein Non-Protein"""


def smd_preparation():
    """
    SMD preparation module
    
    Prepares the system for SMD simulation by:
    - Validating MD results
    - Copying necessary files
    - Creating index files
    - Generating SMD MDP file
    
    Returns:
    --------
    Dict : Preparation results
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output", "./md_results")
    >>> results = pmf.smd_preparation()
    >>> print(f"SMD preparation completed: {results['status']}")
    """
    _check_setup()
    logger.info("=== SMD Preparation Module ===")
    
    if _MD_RESULTS_DIR is None:
        raise ValueError("md_results_dir is required for SMD preparation. Please provide it in pmf.setup()")
    
    smd_dir = _PMF_DIR / "smd"
    smd_dir.mkdir(exist_ok=True)
    
    # Step 1: Validate MD results
    logger.info("Validating MD results...")
    required_files = {
        'gro': 'solv_ions.gro',
        'top': 'topol.top'
    }
    
    found_files = {}
    for file_type, filename in required_files.items():
        file_path = _MD_RESULTS_DIR / filename
        if file_path.exists():
            found_files[file_type] = str(file_path)
            logger.info(f"✓ Found {file_type}: {filename}")
        else:
            raise FileNotFoundError(f"Required file not found: {file_path}")
    
    # Step 2: Copy MD files
    logger.info("Copying MD files...")
    for file_type, file_path in found_files.items():
        dest_path = _PMF_DIR / Path(file_path).name
        shutil.copy2(file_path, dest_path)
        logger.info(f"✓ Copied {file_type}: {Path(file_path).name}")
    
    # Step 3: Create index file script
    logger.info("Creating index file script...")
    moving_group = _get_config('moving_group', default='LIG')
    
    index_script = _PMF_DIR / "create_index.sh"
    script_content = f"""#!/bin/bash
# Create index file for SMD groups
echo "Creating index file for groups..."

cd {_PMF_DIR}

# Create index file
echo "r {moving_group}\\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx

echo "✓ Index file created: index.ndx"
"""
    
    with open(index_script, 'w') as f:
        f.write(script_content)
    os.chmod(index_script, 0o755)
    
    # Step 4: Generate SMD MDP file
    logger.info("Generating SMD MDP file...")
    
    # Get configuration
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    pull_k = _get_config('smd', 'pull_k', 1000.0)
    pull_direction = _get_config('smd', 'pull_direction', [-1, 0, 0])
    dt = _get_config('simulation', 'dt', 0.002)
    
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    
    # Calculate simulation parameters
    distance_range = distance_end - distance_start
    smd_time_ps = distance_range / pull_rate
    nsteps = int(smd_time_ps / dt)
    output_interval = int(5.0 / dt)
    
    smd_content = f"""; SMD (Steered Molecular Dynamics) Parameters
title               = SMD simulation ({distance_range:.2f} nm in {smd_time_ps:.1f} ps)
integrator          = md
nsteps              = {nsteps}
dt                  = {dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}
{_create_common_mdp_sections()}

; Pull code for SMD
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y N N
pull-coord1-vec     = {' '.join(map(str, pull_direction))}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = {pull_rate}
pull_coord1_k       = {pull_k}
pull-pbc-ref-prev-step-com = yes
"""
    
    smd_mdp_file = _MDP_DIR / "smd.mdp"
    with open(smd_mdp_file, 'w') as f:
        f.write(smd_content)
    
    logger.info(f"✓ SMD MDP file: {smd_mdp_file}")
    
    # Step 5: Generate preparation script
    prep_script = _PMF_DIR / "smd_preparation.sh"
    script_content = f"""#!/bin/bash
######################################################
# SMD PREPARATION SCRIPT
######################################################

set -e
echo "=== SMD Preparation ==="

cd {_PMF_DIR}

# Step 1: Create index file
if [ ! -f index.ndx ]; then
    echo "Creating index file..."
    echo "r {moving_group}\\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx
    echo "✓ Index file created"
else
    echo "✓ Index file already exists"
fi

# Step 2: Validate files
echo "Validating required files..."
for file in solv_ions.gro topol.top ../mdps/smd.mdp index.ndx; do
    if [ -f "$file" ]; then
        echo "✓ Found: $file"
    else
        echo "✗ Missing: $file"
        exit 1
    fi
done

echo "✓ SMD preparation completed successfully!"
echo "Next step: Run pmf.smd()"
"""
    
    with open(prep_script, 'w') as f:
        f.write(script_content)
    os.chmod(prep_script, 0o755)
    
    results = {
        'smd_dir': str(smd_dir),
        'smd_mdp': str(smd_mdp_file),
        'preparation_script': str(prep_script),
        'index_script': str(index_script),
        'status': 'ready_for_smd'
    }
    
    logger.info("✓ SMD preparation completed successfully!")
    return results


def smd(force_restart=False):
    """
    SMD simulation module
    
    Runs steered molecular dynamics simulation including:
    - SMD simulation execution
    - Result analysis and plotting
    - Trajectory frame extraction
    
    Parameters:
    -----------
    force_restart : bool
        Whether to restart SMD even if results exist
    
    Returns:
    --------
    Dict : SMD simulation results
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> results = pmf.smd()
    >>> print(f"SMD completed: {results['smd_trajectory']}")
    """
    _check_setup()
    logger.info("=== SMD Simulation Module ===")
    
    smd_dir = _PMF_DIR / "smd"
    smd_dir.mkdir(exist_ok=True)
    
    # Step 1: Validate preparation
    logger.info("Validating SMD preparation...")
    required_files = [
        _PMF_DIR / "solv_ions.gro",
        _PMF_DIR / "topol.top",
        _PMF_DIR / "index.ndx",
        _MDP_DIR / "smd.mdp"
    ]
    
    for file_path in required_files:
        if not file_path.exists():
            raise FileNotFoundError(f"SMD preparation incomplete: {file_path} not found. Run pmf.smd_preparation() first.")
    
    logger.info("✓ SMD preparation validated")
    
    # Step 2: Check if SMD already completed
    smd_gro = smd_dir / "smd.gro"
    smd_xtc = smd_dir / "smd.xtc"
    
    if not force_restart and smd_gro.exists() and smd_xtc.exists():
        logger.info("✓ SMD already completed")
        return _get_smd_results()
    
    # Step 3: Run SMD simulation
    logger.info("Starting SMD simulation...")
    
    # Generate TPR file
    tpr_file = smd_dir / "smd.tpr"
    grompp_cmd = [
        "gmx", "grompp",
        "-f", str(_MDP_DIR / "smd.mdp"),
        "-c", str(_PMF_DIR / "solv_ions.gro"),
        "-n", str(_PMF_DIR / "index.ndx"),
        "-p", str(_PMF_DIR / "topol.top"),
        "-o", str(tpr_file),
        "-maxwarn", "10"
    ]
    
    try:
        subprocess.run(grompp_cmd, check=True, capture_output=True, text=True)
        logger.info("✓ SMD TPR file generated")
    except subprocess.CalledProcessError as e:
        logger.error(f"grompp failed: {e.stderr}")
        raise
    
    # Run MD simulation
    mdrun_cmd = [
        "gmx", "mdrun",
        "-s", str(tpr_file),
        "-deffnm", str(smd_dir / "smd"),
        "-ntmpi", "1",
        "-ntomp", "10",
        "-gpu_id", "0",
        "-v"
    ]
    
    try:
        subprocess.run(mdrun_cmd, check=True, capture_output=True, text=True)
        logger.info("✓ SMD simulation completed")
    except subprocess.CalledProcessError as e:
        logger.error(f"mdrun failed: {e.stderr}")
        raise
    
    # Step 4: Analyze SMD results
    logger.info("Analyzing SMD results...")
    _analyze_smd_results()
    
    # Step 5: Extract trajectory frames
    frames_info = _extract_trajectory_frames()
    
    results = {
        'smd_dir': str(smd_dir),
        'smd_trajectory': str(smd_dir / "smd.xtc"),
        'smd_structure': str(smd_dir / "smd.gro"),
        'smd_tpr': str(tpr_file),
        'trajectory_frames': frames_info,
        'status': 'smd_completed'
    }
    
    logger.info("✓ SMD simulation completed successfully!")
    return results


def _analyze_smd_results():
    """Analyze SMD results and create plots"""
    smd_dir = _PMF_DIR / "smd"
    
    # Plot force profile
    pullf_file = smd_dir / "smd_pullf.xvg"
    if pullf_file.exists():
        try:
            data = np.loadtxt(pullf_file, comments=['#', '@'])
            time = data[:, 0]
            force = data[:, 1]
            
            plt.figure(figsize=(12, 8))
            plt.plot(time, force, 'b-', linewidth=2)
            plt.xlabel('Time (ps)', fontsize=14)
            plt.ylabel('Force (kJ/mol/nm)', fontsize=14)
            plt.title('SMD Pull Force vs Time', fontsize=16)
            plt.grid(True, alpha=0.3)
            
            # Add statistics
            max_force = np.max(force)
            avg_force = np.mean(force)
            plt.text(0.05, 0.95, f'Max Force: {max_force:.2f} kJ/mol/nm\nAvg Force: {avg_force:.2f} kJ/mol/nm', 
                    transform=plt.gca().transAxes, fontsize=12, 
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            plt.tight_layout()
            plot_file = smd_dir / "force_profile.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"✓ Force profile plot: {plot_file}")
        except Exception as e:
            logger.warning(f"Could not plot force profile: {e}")
    
    # Plot distance profile
    pullx_file = smd_dir / "smd_pullx.xvg"
    if pullx_file.exists():
        try:
            data = np.loadtxt(pullx_file, comments=['#', '@'])
            time = data[:, 0]
            distance = data[:, 1]
            
            plt.figure(figsize=(12, 8))
            plt.plot(time, distance, 'r-', linewidth=2)
            plt.xlabel('Time (ps)', fontsize=14)
            plt.ylabel('Distance (nm)', fontsize=14)
            plt.title('SMD Distance vs Time', fontsize=16)
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plot_file = smd_dir / "distance_profile.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"✓ Distance profile plot: {plot_file}")
        except Exception as e:
            logger.warning(f"Could not plot distance profile: {e}")


def _extract_trajectory_frames(n_frames=300):
    """Extract trajectory frames for umbrella sampling"""
    logger.info(f"Extracting {n_frames} trajectory frames...")
    
    frames_dir = _PMF_DIR / "trajectory_frames"
    frames_dir.mkdir(exist_ok=True)
    
    # Calculate frame interval
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    
    distance_range = distance_end - distance_start
    smd_time_ps = distance_range / pull_rate
    frame_interval = int(smd_time_ps / n_frames)
    
    # Extract frames
    trjconv_cmd = [
        "gmx", "trjconv",
        "-f", str(_PMF_DIR / "smd" / "smd.xtc"),
        "-s", str(_PMF_DIR / "smd" / "smd.tpr"),
        "-o", str(frames_dir / "frame.gro"),
        "-sep",
        "-dt", str(frame_interval)
    ]
    
    try:
        process = subprocess.Popen(trjconv_cmd, stdin=subprocess.PIPE, 
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input="0\n")  # Select System
        
        if process.returncode == 0:
            logger.info(f"✓ Extracted trajectory frames to: {frames_dir}")
        else:
            logger.error(f"Frame extraction failed: {stderr}")
            raise RuntimeError("Frame extraction failed")
    except Exception as e:
        logger.error(f"Frame extraction error: {e}")
        raise
    
    return {
        'frames_dir': str(frames_dir),
        'n_frames': n_frames,
        'frame_interval': frame_interval
    }


def _get_smd_results():
    """Get SMD results when already completed"""
    smd_dir = _PMF_DIR / "smd"
    return {
        'smd_dir': str(smd_dir),
        'smd_trajectory': str(smd_dir / "smd.xtc"),
        'smd_structure': str(smd_dir / "smd.gro"),
        'smd_tpr': str(smd_dir / "smd.tpr"),
        'trajectory_frames': {
            'frames_dir': str(_PMF_DIR / "trajectory_frames"),
            'status': 'existing'
        },
        'status': 'smd_completed'
    }


def umbrella_sampling(auto_run=False):
    """
    Umbrella sampling module
    
    Generates umbrella sampling windows and optionally runs them:
    - Distance calculation for all frames
    - Adaptive window generation
    - Umbrella MDP file creation
    - Window setup and run scripts
    
    Parameters:
    -----------
    auto_run : bool
        Whether to automatically run all umbrella windows
    
    Returns:
    --------
    Dict : Umbrella sampling results
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> results = pmf.umbrella_sampling(auto_run=False)
    >>> print(f"Generated {results['n_windows']} umbrella windows")
    """
    _check_setup()
    logger.info("=== Umbrella Sampling Module ===")
    
    umbrella_dir = _PMF_DIR / "umbrella"
    umbrella_dir.mkdir(exist_ok=True)
    
    # Step 1: Validate SMD results
    logger.info("Validating SMD results...")
    required_files = [
        _PMF_DIR / "smd" / "smd.xtc",
        _PMF_DIR / "smd" / "smd.tpr",
        _PMF_DIR / "trajectory_frames"
    ]
    
    for file_path in required_files:
        if not file_path.exists():
            raise FileNotFoundError(f"SMD results incomplete: {file_path} not found. Run pmf.smd() first.")
    
    logger.info("✓ SMD results validated")
    
    # Step 2: Calculate distances for all frames
    distance_data = _calculate_frame_distances()
    
    # Step 3: Generate umbrella windows
    windows = _generate_umbrella_windows(distance_data)
    
    # Step 4: Generate umbrella MDP template
    umbrella_mdp = _generate_umbrella_mdp()
    
    # Step 5: Setup all umbrella windows
    _setup_umbrella_windows(windows)
    
    # Step 6: Generate run scripts
    run_scripts = _generate_umbrella_run_scripts(windows)
    
    results = {
        'umbrella_dir': str(umbrella_dir),
        'windows': windows,
        'umbrella_mdp': umbrella_mdp,
        'run_scripts': run_scripts,
        'n_windows': len(windows),
        'status': 'ready_for_umbrella'
    }
    
    # Step 7: Optionally run all windows
    if auto_run:
        logger.info("Running all umbrella sampling windows...")
        execution_results = {}
        for window in windows:
            try:
                _run_single_umbrella_window(window)
                execution_results[window['window_id']] = 'completed'
                logger.info(f"✓ Window {window['window_id']:03d} completed")
            except Exception as e:
                execution_results[window['window_id']] = f'failed: {e}'
                logger.error(f"✗ Window {window['window_id']:03d} failed: {e}")
        
        results['execution_results'] = execution_results
        results['status'] = 'umbrella_completed'
    
    logger.info(f"✓ Umbrella sampling setup completed ({len(windows)} windows)")
    return results


def _calculate_frame_distances():
    """Calculate distances for all trajectory frames"""
    logger.info("Calculating distances for trajectory frames...")
    
    frames_dir = _PMF_DIR / "trajectory_frames"
    distance_dir = _PMF_DIR / "distance_analysis"
    distance_dir.mkdir(exist_ok=True)
    
    summary_file = distance_dir / "summary_distances.dat"
    
    if summary_file.exists():
        logger.info("✓ Distance data already exists")
        return _read_distance_data(summary_file)
    
    # Get configuration
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    
    # Find all frame files
    frame_files = sorted(frames_dir.glob("frame_*.gro"))
    
    with open(summary_file, 'w') as f:
        f.write("# Frame Distance(nm)\n")
        
        for i, frame_file in enumerate(frame_files):
            # Calculate distance using gmx distance
            distance_cmd = [
                "gmx", "distance",
                "-s", str(_PMF_DIR / "smd" / "smd.tpr"),
                "-f", str(frame_file),
                "-n", str(_PMF_DIR / "index.ndx"),
                "-select", f"com of group {reference_group} plus com of group {moving_group}",
                "-oall", str(distance_dir / f"dist_{i}.xvg")
            ]
            
            try:
                subprocess.run(distance_cmd, check=True, capture_output=True, text=True)
                
                # Read distance value
                dist_file = distance_dir / f"dist_{i}.xvg"
                with open(dist_file, 'r') as df:
                    lines = df.readlines()
                    for line in reversed(lines):
                        if not line.startswith('#') and not line.startswith('@'):
                            distance = float(line.split()[1])
                            f.write(f"{i} {distance:.6f}\n")
                            break
                
                # Clean up temporary file
                dist_file.unlink()
                
            except subprocess.CalledProcessError as e:
                logger.warning(f"Could not calculate distance for frame {i}: {e}")
                continue
    
    logger.info(f"✓ Distance calculation completed: {summary_file}")
    return _read_distance_data(summary_file)


def _read_distance_data(summary_file):
    """Read distance data from file"""
    data = []
    with open(summary_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    frame = int(parts[0])
                    distance = float(parts[1])
                    data.append((frame, distance))
    
    logger.info(f"✓ Read {len(data)} distance points")
    return data


def _generate_umbrella_windows(distance_data):
    """Generate adaptive umbrella sampling windows"""
    logger.info("Generating umbrella sampling windows...")
    
    if not distance_data:
        raise RuntimeError("No distance data available for window generation")
    
    # Get umbrella configuration
    sample_interval_near = _get_config('umbrella', 'sample_interval_near', 0.1)
    sample_interval_far = _get_config('umbrella', 'sample_interval_far', 0.2)
    cutoff_distance = _get_config('umbrella', 'cutoff_distance', 1.5)
    
    # Adaptive sampling
    sampled_indices = []
    current_idx = 0
    sampled_indices.append(current_idx)
    
    distances = [d[1] for d in distance_data]
    
    while current_idx < len(distances):
        current_distance = distances[current_idx]
        
        # Choose interval based on distance
        if current_distance < cutoff_distance:
            target_interval = sample_interval_near
        else:
            target_interval = sample_interval_far
        
        target_distance = current_distance + target_interval
        
        # Find next frame closest to target distance
        remaining_distances = distances[current_idx:]
        diffs = [abs(target_distance - d) for d in remaining_distances]
        next_relative_idx = diffs.index(min(diffs))
        next_idx = current_idx + next_relative_idx
        
        if next_idx == current_idx:
            break
        
        sampled_indices.append(next_idx)
        current_idx = next_idx
    
    # Create window information
    windows = []
    umbrella_dir = _PMF_DIR / "umbrella"
    
    for i, idx in enumerate(sampled_indices):
        frame, distance = distance_data[idx]
        windows.append({
            'window_id': i,
            'frame': frame,
            'distance': distance,
            'directory': umbrella_dir / f"window_{i:03d}"
        })
    
    logger.info(f"✓ Generated {len(windows)} umbrella windows")
    return windows


def _generate_umbrella_mdp():
    """Generate umbrella sampling MDP template"""
    # Get configuration
    dt = _get_config('simulation', 'dt', 0.002)
    production_time = _get_config('umbrella', 'production_time_ps', 22000)
    sampling_interval = _get_config('umbrella', 'sampling_interval_ps', 10)
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    pull_k = _get_config('smd', 'pull_k', 1000.0)
    pull_direction = _get_config('smd', 'pull_direction', [-1, 0, 0])
    
    nsteps = int(production_time / dt)
    output_interval = int(sampling_interval / dt)
    
    umbrella_content = f"""; Umbrella Sampling Parameters
title               = Umbrella sampling (22.0 ns)
integrator          = md
nsteps              = {nsteps}
dt                  = {dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}
{_create_common_mdp_sections()}

; Pull code for umbrella sampling (static constraint)
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y N N
pull-coord1-vec     = {' '.join(map(str, pull_direction))}
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = 0.0
pull_coord1_k       = {pull_k}
pull_coord1_r0      = DISTANCE_PLACEHOLDER
pull-pbc-ref-prev-step-com = yes
"""
    
    umbrella_mdp_file = _MDP_DIR / "umbrella.mdp"
    with open(umbrella_mdp_file, 'w') as f:
        f.write(umbrella_content)
    
    logger.info(f"✓ Umbrella MDP template: {umbrella_mdp_file}")
    return str(umbrella_mdp_file)


def _setup_umbrella_windows(windows):
    """Setup all umbrella sampling windows"""
    logger.info("Setting up umbrella sampling windows...")
    
    for window in windows:
        window_dir = window['directory']
        window_dir.mkdir(exist_ok=True)
        
        # Create window-specific MDP file
        template_file = _MDP_DIR / "umbrella.mdp"
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        content = content.replace("DISTANCE_PLACEHOLDER", f"{window['distance']:.6f}")
        content = content.replace(
            "title               = Umbrella sampling (22.0 ns)",
            f"title               = Umbrella sampling at {window['distance']:.3f} nm (window {window['window_id']:03d})"
        )
        
        # Write window-specific MDP
        window_mdp = window_dir / "umbrella.mdp"
        with open(window_mdp, 'w') as f:
            f.write(content)
        
        # Copy starting structure
        frame_file = _PMF_DIR / "trajectory_frames" / f"frame_{window['frame']}.gro"
        dest_file = window_dir / "start.gro"
        
        if frame_file.exists():
            shutil.copy2(frame_file, dest_file)
        else:
            logger.warning(f"Frame file not found: {frame_file}")
    
    logger.info("✓ All umbrella windows setup completed")


def _generate_umbrella_run_scripts(windows):
    """Generate run scripts for umbrella sampling"""
    umbrella_dir = _PMF_DIR / "umbrella"
    
    # Individual window scripts
    for window in windows:
        script_content = f"""#!/bin/bash
set -e

echo "=== Running Umbrella Window {window['window_id']:03d} ==="
echo "Distance: {window['distance']:.3f} nm (Frame: {window['frame']})"

cd {window['directory']}

# Generate TPR file
echo "Generating TPR file..."
gmx grompp -f umbrella.mdp -c start.gro -n ../../index.ndx -p ../../topol.top -o umbrella.tpr -maxwarn 10

# Run umbrella sampling
echo "Running umbrella sampling..."
gmx mdrun -s umbrella.tpr -deffnm umbrella -ntmpi 1 -ntomp 10 -gpu_id 0 -v

echo "✓ Window {window['window_id']:03d} completed!"
"""
        
        script_file = window['directory'] / "run_umbrella.sh"
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_file, 0o755)
    
    # Master run script
    script_content = f"""#!/bin/bash
######################################################
# MASTER UMBRELLA SAMPLING SCRIPT
######################################################

set -e

echo "=== Running All Umbrella Sampling Windows ==="
echo "Total windows: {len(windows)}"

cd {umbrella_dir}

# Run each window
"""
    
    for window in windows:
        script_content += f"""
echo "Starting window {window['window_id']:03d}/{len(windows)-1}: {window['distance']:.3f} nm"
cd window_{window['window_id']:03d}
bash run_umbrella.sh
cd ..
"""
    
    script_content += """
echo "✓ All umbrella sampling windows completed!"
echo "Results ready for PMF analysis"
"""
    
    master_script = umbrella_dir / "run_all_umbrella.sh"
    with open(master_script, 'w') as f:
        f.write(script_content)
    
    os.chmod(master_script, 0o755)
    
    logger.info(f"✓ Umbrella run scripts generated")
    
    return {
        'master_script': str(master_script),
        'individual_scripts': [str(w['directory'] / "run_umbrella.sh") for w in windows]
    }


def _run_single_umbrella_window(window):
    """Run single umbrella sampling window"""
    script_file = window['directory'] / "run_umbrella.sh"
    subprocess.run(["bash", str(script_file)], check=True)


def pmf_analysis(umbrella_dir=None):
    """
    PMF analysis module
    
    Performs WHAM analysis and generates visualizations:
    - Validates umbrella sampling results
    - Runs WHAM analysis
    - Calculates binding energy
    - Generates PMF plots and reports
    
    Parameters:
    -----------
    umbrella_dir : str, optional
        Directory containing umbrella sampling results
    
    Returns:
    --------
    Dict : PMF analysis results
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> results = pmf.pmf_analysis()
    >>> print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
    """
    _check_setup()
    logger.info("=== PMF Analysis Module ===")
    
    if umbrella_dir is None:
        umbrella_dir = _PMF_DIR / "umbrella"
    else:
        umbrella_dir = Path(umbrella_dir)
    
    analysis_dir = _PMF_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    
    # Step 1: Validate umbrella results
    window_files = _validate_umbrella_results(umbrella_dir)
    
    # Step 2: Run WHAM analysis
    wham_results = _run_wham_analysis(window_files, analysis_dir)
    
    # Step 3: Generate visualizations
    plots = _generate_pmf_visualizations(wham_results)
    
    # Step 4: Calculate binding energy
    binding_energy = _calculate_binding_energy(wham_results)
    
    # Step 5: Generate analysis report
    report = _generate_pmf_report(wham_results, binding_energy, plots, analysis_dir)
    
    results = {
        'analysis_dir': str(analysis_dir),
        'wham_results': wham_results,
        'binding_energy': binding_energy,
        'plots': plots,
        'report': report,
        'status': 'analysis_completed'
    }
    
    logger.info("✓ PMF analysis completed successfully!")
    return results


def _validate_umbrella_results(umbrella_dir):
    """Validate umbrella sampling results"""
    logger.info("Validating umbrella sampling results...")
    
    if not umbrella_dir.exists():
        raise FileNotFoundError(f"Umbrella directory not found: {umbrella_dir}")
    
    window_dirs = sorted(umbrella_dir.glob("window_*"))
    valid_windows = []
    
    for window_dir in window_dirs:
        tpr_file = window_dir / "umbrella.tpr"
        pullf_file = window_dir / "umbrella_pullf.xvg"
        
        if tpr_file.exists() and pullf_file.exists():
            valid_windows.append({
                'window_dir': window_dir,
                'tpr_file': tpr_file,
                'pullf_file': pullf_file
            })
        else:
            logger.warning(f"Incomplete window: {window_dir}")
    
    if len(valid_windows) < 5:
        raise RuntimeError(f"Insufficient umbrella windows: {len(valid_windows)} (need ≥5)")
    
    logger.info(f"✓ Found {len(valid_windows)} valid umbrella windows")
    return valid_windows


def _run_wham_analysis(window_files, analysis_dir):
    """Run WHAM analysis"""
    logger.info("Running WHAM analysis...")
    
    wham_dir = analysis_dir / "wham"
    wham_dir.mkdir(exist_ok=True)
    
    # Get analysis configuration
    begin_time = _get_config('analysis', 'begin_time_ps', 2000)
    bootstrap_iterations = _get_config('analysis', 'bootstrap_iterations', 50)
    energy_unit = _get_config('analysis', 'energy_unit', 'kCal')
    
    # Generate file lists
    tpr_list = wham_dir / "tpr-files.dat"
    pullf_list = wham_dir / "pullf-files.dat"
    
    with open(tpr_list, 'w') as f:
        for window in window_files:
            relative_path = os.path.relpath(window['tpr_file'], wham_dir)
            f.write(f"{relative_path}\n")
    
    with open(pullf_list, 'w') as f:
        for window in window_files:
            relative_path = os.path.relpath(window['pullf_file'], wham_dir)
            f.write(f"{relative_path}\n")
    
    # Run gmx wham
    wham_cmd = [
        "gmx", "wham",
        "-it", "tpr-files.dat",
        "-if", "pullf-files.dat",
        "-b", str(begin_time),
        "-o", "pmf.xvg",
        "-unit", energy_unit,
        "-bsres", "pmferror.xvg",
        "-bsprof", "zerrorprofile.xvg",
        "-nBootstrap", str(bootstrap_iterations),
        "-bs-method", "b-hist"
    ]
    
    try:
        result = subprocess.run(wham_cmd, cwd=wham_dir, capture_output=True, text=True, check=True)
        logger.info("✓ WHAM analysis completed successfully")
        
        return {
            'wham_dir': str(wham_dir),
            'pmf_file': str(wham_dir / "pmf.xvg"),
            'error_file': str(wham_dir / "pmferror.xvg"),
            'profile_file': str(wham_dir / "zerrorprofile.xvg")
        }
    except subprocess.CalledProcessError as e:
        logger.error(f"WHAM analysis failed: {e.stderr}")
        raise


def _generate_pmf_visualizations(wham_results):
    """Generate PMF visualizations"""
    logger.info("Generating PMF visualizations...")
    
    pmf_file = Path(wham_results['pmf_file'])
    error_file = Path(wham_results['error_file'])
    
    # Read PMF data
    pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
    distances = pmf_data[:, 0]
    pmf_values = pmf_data[:, 1]
    
    # Read error data
    errors = None
    if error_file.exists():
        error_data = np.loadtxt(error_file, comments=['#', '@'])
        errors = error_data[:, 1]
    
    # Get energy unit
    energy_unit = _get_config('analysis', 'energy_unit', 'kCal')
    
    # Create PMF plot
    plt.figure(figsize=(12, 8))
    
    if errors is not None:
        plt.errorbar(distances, pmf_values, yerr=errors, fmt='o-', capsize=3, linewidth=2, markersize=4)
    else:
        plt.plot(distances, pmf_values, 'o-', linewidth=2, markersize=4)
    
    plt.xlabel('Distance (nm)', fontsize=14)
    plt.ylabel(f'PMF ({energy_unit.lower()}/mol)', fontsize=14)
    plt.title('Potential of Mean Force', fontsize=16)
    plt.grid(True, alpha=0.3)
    
    # Add binding energy annotation
    min_pmf = np.min(pmf_values)
    max_pmf = np.max(pmf_values)
    binding_energy = max_pmf - min_pmf
    
    plt.text(0.05, 0.95, f'Binding Energy: {binding_energy:.2f} {energy_unit.lower()}/mol', 
            transform=plt.gca().transAxes, fontsize=12, 
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    analysis_dir = Path(wham_results['wham_dir']).parent
    plot_file = analysis_dir / "pmf_curve.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ PMF plot saved: {plot_file}")
    
    return {
        'pmf_plot': str(plot_file),
        'pmf_data': {
            'distances': distances.tolist(),
            'pmf_values': pmf_values.tolist(),
            'errors': errors.tolist() if errors is not None else None
        }
    }


def _calculate_binding_energy(wham_results):
    """Calculate binding energy from PMF"""
    pmf_file = Path(wham_results['pmf_file'])
    pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
    pmf_values = pmf_data[:, 1]
    
    min_pmf = np.min(pmf_values)
    max_pmf = np.max(pmf_values)
    binding_energy = max_pmf - min_pmf
    
    energy_unit = _get_config('analysis', 'energy_unit', 'kCal')
    
    logger.info(f"✓ Binding energy: {binding_energy:.2f} {energy_unit.lower()}/mol")
    
    return {
        'value': binding_energy,
        'unit': f'{energy_unit.lower()}/mol',
        'min_pmf': min_pmf,
        'max_pmf': max_pmf
    }


def _generate_pmf_report(wham_results, binding_energy, plots, analysis_dir):
    """Generate comprehensive analysis report"""
    report_file = analysis_dir / "pmf_analysis_report.txt"
    
    # Get configuration for report
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    begin_time = _get_config('analysis', 'begin_time_ps', 2000)
    bootstrap_iterations = _get_config('analysis', 'bootstrap_iterations', 50)
    energy_unit = _get_config('analysis', 'energy_unit', 'kCal')
    
    report_content = f"""
PMF Analysis Report
==================

Analysis Date: {os.popen('date').read().strip()}
Analysis Directory: {analysis_dir}

WHAM Analysis Results:
---------------------
- PMF file: {wham_results['pmf_file']}
- Error file: {wham_results['error_file']}
- Profile file: {wham_results['profile_file']}

Binding Energy:
--------------
- Binding Energy: {binding_energy['value']:.2f} {binding_energy['unit']}
- Minimum PMF: {binding_energy['min_pmf']:.2f} {binding_energy['unit']}
- Maximum PMF: {binding_energy['max_pmf']:.2f} {binding_energy['unit']}

Configuration Parameters:
------------------------
- Reference Group: {reference_group}
- Moving Group: {moving_group}
- WHAM Begin Time: {begin_time} ps
- Bootstrap Iterations: {bootstrap_iterations}
- Energy Unit: {energy_unit}

Generated Files:
---------------
- PMF curve plot: {plots['pmf_plot']}
- Analysis report: {report_file}

Analysis completed successfully!
"""
    
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    logger.info(f"✓ Analysis report: {report_file}")
    
    return str(report_file)


# Utility functions
def get_status():
    """
    Get current PMF workflow status
    
    Returns:
    --------
    Dict : Status of each module
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> status = pmf.get_status()
    >>> print(status)
    """
    _check_setup()
    
    status = {}
    
    # Check SMD preparation
    smd_prep_files = [
        _PMF_DIR / "solv_ions.gro",
        _PMF_DIR / "topol.top",
        _PMF_DIR / "index.ndx",
        _MDP_DIR / "smd.mdp"
    ]
    status['smd_preparation'] = all(f.exists() for f in smd_prep_files)
    
    # Check SMD
    smd_files = [
        _PMF_DIR / "smd" / "smd.gro",
        _PMF_DIR / "smd" / "smd.xtc",
        _PMF_DIR / "trajectory_frames"
    ]
    status['smd'] = all(f.exists() for f in smd_files)
    
    # Check umbrella sampling
    umbrella_dir = _PMF_DIR / "umbrella"
    if umbrella_dir.exists():
        window_dirs = list(umbrella_dir.glob("window_*"))
        status['umbrella_sampling'] = {
            'setup': len(window_dirs) > 0,
            'n_windows': len(window_dirs),
            'completed_windows': len([w for w in window_dirs if (w / "umbrella.gro").exists()])
        }
    else:
        status['umbrella_sampling'] = {'setup': False, 'n_windows': 0, 'completed_windows': 0}
    
    # Check PMF analysis
    analysis_files = [
        _PMF_DIR / "analysis" / "wham" / "pmf.xvg",
        _PMF_DIR / "analysis" / "pmf_curve.png"
    ]
    status['pmf_analysis'] = all(f.exists() for f in analysis_files)
    
    return status


def clean(modules=None):
    """
    Clean PMF results
    
    Parameters:
    -----------
    modules : list, optional
        List of modules to clean ['smd_preparation', 'smd', 'umbrella_sampling', 'pmf_analysis']
        If None, cleans all modules
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> pmf.clean(['smd'])  # Only clean SMD results
    >>> pmf.clean()  # Clean all results
    """
    _check_setup()
    
    if modules is None:
        modules = ['smd_preparation', 'smd', 'umbrella_sampling', 'pmf_analysis']
    
    logger.info(f"Cleaning PMF modules: {modules}")
    
    if 'smd_preparation' in modules:
        files_to_remove = [
            _PMF_DIR / "solv_ions.gro",
            _PMF_DIR / "topol.top", 
            _PMF_DIR / "index.ndx",
            _MDP_DIR / "smd.mdp"
        ]
        for file_path in files_to_remove:
            if file_path.exists():
                file_path.unlink()
                logger.info(f"✓ Removed: {file_path}")
    
    if 'smd' in modules:
        smd_dir = _PMF_DIR / "smd"
        frames_dir = _PMF_DIR / "trajectory_frames"
        distance_dir = _PMF_DIR / "distance_analysis"
        
        for dir_path in [smd_dir, frames_dir, distance_dir]:
            if dir_path.exists():
                shutil.rmtree(dir_path)
                logger.info(f"✓ Removed: {dir_path}")
    
    if 'umbrella_sampling' in modules:
        umbrella_dir = _PMF_DIR / "umbrella"
        umbrella_mdp = _MDP_DIR / "umbrella.mdp"
        
        if umbrella_dir.exists():
            shutil.rmtree(umbrella_dir)
            logger.info(f"✓ Removed: {umbrella_dir}")
        
        if umbrella_mdp.exists():
            umbrella_mdp.unlink()
            logger.info(f"✓ Removed: {umbrella_mdp}")
    
    if 'pmf_analysis' in modules:
        analysis_dir = _PMF_DIR / "analysis"
        if analysis_dir.exists():
            shutil.rmtree(analysis_dir)
            logger.info(f"✓ Removed: {analysis_dir}")
    
    logger.info("✓ Cleaning completed")


# Example usage and documentation
if __name__ == "__main__":
    # Example configuration
    example_config = {
        'reference_group': 'Protein',
        'moving_group': 'LIG',
        'smd': {
            'pull_rate': 0.005,
            'pull_k': 1000.0,
            'pull_direction': [-1, 0, 0]
        },
        'distance': {
            'start': 0.3,
            'end': 2.0
        },
        'simulation': {
            'dt': 0.002,
            'temperature': 310.0,
            'pressure': 1.0
        },
        'umbrella': {
            'production_time_ps': 22000,
            'sampling_interval_ps': 10,
            'sample_interval_near': 0.1,
            'sample_interval_far': 0.2,
            'cutoff_distance': 1.5
        },
        'analysis': {
            'begin_time_ps': 2000,
            'bootstrap_iterations': 50,
            'energy_unit': 'kCal'
        }
    }
    
    print("PMF Simple API - Usage Examples:")
    print("=" * 50)
    print("""
# Setup PMF system
import pmf
pmf.setup(config, "./pmf_output", "./md_results")

# Run individual modules
pmf.smd_preparation()  # Prepare system for SMD
pmf.smd()              # Run SMD simulation
pmf.umbrella_sampling() # Setup umbrella sampling windows
pmf.pmf_analysis()     # Run WHAM analysis

# Check status
status = pmf.get_status()
print(status)

# Clean specific modules
pmf.clean(['smd'])     # Clean only SMD results
pmf.clean()            # Clean all results
""")
    print("PMF Simple API ready for use!")