"""
Step-by-Step PMF (Potential of Mean Force) calculation system

This system is designed for non-continuous execution where users manually
run GROMACS simulations between each step.

Workflow:
1. pmf_step1_smd_preparation() - Prepare SMD setup and generate run scripts
2. [USER RUNS GROMACS SMD MANUALLY]
3. pmf_step2_umbrella_preparation() - Analyze SMD results and setup umbrella windows
4. [USER RUNS GROMACS UMBRELLA SAMPLING MANUALLY]
5. pmf_step3_wham_analysis() - Perform WHAM analysis and generate PMF

Analysis modules:
- pmf_analyze_smd() - Analyze SMD force profiles
- pmf_analyze_pmf() - Analyze final PMF results

Author: PRISM Team (Step-by-Step API)
Version: 5.1 (Optimized)
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


def setup(config: Dict, output_dir: str, modeling_task_dir: str = None):
    """
    Setup global PMF configuration with streamlined validation
    
    Parameters:
    -----------
    config : Dict
        PMF configuration dictionary
    output_dir : str
        Output directory for PMF results
    modeling_task_dir : str, optional
        Directory containing modeling task results
        If not specified, defaults to "gaff_model"
    
    Example:
    --------
    >>> import pmf
    >>> config = {
    ...     'reference_group': 'Protein',
    ...     'moving_group': 'LIG',
    ...     'smd': {'pull_rate': 0.005, 'pull_k': 1000.0}
    ... }
    >>> pmf.setup(config, "./GMX_PROLIG_PMF", "./gaff_model")
    """
    global _CONFIG, _OUTPUT_DIR, _MD_RESULTS_DIR
    
    _CONFIG = config
    _OUTPUT_DIR = Path(output_dir).resolve()
    
    # Set default modeling_task_dir if not provided
    if modeling_task_dir is None:
        modeling_task_dir = "gaff_model"
        logger.info(f"No modeling_task_dir specified, using default: {modeling_task_dir}")
    
    _MD_RESULTS_DIR = Path(modeling_task_dir).resolve()
    
    # Create base output directory
    _OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Basic validation - only check if directory exists
    if not _MD_RESULTS_DIR.exists():
        raise FileNotFoundError(f"Modeling task directory not found: {_MD_RESULTS_DIR}")
    
    logger.info("PMF system configured successfully")
    logger.info(f"Output directory: {_OUTPUT_DIR}")
    logger.info(f"Modeling task directory: {_MD_RESULTS_DIR}")


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


def pmf_step1_smd_preparation():
    """
    STEP 1: SMD Preparation
    
    Prepares the system for SMD simulation:
    - Copies necessary files from modeling task results
    - Creates index files automatically
    - Generates SMD MDP file
    - Creates shell script to run SMD simulation
    - Creates analysis script for SMD results
    
    Returns:
    --------
    Dict : Preparation results with paths to generated files
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output", "./gaff_model")
    >>> results = pmf_step1_smd_preparation()
    >>> print(f"Run SMD using: {results['run_script']}")
    """
    _check_setup()
    logger.info("=== STEP 1: SMD Preparation ===")
    
    # Create SMD directory
    smd_dir = _OUTPUT_DIR / "smd"
    smd_dir.mkdir(exist_ok=True)
    
    # Copy modeling task results
    logger.info("Copying modeling task results...")
    _copy_modeling_results(smd_dir)
    
    # Create index file automatically
    logger.info("Creating index file...")
    _create_index_file(smd_dir)
    
    # Generate SMD MDP file
    logger.info("Generating SMD MDP file...")
    smd_mdp_file = _generate_smd_mdp(smd_dir)
    
    # Generate SMD run script
    logger.info("Generating SMD run script...")
    run_script = _generate_smd_run_script(smd_dir)
    
    # Generate SMD analysis script
    logger.info("Generating SMD analysis script...")
    analysis_script = _generate_smd_analysis_script(smd_dir)
    
    results = {
        'smd_dir': str(smd_dir),
        'smd_mdp': str(smd_mdp_file),
        'run_script': str(run_script),
        'analysis_script': str(analysis_script),
        'status': 'ready_for_smd',
        'next_step': 'Run SMD simulation using the generated script, then call pmf_step2_umbrella_preparation()'
    }
    
    logger.info("✓ STEP 1 completed successfully!")
    logger.info(f"Next: Run SMD simulation using: {run_script}")
    
    return results


def _copy_modeling_results(smd_dir):
    """Copy required files from modeling task directory"""
    # Define required files with their relative paths
    required_files = {
        'gro': {
            'source': 'GMX_PROLIG_MD/prod/md.gro',
            'dest': 'md.gro',
            'description': 'Final MD structure file'
        },
        'top': {
            'source': 'GMX_PROLIG_MD/topol.top',
            'dest': 'topol.top',
            'description': 'Topology file'
        },
        'posre': {
            'source': 'GMX_PROLIG_MD/posre.itp',
            'dest': 'posre.itp',
            'description': 'Position restraint file'
        }
    }
    
    # Copy individual files
    for file_type, file_info in required_files.items():
        source_file = _MD_RESULTS_DIR / file_info['source']
        dest_file = smd_dir / file_info['dest']
        
        if not source_file.exists():
            raise FileNotFoundError(f"Required file not found: {source_file}")
        
        shutil.copy2(source_file, dest_file)
        logger.info(f"✓ Copied {file_type}: {file_info['source']} -> {file_info['dest']}")
    
    # Copy LIG.amb2gmx directory
    lig_amb2gmx_source = _MD_RESULTS_DIR / "LIG.amb2gmx"
    lig_amb2gmx_dest = _OUTPUT_DIR / "LIG.amb2gmx"
    
    if not lig_amb2gmx_source.exists():
        raise FileNotFoundError(f"Required directory not found: {lig_amb2gmx_source}")
    
    # Remove destination directory if it exists, then copy
    if lig_amb2gmx_dest.exists():
        shutil.rmtree(lig_amb2gmx_dest)
    
    shutil.copytree(lig_amb2gmx_source, lig_amb2gmx_dest)
    file_count = len(list(lig_amb2gmx_dest.glob('*')))
    logger.info(f"✓ Copied LIG.amb2gmx directory: {file_count} files")


def _create_index_file(smd_dir):
    """Create index file for group selections"""
    moving_group = _get_config('moving_group', default='LIG')
    index_file = smd_dir / "index.ndx"
    
    if index_file.exists():
        logger.info("✓ Index file already exists")
        return
    
    try:
        make_ndx_cmd = [
            "gmx", "make_ndx",
            "-f", "md.gro",
            "-o", "index.ndx"
        ]
        
        process = subprocess.Popen(
            make_ndx_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(smd_dir)
        )
        
        # Create the moving group and quit
        input_commands = f"r {moving_group}\nq\n"
        stdout, stderr = process.communicate(input=input_commands)
        
        if process.returncode == 0:
            logger.info("✓ Index file created successfully")
        else:
            logger.warning(f"Index creation warning: {stderr}")
            # Create basic index file
            basic_cmd = ["gmx", "make_ndx", "-f", "md.gro", "-o", "index.ndx"]
            subprocess.run(basic_cmd, input="q\n", text=True, cwd=str(smd_dir), check=True)
            logger.info("✓ Basic index file created")
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create index file: {e}")
        raise RuntimeError(f"Index file creation failed: {e}")


def _generate_smd_mdp(smd_dir):
    """Generate SMD MDP file"""
    # Get configuration
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    pull_k = _get_config('smd', 'pull_k', 1000.0)
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
define              = -DPOSRE
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
pull_coord1_dim     = Y Y Y
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = {pull_rate}
pull_coord1_k       = {pull_k}
pull-pbc-ref-prev-step-com = yes
"""
    
    smd_mdp_file = smd_dir / "smd.mdp"
    with open(smd_mdp_file, 'w') as f:
        f.write(smd_content)
    
    logger.info(f"✓ SMD MDP file created: smd.mdp")
    return smd_mdp_file


def _generate_smd_run_script(smd_dir):
    """Generate shell script to run SMD simulation"""
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    
    distance_range = distance_end - distance_start
    smd_time_ps = distance_range / pull_rate
    
    script_content = f"""#!/bin/bash
######################################################
# SMD SIMULATION SCRIPT
# Distance: {distance_start:.2f} -> {distance_end:.2f} nm
# Pull rate: {pull_rate} nm/ps
# Simulation time: {smd_time_ps:.1f} ps
######################################################

set -e
echo "=== SMD Simulation ==="

# Create necessary directories
mkdir -p results
mkdir -p analysis

echo "Step 1: Generating TPR file..."
gmx grompp -f smd.mdp -c md.gro -n index.ndx -p topol.top -r md.gro -o smd.tpr -maxwarn 10

echo "Step 2: Running SMD simulation..."
echo "This may take a while depending on your system size and simulation time..."
gmx mdrun -s smd.tpr -deffnm smd -ntmpi 1 -ntomp 10 -v

echo "Step 3: Moving results to results directory..."
mv smd.* results/

echo "Step 4: Basic result check..."
if [ -f "results/smd.gro" ] && [ -f "results/smd.xtc" ] && [ -f "results/smd_pullf.xvg" ]; then
    echo "✓ SMD simulation completed successfully!"
    echo "✓ Generated files:"
    echo "  - results/smd.gro (final structure)"
    echo "  - results/smd.xtc (trajectory)"
    echo "  - results/smd_pullf.xvg (pull force)"
    echo "  - results/smd_pullx.xvg (pull distance)"
    echo ""
    echo "Next steps:"
    echo "1. Run analysis: bash analyze_smd.sh"
    echo "2. Then call: pmf_step2_umbrella_preparation()"
else
    echo "✗ SMD simulation failed or incomplete!"
    echo "Check the mdrun output for errors."
    exit 1
fi
"""
    
    run_script = smd_dir / "run_smd.sh"
    with open(run_script, 'w') as f:
        f.write(script_content)
    
    os.chmod(run_script, 0o755)
    logger.info(f"✓ SMD run script created: run_smd.sh")
    
    return run_script


def _generate_smd_analysis_script(smd_dir):
    """Generate script to analyze SMD results"""
    script_content = """#!/bin/bash
######################################################
# SMD ANALYSIS SCRIPT
######################################################

set -e
echo "=== SMD Analysis ==="

# Check if results exist
if [ ! -f "results/smd_pullf.xvg" ] || [ ! -f "results/smd_pullx.xvg" ]; then
    echo "✗ SMD results not found! Run SMD simulation first."
    exit 1
fi

echo "Analyzing SMD results..."

# Copy analysis files to analysis directory
cp results/smd_pullf.xvg analysis/
cp results/smd_pullx.xvg analysis/

echo "✓ SMD analysis files prepared"
echo "✓ Files available for plotting:"
echo "  - analysis/smd_pullf.xvg (force vs time)"
echo "  - analysis/smd_pullx.xvg (distance vs time)"
echo ""
echo "Use pmf_analyze_smd() in Python to generate plots"
"""
    
    analysis_script = smd_dir / "analyze_smd.sh"
    with open(analysis_script, 'w') as f:
        f.write(script_content)
    
    os.chmod(analysis_script, 0o755)
    logger.info(f"✓ SMD analysis script created: analyze_smd.sh")
    
    return analysis_script


def pmf_step2_umbrella_preparation():
    """
    STEP 2: Umbrella Sampling Preparation
    
    Analyzes SMD results and prepares umbrella sampling:
    - Extracts trajectory frames
    - Calculates distances and generates windows
    - Creates umbrella MDP files
    - Generates shell scripts to run umbrella sampling
    
    Returns:
    --------
    Dict : Umbrella preparation results
    
    Example:
    --------
    >>> import pmf
    >>> results = pmf_step2_umbrella_preparation()
    >>> print(f"Generated {results['n_windows']} windows")
    >>> print(f"Run umbrella sampling: {results['run_script']}")
    """
    _check_setup()
    logger.info("=== STEP 2: Umbrella Sampling Preparation ===")
    
    smd_dir = _OUTPUT_DIR / "smd"
    
    # Check SMD completion
    logger.info("Checking SMD results...")
    _check_smd_completion(smd_dir)
    
    # Extract trajectory frames
    logger.info("Extracting trajectory frames...")
    frames_info = _extract_trajectory_frames(smd_dir)
    
    # Calculate distances and generate windows
    logger.info("Generating umbrella windows...")
    windows = _generate_umbrella_windows(smd_dir, frames_info)
    
    # Create umbrella directory structure and MDP files
    logger.info("Setting up umbrella directories...")
    umbrella_dir = _setup_umbrella_directories(windows)
    
    # Generate umbrella run scripts
    logger.info("Generating umbrella run scripts...")
    run_scripts = _generate_umbrella_run_scripts(umbrella_dir, windows)
    
    # Clean up temporary trajectory frames
    logger.info("Cleaning up temporary trajectory frames...")
    frames_dir = smd_dir / "trajectory_frames"
    if frames_dir.exists():
        shutil.rmtree(frames_dir)
        logger.info("✓ Temporary frames cleaned up")
    
    results = {
        'umbrella_dir': str(umbrella_dir),
        'windows': [{'id': w['window_id'], 'distance': w['distance']} for w in windows],
        'n_windows': len(windows),
        'run_script': str(run_scripts['master_script']),
        'individual_scripts': run_scripts['individual_scripts'],
        'status': 'ready_for_umbrella',
        'next_step': 'Run umbrella sampling using the generated scripts, then call pmf_step3_wham_analysis()'
    }
    
    logger.info("✓ STEP 2 completed successfully!")
    logger.info(f"Generated {len(windows)} umbrella windows")
    logger.info(f"Run umbrella sampling using: {run_scripts['master_script']}")
    
    return results


def _check_smd_completion(smd_dir):
    """Check that SMD simulation has been completed"""
    required_files = [
        smd_dir / "results" / "smd.gro",
        smd_dir / "results" / "smd.xtc",
        smd_dir / "results" / "smd_pullf.xvg",
        smd_dir / "results" / "smd_pullx.xvg",
        smd_dir / "results" / "smd.tpr"
    ]
    
    missing_files = [str(f) for f in required_files if not f.exists()]
    
    if missing_files:
        raise FileNotFoundError(
            f"SMD simulation incomplete. Missing files:\n" + 
            "\n".join(f"  - {f}" for f in missing_files) +
            f"\n\nPlease run the SMD simulation first using: {smd_dir}/run_smd.sh"
        )
    
    logger.info("✓ SMD results validated")


def _extract_trajectory_frames(smd_dir, n_frames=300):
    """Extract frames from SMD trajectory"""
    logger.info(f"Extracting {n_frames} trajectory frames...")
    
    frames_dir = smd_dir / "trajectory_frames"
    frames_dir.mkdir(exist_ok=True)
    
    # Calculate frame extraction parameters
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    
    distance_range = distance_end - distance_start
    smd_time_ps = distance_range / pull_rate
    frame_interval_ps = smd_time_ps / n_frames
    
    # Extract frames using trjconv
    trjconv_cmd = [
        "gmx", "trjconv",
        "-f", str(smd_dir / "results" / "smd.xtc"),
        "-s", str(smd_dir / "results" / "smd.tpr"),
        "-o", str(frames_dir / "frame.gro"),
        "-sep",
        "-dt", str(frame_interval_ps)
    ]
    
    try:
        process = subprocess.Popen(
            trjconv_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(smd_dir)
        )
        
        stdout, stderr = process.communicate(input="0\n")  # Select System
        
        if process.returncode != 0:
            raise RuntimeError(f"Frame extraction failed: {stderr}")
        
        # Count extracted frames
        extracted_frames = len(list(frames_dir.glob("frame*.gro")))
        logger.info(f"✓ Extracted {extracted_frames} trajectory frames")
        
        return {
            'frames_dir': str(frames_dir),
            'n_frames': extracted_frames,
            'frame_interval_ps': frame_interval_ps
        }
        
    except Exception as e:
        logger.error(f"Frame extraction failed: {e}")
        raise


def _generate_umbrella_windows(smd_dir, frames_info):
    """Generate umbrella sampling windows based on trajectory frames"""
    frames_dir = Path(frames_info['frames_dir'])
    
    # Calculate distances for all frames
    logger.info("Calculating distances for frames...")
    distance_data = _calculate_frame_distances(smd_dir, frames_dir)
    
    # Generate adaptive windows
    sample_interval_near = _get_config('umbrella', 'sample_interval_near', 0.1)
    sample_interval_far = _get_config('umbrella', 'sample_interval_far', 0.2)
    cutoff_distance = _get_config('umbrella', 'cutoff_distance', 1.5)
    
    # Adaptive sampling logic
    windows = []
    selected_indices = []
    
    if not distance_data:
        raise RuntimeError("No distance data available for window generation")
    
    current_idx = 0
    selected_indices.append(current_idx)
    
    distances = [d[1] for d in distance_data]
    
    while current_idx < len(distances) - 1:
        current_distance = distances[current_idx]
        
        # Choose interval based on distance
        target_interval = sample_interval_near if current_distance < cutoff_distance else sample_interval_far
        target_distance = current_distance + target_interval
        
        # Find next frame closest to target distance
        best_idx = current_idx
        best_diff = float('inf')
        
        for idx in range(current_idx + 1, len(distances)):
            diff = abs(distances[idx] - target_distance)
            if diff < best_diff:
                best_diff = diff
                best_idx = idx
        
        if best_idx > current_idx:
            selected_indices.append(best_idx)
            current_idx = best_idx
        else:
            break
    
    # Create window information
    for i, idx in enumerate(selected_indices):
        frame, distance = distance_data[idx]
        windows.append({
            'window_id': i,
            'frame': frame,
            'distance': distance,
            'frame_file': frames_dir / f"frame{frame}.gro"
        })
    
    logger.info(f"✓ Generated {len(windows)} umbrella windows")
    return windows


def _calculate_frame_distances(smd_dir, frames_dir):
    """Calculate distances for trajectory frames"""
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    
    frame_files = sorted(frames_dir.glob("frame*.gro"))
    distance_data = []
    
    temp_dir = smd_dir / "temp_distance"
    temp_dir.mkdir(exist_ok=True)
    
    try:
        for frame_file in frame_files:
            # Extract frame number from filename
            frame_num = int(frame_file.stem.replace('frame', ''))
            
            # Calculate distance using gmx distance
            distance_cmd = [
                "gmx", "distance",
                "-s", str(smd_dir / "results" / "smd.tpr"),
                "-f", str(frame_file),
                "-n", str(smd_dir / "index.ndx"),
                "-select", f"com of group {reference_group} plus com of group {moving_group}",
                "-oall", str(temp_dir / "dist_temp.xvg")
            ]
            
            try:
                subprocess.run(distance_cmd, check=True, capture_output=True, text=True)
                
                # Read distance value
                dist_file = temp_dir / "dist_temp.xvg"
                with open(dist_file, 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        if not line.startswith('#') and not line.startswith('@') and line.strip():
                            distance = float(line.split()[1])
                            distance_data.append((frame_num, distance))
                            break
                
            except subprocess.CalledProcessError as e:
                logger.warning(f"Could not calculate distance for frame {frame_num}: {e}")
                continue
    
    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    logger.info(f"✓ Calculated distances for {len(distance_data)} frames")
    return distance_data


def _setup_umbrella_directories(windows):
    """Setup umbrella sampling directory structure"""
    umbrella_dir = _OUTPUT_DIR / "umbrella"
    umbrella_dir.mkdir(exist_ok=True)
    
    # Generate umbrella MDP template
    umbrella_mdp_content = _generate_umbrella_mdp_content()
    
    for window in windows:
        window_dir = umbrella_dir / f"window_{window['window_id']:03d}"
        window_dir.mkdir(exist_ok=True)
        
        # Create window-specific MDP file
        window_mdp_content = umbrella_mdp_content.replace(
            "DISTANCE_PLACEHOLDER", f"{window['distance']:.6f}"
        ).replace(
            "WINDOW_TITLE", f"Umbrella sampling at {window['distance']:.3f} nm (window {window['window_id']:03d})"
        )
        
        window_mdp = window_dir / "umbrella.mdp"
        with open(window_mdp, 'w') as f:
            f.write(window_mdp_content)
        
        # Copy starting structure
        if window['frame_file'].exists():
            dest_file = window_dir / "start.gro"
            shutil.copy2(window['frame_file'], dest_file)
        
        # Update window info with directory
        window['directory'] = window_dir
    
    logger.info(f"✓ Setup {len(windows)} umbrella directories")
    return umbrella_dir


def _generate_umbrella_mdp_content():
    """Generate umbrella sampling MDP template content"""
    dt = _get_config('simulation', 'dt', 0.002)
    production_time = _get_config('umbrella', 'production_time_ps', 22000)
    sampling_interval = _get_config('umbrella', 'sampling_interval_ps', 10)
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    pull_k = _get_config('smd', 'pull_k', 1000.0)
    
    nsteps = int(production_time / dt)
    output_interval = int(sampling_interval / dt)
    
    return f"""; WINDOW_TITLE
title               = WINDOW_TITLE
define              = -DPOSRE
integrator          = md
nsteps              = {nsteps}
dt                  = {dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}
{_create_common_mdp_sections()}

; Pull code for umbrella sampling
pull                = yes
pull_ncoords        = 1
pull_ngroups        = 2
pull_group1_name    = {reference_group}
pull_group2_name    = {moving_group}
pull_coord1_type    = umbrella
pull_coord1_geometry = distance
pull_coord1_dim     = Y Y Y
pull_coord1_groups  = 1 2
pull_coord1_start   = yes
pull_coord1_rate    = 0.0
pull_coord1_k       = {pull_k}
pull_coord1_r0      = DISTANCE_PLACEHOLDER
pull-pbc-ref-prev-step-com = yes
"""


def _generate_umbrella_run_scripts(umbrella_dir, windows):
    """Generate shell scripts to run umbrella sampling"""
    smd_dir = _OUTPUT_DIR / "smd"
    
    # Copy shared files to umbrella directory
    shared_files = ['topol.top', 'index.ndx']
    for filename in shared_files:
        source_file = smd_dir / filename
        dest_file = umbrella_dir / filename
        if source_file.exists():
            shutil.copy2(source_file, dest_file)
    
    # Copy LIG.amb2gmx directory
    lig_source = _OUTPUT_DIR / "LIG.amb2gmx"
    lig_dest = umbrella_dir / "LIG.amb2gmx"
    if lig_source.exists() and not lig_dest.exists():
        shutil.copytree(lig_source, lig_dest)
    
    # Individual window scripts
    individual_scripts = []
    for window in windows:
        window_dir = window['directory']
        
        script_content = f"""#!/bin/bash
######################################################
# UMBRELLA WINDOW {window['window_id']:03d} SCRIPT
# Distance: {window['distance']:.3f} nm
######################################################

set -e
echo "=== Running Umbrella Window {window['window_id']:03d} ==="
echo "Target distance: {window['distance']:.3f} nm"

# Create results directory
mkdir -p results

echo "Step 1: Generating TPR file..."
gmx grompp -f umbrella.mdp -c start.gro -n ../index.ndx -p ../topol.top -o umbrella.tpr -maxwarn 10

echo "Step 2: Running umbrella sampling..."
echo "This will take approximately 22 ns simulation time..."
gmx mdrun -s umbrella.tpr -deffnm umbrella -ntmpi 1 -ntomp 10 -v

echo "Step 3: Moving results..."
mv umbrella.* results/

echo "Step 4: Validating results..."
if [ -f "results/umbrella.gro" ] && [ -f "results/umbrella_pullf.xvg" ]; then
    echo "✓ Window {window['window_id']:03d} completed successfully!"
else
    echo "✗ Window {window['window_id']:03d} failed!"
    exit 1
fi
"""

        script_file = window_dir / "run_window.sh"
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_file, 0o755)
        individual_scripts.append(str(script_file))
    
    # Master run script
    master_script_content = f"""#!/bin/bash
######################################################
# MASTER UMBRELLA SAMPLING SCRIPT
# Total windows: {len(windows)}
######################################################

set -e
echo "=== Running All Umbrella Sampling Windows ==="
echo "Total windows: {len(windows)}"
echo "Estimated total time: Very long! Consider running in parallel."
echo ""

# Function to run a single window
run_window() {{
    local window_id=$1
    local window_dir="window_$(printf "%03d" $window_id)"
    
    echo "Starting window $window_id..."
    cd "$window_dir"
    
    if bash run_window.sh; then
        echo "✓ Window $window_id completed"
        cd ..
        return 0
    else
        echo "✗ Window $window_id failed"
        cd ..
        return 1
    fi
}}

# Option 1: Sequential execution (very slow)
if [ "$1" = "sequential" ]; then
    echo "Running windows sequentially..."
    failed_windows=()
    
    for i in {{0..{len(windows)-1}}}; do
        if ! run_window $i; then
            failed_windows+=($i)
        fi
    done
    
    if [ ${{#failed_windows[@]}} -eq 0 ]; then
        echo "✓ All windows completed successfully!"
    else
        echo "✗ Failed windows: ${{failed_windows[*]}}"
        exit 1
    fi

# Option 2: Parallel execution (recommended)
elif [ "$1" = "parallel" ]; then
    echo "Running windows in parallel..."
    echo "Make sure your system can handle multiple GROMACS jobs!"
    
    # Run windows in parallel (adjust -P based on your system)
    seq 0 {len(windows)-1} | xargs -n 1 -P 4 -I {{}} bash -c 'run_window {{}}'
    
    echo "✓ Parallel execution completed!"

# Option 3: Generate job scripts for cluster
elif [ "$1" = "cluster" ]; then
    echo "Generating cluster job scripts..."
    
    for i in {{0..{len(windows)-1}}}; do
        cat > "job_window_$(printf "%03d" $i).sh" << EOF
#!/bin/bash
#SBATCH --job-name=umbrella_$i
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

cd window_$(printf "%03d" $i)
bash run_window.sh
EOF
        echo "Created job_window_$(printf "%03d" $i).sh"
    done
    
    echo "Submit jobs using: sbatch job_window_*.sh"

else
    echo "Usage: $0 [sequential|parallel|cluster]"
    echo ""
    echo "Options:"
    echo "  sequential - Run windows one by one (very slow)"
    echo "  parallel   - Run multiple windows simultaneously"
    echo "  cluster    - Generate SLURM job scripts"
    echo ""
    echo "Alternatively, run individual windows manually:"
"""

    for window in windows:
        master_script_content += f"""    echo "cd window_{window['window_id']:03d} && bash run_window.sh"
"""

    master_script_content += """
fi
"""

    master_script = umbrella_dir / "run_all_umbrella.sh"
    with open(master_script, 'w') as f:
        f.write(master_script_content)
    
    os.chmod(master_script, 0o755)
    
    logger.info("✓ Umbrella run scripts generated")
    
    return {
        'master_script': master_script,
        'individual_scripts': individual_scripts
    }


def pmf_step3_wham_analysis():
    """
    STEP 3: WHAM Analysis and PMF Calculation
    
    Performs WHAM analysis on umbrella sampling results:
    - Validates umbrella sampling completion
    - Runs WHAM analysis with bootstrap error estimation
    - Calculates binding energy
    - Generates PMF plots and comprehensive report
    
    Returns:
    --------
    Dict : WHAM analysis results
    
    Example:
    --------
    >>> import pmf
    >>> results = pmf_step3_wham_analysis()
    >>> print(f"Binding energy: {results['binding_energy']['value']:.2f} kcal/mol")
    """
    _check_setup()
    logger.info("=== STEP 3: WHAM Analysis ===")
    
    umbrella_dir = _OUTPUT_DIR / "umbrella"
    
    # Check umbrella sampling completion
    logger.info("Checking umbrella sampling results...")
    valid_windows = _check_umbrella_completion(umbrella_dir)
    
    # Create analysis directory
    analysis_dir = _OUTPUT_DIR / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    
    # Run WHAM analysis
    logger.info("Running WHAM analysis...")
    wham_results = _run_wham_analysis(valid_windows, analysis_dir)
    
    # Generate visualizations
    logger.info("Generating PMF visualizations...")
    plots = _generate_pmf_plots(wham_results, analysis_dir)
    
    # Calculate binding energy
    logger.info("Calculating binding energy...")
    binding_energy = _calculate_binding_energy(wham_results)
    
    # Generate comprehensive report
    logger.info("Generating analysis report...")
    report = _generate_final_report(wham_results, binding_energy, plots, analysis_dir, len(valid_windows))
    
    results = {
        'analysis_dir': str(analysis_dir),
        'wham_results': wham_results,
        'binding_energy': binding_energy,
        'plots': plots,
        'report': str(report),
        'n_windows_analyzed': len(valid_windows),
        'status': 'pmf_completed'
    }
    
    logger.info("✓ STEP 3 completed successfully!")
    logger.info(f"Binding energy: {binding_energy['value']:.2f} {binding_energy['unit']}")
    logger.info(f"Full report: {report}")
    
    return results


def _check_umbrella_completion(umbrella_dir):
    """Check umbrella sampling completion and return valid windows"""
    if not umbrella_dir.exists():
        raise FileNotFoundError(f"Umbrella directory not found: {umbrella_dir}")
    
    window_dirs = sorted(umbrella_dir.glob("window_*"))
    valid_windows = []
    
    for window_dir in window_dirs:
        tpr_file = window_dir / "results" / "umbrella.tpr"
        pullf_file = window_dir / "results" / "umbrella_pullf.xvg"
        
        if tpr_file.exists() and pullf_file.exists():
            valid_windows.append({
                'window_dir': window_dir,
                'tpr_file': tpr_file,
                'pullf_file': pullf_file,
                'window_id': int(window_dir.name.split('_')[1])
            })
        else:
            logger.warning(f"Incomplete umbrella window: {window_dir}")
    
    if len(valid_windows) < 5:
        raise RuntimeError(
            f"Insufficient completed umbrella windows: {len(valid_windows)} (need ≥5)\n"
            f"Please ensure umbrella sampling has been completed for more windows."
        )
    
    logger.info(f"✓ Found {len(valid_windows)} completed umbrella windows")
    return valid_windows


def _run_wham_analysis(valid_windows, analysis_dir):
    """Run WHAM analysis on umbrella sampling results"""
    wham_dir = analysis_dir / "wham"
    wham_dir.mkdir(exist_ok=True)
    
    # Get analysis configuration
    begin_time = _get_config('analysis', 'begin_time_ps', 2000)
    bootstrap_iterations = _get_config('analysis', 'bootstrap_iterations', 50)
    energy_unit = _get_config('analysis', 'energy_unit', 'kCal')
    
    # Generate file lists for WHAM
    tpr_list = wham_dir / "tpr-files.dat"
    pullf_list = wham_dir / "pullf-files.dat"
    
    with open(tpr_list, 'w') as f:
        for window in sorted(valid_windows, key=lambda x: x['window_id']):
            # Use relative path from wham directory
            rel_path = os.path.relpath(window['tpr_file'], wham_dir)
            f.write(f"{rel_path}\n")
    
    with open(pullf_list, 'w') as f:
        for window in sorted(valid_windows, key=lambda x: x['window_id']):
            # Use relative path from wham directory
            rel_path = os.path.relpath(window['pullf_file'], wham_dir)
            f.write(f"{rel_path}\n")
    
    # Run WHAM analysis
    wham_cmd = [
        "gmx", "wham",
        "-it", "tpr-files.dat",
        "-if", "pullf-files.dat",
        "-b", str(begin_time),
        "-o", "pmf.xvg",
        "-unit", energy_unit,
        "-bsres", "pmferror.xvg",
        "-bsprof", "bsprofile.xvg",
        "-nBootstrap", str(bootstrap_iterations),
        "-bs-method", "b-hist",
        "-v"
    ]
    
    try:
        result = subprocess.run(
            wham_cmd, 
            cwd=wham_dir, 
            capture_output=True, 
            text=True, 
            check=True,
            timeout=3600  # 1 hour timeout
        )
        
        logger.info("✓ WHAM analysis completed successfully")
        
        return {
            'wham_dir': str(wham_dir),
            'pmf_file': str(wham_dir / "pmf.xvg"),
            'error_file': str(wham_dir / "pmferror.xvg"),
            'profile_file': str(wham_dir / "bsprofile.xvg"),
            'begin_time': begin_time,
            'bootstrap_iterations': bootstrap_iterations,
            'energy_unit': energy_unit
        }
        
    except subprocess.CalledProcessError as e:
        logger.error(f"WHAM analysis failed: {e.stderr}")
        raise RuntimeError(f"WHAM analysis failed. Check input files and parameters.")
    except subprocess.TimeoutExpired:
        logger.error("WHAM analysis timed out")
        raise RuntimeError("WHAM analysis timed out. Consider reducing bootstrap iterations.")


def _generate_pmf_plots(wham_results, analysis_dir):
    """Generate PMF visualizations and plots"""
    pmf_file = Path(wham_results['pmf_file'])
    error_file = Path(wham_results['error_file'])
    energy_unit = wham_results['energy_unit']
    
    # Read PMF data
    try:
        pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
        distances = pmf_data[:, 0]
        pmf_values = pmf_data[:, 1]
    except Exception as e:
        raise RuntimeError(f"Failed to read PMF data: {e}")
    
    # Read error data if available
    errors = None
    if error_file.exists():
        try:
            error_data = np.loadtxt(error_file, comments=['#', '@'])
            if error_data.shape[1] >= 2:
                errors = error_data[:, 1]
        except Exception as e:
            logger.warning(f"Could not read error data: {e}")
    
    # Create main PMF plot
    plt.figure(figsize=(12, 8))
    
    if errors is not None:
        plt.errorbar(distances, pmf_values, yerr=errors, 
                    fmt='o-', capsize=3, linewidth=2, markersize=4,
                    color='blue', ecolor='lightblue', label='PMF with errors')
    else:
        plt.plot(distances, pmf_values, 'o-', linewidth=2, markersize=4,
                color='blue', label='PMF')
    
    plt.xlabel('Distance (nm)', fontsize=14)
    plt.ylabel(f'PMF ({energy_unit.lower()}/mol)', fontsize=14)
    plt.title('Potential of Mean Force', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Add statistical information
    min_pmf = np.min(pmf_values)
    max_pmf = np.max(pmf_values)
    binding_energy = max_pmf - min_pmf
    min_distance = distances[np.argmin(pmf_values)]
    
    stats_text = f'Binding Energy: {binding_energy:.2f} {energy_unit.lower()}/mol\n'
    stats_text += f'Minimum at: {min_distance:.3f} nm\n'
    stats_text += f'PMF range: {min_pmf:.2f} to {max_pmf:.2f}'
    
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save main plot
    pmf_plot = analysis_dir / "pmf_curve.png"
    plt.savefig(pmf_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create convergence plot if bootstrap data available
    convergence_plot = None
    profile_file = Path(wham_results['profile_file'])
    if profile_file.exists():
        try:
            convergence_plot = _create_convergence_plot(profile_file, analysis_dir, energy_unit)
        except Exception as e:
            logger.warning(f"Could not create convergence plot: {e}")
    
    logger.info(f"✓ PMF plots saved to: {analysis_dir}")
    
    return {
        'pmf_plot': str(pmf_plot),
        'convergence_plot': str(convergence_plot) if convergence_plot else None,
        'pmf_data': {
            'distances': distances.tolist(),
            'pmf_values': pmf_values.tolist(),
            'errors': errors.tolist() if errors is not None else None
        }
    }


def _create_convergence_plot(profile_file, analysis_dir, energy_unit):
    """Create bootstrap convergence plot"""
    try:
        profile_data = np.loadtxt(profile_file, comments=['#', '@'])
        
        plt.figure(figsize=(10, 6))
        
        # Plot bootstrap profiles (sample a few for clarity)
        n_profiles = profile_data.shape[1] - 1  # First column is distance
        distances = profile_data[:, 0]
        
        # Plot every 10th profile to avoid clutter
        for i in range(1, min(n_profiles, 50), 5):
            plt.plot(distances, profile_data[:, i], 'gray', alpha=0.3, linewidth=0.5)
        
        # Plot average (if available)
        if n_profiles > 1:
            avg_profile = np.mean(profile_data[:, 1:], axis=1)
            plt.plot(distances, avg_profile, 'red', linewidth=2, label='Average')
        
        plt.xlabel('Distance (nm)', fontsize=12)
        plt.ylabel(f'PMF ({energy_unit.lower()}/mol)', fontsize=12)
        plt.title('Bootstrap Convergence Analysis', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        convergence_plot = analysis_dir / "pmf_convergence.png"
        plt.savefig(convergence_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        return convergence_plot
        
    except Exception as e:
        logger.warning(f"Convergence plot creation failed: {e}")
        return None


def _calculate_binding_energy(wham_results):
    """Calculate binding energy from PMF data"""
    pmf_file = Path(wham_results['pmf_file'])
    energy_unit = wham_results['energy_unit']
    
    try:
        pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
        pmf_values = pmf_data[:, 1]
        distances = pmf_data[:, 0]
        
        min_pmf = np.min(pmf_values)
        max_pmf = np.max(pmf_values)
        binding_energy = max_pmf - min_pmf
        
        min_distance = distances[np.argmin(pmf_values)]
        max_distance = distances[np.argmax(pmf_values)]
        
        logger.info(f"✓ Binding energy: {binding_energy:.2f} {energy_unit.lower()}/mol")
        
        return {
            'value': binding_energy,
            'unit': f'{energy_unit.lower()}/mol',
            'min_pmf': min_pmf,
            'max_pmf': max_pmf,
            'min_distance': min_distance,
            'max_distance': max_distance,
            'pmf_range': max_pmf - min_pmf
        }
        
    except Exception as e:
        logger.error(f"Binding energy calculation failed: {e}")
        raise RuntimeError(f"Failed to calculate binding energy: {e}")


def _generate_final_report(wham_results, binding_energy, plots, analysis_dir, n_windows):
    """Generate comprehensive final analysis report"""
    from datetime import datetime
    
    report_file = analysis_dir / "pmf_analysis_report.txt"
    
    # Get configuration for report
    reference_group = _get_config('reference_group', default='Protein')
    moving_group = _get_config('moving_group', default='LIG')
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    
    report_content = f"""
PMF Analysis Report
==================

Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Analysis Directory: {analysis_dir}

SYSTEM INFORMATION
-----------------
Reference Group: {reference_group}
Moving Group: {moving_group}
Distance Range: {distance_start:.2f} - {distance_end:.2f} nm
SMD Pull Rate: {pull_rate} nm/ps

SIMULATION SUMMARY
-----------------
Number of Umbrella Windows: {n_windows}
WHAM Begin Time: {wham_results['begin_time']} ps
Bootstrap Iterations: {wham_results['bootstrap_iterations']}
Energy Unit: {wham_results['energy_unit']}

BINDING ENERGY RESULTS
---------------------
Binding Energy: {binding_energy['value']:.2f} {binding_energy['unit']}
Minimum PMF: {binding_energy['min_pmf']:.2f} {binding_energy['unit']}
Maximum PMF: {binding_energy['max_pmf']:.2f} {binding_energy['unit']}
PMF Range: {binding_energy['pmf_range']:.2f} {binding_energy['unit']}

Distance at Minimum: {binding_energy['min_distance']:.3f} nm
Distance at Maximum: {binding_energy['max_distance']:.3f} nm

OUTPUT FILES
-----------
PMF Data: {wham_results['pmf_file']}
Error Data: {wham_results['error_file']}
Bootstrap Profile: {wham_results['profile_file']}
PMF Plot: {plots['pmf_plot']}"""

    if plots['convergence_plot']:
        report_content += f"""
Convergence Plot: {plots['convergence_plot']}"""

    report_content += f"""

INTERPRETATION
-------------
The binding energy of {binding_energy['value']:.2f} {binding_energy['unit']} represents
the free energy difference between the bound and unbound states.

- Negative binding energy indicates favorable binding
- Positive binding energy indicates unfavorable binding
- The minimum at {binding_energy['min_distance']:.3f} nm represents the most stable configuration

QUALITY ASSESSMENT
-----------------
Number of windows used: {n_windows}
Window coverage: {'Good' if n_windows >= 20 else 'Adequate' if n_windows >= 10 else 'Limited'}
Bootstrap analysis: {'Completed' if wham_results['bootstrap_iterations'] > 0 else 'Not performed'}

RECOMMENDATIONS
--------------
"""

    if n_windows < 10:
        report_content += "- Consider using more umbrella windows for better accuracy\n"
    if binding_energy['value'] > 0:
        report_content += "- Positive binding energy suggests weak or unfavorable binding\n"
    if binding_energy['value'] > 10:
        report_content += "- Very high binding energy may indicate sampling issues\n"

    report_content += """
PMF analysis completed successfully!
For questions about interpretation, consult molecular dynamics literature.
"""

    with open(report_file, 'w') as f:
        f.write(report_content)
    
    logger.info(f"✓ Final report generated: {report_file}")
    return report_file


def pmf_analyze_smd():
    """
    Analyze SMD results and generate force/distance plots
    
    Returns:
    --------
    Dict : SMD analysis results with plot paths
    
    Example:
    --------
    >>> import pmf
    >>> results = pmf_analyze_smd()
    >>> print(f"Force plot: {results['force_plot']}")
    """
    _check_setup()
    logger.info("=== SMD Analysis ===")
    
    smd_dir = _OUTPUT_DIR / "smd"
    
    # Check SMD results exist
    pullf_file = smd_dir / "results" / "smd_pullf.xvg"
    pullx_file = smd_dir / "results" / "smd_pullx.xvg"
    
    if not pullf_file.exists() or not pullx_file.exists():
        raise FileNotFoundError(
            f"SMD results not found. Expected files:\n"
            f"  - {pullf_file}\n"
            f"  - {pullx_file}\n"
            f"Run SMD simulation first using the generated scripts."
        )
    
    analysis_dir = smd_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    
    plots = {}
    
    # Plot force profile
    logger.info("Analyzing force profile...")
    try:
        force_data = np.loadtxt(pullf_file, comments=['#', '@'])
        time = force_data[:, 0]
        force = force_data[:, 1]
        
        plt.figure(figsize=(12, 8))
        plt.plot(time, force, 'b-', linewidth=1.5)
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('Force (kJ/mol/nm)', fontsize=14)
        plt.title('SMD Pull Force vs Time', fontsize=16)
        plt.grid(True, alpha=0.3)
        
        # Add statistics
        max_force = np.max(force)
        avg_force = np.mean(force)
        final_force = force[-1]
        
        stats_text = f'Max Force: {max_force:.2f} kJ/mol/nm\n'
        stats_text += f'Avg Force: {avg_force:.2f} kJ/mol/nm\n'
        stats_text += f'Final Force: {final_force:.2f} kJ/mol/nm'
        
        plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        force_plot = analysis_dir / "smd_force_profile.png"
        plt.savefig(force_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['force_plot'] = str(force_plot)
        logger.info(f"✓ Force profile plot: {force_plot}")
        
    except Exception as e:
        logger.warning(f"Could not create force profile plot: {e}")
    
    # Plot distance profile
    logger.info("Analyzing distance profile...")
    try:
        distance_data = np.loadtxt(pullx_file, comments=['#', '@'])
        time = distance_data[:, 0]
        distance = distance_data[:, 1]
        
        plt.figure(figsize=(12, 8))
        plt.plot(time, distance, 'r-', linewidth=1.5)
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('Distance (nm)', fontsize=14)
        plt.title('SMD Distance vs Time', fontsize=16)
        plt.grid(True, alpha=0.3)
        
        # Add statistics
        initial_dist = distance[0]
        final_dist = distance[-1]
        total_change = final_dist - initial_dist
        
        stats_text = f'Initial Distance: {initial_dist:.3f} nm\n'
        stats_text += f'Final Distance: {final_dist:.3f} nm\n'
        stats_text += f'Total Change: {total_change:.3f} nm'
        
        plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        plt.tight_layout()
        distance_plot = analysis_dir / "smd_distance_profile.png"
        plt.savefig(distance_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['distance_plot'] = str(distance_plot)
        logger.info(f"✓ Distance profile plot: {distance_plot}")
        
    except Exception as e:
        logger.warning(f"Could not create distance profile plot: {e}")
    
    # Create combined plot
    if 'force_plot' in plots and 'distance_plot' in plots:
        logger.info("Creating combined analysis plot...")
        try:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
            
            # Force subplot
            force_data = np.loadtxt(pullf_file, comments=['#', '@'])
            ax1.plot(force_data[:, 0], force_data[:, 1], 'b-', linewidth=1.5)
            ax1.set_ylabel('Force (kJ/mol/nm)', fontsize=12)
            ax1.set_title('SMD Analysis: Force and Distance Profiles', fontsize=14)
            ax1.grid(True, alpha=0.3)
            
            # Distance subplot
            distance_data = np.loadtxt(pullx_file, comments=['#', '@'])
            ax2.plot(distance_data[:, 0], distance_data[:, 1], 'r-', linewidth=1.5)
            ax2.set_xlabel('Time (ps)', fontsize=12)
            ax2.set_ylabel('Distance (nm)', fontsize=12)
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            combined_plot = analysis_dir / "smd_combined_analysis.png"
            plt.savefig(combined_plot, dpi=300, bbox_inches='tight')
            plt.close()
            
            plots['combined_plot'] = str(combined_plot)
            logger.info(f"✓ Combined analysis plot: {combined_plot}")
            
        except Exception as e:
            logger.warning(f"Could not create combined plot: {e}")
    
    results = {
        'analysis_dir': str(analysis_dir),
        'plots': plots,
        'status': 'smd_analysis_completed'
    }
    
    logger.info("✓ SMD analysis completed!")
    return results


def pmf_analyze_pmf():
    """
    Analyze final PMF results and generate detailed plots
    
    Returns:
    --------
    Dict : PMF analysis results
    
    Example:
    --------
    >>> import pmf
    >>> results = pmf_analyze_pmf()
    >>> print(f"PMF analysis: {results['plots']}")
    """
    _check_setup()
    logger.info("=== PMF Results Analysis ===")
    
    analysis_dir = _OUTPUT_DIR / "analysis"
    if not analysis_dir.exists():
        raise FileNotFoundError(
            f"PMF analysis not found. Run pmf_step3_wham_analysis() first."
        )
    
    wham_dir = analysis_dir / "wham"
    pmf_file = wham_dir / "pmf.xvg"
    error_file = wham_dir / "pmferror.xvg"
    
    if not pmf_file.exists():
        raise FileNotFoundError(f"PMF data file not found: {pmf_file}")
    
    # Read PMF data
    try:
        pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
        distances = pmf_data[:, 0]
        pmf_values = pmf_data[:, 1]
    except Exception as e:
        raise RuntimeError(f"Failed to read PMF data: {e}")
    
    # Read error data if available
    errors = None
    if error_file.exists():
        try:
            error_data = np.loadtxt(error_file, comments=['#', '@'])
            if error_data.shape[1] >= 2:
                errors = error_data[:, 1]
        except Exception as e:
            logger.warning(f"Could not read error data: {e}")
    
    plots_dir = analysis_dir / "detailed_plots"
    plots_dir.mkdir(exist_ok=True)
    
    plots = {}
    
    # Create detailed PMF plot with enhanced features
    logger.info("Creating detailed PMF plot...")
    try:
        fig, ax = plt.subplots(figsize=(14, 10))
        
        if errors is not None:
            ax.errorbar(distances, pmf_values, yerr=errors, 
                       fmt='o-', capsize=3, linewidth=2, markersize=6,
                       color='blue', ecolor='lightblue', label='PMF ± Bootstrap Error')
        else:
            ax.plot(distances, pmf_values, 'o-', linewidth=2, markersize=6,
                   color='blue', label='PMF')
        
        # Add key points
        min_idx = np.argmin(pmf_values)
        max_idx = np.argmax(pmf_values)
        
        ax.plot(distances[min_idx], pmf_values[min_idx], 'go', markersize=10, 
               label=f'Minimum ({distances[min_idx]:.3f} nm)')
        ax.plot(distances[max_idx], pmf_values[max_idx], 'ro', markersize=10, 
               label=f'Maximum ({distances[max_idx]:.3f} nm)')
        
        # Add binding energy annotation
        binding_energy = pmf_values[max_idx] - pmf_values[min_idx]
        ax.annotate(f'ΔG = {binding_energy:.2f} kcal/mol', 
                   xy=(distances[max_idx], pmf_values[max_idx]),
                   xytext=(distances[max_idx] + 0.2, pmf_values[max_idx] + 2),
                   arrowprops=dict(arrowstyle='->', color='red', lw=2),
                   fontsize=14, fontweight='bold')
        
        ax.set_xlabel('Distance (nm)', fontsize=16)
        ax.set_ylabel('PMF (kcal/mol)', fontsize=16)
        ax.set_title('Potential of Mean Force - Detailed Analysis', fontsize=18)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12)
        
        # Add statistical info box
        stats_text = f'Binding Energy: {binding_energy:.2f} kcal/mol\n'
        stats_text += f'Min PMF: {pmf_values[min_idx]:.2f} kcal/mol\n'
        stats_text += f'Max PMF: {pmf_values[max_idx]:.2f} kcal/mol\n'
        stats_text += f'Distance Range: {distances[0]:.2f} - {distances[-1]:.2f} nm'
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               fontsize=12, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        detailed_plot = plots_dir / "pmf_detailed.png"
        plt.savefig(detailed_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['detailed_plot'] = str(detailed_plot)
        logger.info(f"✓ Detailed PMF plot: {detailed_plot}")
        
    except Exception as e:
        logger.warning(f"Could not create detailed PMF plot: {e}")
    
    # Create energy landscape plot
    logger.info("Creating energy landscape plot...")
    try:
        plt.figure(figsize=(12, 8))
        
        # Fill area under curve
        plt.fill_between(distances, pmf_values, alpha=0.3, color='lightblue', label='Energy Surface')
        plt.plot(distances, pmf_values, 'b-', linewidth=3, label='PMF')
        
        # Mark important regions
        min_idx = np.argmin(pmf_values)
        threshold = pmf_values[min_idx] + 2.0  # 2 kcal/mol above minimum
        bound_region = distances[pmf_values <= threshold]
        if len(bound_region) > 0:
            plt.axvspan(bound_region[0], bound_region[-1], alpha=0.2, color='green', 
                       label='Favorable Binding Region (< 2 kcal/mol)')
        
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.ylabel('PMF (kcal/mol)', fontsize=14)
        plt.title('Energy Landscape for Protein-Ligand Binding', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        landscape_plot = plots_dir / "energy_landscape.png"
        plt.savefig(landscape_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['landscape_plot'] = str(landscape_plot)
        logger.info(f"✓ Energy landscape plot: {landscape_plot}")
        
    except Exception as e:
        logger.warning(f"Could not create energy landscape plot: {e}")
    
    # Create derivative plot (force profile)
    logger.info("Creating force profile from PMF...")
    try:
        # Calculate derivative (negative gradient = force)
        pmf_gradient = np.gradient(pmf_values, distances)
        force_profile = -pmf_gradient  # Force is negative gradient
        
        plt.figure(figsize=(12, 8))
        plt.plot(distances, force_profile, 'r-', linewidth=2, label='Force Profile')
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='Zero Force')
        
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.ylabel('Force (kcal/mol/nm)', fontsize=14)
        plt.title('Force Profile Derived from PMF', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        
        force_plot = plots_dir / "pmf_force_profile.png"
        plt.savefig(force_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots['force_profile'] = str(force_plot)
        logger.info(f"✓ Force profile plot: {force_plot}")
        
    except Exception as e:
        logger.warning(f"Could not create force profile plot: {e}")
    
    # Generate detailed analysis report
    logger.info("Generating detailed PMF analysis report...")
    detailed_report = _generate_detailed_pmf_report(
        distances, pmf_values, errors, plots_dir, plots
    )
    
    results = {
        'plots_dir': str(plots_dir),
        'plots': plots,
        'detailed_report': str(detailed_report),
        'pmf_summary': {
            'binding_energy': binding_energy,
            'min_distance': distances[min_idx],
            'max_distance': distances[max_idx],
            'distance_range': [distances[0], distances[-1]],
            'pmf_range': [pmf_values[min_idx], pmf_values[max_idx]]
        },
        'status': 'detailed_analysis_completed'
    }
    
    logger.info("✓ Detailed PMF analysis completed!")
    return results


def _generate_detailed_pmf_report(distances, pmf_values, errors, plots_dir, plots):
    """Generate detailed PMF analysis report"""
    from datetime import datetime
    
    report_file = plots_dir / "detailed_pmf_analysis.txt"
    
    # Calculate statistics
    min_idx = np.argmin(pmf_values)
    max_idx = np.argmax(pmf_values)
    binding_energy = pmf_values[max_idx] - pmf_values[min_idx]
    
    # Calculate additional statistics
    pmf_std = np.std(pmf_values)
    distance_span = distances[-1] - distances[0]
    
    # Estimate equilibrium constant (rough approximation)
    RT = 0.593  # kcal/mol at 298K
    if binding_energy < 0:
        keq = np.exp(-binding_energy / RT)
    else:
        keq = np.exp(-binding_energy / RT)
    
    report_content = f"""
DETAILED PMF ANALYSIS REPORT
===========================

Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Analysis Directory: {plots_dir}

PMF STATISTICS
-------------
Binding Free Energy: {binding_energy:.3f} kcal/mol
Minimum PMF: {pmf_values[min_idx]:.3f} kcal/mol at {distances[min_idx]:.3f} nm
Maximum PMF: {pmf_values[max_idx]:.3f} kcal/mol at {distances[max_idx]:.3f} nm
PMF Standard Deviation: {pmf_std:.3f} kcal/mol
Distance Range: {distances[0]:.3f} - {distances[-1]:.3f} nm ({distance_span:.3f} nm span)

THERMODYNAMIC INTERPRETATION
---------------------------
Estimated Equilibrium Constant (Keq): {keq:.2e}
Binding Affinity: {'Favorable' if binding_energy < 0 else 'Unfavorable'}
Binding Strength: {_classify_binding_strength(binding_energy)}

ERROR ANALYSIS
-------------"""

    if errors is not None:
        avg_error = np.mean(errors)
        max_error = np.max(errors)
        report_content += f"""
Average Bootstrap Error: {avg_error:.3f} kcal/mol
Maximum Bootstrap Error: {max_error:.3f} kcal/mol
Error/Signal Ratio: {avg_error/abs(binding_energy)*100:.1f}%"""
    else:
        report_content += """
Bootstrap errors not available."""

    report_content += f"""

GENERATED PLOTS
--------------"""
    
    for plot_name, plot_path in plots.items():
        report_content += f"""
{plot_name}: {plot_path}"""

    report_content += f"""

RECOMMENDATIONS
--------------
{_generate_pmf_recommendations(binding_energy, pmf_std, errors)}

TECHNICAL NOTES
--------------
- PMF calculated using WHAM method
- Errors estimated using bootstrap resampling
- Force profile calculated as negative gradient of PMF
- All energies relative to PMF minimum

For detailed methodology, refer to:
- Kumar et al. J. Comput. Chem. 13:1011-1021 (1992)
- Hub et al. J. Chem. Theory Comput. 6:3713-3720 (2010)
"""

    with open(report_file, 'w') as f:
        f.write(report_content)
    
    logger.info(f"✓ Detailed PMF report: {report_file}")
    return report_file


def _classify_binding_strength(binding_energy):
    """Classify binding strength based on binding energy"""
    if binding_energy < -10:
        return "Very Strong"
    elif binding_energy < -5:
        return "Strong"
    elif binding_energy < -2:
        return "Moderate"
    elif binding_energy < 0:
        return "Weak"
    elif binding_energy < 5:
        return "Very Weak/Unfavorable"
    else:
        return "Highly Unfavorable"


def _generate_pmf_recommendations(binding_energy, pmf_std, errors):
    """Generate recommendations based on PMF analysis"""
    recommendations = []
    
    if binding_energy > 0:
        recommendations.append("- Positive binding energy suggests unfavorable binding")
        recommendations.append("- Consider checking system setup or extending sampling time")
    
    if abs(binding_energy) < 1:
        recommendations.append("- Very small binding energy - results may be within error margins")
    
    if errors is not None:
        avg_error = np.mean(errors)
        if avg_error > abs(binding_energy) * 0.5:
            recommendations.append("- Large bootstrap errors relative to binding energy")
            recommendations.append("- Consider extending simulation time or using more windows")
    
    if pmf_std > 5:
        recommendations.append("- High PMF fluctuations suggest possible convergence issues")
    
    if len(recommendations) == 0:
        recommendations.append("- PMF results appear reasonable")
        recommendations.append("- Consider validation with experimental data if available")
    
    return "\n".join(recommendations)


def get_workflow_status():
    """
    Get comprehensive workflow status
    
    Returns:
    --------
    Dict : Complete status of all workflow steps
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> status = pmf.get_workflow_status()
    >>> print(status)
    """
    _check_setup()
    
    status = {
        'step1_smd_preparation': False,
        'step1_details': {},
        'step2_umbrella_preparation': False,
        'step2_details': {},
        'step3_wham_analysis': False,
        'step3_details': {},
        'analysis_modules': {}
    }
    
    # Check Step 1: SMD Preparation
    smd_dir = _OUTPUT_DIR / "smd"
    if smd_dir.exists():
        smd_files = {
            'md.gro': smd_dir / "md.gro",
            'topol.top': smd_dir / "topol.top",
            'index.ndx': smd_dir / "index.ndx",
            'smd.mdp': smd_dir / "smd.mdp",
            'run_script': smd_dir / "run_smd.sh",
            'analysis_script': smd_dir / "analyze_smd.sh"
        }
        
        status['step1_smd_preparation'] = all(f.exists() for f in smd_files.values())
        status['step1_details'] = {k: v.exists() for k, v in smd_files.items()}
        
        # Check if SMD has been run
        smd_results = {
            'smd.gro': smd_dir / "results" / "smd.gro",
            'smd.xtc': smd_dir / "results" / "smd.xtc",
            'smd_pullf.xvg': smd_dir / "results" / "smd_pullf.xvg",
            'smd_pullx.xvg': smd_dir / "results" / "smd_pullx.xvg"
        }
        status['step1_details']['smd_completed'] = all(f.exists() for f in smd_results.values())
        status['step1_details']['smd_results'] = {k: v.exists() for k, v in smd_results.items()}
    
    # Check Step 2: Umbrella Preparation
    umbrella_dir = _OUTPUT_DIR / "umbrella"
    if umbrella_dir.exists():
        window_dirs = list(umbrella_dir.glob("window_*"))
        status['step2_details']['n_windows_setup'] = len(window_dirs)
        
        if len(window_dirs) > 0:
            status['step2_umbrella_preparation'] = True
            
            # Check window completion
            completed_windows = 0
            for window_dir in window_dirs:
                if (window_dir / "results" / "umbrella.gro").exists():
                    completed_windows += 1
            
            status['step2_details']['n_windows_completed'] = completed_windows
            status['step2_details']['umbrella_completed'] = completed_windows == len(window_dirs)
    
    # Check Step 3: WHAM Analysis
    analysis_dir = _OUTPUT_DIR / "analysis"
    if analysis_dir.exists():
        wham_files = {
            'pmf.xvg': analysis_dir / "wham" / "pmf.xvg",
            'pmferror.xvg': analysis_dir / "wham" / "pmferror.xvg",
            'pmf_plot': analysis_dir / "pmf_curve.png",
            'report': analysis_dir / "pmf_analysis_report.txt"
        }
        
        status['step3_wham_analysis'] = all(f.exists() for f in wham_files.values())
        status['step3_details'] = {k: v.exists() for k, v in wham_files.items()}
    
    # Check analysis modules
    status['analysis_modules'] = {
        'smd_analysis': (smd_dir / "analysis").exists() if smd_dir.exists() else False,
        'detailed_pmf_analysis': (analysis_dir / "detailed_plots").exists() if analysis_dir.exists() else False
    }
    
    return status


def clean_workflow(steps=None):
    """
    Clean PMF workflow results
    
    Parameters:
    -----------
    steps : list, optional
        List of steps to clean ['step1', 'step2', 'step3', 'analysis']
        If None, prompts user for confirmation to clean all
    
    Example:
    --------
    >>> import pmf
    >>> pmf.setup(config, "./output")
    >>> pmf.clean_workflow(['step1'])  # Only clean step 1
    >>> pmf.clean_workflow()  # Clean all (with confirmation)
    """
    _check_setup()
    
    if steps is None:
        print("This will clean ALL PMF workflow results.")
        confirm = input("Are you sure? (yes/no): ").lower().strip()
        if confirm not in ['yes', 'y']:
            print("Cleaning cancelled.")
            return
        steps = ['step1', 'step2', 'step3', 'analysis']
    
    logger.info(f"Cleaning PMF workflow steps: {steps}")
    
    if 'step1' in steps:
        smd_dir = _OUTPUT_DIR / "smd"
        if smd_dir.exists():
            shutil.rmtree(smd_dir)
            logger.info(f"✓ Cleaned step 1: {smd_dir}")
    
    if 'step2' in steps:
        umbrella_dir = _OUTPUT_DIR / "umbrella"
        if umbrella_dir.exists():
            shutil.rmtree(umbrella_dir)
            logger.info(f"✓ Cleaned step 2: {umbrella_dir}")
    
    if 'step3' in steps or 'analysis' in steps:
        analysis_dir = _OUTPUT_DIR / "analysis"
        if analysis_dir.exists():
            shutil.rmtree(analysis_dir)
            logger.info(f"✓ Cleaned step 3/analysis: {analysis_dir}")
    
    logger.info("✓ Workflow cleaning completed")


def get_expected_runtime_estimates():
    """
    Get estimated runtime for different steps
    
    Returns:
    --------
    Dict : Runtime estimates for planning
    
    Example:
    --------
    >>> import pmf
    >>> estimates = pmf.get_expected_runtime_estimates()
    >>> print(f"SMD estimated time: {estimates['smd_simulation']}")
    """
    _check_setup()
    
    # Calculate based on configuration
    distance_start = _get_config('distance', 'start', 0.3)
    distance_end = _get_config('distance', 'end', 2.0)
    pull_rate = _get_config('smd', 'pull_rate', 0.005)
    production_time = _get_config('umbrella', 'production_time_ps', 22000)
    
    distance_range = distance_end - distance_start
    smd_time_ps = distance_range / pull_rate
    
    # Runtime estimates (highly system-dependent)
    estimates = {
        'smd_simulation': {
            'time_ps': smd_time_ps,
            'estimated_walltime': f"{smd_time_ps/1000:.1f}ns simulation",
            'typical_runtime': "1-6 hours (depends on system size)"
        },
        'umbrella_per_window': {
            'time_ps': production_time,
            'estimated_walltime': f"{production_time/1000:.1f}ns per window",
            'typical_runtime': "2-12 hours per window"
        },
        'total_umbrella_sequential': {
            'note': "If running all windows sequentially (NOT recommended)",
            'typical_runtime': "Days to weeks"
        },
        'wham_analysis': {
            'typical_runtime': "Minutes to hours (depends on bootstrap iterations)"
        },
        'recommendations': [
            "Run umbrella windows in parallel if possible",
            "Use cluster/HPC resources for umbrella sampling",
            "SMD can typically run on single GPU",
            "WHAM analysis is CPU-only but relatively fast"
        ]
    }
    
    return estimates


def print_workflow_summary():
    """
    Print a comprehensive workflow summary
    
    Example:
    --------
    >>> import pmf  
    >>> pmf.setup(config, "./output")
    >>> pmf.print_workflow_summary()
    """
    _check_setup()
    
    print("\n" + "="*60)
    print("PMF WORKFLOW SUMMARY")
    print("="*60)
    
    # Configuration summary
    print("\n🔧 CONFIGURATION:")
    print(f"   Reference Group: {_get_config('reference_group', default='Protein')}")
    print(f"   Moving Group: {_get_config('moving_group', default='LIG')}")
    print(f"   Distance Range: {_get_config('distance', 'start', 0.3):.2f} - {_get_config('distance', 'end', 2.0):.2f} nm")
    print(f"   SMD Pull Rate: {_get_config('smd', 'pull_rate', 0.005)} nm/ps")
    print(f"   Output Directory: {_OUTPUT_DIR}")
    
    # Runtime estimates
    estimates = get_expected_runtime_estimates()
    print(f"\n⏱️  ESTIMATED RUNTIMES:")
    print(f"   SMD Simulation: {estimates['smd_simulation']['typical_runtime']}")
    print(f"   Per Umbrella Window: {estimates['umbrella_per_window']['typical_runtime']}")
    print(f"   WHAM Analysis: {estimates['wham_analysis']['typical_runtime']}")
    
    # Workflow steps
    print(f"\n📋 WORKFLOW STEPS:")
    print(f"   1. pmf_step1_smd_preparation()     → Generate SMD setup")
    print(f"      → Run: bash GMX_PROLIG_PMF/smd/run_smd.sh")
    print(f"      → Analyze: pmf_analyze_smd()")
    print(f"   ")
    print(f"   2. pmf_step2_umbrella_preparation() → Setup umbrella windows")
    print(f"      → Run: bash GMX_PROLIG_PMF/umbrella/run_all_umbrella.sh parallel")
    print(f"   ")
    print(f"   3. pmf_step3_wham_analysis()       → Calculate PMF")
    print(f"      → Analyze: pmf_analyze_pmf()")
    
    # Status check
    status = get_workflow_status()
    print(f"\n📊 CURRENT STATUS:")
    print(f"   Step 1 Prepared: {'✅' if status['step1_smd_preparation'] else '❌'}")
    if status['step1_details'].get('smd_completed'):
        print(f"   SMD Completed: ✅")
    print(f"   Step 2 Prepared: {'✅' if status['step2_umbrella_preparation'] else '❌'}")
    if status['step2_details'].get('umbrella_completed'):
        print(f"   Umbrella Completed: ✅")
    print(f"   Step 3 Completed: {'✅' if status['step3_wham_analysis'] else '❌'}")
    
    print(f"\n💡 RECOMMENDATIONS:")
    for rec in estimates['recommendations']:
        print(f"   • {rec}")
    
    print("="*60)


# Example configuration for reference
def get_example_config():
    """
    Get example configuration dictionary
    
    Returns:
    --------
    Dict : Example configuration with all available options
    
    Example:
    --------
    >>> import pmf
    >>> config = pmf.get_example_config()
    >>> pmf.setup(config, "./output")
    """
    return {
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


if __name__ == "__main__":
    print("PMF Step-by-Step API - Optimized Version")
    print("=" * 50)