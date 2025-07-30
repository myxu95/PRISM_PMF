#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modular PMF (Potential of Mean Force) calculation system

This module provides four independent components for PMF calculations:
1. SMDPreparation - System preparation for SMD
2. SMD - Steered Molecular Dynamics execution  
3. UmbrellaSampling - Umbrella sampling window generation and execution
4. PMFAnalysis - WHAM analysis and visualization

Author: PRISM Team (Modular)
Version: 3.0
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
from abc import ABC, abstractmethod

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PMFBase(ABC):
    """Base class for all PMF modules"""
    
    def __init__(self, config: Dict, output_dir: str):
        """
        Initialize PMF base module
        
        Parameters:
        -----------
        config : Dict
            Module configuration dictionary
        output_dir : str
            Output directory for results
        """
        self.config = config
        self.output_dir = Path(output_dir).resolve()
        self._extract_common_config()
        self._setup_base_directories()
    
    def _extract_common_config(self):
        """Extract common configuration parameters"""
        # Basic settings
        self.reference_group = self.config.get('reference_group', 'Protein')
        self.moving_group = self.config.get('moving_group', 'LIG')
        
        # Simulation parameters
        sim = self.config.get('simulation', {})
        self.dt = sim.get('dt', 0.002)
        self.temperature = sim.get('temperature', 310.0)
        self.pressure = sim.get('pressure', 1.0)
        
        # SMD parameters
        smd = self.config.get('smd', {})
        self.smd_pull_rate = smd.get('pull_rate', 0.005)
        self.smd_pull_k = smd.get('pull_k', 1000.0)
        self.smd_pull_direction = smd.get('pull_direction', [-1, 0, 0])
        
        # Distance parameters
        dist = self.config.get('distance', {})
        self.distance_start = dist.get('start', 0.3)
        self.distance_end = dist.get('end', 2.0)
    
    def _setup_base_directories(self):
        """Setup base directory structure"""
        self.pmf_dir = self.output_dir / "pmf"
        self.mdp_dir = self.output_dir / "mdps"
        
        for directory in [self.pmf_dir, self.mdp_dir]:
            directory.mkdir(parents=True, exist_ok=True)
    
    def _create_common_mdp_sections(self):
        """Create common MDP file sections"""
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
comm-grps           = Protein Non-Protein"""
    
    @abstractmethod
    def run(self):
        """Abstract method to be implemented by each module"""
        pass


class SMDPreparation(PMFBase):
    """Module 1: SMD System Preparation"""
    
    def __init__(self, config: Dict, output_dir: str, md_results_dir: str):
        """
        Initialize SMD preparation module
        
        Parameters:
        -----------
        config : Dict
            SMD preparation configuration
        output_dir : str
            Output directory
        md_results_dir : str
            Directory containing MD simulation results
        """
        super().__init__(config, output_dir)
        self.md_results_dir = Path(md_results_dir).resolve()
        self.smd_dir = self.pmf_dir / "smd"
        self.smd_dir.mkdir(exist_ok=True)
    
    def run(self):
        """Run SMD preparation workflow"""
        logger.info("=== SMD Preparation Module ===")
        
        # Step 1: Validate MD results
        md_files = self._validate_md_files()
        
        # Step 2: Copy necessary files
        self._copy_md_files(md_files)
        
        # Step 3: Create index file
        self._create_index_file()
        
        # Step 4: Generate SMD MDP file
        smd_mdp = self._generate_smd_mdp()
        
        # Step 5: Generate preparation script
        prep_script = self._generate_preparation_script()
        
        results = {
            'smd_dir': str(self.smd_dir),
            'smd_mdp': smd_mdp,
            'preparation_script': prep_script,
            'status': 'ready_for_smd'
        }
        
        logger.info("✓ SMD preparation completed")
        return results
    
    def _validate_md_files(self):
        """Validate required MD files exist"""
        required_files = {
            'gro': 'solv_ions.gro',
            'top': 'topol.top'
        }
        
        found_files = {}
        for file_type, filename in required_files.items():
            file_path = self.md_results_dir / filename
            if file_path.exists():
                found_files[file_type] = str(file_path)
                logger.info(f"✓ Found {file_type}: {filename}")
            else:
                raise FileNotFoundError(f"Required file not found: {file_path}")
        
        return found_files
    
    def _copy_md_files(self, md_files):
        """Copy MD files to PMF directory"""
        for file_type, file_path in md_files.items():
            dest_path = self.pmf_dir / Path(file_path).name
            shutil.copy2(file_path, dest_path)
            logger.info(f"✓ Copied {file_type}: {Path(file_path).name}")
    
    def _create_index_file(self):
        """Create index file for group selection"""
        index_script = self.pmf_dir / "create_index.sh"
        
        script_content = f"""#!/bin/bash
# Create index file for SMD groups
echo "Creating index file for {self.reference_group} and {self.moving_group}..."

cd {self.pmf_dir}

# Create index file
echo "r {self.moving_group}\\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx

echo "✓ Index file created: index.ndx"
"""
        
        with open(index_script, 'w') as f:
            f.write(script_content)
        
        os.chmod(index_script, 0o755)
        logger.info(f"✓ Index creation script: {index_script}")
        
        return str(index_script)
    
    def _generate_smd_mdp(self):
        """Generate SMD MDP file"""
        distance_range = self.distance_end - self.distance_start
        smd_time_ps = distance_range / self.smd_pull_rate
        nsteps = int(smd_time_ps / self.dt)
        output_interval = int(5.0 / self.dt)
        
        smd_content = f"""; SMD (Steered Molecular Dynamics) Parameters
title               = SMD simulation ({distance_range:.2f} nm in {smd_time_ps:.1f} ps)
integrator          = md
nsteps              = {nsteps}
dt                  = {self.dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}
{self._create_common_mdp_sections()}

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
        
        smd_mdp_file = self.mdp_dir / "smd.mdp"
        with open(smd_mdp_file, 'w') as f:
            f.write(smd_content)
        
        logger.info(f"✓ SMD MDP file: {smd_mdp_file}")
        return str(smd_mdp_file)
    
    def _generate_preparation_script(self):
        """Generate complete preparation script"""
        script_path = self.pmf_dir / "smd_preparation.sh"
        
        script_content = f"""#!/bin/bash
######################################################
# SMD PREPARATION SCRIPT
######################################################

set -e
echo "=== SMD Preparation ==="

cd {self.pmf_dir}

# Step 1: Create index file
if [ ! -f index.ndx ]; then
    echo "Creating index file..."
    echo "r {self.moving_group}\\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx
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
echo "Next step: Run SMD module"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_path, 0o755)
        logger.info(f"✓ Preparation script: {script_path}")
        
        return str(script_path)


class SMD(PMFBase):
    """Module 2: Steered Molecular Dynamics"""
    
    def __init__(self, config: Dict, output_dir: str):
        """
        Initialize SMD module
        
        Parameters:
        -----------
        config : Dict
            SMD configuration
        output_dir : str
            Output directory (should contain prepared files)
        """
        super().__init__(config, output_dir)
        self.smd_dir = self.pmf_dir / "smd"
        self.smd_dir.mkdir(exist_ok=True)
    
    def run(self, force_restart=False):
        """
        Run SMD simulation
        
        Parameters:
        -----------
        force_restart : bool
            Whether to restart SMD even if results exist
        """
        logger.info("=== SMD Simulation Module ===")
        
        # Step 1: Validate preparation
        self._validate_preparation()
        
        # Step 2: Check if SMD already completed
        if not force_restart and self._check_smd_completed():
            logger.info("✓ SMD already completed")
            return self._get_smd_results()
        
        # Step 3: Run SMD simulation
        smd_results = self._run_smd_simulation()
        
        # Step 4: Generate SMD analysis
        self._analyze_smd_results()
        
        # Step 5: Extract trajectory frames
        frames_info = self._extract_trajectory_frames()
        
        results = {
            'smd_dir': str(self.smd_dir),
            'smd_trajectory': smd_results['trajectory'],
            'smd_structure': smd_results['structure'],
            'trajectory_frames': frames_info,
            'status': 'smd_completed'
        }
        
        logger.info("✓ SMD simulation completed")
        return results
    
    def _validate_preparation(self):
        """Validate SMD preparation is complete"""
        required_files = [
            self.pmf_dir / "solv_ions.gro",
            self.pmf_dir / "topol.top",
            self.pmf_dir / "index.ndx",
            self.mdp_dir / "smd.mdp"
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"SMD preparation incomplete: {file_path} not found")
        
        logger.info("✓ SMD preparation validated")
    
    def _check_smd_completed(self):
        """Check if SMD simulation is already completed"""
        smd_gro = self.smd_dir / "smd.gro"
        smd_xtc = self.smd_dir / "smd.xtc"
        return smd_gro.exists() and smd_xtc.exists()
    
    def _run_smd_simulation(self):
        """Run the actual SMD simulation"""
        logger.info("Starting SMD simulation...")
        
        # Step 1: Generate TPR file
        tpr_file = self.smd_dir / "smd.tpr"
        grompp_cmd = [
            "gmx", "grompp",
            "-f", str(self.mdp_dir / "smd.mdp"),
            "-c", str(self.pmf_dir / "solv_ions.gro"),
            "-n", str(self.pmf_dir / "index.ndx"),
            "-p", str(self.pmf_dir / "topol.top"),
            "-o", str(tpr_file),
            "-maxwarn", "10"
        ]
        
        try:
            subprocess.run(grompp_cmd, check=True, capture_output=True, text=True)
            logger.info("✓ SMD TPR file generated")
        except subprocess.CalledProcessError as e:
            logger.error(f"grompp failed: {e.stderr}")
            raise
        
        # Step 2: Run MD simulation
        mdrun_cmd = [
            "gmx", "mdrun",
            "-s", str(tpr_file),
            "-deffnm", str(self.smd_dir / "smd"),
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
        
        return {
            'trajectory': str(self.smd_dir / "smd.xtc"),
            'structure': str(self.smd_dir / "smd.gro"),
            'tpr': str(tpr_file)
        }
    
    def _analyze_smd_results(self):
        """Analyze SMD results and create plots"""
        logger.info("Analyzing SMD results...")
        
        # Plot force profile
        pullf_file = self.smd_dir / "smd_pullf.xvg"
        if pullf_file.exists():
            self._plot_force_profile(pullf_file)
        
        # Plot distance profile  
        pullx_file = self.smd_dir / "smd_pullx.xvg"
        if pullx_file.exists():
            self._plot_distance_profile(pullx_file)
    
    def _plot_force_profile(self, pullf_file):
        """Plot SMD force profile"""
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
            plot_file = self.smd_dir / "force_profile.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"✓ Force profile plot: {plot_file}")
        except Exception as e:
            logger.warning(f"Could not plot force profile: {e}")
    
    def _plot_distance_profile(self, pullx_file):
        """Plot SMD distance profile"""
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
            plot_file = self.smd_dir / "distance_profile.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"✓ Distance profile plot: {plot_file}")
        except Exception as e:
            logger.warning(f"Could not plot distance profile: {e}")
    
    def _extract_trajectory_frames(self, n_frames=300):
        """Extract trajectory frames for umbrella sampling"""
        logger.info(f"Extracting {n_frames} trajectory frames...")
        
        frames_dir = self.pmf_dir / "trajectory_frames"
        frames_dir.mkdir(exist_ok=True)
        
        # Calculate frame interval
        distance_range = self.distance_end - self.distance_start
        smd_time_ps = distance_range / self.smd_pull_rate
        frame_interval = int(smd_time_ps / n_frames)
        
        # Extract frames
        trjconv_cmd = [
            "gmx", "trjconv",
            "-f", str(self.smd_dir / "smd.xtc"),
            "-s", str(self.smd_dir / "smd.tpr"),
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
    
    def _get_smd_results(self):
        """Get SMD results when already completed"""
        return {
            'smd_dir': str(self.smd_dir),
            'smd_trajectory': str(self.smd_dir / "smd.xtc"),
            'smd_structure': str(self.smd_dir / "smd.gro"),
            'trajectory_frames': {
                'frames_dir': str(self.pmf_dir / "trajectory_frames"),
                'status': 'existing'
            },
            'status': 'smd_completed'
        }


class UmbrellaSampling(PMFBase):
    """Module 3: Umbrella Sampling"""
    
    def __init__(self, config: Dict, output_dir: str):
        """
        Initialize Umbrella Sampling module
        
        Parameters:
        -----------
        config : Dict
            Umbrella sampling configuration
        output_dir : str
            Output directory (should contain SMD results)
        """
        super().__init__(config, output_dir)
        
        # Umbrella specific config
        umbrella = self.config.get('umbrella', {})
        self.production_time = umbrella.get('production_time_ps', 22000)
        self.sampling_interval = umbrella.get('sampling_interval_ps', 10)
        self.sample_interval_near = umbrella.get('sample_interval_near', 0.1)
        self.sample_interval_far = umbrella.get('sample_interval_far', 0.2)
        self.cutoff_distance = umbrella.get('cutoff_distance', 1.5)
        
        self.umbrella_dir = self.pmf_dir / "umbrella"
        self.umbrella_dir.mkdir(exist_ok=True)
    
    def run(self, auto_run=False):
        """
        Run umbrella sampling preparation and optionally execution
        
        Parameters:
        -----------
        auto_run : bool
            Whether to automatically run all umbrella windows
        """
        logger.info("=== Umbrella Sampling Module ===")
        
        # Step 1: Validate SMD results
        self._validate_smd_results()
        
        # Step 2: Calculate distances for all frames
        distance_data = self._calculate_frame_distances()
        
        # Step 3: Generate umbrella windows
        windows = self._generate_umbrella_windows(distance_data)
        
        # Step 4: Generate umbrella MDP template
        umbrella_mdp = self._generate_umbrella_mdp()
        
        # Step 5: Setup all umbrella windows
        self._setup_umbrella_windows(windows)
        
        # Step 6: Generate run scripts
        run_scripts = self._generate_run_scripts(windows)
        
        results = {
            'umbrella_dir': str(self.umbrella_dir),
            'windows': windows,
            'umbrella_mdp': umbrella_mdp,
            'run_scripts': run_scripts,
            'n_windows': len(windows),
            'status': 'ready_for_umbrella'
        }
        
        # Step 7: Optionally run all windows
        if auto_run:
            results['execution_results'] = self._run_all_windows(windows)
            results['status'] = 'umbrella_completed'
        
        logger.info(f"✓ Umbrella sampling setup completed ({len(windows)} windows)")
        return results
    
    def _validate_smd_results(self):
        """Validate SMD results are available"""
        required_files = [
            self.pmf_dir / "smd" / "smd.xtc",
            self.pmf_dir / "smd" / "smd.tpr",
            self.pmf_dir / "trajectory_frames"
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"SMD results incomplete: {file_path} not found")
        
        logger.info("✓ SMD results validated")
    
    def _calculate_frame_distances(self):
        """Calculate distances for all trajectory frames"""
        logger.info("Calculating distances for trajectory frames...")
        
        frames_dir = self.pmf_dir / "trajectory_frames"
        distance_dir = self.pmf_dir / "distance_analysis"
        distance_dir.mkdir(exist_ok=True)
        
        summary_file = distance_dir / "summary_distances.dat"
        
        if summary_file.exists():
            logger.info("✓ Distance data already exists")
            return self._read_distance_data(summary_file)
        
        # Find all frame files
        frame_files = sorted(frames_dir.glob("frame_*.gro"))
        
        with open(summary_file, 'w') as f:
            f.write("# Frame Distance(nm)\n")
            
            for i, frame_file in enumerate(frame_files):
                # Calculate distance using gmx distance
                distance_cmd = [
                    "gmx", "distance",
                    "-s", str(self.pmf_dir / "smd" / "smd.tpr"),
                    "-f", str(frame_file),
                    "-n", str(self.pmf_dir / "index.ndx"),
                    "-select", f"com of group {self.reference_group} plus com of group {self.moving_group}",
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
        return self._read_distance_data(summary_file)
    
    def _read_distance_data(self, summary_file):
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
    
    def _generate_umbrella_windows(self, distance_data):
        """Generate adaptive umbrella sampling windows"""
        logger.info("Generating umbrella sampling windows...")
        
        if not distance_data:
            raise RuntimeError("No distance data available for window generation")
        
        # Adaptive sampling
        sampled_indices = []
        current_idx = 0
        sampled_indices.append(current_idx)
        
        distances = [d[1] for d in distance_data]
        
        while current_idx < len(distances):
            current_distance = distances[current_idx]
            
            # Choose interval based on distance
            if current_distance < self.cutoff_distance:
                target_interval = self.sample_interval_near
            else:
                target_interval = self.sample_interval_far
            
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
        for i, idx in enumerate(sampled_indices):
            frame, distance = distance_data[idx]
            windows.append({
                'window_id': i,
                'frame': frame,
                'distance': distance,
                'directory': self.umbrella_dir / f"window_{i:03d}"
            })
        
        logger.info(f"✓ Generated {len(windows)} umbrella windows")
        return windows
    
    def _generate_umbrella_mdp(self):
        """Generate umbrella sampling MDP template"""
        nsteps = int(self.production_time / self.dt)
        output_interval = int(self.sampling_interval / self.dt)
        
        umbrella_content = f"""; Umbrella Sampling Parameters
title               = Umbrella sampling (22.0 ns)
integrator          = md
nsteps              = {nsteps}
dt                  = {self.dt}

; Output control
nstxtcout           = {output_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {output_interval}
{self._create_common_mdp_sections()}

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
        
        umbrella_mdp_file = self.mdp_dir / "umbrella.mdp"
        with open(umbrella_mdp_file, 'w') as f:
            f.write(umbrella_content)
        
        logger.info(f"✓ Umbrella MDP template: {umbrella_mdp_file}")
        return str(umbrella_mdp_file)
    
    def _setup_umbrella_windows(self, windows):
        """Setup all umbrella sampling windows"""
        logger.info("Setting up umbrella sampling windows...")
        
        for window in windows:
            window_dir = window['directory']
            window_dir.mkdir(exist_ok=True)
            
            # Create window-specific MDP file
            self._create_window_mdp(window)
            
            # Copy starting structure
            self._copy_window_structure(window)
        
        logger.info("✓ All umbrella windows setup completed")
    
    def _create_window_mdp(self, window):
        """Create window-specific MDP file"""
        # Read template and replace placeholder
        template_file = self.mdp_dir / "umbrella.mdp"
        with open(template_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        content = content.replace("DISTANCE_PLACEHOLDER", f"{window['distance']:.6f}")
        content = content.replace(
            "title               = Umbrella sampling (22.0 ns)",
            f"title               = Umbrella sampling at {window['distance']:.3f} nm (window {window['window_id']:03d})"
        )
        
        # Write window-specific MDP
        window_mdp = window['directory'] / "umbrella.mdp"
        with open(window_mdp, 'w') as f:
            f.write(content)
    
    def _copy_window_structure(self, window):
        """Copy starting structure for window"""
        frame_file = self.pmf_dir / "trajectory_frames" / f"frame_{window['frame']}.gro"
        dest_file = window['directory'] / "start.gro"
        
        if frame_file.exists():
            shutil.copy2(frame_file, dest_file)
        else:
            logger.warning(f"Frame file not found: {frame_file}")
    
    def _generate_run_scripts(self, windows):
        """Generate run scripts for umbrella sampling"""
        # Individual window scripts
        for window in windows:
            self._create_window_script(window)
        
        # Master run script
        master_script = self._create_master_script(windows)
        
        return {
            'master_script': master_script,
            'individual_scripts': [str(w['directory'] / "run_umbrella.sh") for w in windows]
        }
    
    def _create_window_script(self, window):
        """Create run script for individual window"""
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
    
    def _create_master_script(self, windows):
        """Create master script to run all windows"""
        script_content = f"""#!/bin/bash
######################################################
# MASTER UMBRELLA SAMPLING SCRIPT
######################################################

set -e

echo "=== Running All Umbrella Sampling Windows ==="
echo "Total windows: {len(windows)}"

cd {self.umbrella_dir}

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
        
        script_file = self.umbrella_dir / "run_all_umbrella.sh"
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        os.chmod(script_file, 0o755)
        
        logger.info(f"✓ Master umbrella script: {script_file}")
        return str(script_file)
    
    def _run_all_windows(self, windows):
        """Run all umbrella sampling windows"""
        logger.info("Running all umbrella sampling windows...")
        
        results = {}
        for window in windows:
            try:
                self._run_single_window(window)
                results[window['window_id']] = 'completed'
                logger.info(f"✓ Window {window['window_id']:03d} completed")
            except Exception as e:
                results[window['window_id']] = f'failed: {e}'
                logger.error(f"✗ Window {window['window_id']:03d} failed: {e}")
        
        return results
    
    def _run_single_window(self, window):
        """Run single umbrella sampling window"""
        script_file = window['directory'] / "run_umbrella.sh"
        subprocess.run(["bash", str(script_file)], check=True)


class PMFAnalysis(PMFBase):
    """Module 4: PMF Analysis and Visualization"""
    
    def __init__(self, config: Dict, output_dir: str):
        """
        Initialize PMF Analysis module
        
        Parameters:
        -----------
        config : Dict
            PMF analysis configuration
        output_dir : str
            Output directory (should contain umbrella results)
        """
        super().__init__(config, output_dir)
        
        # Analysis specific config
        analysis = self.config.get('analysis', {})
        self.wham_begin_time = analysis.get('begin_time_ps', 2000)
        self.bootstrap_iterations = analysis.get('bootstrap_iterations', 50)
        self.energy_unit = analysis.get('energy_unit', 'kCal')
        
        self.analysis_dir = self.pmf_dir / "analysis"
        self.analysis_dir.mkdir(exist_ok=True)
    
    def run(self, umbrella_dir=None):
        """
        Run complete PMF analysis
        
        Parameters:
        -----------
        umbrella_dir : str, optional
            Directory containing umbrella sampling results
        """
        logger.info("=== PMF Analysis Module ===")
        
        if umbrella_dir is None:
            umbrella_dir = self.pmf_dir / "umbrella"
        else:
            umbrella_dir = Path(umbrella_dir)
        
        # Step 1: Validate umbrella results
        window_files = self._validate_umbrella_results(umbrella_dir)
        
        # Step 2: Run WHAM analysis
        wham_results = self._run_wham_analysis(window_files)
        
        # Step 3: Generate visualizations
        plots = self._generate_visualizations(wham_results)
        
        # Step 4: Calculate binding energy
        binding_energy = self._calculate_binding_energy(wham_results)
        
        # Step 5: Generate analysis report
        report = self._generate_analysis_report(wham_results, binding_energy, plots)
        
        results = {
            'analysis_dir': str(self.analysis_dir),
            'wham_results': wham_results,
            'binding_energy': binding_energy,
            'plots': plots,
            'report': report,
            'status': 'analysis_completed'
        }
        
        logger.info("✓ PMF analysis completed")
        return results
    
    def _validate_umbrella_results(self, umbrella_dir):
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
    
    def _run_wham_analysis(self, window_files):
        """Run WHAM analysis"""
        logger.info("Running WHAM analysis...")
        
        wham_dir = self.analysis_dir / "wham"
        wham_dir.mkdir(exist_ok=True)
        
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
            "-b", str(self.wham_begin_time),
            "-o", "pmf.xvg",
            "-unit", self.energy_unit,
            "-bsres", "pmferror.xvg",
            "-bsprof", "zerrorprofile.xvg",
            "-nBootstrap", str(self.bootstrap_iterations),
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
    
    def _generate_visualizations(self, wham_results):
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
        
        # Create PMF plot
        plt.figure(figsize=(12, 8))
        
        if errors is not None:
            plt.errorbar(distances, pmf_values, yerr=errors, fmt='o-', capsize=3, linewidth=2, markersize=4)
        else:
            plt.plot(distances, pmf_values, 'o-', linewidth=2, markersize=4)
        
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.ylabel(f'PMF ({self.energy_unit.lower()}/mol)', fontsize=14)
        plt.title('Potential of Mean Force', fontsize=16)
        plt.grid(True, alpha=0.3)
        
        # Add binding energy annotation
        min_pmf = np.min(pmf_values)
        max_pmf = np.max(pmf_values)
        binding_energy = max_pmf - min_pmf
        
        plt.text(0.05, 0.95, f'Binding Energy: {binding_energy:.2f} {self.energy_unit.lower()}/mol', 
                transform=plt.gca().transAxes, fontsize=12, 
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        
        # Save plot
        plot_file = self.analysis_dir / "pmf_curve.png"
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
    
    def _calculate_binding_energy(self, wham_results):
        """Calculate binding energy from PMF"""
        pmf_file = Path(wham_results['pmf_file'])
        pmf_data = np.loadtxt(pmf_file, comments=['#', '@'])
        pmf_values = pmf_data[:, 1]
        
        min_pmf = np.min(pmf_values)
        max_pmf = np.max(pmf_values)
        binding_energy = max_pmf - min_pmf
        
        logger.info(f"✓ Binding energy: {binding_energy:.2f} {self.energy_unit.lower()}/mol")
        
        return {
            'value': binding_energy,
            'unit': f'{self.energy_unit.lower()}/mol',
            'min_pmf': min_pmf,
            'max_pmf': max_pmf
        }
    
    def _generate_analysis_report(self, wham_results, binding_energy, plots):
        """Generate comprehensive analysis report"""
        report_file = self.analysis_dir / "pmf_analysis_report.txt"
        
        report_content = f"""
PMF Analysis Report
==================

Analysis Date: {os.popen('date').read().strip()}
Analysis Directory: {self.analysis_dir}

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
- Reference Group: {self.reference_group}
- Moving Group: {self.moving_group}
- WHAM Begin Time: {self.wham_begin_time} ps
- Bootstrap Iterations: {self.bootstrap_iterations}
- Energy Unit: {self.energy_unit}

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


# Unified PMF workflow manager
class PMFWorkflow:
    """Unified PMF workflow manager for running all modules"""
    
    def __init__(self, config: Dict, output_dir: str, md_results_dir: str):
        """
        Initialize PMF workflow manager
        
        Parameters:
        -----------
        config : Dict
            Complete PMF configuration
        output_dir : str
            Output directory
        md_results_dir : str
            MD results directory
        """
        self.config = config
        self.output_dir = output_dir
        self.md_results_dir = md_results_dir
        
        # Initialize modules
        self.smd_prep = SMDPreparation(config, output_dir, md_results_dir)
        self.smd = SMD(config, output_dir)
        self.umbrella = UmbrellaSampling(config, output_dir)
        self.analysis = PMFAnalysis(config, output_dir)
    
    def run_complete_workflow(self, auto_run_umbrella=False):
        """Run complete PMF workflow"""
        logger.info("=" * 60)
        logger.info("STARTING COMPLETE PMF WORKFLOW")
        logger.info("=" * 60)
        
        results = {}
        
        # Step 1: SMD Preparation
        logger.info("Step 1: SMD Preparation")
        results['smd_preparation'] = self.smd_prep.run()
        
        # Step 2: SMD Simulation
        logger.info("Step 2: SMD Simulation")
        results['smd'] = self.smd.run()
        
        # Step 3: Umbrella Sampling
        logger.info("Step 3: Umbrella Sampling")
        results['umbrella_sampling'] = self.umbrella.run(auto_run=auto_run_umbrella)
        
        # Step 4: PMF Analysis (only if umbrella sampling was run)
        if auto_run_umbrella or results['umbrella_sampling']['status'] == 'umbrella_completed':
            logger.info("Step 4: PMF Analysis")
            results['pmf_analysis'] = self.analysis.run()
        else:
            logger.info("Step 4: PMF Analysis (skipped - umbrella sampling not run)")
            results['pmf_analysis'] = {'status': 'skipped'}
        
        logger.info("=" * 60)
        logger.info("PMF WORKFLOW COMPLETED")
        logger.info("=" * 60)
        
        return results
    
    def run_module(self, module_name: str, **kwargs):
        """Run individual module"""
        modules = {
            'smd_preparation': self.smd_prep,
            'smd': self.smd,
            'umbrella_sampling': self.umbrella,
            'pmf_analysis': self.analysis
        }
        
        if module_name not in modules:
            raise ValueError(f"Unknown module: {module_name}. Available: {list(modules.keys())}")
        
        logger.info(f"Running individual module: {module_name}")
        return modules[module_name].run(**kwargs)


# Main entry points
def run_pmf_analysis(config: Dict, output_dir: str, md_results_dir: str, auto_run_umbrella: bool = False) -> Dict:
    """
    Main entry point for complete PMF analysis
    
    Parameters:
    -----------
    config : Dict
        PMF configuration dictionary
    output_dir : str
        Output directory
    md_results_dir : str
        MD results directory
    auto_run_umbrella : bool
        Whether to automatically run umbrella sampling
    
    Returns:
    --------
    Dict : Complete workflow results
    """
    if not config.get('enabled', False):
        logger.info("PMF analysis is disabled in configuration")
        return None
    
    workflow = PMFWorkflow(config, output_dir, md_results_dir)
    return workflow.run_complete_workflow(auto_run_umbrella=auto_run_umbrella)


def run_pmf_module(module_name: str, config: Dict, output_dir: str, md_results_dir: str = None, **kwargs) -> Dict:
    """
    Run individual PMF module
    
    Parameters:
    -----------
    module_name : str
        Module name ('smd_preparation', 'smd', 'umbrella_sampling', 'pmf_analysis')
    config : Dict
        PMF configuration dictionary
    output_dir : str
        Output directory
    md_results_dir : str, optional
        MD results directory (required for smd_preparation)
    **kwargs : additional arguments for the module
    
    Returns:
    --------
    Dict : Module execution results
    """
    if module_name == 'smd_preparation':
        if md_results_dir is None:
            raise ValueError("md_results_dir is required for smd_preparation module")
        module = SMDPreparation(config, output_dir, md_results_dir)
    elif module_name == 'smd':
        module = SMD(config, output_dir)
    elif module_name == 'umbrella_sampling':
        module = UmbrellaSampling(config, output_dir)
    elif module_name == 'pmf_analysis':
        module = PMFAnalysis(config, output_dir)
    else:
        raise ValueError(f"Unknown module: {module_name}")
    
    return module.run(**kwargs)


# Usage examples
if __name__ == "__main__":
    # Example configuration
    config = {
        'enabled': True,
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
    
    # Example 1: Run complete workflow
    print("Example 1: Complete workflow")
    # results = run_pmf_analysis(config, "./pmf_output", "./md_results", auto_run_umbrella=False)
    
    # Example 2: Run individual modules
    print("\nExample 2: Individual modules")
    
    # Run SMD preparation only
    # smd_prep_results = run_pmf_module('smd_preparation', config, "./pmf_output", "./md_results")
    
    # Run SMD simulation only
    # smd_results = run_pmf_module('smd', config, "./pmf_output")
    
    # Run umbrella sampling setup only
    # umbrella_results = run_pmf_module('umbrella_sampling', config, "./pmf_output", auto_run=False)
    
    # Run PMF analysis only
    # analysis_results = run_pmf_module('pmf_analysis', config, "./pmf_output")
    
    print("PMF modular system ready for use!")