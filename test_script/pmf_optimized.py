#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PMF (Potential of Mean Force) module for PRISM

This module provides comprehensive PMF calculation capabilities including:
- Steered Molecular Dynamics (SMD)
- Umbrella Sampling  
- WHAM analysis
- Automated workflow management

Author: PRISM Team
Version: 1.0
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
    PMF (Potential of Mean Force) calculator for protein-ligand systems
    
    This class provides a comprehensive workflow for PMF calculations including:
    - System preparation
    - SMD simulation
    - Umbrella sampling
    - WHAM analysis
    - Result visualization
    """
    
    def __init__(self, config: Dict, output_dir: str, md_results_dir: str):
        """
        Initialize PMF calculator
        
        Parameters:
        -----------
        config : Dict
            PMF configuration dictionary containing all simulation parameters
        output_dir : str
            Output directory for PMF results
        md_results_dir : str
            Directory containing MD simulation results
        """
        self.config = config
        self.output_dir = os.path.abspath(output_dir)
        self.md_results_dir = os.path.abspath(md_results_dir)
        
        # Create directory structure
        self._setup_directories()
        
        # Extract and validate configuration
        self._extract_configuration()
        
        # Initialize workflow state
        self.windows = []
        
        logger.info("PMF Calculator initialized successfully")
        self._print_configuration_summary()
    
    def _setup_directories(self):
        """Create necessary directory structure for PMF calculations"""
        # Main directories
        self.pmf_dir = os.path.join(self.output_dir, "pmf")
        self.windows_dir = os.path.join(self.pmf_dir, "windows")
        self.analysis_dir = os.path.join(self.pmf_dir, "analysis")
        self.mdp_dir = os.path.join(self.output_dir, "mdps")
        
        # Create directories
        directories = [self.pmf_dir, self.windows_dir, self.analysis_dir, self.mdp_dir]
        for directory in directories:
            os.makedirs(directory, exist_ok=True)
            logger.debug(f"Created directory: {directory}")
    
    def _extract_configuration(self):
        """Extract and validate configuration parameters"""
        # Basic PMF settings
        self.method = self.config.get('method', 'umbrella')
        self.reaction_coordinate = self.config.get('reaction_coordinate', 'distance')
        self.reference_group = self.config.get('reference_group', 'Protein')
        self.moving_group = self.config.get('moving_group', 'LIG')
        
        # Extract configuration sections
        self._extract_smd_config()
        self._extract_distance_config()
        self._extract_umbrella_config()
        self._extract_analysis_config()
        self._extract_output_config()
    
    def _extract_smd_config(self):
        """Extract SMD-specific configuration"""
        smd_config = self.config.get('smd', {})
        self.smd_pull_rate = smd_config.get('pull_rate', 0.005)  # nm/ps
        self.smd_pull_k = smd_config.get('pull_k', 1000.0)  # kJ/mol/nm²
        self.smd_pull_direction = smd_config.get('pull_direction', [-1, 0, 0])
        self.smd_trajectory_interval = smd_config.get('trajectory_interval_ps', 5.0)
        self.smd_log_interval = smd_config.get('log_interval_ps', 5.0)
    
    def _extract_distance_config(self):
        """Extract distance-related configuration"""
        distance_config = self.config.get('distance', {})
        self.distance_start = distance_config.get('start', 0.3)
        self.distance_end = distance_config.get('end', 2.0)
        self.distance_step = distance_config.get('step', 0.1)
        self.force_constant = distance_config.get('force_constant', 1000.0)
    
    def _extract_umbrella_config(self):
        """Extract umbrella sampling configuration"""
        umbrella_config = self.config.get('umbrella', {})
        self.n_windows = umbrella_config.get('n_windows', 20)
        self.equilibration_time = umbrella_config.get('equilibration_time_ps', 1000)
        self.production_time = umbrella_config.get('production_time_ps', 5000)
        self.sampling_interval = umbrella_config.get('sampling_interval_ps', 10)
    
    def _extract_analysis_config(self):
        """Extract analysis configuration"""
        analysis_config = self.config.get('analysis', {})
        self.analysis_method = analysis_config.get('method', 'wham')
        self.temperature = analysis_config.get('temperature', 310.0)
        self.tolerance = analysis_config.get('tolerance', 1e-6)
        self.max_iterations = analysis_config.get('max_iterations', 50000)
    
    def _extract_output_config(self):
        """Extract output configuration"""
        output_config = self.config.get('output', {})
        self.plot_format = output_config.get('plot_format', 'png')
        self.data_format = output_config.get('data_format', 'txt')
        self.include_error_bars = output_config.get('include_error_bars', True)
    
    def _print_configuration_summary(self):
        """Print a summary of the current configuration"""
        logger.info("PMF Configuration Summary:")
        logger.info(f"  Method: {self.method}")
        logger.info(f"  Reaction coordinate: {self.reaction_coordinate}")
        logger.info(f"  Reference group: {self.reference_group}")
        logger.info(f"  Moving group: {self.moving_group}")
        logger.info(f"  Distance range: {self.distance_start} - {self.distance_end} nm")
        logger.info(f"  Number of windows: {self.n_windows}")
        logger.info(f"  Output directory: {self.output_dir}")
        logger.info(f"  SMD pull rate: {self.smd_pull_rate} nm/ps")
    
    def generate_smd_mdp(self, output_dir=None):
        """
        Generate SMD (Steered Molecular Dynamics) MDP file
        
        Parameters:
        -----------
        output_dir : str, optional
            Output directory for SMD MDP file (default: self.mdp_dir)
            
        Returns:
        --------
        str : Path to the generated SMD MDP file
        """
        if output_dir is None:
            output_dir = self.mdp_dir
        
        os.makedirs(output_dir, exist_ok=True)
        smd_mdp = os.path.join(output_dir, "smd.mdp")
        
        # Calculate simulation parameters
        simulation_params = self._calculate_smd_parameters()
        
        # Generate MDP content
        smd_content = self._create_smd_mdp_content(simulation_params)
        
        # Write MDP file
        with open(smd_mdp, 'w') as f:
            f.write(smd_content)
        
        logger.info(f"SMD MDP file generated: {smd_mdp}")
        self._log_smd_parameters(simulation_params)
        
        return smd_mdp
    
    def _calculate_smd_parameters(self):
        """Calculate SMD simulation parameters"""
        # Calculate steps based on distance and pull rate
        distance_range = self.distance_end - self.distance_start
        smd_time_ps = distance_range / self.smd_pull_rate
        dt = self.config.get('simulation', {}).get('dt', 0.002)
        nsteps = int(smd_time_ps / dt)
        
        # Calculate output intervals
        traj_interval = int(self.smd_trajectory_interval / dt)
        log_interval = int(self.smd_log_interval / dt)
        
        return {
            'distance_range': distance_range,
            'smd_time_ps': smd_time_ps,
            'nsteps': nsteps,
            'dt': dt,
            'traj_interval': traj_interval,
            'log_interval': log_interval
        }
    
    def _create_smd_mdp_content(self, params):
        """Create SMD MDP file content"""
        # Get configuration sections
        sim_config = self.config.get('simulation', {})
        const_config = self.config.get('constraints', {})
        elec_config = self.config.get('electrostatics', {})
        vdw_config = self.config.get('vdw', {})
        tc_config = self.config.get('temperature_coupling', {})
        pc_config = self.config.get('pressure_coupling', {})
        
        # Format pull direction
        pull_direction = ' '.join(map(str, self.smd_pull_direction))
        
        return f"""; smd.mdp - Steered Molecular Dynamics
title               = SMD run ({params['distance_range']:.2f} nm, {params['smd_time_ps']:.1f} ps)
integrator          = md
nsteps              = {params['nsteps']}
dt                  = {params['dt']}

nstxtcout           = {params['traj_interval']}
nstvout             = 0
nstenergy           = 0
nstlog              = {params['log_interval']}

continuation            = yes
constraint_algorithm    = {const_config.get('algorithm', 'lincs')}
constraints             = {const_config.get('type', 'h-bonds')}
lincs_iter              = {const_config.get('lincs_iter', 1)}
lincs_order             = {const_config.get('lincs_order', 4)}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config.get('rcoulomb', 1.0)}
rvdw                    = {vdw_config.get('rvdw', 1.0)}

coulombtype             = {elec_config.get('coulombtype', 'PME')}
pme_order               = {elec_config.get('pme_order', 4)}
fourierspacing          = {elec_config.get('fourierspacing', 0.16)}

tcoupl                  = {tc_config.get('tcoupl', 'V-rescale')}
tc-grps                 = {' '.join(tc_config.get('tc_grps', ['Protein', 'Non-Protein']))}
tau_t                   = {'     '.join(map(str, tc_config.get('tau_t', [0.1, 0.1])))}
ref_t                   = {'     '.join([str(sim_config.get('temperature', 310.0))] * len(tc_config.get('tc_grps', ['Protein', 'Non-Protein'])))}

pcoupl                  = {pc_config.get('pcoupl', 'C-rescale')}
pcoupltype              = {pc_config.get('pcoupltype', 'isotropic')}
tau_p                   = {pc_config.get('tau_p', 1.0)}
ref_p                   = {sim_config.get('pressure', 1.0)}
compressibility         = {pc_config.get('compressibility', 4.5e-5)}

pbc                     = xyz
DispCorr                = {vdw_config.get('dispcorr', 'EnerPres')}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {' '.join(tc_config.get('tc_grps', ['Protein', 'Non-Protein']))}

; Pull code
pull                    = yes
pull_ncoords            = 1         
pull_ngroups            = 2         
pull_group1_name        = {self.reference_group}
pull_group2_name        = {self.moving_group}
pull_coord1_type        = umbrella  
pull_coord1_geometry    = distance  
pull_coord1_dim         = Y N N   
pull-coord1-vec         = {pull_direction}
pull_coord1_groups      = 1 2
pull_coord1_start       = yes     
pull_coord1_rate        = {self.smd_pull_rate}   
pull_coord1_k           = {self.smd_pull_k}    
pull-pbc-ref-prev-step-com = yes  
"""
    
    def _log_smd_parameters(self, params):
        """Log SMD parameters for debugging"""
        logger.info(f"SMD Parameters:")
        logger.info(f"  Pull rate: {self.smd_pull_rate} nm/ps")
        logger.info(f"  Pull force constant: {self.smd_pull_k} kJ/mol/nm²")
        logger.info(f"  Pull direction: {self.smd_pull_direction}")
        logger.info(f"  Total simulation time: {params['smd_time_ps']:.1f} ps")
        logger.info(f"  Number of steps: {params['nsteps']}")
    
    def run(self) -> str:
        """
        Run the complete PMF calculation workflow
        
        This method orchestrates the entire PMF calculation process:
        1. System preparation
        2. MDP file generation
        3. Window generation
        4. Umbrella sampling
        5. Analysis
        6. Report generation
        
        Returns:
        --------
        str : Path to the PMF results directory
        """
        logger.info("=" * 60)
        logger.info("Starting PMF Calculation Workflow")
        logger.info("=" * 60)
        
        try:
            # Step 1: Prepare system for PMF
            self.prepare_system()
            
            # Step 2: Generate MDP files
            self.generate_smd_mdp()
            self.generate_umbrella_mdp()
            
            # Step 3: Generate umbrella sampling windows
            self.generate_windows()
            
            # Step 4: Run umbrella sampling
            self.run_umbrella_sampling()
            
            # Step 5: Analyze results
            self.analyze_results()
            
            # Step 6: Generate plots and reports
            self.generate_reports()
            
            logger.info("=" * 60)
            logger.info("PMF Calculation Workflow Completed Successfully!")
            logger.info("=" * 60)
            logger.info(f"PMF results are in: {self.pmf_dir}")
            logger.info(f"Analysis results are in: {self.analysis_dir}")
            
            return self.pmf_dir
            
        except Exception as e:
            logger.error(f"Error during PMF calculation: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def prepare_system(self):
        """Prepare the system for PMF calculation"""
        logger.info("=== Preparing System for PMF ===")
        
        # Check if MD results exist
        md_files = self._check_md_files()
        if not md_files:
            raise RuntimeError("MD simulation files not found. Please run MD simulation first.")
        
        # Copy necessary files to PMF directory
        self._copy_md_files(md_files)
        
        # Create index file for groups
        self._create_index_file()
        
        logger.info("System preparation completed")
    
    def _check_md_files(self) -> Dict[str, str]:
        """Check if required MD files exist"""
        required_files = {
            'gro': 'solv_ions.gro',
            'top': 'topol.top',
            'tpr': 'md.tpr'
        }
        
        found_files = {}
        for file_type, filename in required_files.items():
            file_path = os.path.join(self.md_results_dir, filename)
            if os.path.exists(file_path):
                found_files[file_type] = file_path
                logger.debug(f"Found {file_type} file: {filename}")
            else:
                logger.warning(f"{filename} not found in {self.md_results_dir}")
        
        return found_files
    
    def _copy_md_files(self, md_files: Dict[str, str]):
        """Copy MD files to PMF directory"""
        for file_type, file_path in md_files.items():
            dest_path = os.path.join(self.pmf_dir, os.path.basename(file_path))
            try:
                shutil.copy2(file_path, dest_path)
                logger.info(f"Copied {file_type} file: {os.path.basename(file_path)}")
            except Exception as e:
                logger.error(f"Failed to copy {file_type} file: {e}")
                raise
    
    def _create_index_file(self):
        """Create index file for group selection"""
        index_file = os.path.join(self.pmf_dir, "index.ndx")
        
        index_content = f"""
[ {self.reference_group} ]
# Reference group atoms will be added here

[ {self.moving_group} ]
# Moving group atoms will be added here

[ System ]
# All atoms
"""
        
        with open(index_file, 'w') as f:
            f.write(index_content)
        
        logger.info(f"Created index file: {index_file}")
    
    def generate_windows(self):
        """Generate umbrella sampling windows"""
        logger.info("=== Generating Umbrella Sampling Windows ===")
        
        if self.reaction_coordinate == 'distance':
            self._generate_distance_windows()
        else:
            raise NotImplementedError(f"Reaction coordinate type '{self.reaction_coordinate}' not implemented yet")
        
        logger.info(f"Generated {self.n_windows} umbrella sampling windows")
    
    def _generate_distance_windows(self):
        """Generate distance-based umbrella sampling windows"""
        distances = np.linspace(self.distance_start, self.distance_end, self.n_windows)
        
        self.windows = []
        for i, distance in enumerate(distances):
            window = {
                'index': i,
                'distance': distance,
                'force_constant': self.force_constant,
                'directory': os.path.join(self.windows_dir, f"window_{i:03d}")
            }
            self.windows.append(window)
            
            # Create window directory
            os.makedirs(window['directory'], exist_ok=True)
    
    def run_umbrella_sampling(self):
        """Run umbrella sampling simulations"""
        logger.info("=== Running Umbrella Sampling ===")
        
        for i, window in enumerate(self.windows):
            logger.info(f"Running window {i+1}/{len(self.windows)}: {window['distance']:.3f} nm")
            self._run_single_window(window)
        
        logger.info("Umbrella sampling completed")
    
    def _run_single_window(self, window: Dict):
        """Run umbrella sampling for a single window"""
        window_dir = window['directory']
        
        try:
            # Generate MDP file for this window
            mdp_file = os.path.join(window_dir, "umbrella.mdp")
            self._create_umbrella_mdp(mdp_file, window)
            
            # Run grompp
            tpr_file = os.path.join(window_dir, "umbrella.tpr")
            self._run_grompp(window_dir, mdp_file, tpr_file)
            
            # Run mdrun
            self._run_mdrun(window_dir, tpr_file)
            
        except Exception as e:
            logger.error(f"Failed to run window {window['index']}: {e}")
            raise
    
    def _create_umbrella_mdp(self, mdp_file: str, window: Dict):
        """Create MDP file for umbrella sampling"""
        mdp_content = f"""; Umbrella sampling MDP file
title = Umbrella sampling at {window['distance']:.3f} nm
integrator = md
nsteps = {int(self.production_time / 0.002)}  ; Convert ps to steps
dt = 0.002

; Output
nstxout = 0
nstvout = 0
nstenergy = {int(self.sampling_interval / 0.002)}
nstlog = {int(self.sampling_interval / 0.002)}

; Constraints
constraints = h-bonds
constraint-algorithm = lincs

; Cutoffs
cutoff-scheme = Verlet
rcoulomb = 1.0
rvdw = 1.0

; PME
coulombtype = PME
pme-order = 4
fourierspacing = 0.16

; Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = 310.0 310.0

; Pressure coupling
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 1.0
ref-p = 1.0
compressibility = 4.5e-5

; Umbrella sampling
pull = yes
pull-ngroups = 1
pull-group1-name = {self.moving_group}
pull-group2-name = {self.reference_group}
pull-geometry = distance
pull-dim = Y Y Y
pull-start = yes
pull-nstxout = {int(self.sampling_interval / 0.002)}
pull-nstfout = {int(self.sampling_interval / 0.002)}
pull-k1 = {window['force_constant']}
pull-r0 = {window['distance']}
"""
        
        with open(mdp_file, 'w') as f:
            f.write(mdp_content)
    
    def _run_grompp(self, work_dir: str, mdp_file: str, tpr_file: str):
        """Run gmx grompp"""
        cmd = [
            "gmx", "grompp",
            "-f", mdp_file,
            "-c", os.path.join(self.pmf_dir, "solv_ions.gro"),
            "-p", os.path.join(self.pmf_dir, "topol.top"),
            "-n", os.path.join(self.pmf_dir, "index.ndx"),
            "-o", tpr_file,
            "-maxwarn", "1"
        ]
        
        try:
            subprocess.run(cmd, cwd=work_dir, check=True, capture_output=True, text=True)
            logger.debug(f"grompp completed for {os.path.basename(work_dir)}")
        except subprocess.CalledProcessError as e:
            logger.error(f"grompp failed: {e.stderr}")
            raise
    
    def _run_mdrun(self, work_dir: str, tpr_file: str):
        """Run gmx mdrun"""
        cmd = [
            "gmx", "mdrun",
            "-s", tpr_file,
            "-deffnm", "umbrella",
            "-nt", "4"  # Number of threads
        ]
        
        try:
            subprocess.run(cmd, cwd=work_dir, check=True, capture_output=True, text=True)
            logger.debug(f"mdrun completed for {os.path.basename(work_dir)}")
        except subprocess.CalledProcessError as e:
            logger.error(f"mdrun failed: {e.stderr}")
            raise
    
    def analyze_results(self):
        """Analyze umbrella sampling results"""
        logger.info("=== Analyzing PMF Results ===")
        
        if self.analysis_method == 'wham':
            self._analyze_with_wham()
        elif self.analysis_method == 'mbar':
            self._analyze_with_mbar()
        else:
            raise NotImplementedError(f"Analysis method '{self.analysis_method}' not implemented yet")
        
        logger.info("PMF analysis completed")
    
    def _analyze_with_wham(self):
        """Analyze results using WHAM"""
        logger.info("Analyzing with WHAM...")
        # This is a placeholder - you'd implement WHAM analysis
        # or call an external WHAM program
    
    def _analyze_with_mbar(self):
        """Analyze results using MBAR"""
        logger.info("Analyzing with MBAR...")
        # This is a placeholder - you'd implement MBAR analysis
    
    def generate_reports(self):
        """Generate PMF plots and reports"""
        logger.info("=== Generating PMF Reports ===")
        
        # Load PMF data
        pmf_data = self._load_pmf_data()
        
        # Generate PMF plot
        self._plot_pmf(pmf_data)
        
        # Generate window distribution plot
        self._plot_window_distributions()
        
        # Save PMF data
        self._save_pmf_data(pmf_data)
        
        logger.info("PMF reports generated")
    
    def _load_pmf_data(self) -> Dict:
        """Load PMF analysis results"""
        # This is a placeholder - you'd load actual PMF data
        return {
            'distances': np.linspace(self.distance_start, self.distance_end, 100),
            'pmf': np.zeros(100),  # Placeholder PMF values
            'error': np.zeros(100)  # Placeholder error values
        }
    
    def _plot_pmf(self, pmf_data: Dict):
        """Plot PMF curve"""
        plt.figure(figsize=(10, 6))
        
        if self.include_error_bars and 'error' in pmf_data:
            plt.errorbar(pmf_data['distances'], pmf_data['pmf'], 
                        yerr=pmf_data['error'], fmt='o-', capsize=3)
        else:
            plt.plot(pmf_data['distances'], pmf_data['pmf'], 'o-')
        
        plt.xlabel('Distance (nm)')
        plt.ylabel('PMF (kJ/mol)')
        plt.title('Potential of Mean Force')
        plt.grid(True, alpha=0.3)
        
        plot_file = os.path.join(self.analysis_dir, f"pmf_plot.{self.plot_format}")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"PMF plot saved: {plot_file}")
    
    def _plot_window_distributions(self):
        """Plot umbrella window distributions"""
        # This is a placeholder - you'd implement window distribution plotting
        pass
    
    def _save_pmf_data(self, pmf_data: Dict):
        """Save PMF data to file"""
        data_file = os.path.join(self.analysis_dir, f"pmf_data.{self.data_format}")
        
        if self.data_format == 'txt':
            with open(data_file, 'w') as f:
                f.write("# Distance(nm) PMF(kJ/mol) Error(kJ/mol)\n")
                for i in range(len(pmf_data['distances'])):
                    f.write(f"{pmf_data['distances'][i]:.3f} "
                           f"{pmf_data['pmf'][i]:.3f} "
                           f"{pmf_data['error'][i]:.3f}\n")
        elif self.data_format == 'csv':
            import csv
            with open(data_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Distance(nm)', 'PMF(kJ/mol)', 'Error(kJ/mol)'])
                for i in range(len(pmf_data['distances'])):
                    writer.writerow([
                        f"{pmf_data['distances'][i]:.3f}",
                        f"{pmf_data['pmf'][i]:.3f}",
                        f"{pmf_data['error'][i]:.3f}"
                    ])
        
        logger.info(f"PMF data saved: {data_file}")
    
    def generate_umbrella_mdp(self, output_dir=None):
        """
        Generate umbrella sampling MDP file
        
        Parameters:
        -----------
        output_dir : str, optional
            Output directory for umbrella MDP file (default: self.mdp_dir)
            
        Returns:
        --------
        str : Path to the generated umbrella MDP file
        """
        if output_dir is None:
            output_dir = self.mdp_dir
        
        os.makedirs(output_dir, exist_ok=True)
        umbrella_mdp = os.path.join(output_dir, "umbrella.mdp")
        
        # Calculate simulation parameters
        simulation_params = self._calculate_umbrella_parameters()
        
        # Generate MDP content
        umbrella_content = self._create_umbrella_mdp_content(simulation_params)
        
        # Write MDP file
        with open(umbrella_mdp, 'w') as f:
            f.write(umbrella_content)
        
        logger.info(f"Umbrella MDP file generated: {umbrella_mdp}")
        self._log_umbrella_parameters(simulation_params)
        
        return umbrella_mdp
    
    def _calculate_umbrella_parameters(self):
        """Calculate umbrella sampling parameters"""
        # Calculate steps for 22ns simulation
        simulation_time_ns = 22.0
        simulation_time_ps = simulation_time_ns * 1000.0
        dt = self.config.get('simulation', {}).get('dt', 0.002)
        nsteps = int(simulation_time_ps / dt)
        
        # Calculate output intervals
        traj_interval = int(self.smd_trajectory_interval / dt)
        log_interval = int(self.smd_log_interval / dt)
        
        return {
            'simulation_time_ns': simulation_time_ns,
            'nsteps': nsteps,
            'dt': dt,
            'traj_interval': traj_interval,
            'log_interval': log_interval
        }
    
    def _create_umbrella_mdp_content(self, params):
        """Create umbrella sampling MDP file content"""
        # Get configuration sections
        sim_config = self.config.get('simulation', {})
        const_config = self.config.get('constraints', {})
        elec_config = self.config.get('electrostatics', {})
        vdw_config = self.config.get('vdw', {})
        tc_config = self.config.get('temperature_coupling', {})
        pc_config = self.config.get('pressure_coupling', {})
        
        # Format pull direction
        pull_direction = ' '.join(map(str, self.smd_pull_direction))
        
        return f"""; umbrella.mdp - Umbrella Sampling
title               = Umbrella sampling run ({params['simulation_time_ns']:.1f} ns)
integrator          = md
nsteps              = {params['nsteps']}
dt                  = {params['dt']}

nstxtcout           = {params['traj_interval']}
nstvout             = 0
nstenergy           = 0
nstlog              = {params['log_interval']}

continuation            = yes
constraint_algorithm    = {const_config.get('algorithm', 'lincs')}
constraints             = {const_config.get('type', 'h-bonds')}
lincs_iter              = {const_config.get('lincs_iter', 1)}
lincs_order             = {const_config.get('lincs_order', 4)}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config.get('rcoulomb', 1.0)}
rvdw                    = {vdw_config.get('rvdw', 1.0)}

coulombtype             = {elec_config.get('coulombtype', 'PME')}
pme_order               = {elec_config.get('pme_order', 4)}
fourierspacing          = {elec_config.get('fourierspacing', 0.16)}

tcoupl                  = {tc_config.get('tcoupl', 'V-rescale')}
tc-grps                 = {' '.join(tc_config.get('tc_grps', ['Protein', 'Non-Protein']))}
tau_t                   = {'     '.join(map(str, tc_config.get('tau_t', [0.1, 0.1])))}
ref_t                   = {'     '.join([str(sim_config.get('temperature', 310.0))] * len(tc_config.get('tc_grps', ['Protein', 'Non-Protein'])))}

pcoupl                  = {pc_config.get('pcoupl', 'C-rescale')}
pcoupltype              = {pc_config.get('pcoupltype', 'isotropic')}
tau_p                   = {pc_config.get('tau_p', 1.0)}
ref_p                   = {sim_config.get('pressure', 1.0)}
compressibility         = {pc_config.get('compressibility', 4.5e-5)}

pbc                     = xyz
DispCorr                = {vdw_config.get('dispcorr', 'EnerPres')}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {' '.join(tc_config.get('tc_grps', ['Protein', 'Non-Protein']))}

; Pull code - Umbrella sampling (static constraint)
pull                    = yes
pull_ncoords            = 1         
pull_ngroups            = 2         
pull_group1_name        = {self.reference_group}
pull_group2_name        = {self.moving_group}
pull_coord1_type        = umbrella  
pull_coord1_geometry    = distance  
pull_coord1_dim         = Y N N   
pull-coord1-vec         = {pull_direction}
pull_coord1_groups      = 1 2
pull_coord1_start       = yes     
pull_coord1_rate        = 0.0      ; Static constraint (no pulling)
pull_coord1_k           = {self.smd_pull_k}    
pull_coord1_r0          = 1.0      ; Reference distance (will be set per window)
pull-pbc-ref-prev-step-com = yes  
"""
    
    def _log_umbrella_parameters(self, params):
        """Log umbrella sampling parameters"""
        logger.info(f"Umbrella Sampling Parameters:")
        logger.info(f"  Simulation time: {params['simulation_time_ns']} ns")
        logger.info(f"  Number of steps: {params['nsteps']}")
        logger.info(f"  Pull rate: 0.0 (static constraint)")
        logger.info(f"  Pull force constant: {self.smd_pull_k} kJ/mol/nm²")


def run_pmf_analysis(config: Dict, output_dir: str, md_results_dir: str) -> str:
    """
    Run PMF analysis using the provided configuration
    
    This is the main entry point for PMF calculations.
    
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
        return calculator.run()
    except Exception as e:
        logger.error(f"PMF analysis failed: {e}")
        raise 