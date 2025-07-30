#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM_Builder - Complete protein-ligand system builder for GROMACS MD simulations

Enhanced version with multiple force field support (GAFF and OpenFF).

Usage:
    python prism_builder.py <protein_pdb> <ligand_file> --output <output_dir> --forcefield <gaff|openff> --config <config.yaml>

Author: Enhanced version with OpenFF support
"""

import os
import sys
import subprocess
import shutil
import re
import argparse
import yaml
from pathlib import Path

# Import force field generators
try:
    from .forcefield.gaff import GAFFForceFieldGenerator
    from .forcefield.openff import OpenFFForceFieldGenerator
    from .pmf import run_pmf_analysis
except ImportError:
    # For standalone usage
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from prism.forcefield.gaff import GAFFForceFieldGenerator
        from prism.forcefield.openff import OpenFFForceFieldGenerator
        from prism.pmf import run_pmf_analysis
    except ImportError:
        print("Error: Cannot import force field generators or PMF module")
        print("Please check your PRISM installation")
        sys.exit(1)


class PRISMBuilder:
    """Complete system builder for protein-ligand MD simulations with multiple force field support"""
    
    def __init__(self, protein_path, ligand_path, output_dir, ligand_forcefield='gaff',
                 config_path=None, forcefield=None, water_model=None, overwrite=None):
        """
        Initialize PRISM Builder with configuration support
        
        Parameters:
        -----------
        protein_path : str
            Path to the protein PDB file
        ligand_path : str
            Path to the ligand file (MOL2/SDF)
        output_dir : str
            Directory where output files will be stored
        ligand_forcefield : str
            Force field for ligand ('gaff' or 'openff')
        config_path : str, optional
            Path to configuration YAML file
        forcefield : int, optional
            Protein force field index (overrides config)
        water_model : int, optional
            Water model index (overrides config)
        overwrite : bool, optional
            Whether to overwrite existing files (overrides config)
        """
        self.protein_path = os.path.abspath(protein_path)
        self.ligand_path = os.path.abspath(ligand_path)
        self.output_dir = os.path.abspath(output_dir)
        self.ligand_forcefield = ligand_forcefield.lower()
        
        # Validate ligand force field
        if self.ligand_forcefield not in ['gaff', 'openff']:
            raise ValueError(f"Unsupported ligand force field: {self.ligand_forcefield}. Use 'gaff' or 'openff'")
        
        # Load configuration
        self.config = self._load_config(config_path)
        
        # Add ligand force field to config
        self.config['ligand_forcefield'] = {
            'type': self.ligand_forcefield,
            'charge': self.config.get('simulation', {}).get('ligand_charge', 0)
        }
        
        # Override config with explicit parameters if provided
        if forcefield is not None:
            self.config['forcefield']['index'] = forcefield
        if water_model is not None:
            self.config['water_model']['index'] = water_model
        if overwrite is not None:
            self.config['general']['overwrite'] = overwrite
        
        # Extract configuration values
        self.overwrite = self.config['general']['overwrite']
        self.forcefield_idx = self.config['forcefield']['index']
        self.water_model_idx = self.config['water_model']['index']
        
        # Get force field and water model names
        self.forcefield = self._get_forcefield_info()
        self.water_model = self._get_water_model_info()
        
        # Extract names
        self.protein_name = Path(protein_path).stem
        self.ligand_name = Path(ligand_path).stem
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Subdirectories
        self.lig_ff_dir = None
        self.mdp_dir = os.path.join(self.output_dir, "mdps")
        
        print(f"\nInitialized PRISM Builder:")
        print(f"  Protein: {self.protein_path}")
        print(f"  Ligand: {self.ligand_path}")
        print(f"  Ligand force field: {self.ligand_forcefield.upper()}")
        print(f"  Output directory: {self.output_dir}")
        print(f"  Protein force field: {self.forcefield['name']} (#{self.forcefield_idx})")
        print(f"  Water model: {self.water_model['name']} (#{self.water_model_idx})")
        print(f"  Box distance: {self.config['box']['distance']} nm")
        print(f"  Temperature: {self.config['simulation']['temperature']} K")
        print(f"  pH: {self.config['simulation']['pH']}")
        print(f"  Production time: {self.config['simulation']['production_time_ns']} ns")
    
    def _load_config(self, config_path):
        """Load configuration from YAML file or use defaults"""
        default_config = self._get_default_config()
        
        if config_path and os.path.exists(config_path):
            print(f"Loading configuration from: {config_path}")
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f)
            
            # Merge user config with defaults
            config = self._merge_configs(default_config, user_config)
        else:
            if config_path:
                print(f"Config file not found: {config_path}. Using default configuration.")
            else:
                print("No config file specified. Using default configuration.")
            config = default_config
        
        return config
    
    def _get_default_config(self):
        """Get default configuration"""
        return {
            'general': {
                'overwrite': False
            },
            'forcefield': {
                'index': 3,
                'custom_forcefields': {
                    1: {"name": "amber99sb", "dir": "amber99sb.ff"},
                    2: {"name": "amber99sb-ildn", "dir": "amber99sb-ildn.ff"},
                    3: {"name": "amber14sb_OL15", "dir": "amber14sb_OL15_cufix_zn.ff"},
                    4: {"name": "amber03", "dir": "amber03.ff"},
                    5: {"name": "amber96", "dir": "amber96.ff"},
                    6: {"name": "amber94", "dir": "amber94.ff"},
                    7: {"name": "charmm27", "dir": "charmm27.ff"},
                    8: {"name": "oplsaa", "dir": "oplsaa.ff"}
                }
            },
            'water_model': {
                'index': 1,
                'custom_water_models': {
                    1: {"name": "tip3p"},
                    2: {"name": "tip4p"},
                    3: {"name": "spc"},
                    4: {"name": "spce"},
                    5: {"name": "none"}
                }
            },
            'box': {
                'distance': 1.5,  # nm
                'shape': 'cubic',  # cubic, dodecahedron, octahedron
                'center': True
            },
            'simulation': {
                'temperature': 310,  # K
                'pressure': 1.0,  # bar
                'pH': 7.0,
                'ligand_charge': 0,  # Default ligand charge
                'production_time_ns': 500,  # ns
                'dt': 0.002,  # ps
                'equilibration_nvt_time_ps': 500,  # ps
                'equilibration_npt_time_ps': 500  # ps
            },
            'ions': {
                'neutral': True,
                'concentration': 0.15,  # M
                'positive_ion': 'NA',
                'negative_ion': 'CL'
            },
            'constraints': {
                'algorithm': 'lincs',
                'type': 'h-bonds',  # none, h-bonds, all-bonds
                'lincs_iter': 1,
                'lincs_order': 4
            },
            'energy_minimization': {
                'integrator': 'steep',
                'emtol': 200.0,  # kJ/mol/nm
                'emstep': 0.01,
                'nsteps': 10000
            },
            'output': {
                'trajectory_interval_ps': 500,  # ps
                'energy_interval_ps': 10,  # ps
                'log_interval_ps': 10,  # ps
                'compressed_trajectory': True
            },
            'electrostatics': {
                'coulombtype': 'PME',
                'rcoulomb': 1.0,  # nm
                'pme_order': 4,
                'fourierspacing': 0.16  # nm
            },
            'vdw': {
                'rvdw': 1.0,  # nm
                'dispcorr': 'EnerPres'
            },
            'temperature_coupling': {
                'tcoupl': 'V-rescale',
                'tc_grps': ['Protein', 'Non-Protein'],
                'tau_t': [0.1, 0.1],  # ps
            },
            'pressure_coupling': {
                'pcoupl': 'C-rescale',
                'pcoupltype': 'isotropic',
                'tau_p': 1.0,  # ps
                'compressibility': 4.5e-5  # bar^-1
            },
            'pmf': {
                'enabled': False,
                'mdp_path': 'md.mdp',
                'tpr_path': 'md.tpr',
                'output_dir': 'pmf_results',
                'temperature': 310,
                'pressure': 1.0,
                'pH': 7.0,
                'dt': 0.002,
                'production_time_ns': 500,
                'trajectory_interval_ps': 500,
                'energy_interval_ps': 10,
                'log_interval_ps': 10,
                'coulombtype': 'PME',
                'rcoulomb': 1.0,
                'pme_order': 4,
                'fourierspacing': 0.16,
                'vdw': {
                    'rvdw': 1.0,
                    'dispcorr': 'EnerPres'
                },
                'tcoupl': 'V-rescale',
                'tc_grps': ['Protein', 'Non-Protein'],
                'tau_t': [0.1, 0.1],
                'pcoupl': 'C-rescale',
                'pcoupltype': 'isotropic',
                'tau_p': 1.0,
                'compressibility': 4.5e-5
            }
        }
    
    def _merge_configs(self, default, user):
        """Recursively merge user config with default config"""
        if isinstance(default, dict):
            if not isinstance(user, dict):
                return default
            merged = default.copy()
            for key, value in user.items():
                if key in merged:
                    merged[key] = self._merge_configs(merged[key], value)
                else:
                    merged[key] = value
            return merged
        else:
            return user
    
    def _get_forcefield_info(self):
        """Get force field information from config"""
        ff_idx = self.config['forcefield']['index']
        custom_ff = self.config['forcefield']['custom_forcefields']
        
        if ff_idx in custom_ff:
            return custom_ff[ff_idx]
        else:
            raise ValueError(f"Force field index {ff_idx} not found in configuration")
    
    def _get_water_model_info(self):
        """Get water model information from config"""
        wm_idx = self.config['water_model']['index']
        custom_wm = self.config['water_model']['custom_water_models']
        
        if wm_idx in custom_wm:
            return custom_wm[wm_idx]
        else:
            raise ValueError(f"Water model index {wm_idx} not found in configuration")
    
    def run_command(self, cmd, input_text=None, cwd=None, check=True, shell=False):
        """Run a shell command with optional input"""
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        print(f"Running: {cmd_str}")
        
        try:
            if input_text is not None:
                if shell:
                    process = subprocess.Popen(
                        cmd_str,
                        cwd=cwd,
                        shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                else:
                    process = subprocess.Popen(
                        cmd,
                        cwd=cwd,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                stdout, stderr = process.communicate(input=input_text)
                
                if process.returncode != 0 and check:
                    print(f"Error: {stderr}")
                    raise subprocess.CalledProcessError(process.returncode, cmd, stdout, stderr)
                
                return stdout
            else:
                if shell:
                    result = subprocess.run(
                        cmd_str,
                        cwd=cwd,
                        shell=True,
                        check=check,
                        capture_output=True,
                        text=True
                    )
                else:
                    result = subprocess.run(
                        cmd,
                        cwd=cwd,
                        check=check,
                        capture_output=True,
                        text=True
                    )
                return result.stdout
                
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {cmd_str}")
            print(f"Error: {e.stderr}")
            if check:
                raise
            return ""
    
    def generate_ligand_forcefield(self):
        """Generate ligand force field using selected force field generator"""
        print(f"\n=== Generating Ligand Force Field ({self.ligand_forcefield.upper()}) ===")
        
        if self.ligand_forcefield == 'gaff':
            # Use GAFF force field generator
            generator = GAFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()
            
        elif self.ligand_forcefield == 'openff':
            # Use OpenFF force field generator
            generator = OpenFFForceFieldGenerator(
                self.ligand_path,
                self.output_dir,
                charge=self.config['ligand_forcefield']['charge'],
                overwrite=self.overwrite
            )
            self.lig_ff_dir = generator.run()
        
        # Verify required files exist
        required_files = ["LIG.itp", "LIG.gro", "atomtypes_LIG.itp"]
        for filename in required_files:
            filepath = os.path.join(self.lig_ff_dir, filename)
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Required file not found: {filepath}")
        
        print(f"Ligand force field files generated in: {self.lig_ff_dir}")
        return self.lig_ff_dir
    
    def clean_protein(self):
        """Clean the protein PDB file"""
        print("\n=== Cleaning Protein ===")
        
        cleaned_pdb = os.path.join(self.output_dir, f"{self.protein_name}_clean.pdb")
        
        if os.path.exists(cleaned_pdb) and not self.overwrite:
            print(f"Using existing cleaned protein: {cleaned_pdb}")
            return cleaned_pdb
        
        with open(self.protein_path, 'r') as f_in, open(cleaned_pdb, 'w') as f_out:
            for line in f_in:
                if not line.startswith('HETATM'):
                    # Fix terminal hydrogen names for AMBER compatibility
                    if 'HN1' in line:
                        line = line.replace('HN1', 'H1 ')
                    if 'HN2' in line:
                        line = line.replace('HN2', 'H2 ')
                    if 'HN3' in line:
                        line = line.replace('HN3', 'H3 ')
                    f_out.write(line)
        
        print(f"Protein cleaned and saved to: {cleaned_pdb}")
        return cleaned_pdb
    
    def build_model(self, cleaned_protein):
        """Build the GROMACS model"""
        print("\n=== Building GROMACS Model ===")
        
        model_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
        
        if os.path.exists(os.path.join(model_dir, "solv_ions.gro")) and not self.overwrite:
            print(f"Final model file solv_ions.gro already exists in {model_dir}, skipping model building")
            return model_dir
            
        os.makedirs(model_dir, exist_ok=True)
        
        # Step 1: Run pdbfixer to fix missing atoms
        fixed_pdb = os.path.join(model_dir, "fixed_clean_protein.pdb")
        if not os.path.exists(fixed_pdb) or self.overwrite:
            print("Running pdbfixer to fix missing atoms...")
            self.run_command([
                "pdbfixer", 
                cleaned_protein, 
                "--output", fixed_pdb, 
                "--add-atoms", "heavy", 
                "--keep-heterogens", "none",
                f"--ph={self.config['simulation']['pH']}"
            ])
        else:
            print(f"Fixed protein PDB already exists at {fixed_pdb}, skipping pdbfixer")
        
        # Step 2: Generate topology using gmx pdb2gmx
        current_dir = os.getcwd()
        os.chdir(model_dir)
        
        pro_gro = os.path.join(model_dir, "pro.gro")
        if not os.path.exists(pro_gro) or self.overwrite:
            histidine_count = self._count_histidines(fixed_pdb)
            
            input_lines = [str(self.forcefield_idx), str(self.water_model_idx)]
            
            if histidine_count > 0:
                print(f"Found {histidine_count} histidine residue(s), defaulting to HIE")
                for _ in range(histidine_count):
                    input_lines.append("1")  # HIE
            
            input_text = '\n'.join(input_lines) + '\n'
            
            print(f"Running pdb2gmx with force field #{self.forcefield_idx} and water model #{self.water_model_idx}")
            self.run_command(
                ["gmx", "pdb2gmx", "-f", "fixed_clean_protein.pdb", "-o", "pro.gro", "-ignh"],
                input_text=input_text
            )
        else:
            print(f"Protein topology already generated at {pro_gro}, skipping pdb2gmx")
        
        # Step 3: Combine protein and ligand coordinates
        pro_lig_gro = os.path.join(model_dir, "pro_lig.gro")
        if not os.path.exists(pro_lig_gro) or self.overwrite:
            print("Combining protein and ligand coordinates...")
            
            lig_gro = os.path.join(self.lig_ff_dir, "LIG.gro")
            with open(lig_gro, 'r') as f:
                lines = f.readlines()
                lig_atom_count = int(lines[1].strip())
                lig_coords = ''.join(lines[2:-1])
            
            with open(pro_gro, 'r') as f:
                lines = f.readlines()
                protein_atom_count = int(lines[1].strip())
                gro_header = lines[0]
                gro_box = lines[-1]
                protein_coords = ''.join(lines[2:-1])
            
            total_atom_count = protein_atom_count + lig_atom_count
            with open(pro_lig_gro, 'w') as f:
                f.write(gro_header)
                f.write(f"{total_atom_count}\n")
                f.write(protein_coords)
                f.write(lig_coords)
                f.write(gro_box)
            
            print("Protein-ligand complex created")
        
        # Step 4: Fix topology file
        topol_path = os.path.join(model_dir, "topol.top")
        if os.path.exists(topol_path):
            print("Fixing topology file...")
            self._fix_topology(topol_path)
        
        # Step 5: Create box
        box_gro = os.path.join(model_dir, "pro_lig_newbox.gro")
        if not os.path.exists(box_gro) or self.overwrite:
            box_cmd = [
                "gmx", "editconf",
                "-f", "pro_lig.gro",
                "-o", "pro_lig_newbox.gro"
            ]
            
            if self.config['box']['center']:
                box_cmd.append("-c")
            
            # Add box shape
            if self.config['box']['shape'] == 'dodecahedron':
                box_cmd.extend(["-bt", "dodecahedron"])
            elif self.config['box']['shape'] == 'octahedron':
                box_cmd.extend(["-bt", "octahedron"])
            
            # Add distance
            box_cmd.extend(["-d", str(self.config['box']['distance'])])
            
            self.run_command(box_cmd)
        
        # Step 6: Solvate
        solv_gro = os.path.join(model_dir, "solv.gro")
        if not os.path.exists(solv_gro) or self.overwrite:
            self.run_command([
                "gmx", "solvate",
                "-cp", "pro_lig_newbox.gro",
                "-p", "topol.top",
                "-o", "solv.gro"
            ])
        
        # Step 7: Add ions
        ions_mdp = os.path.join(model_dir, "ions.mdp")
        if not os.path.exists(ions_mdp) or self.overwrite:
            self._create_ions_mdp(ions_mdp)
        
        ions_tpr = os.path.join(model_dir, "ions.tpr")
        if not os.path.exists(ions_tpr) or self.overwrite:
            self.run_command([
                "gmx", "grompp",
                "-f", ions_mdp,
                "-c", "solv.gro",
                "-p", "topol.top",
                "-o", "ions.tpr",
                "-maxwarn", "99999"
            ])
        
        solv_ions_gro = os.path.join(model_dir, "solv_ions.gro")
        if not os.path.exists(solv_ions_gro) or self.overwrite:
            print("Adding ions...")
            genion_cmd = [
                "gmx", "genion", "-s", "ions.tpr", "-o", "solv_ions.gro", 
                "-p", "topol.top", 
                "-pname", self.config['ions']['positive_ion'], 
                "-nname", self.config['ions']['negative_ion']
            ]
            
            if self.config['ions']['neutral']:
                genion_cmd.append("-neutral")
            
            if self.config['ions']['concentration'] > 0:
                genion_cmd.extend(["-conc", str(self.config['ions']['concentration'])])
            
            # Try different group numbers for SOL
            for group in ["13", "14", "15"]:
                try:
                    self.run_command(genion_cmd, input_text=group)
                    print(f"Successfully added ions using group {group}")
                    break
                except:
                    if group == "15":
                        genion_cmd_str = ' '.join(genion_cmd)
                        self.run_command(f"echo {group} | {genion_cmd_str}", shell=True)
        
        os.chdir(current_dir)
        print(f"\nModel build completed. System files are in {model_dir}")
        return model_dir
    
    def _count_histidines(self, pdb_file):
        """Count histidine residues in PDB file"""
        his_count = 0
        his_residues = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') and line[17:20].strip() in ['HIS', 'HID', 'HIE', 'HIP']:
                    resnum = line[22:26].strip()
                    chain = line[21].strip()
                    res_id = f"{chain}:{resnum}" if chain else resnum
                    
                    if res_id not in his_residues:
                        his_residues.append(res_id)
                        his_count += 1
        
        if his_count > 0:
            print(f"Histidine residues found: {', '.join(his_residues)}")
        
        return his_count
    
    def _fix_topology(self, topol_path):
        """Fix topology file to include ligand parameters in correct order"""
        with open(topol_path, 'r') as f:
            content = f.read()
        
        # Get the correct directory name based on force field
        lig_dir_name = os.path.basename(self.lig_ff_dir)
        
        if f"{lig_dir_name}/LIG.itp" in content:
            print("Topology already includes ligand parameters")
            return
        
        atomtypes_path = os.path.join(self.lig_ff_dir, "atomtypes_LIG.itp")
        with open(atomtypes_path, 'r') as f:
            atomtypes_content = f.read()
        
        new_content = []
        lines = content.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i]
            new_content.append(line)
            
            if '#include' in line and '.ff/forcefield.itp' in line:
                new_content.append("")
                new_content.append("; Ligand atomtypes")
                new_content.extend(atomtypes_content.strip().split('\n'))
                new_content.append("")
                new_content.append("; Ligand topology")
                new_content.append(f'#include "../{lig_dir_name}/LIG.itp"')
                new_content.append("")
            
            if '[ molecules ]' in line:
                j = i + 1
                while j < len(lines) and lines[j].strip() and not lines[j].startswith('['):
                    new_content.append(lines[j])
                    j += 1
                new_content.append("LIG                 1")
                i = j - 1
            
            i += 1
        
        with open(topol_path, 'w') as f:
            f.write('\n'.join(new_content))
        
        print("Topology file updated with ligand parameters")
    
    def _create_ions_mdp(self, mdp_path):
        """Create ions.mdp file using config parameters"""
        em_config = self.config['energy_minimization']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']
        
        content = f"""; ions.mdp - Parameters for adding ions
integrator  = {em_config['integrator']}
emtol       = {em_config['emtol']}
emstep      = {em_config['emstep']}
nsteps      = {em_config['nsteps']}

; Output
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = {elec_config['rcoulomb']}
rvdw            = {vdw_config['rvdw']}

; Electrostatics
coulombtype     = {elec_config['coulombtype']}
pme_order       = {elec_config['pme_order']}
fourierspacing  = {elec_config['fourierspacing']}

; Temperature and pressure coupling are off
tcoupl      = no
pcoupl      = no

; Periodic boundary conditions
pbc         = xyz

; Constraints
constraints             = {self.config['constraints']['type']}
constraint-algorithm    = {self.config['constraints']['algorithm']}
"""
        
        with open(mdp_path, 'w') as f:
            f.write(content)
    
    def generate_mdp_files(self):
        """Generate MDP files for MD simulations using config parameters"""
        print("\n=== Generating MDP Files ===")
        
        os.makedirs(self.mdp_dir, exist_ok=True)
        
        # Get config sections
        sim_config = self.config['simulation']
        em_config = self.config['energy_minimization']
        const_config = self.config['constraints']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']
        tc_config = self.config['temperature_coupling']
        pc_config = self.config['pressure_coupling']
        output_config = self.config['output']
        
        # Calculate steps
        nvt_steps = int(sim_config['equilibration_nvt_time_ps'] / sim_config['dt'])
        npt_steps = int(sim_config['equilibration_npt_time_ps'] / sim_config['dt'])
        prod_steps = int(sim_config['production_time_ns'] * 1000 / sim_config['dt'])
        
        # Calculate output intervals
        energy_interval = int(output_config['energy_interval_ps'] / sim_config['dt'])
        log_interval = int(output_config['log_interval_ps'] / sim_config['dt'])
        traj_interval = int(output_config['trajectory_interval_ps'] / sim_config['dt'])
        
        # Energy minimization
        em_mdp = os.path.join(self.mdp_dir, "em.mdp")
        with open(em_mdp, 'w') as f:
            f.write(f"""; em.mdp - Energy minimization
title           = Energy minimization
define          = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator      = {em_config['integrator']}
emtol           = {em_config['emtol']}
emstep          = {em_config['emstep']}
nsteps          = {em_config['nsteps']}

nstxout         = 0
nstvout         = 0
nstenergy       = 0
nstlog          = 0
nstxout-compressed = 0

cutoff-scheme   = Verlet
coulombtype     = {elec_config['coulombtype']}
rcoulomb        = {elec_config['rcoulomb']}
rvdw            = {vdw_config['rvdw']}
pbc             = xyz
""")
        
        # NVT equilibration
        nvt_mdp = os.path.join(self.mdp_dir, "nvt.mdp")
        with open(nvt_mdp, 'w') as f:
            tc_grps = ' '.join(tc_config['tc_grps'])
            tau_t = '     '.join(map(str, tc_config['tau_t']))
            ref_t = '     '.join([str(sim_config['temperature'])] * len(tc_config['tc_grps']))
            
            f.write(f"""; nvt.mdp - NVT equilibration
title               = NVT equilibration with the complex restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {nvt_steps}
dt                  = {sim_config['dt']}

nstxout             = 0
nstvout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

continuation            = no
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = no
pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}

gen_vel                 = yes
gen_temp                = {sim_config['temperature']}
gen_seed                = -1

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")
        
        # NPT equilibration
        npt_mdp = os.path.join(self.mdp_dir, "npt.mdp")
        with open(npt_mdp, 'w') as f:
            f.write(f"""; npt.mdp - NPT equilibration
title               = NPT equilibration with protein restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {npt_steps}
dt                  = {sim_config['dt']}

nstxout             = 0
nstvout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")
        
        # Production MD
        md_mdp = os.path.join(self.mdp_dir, "md.mdp")
        with open(md_mdp, 'w') as f:
            output_line = "nstxtcout" if output_config['compressed_trajectory'] else "nstxout"
            
            f.write(f"""; md.mdp - Production MD
title               = Production run ({sim_config['production_time_ns']} ns)
integrator          = md
nsteps              = {prod_steps}
dt                  = {sim_config['dt']}

{output_line}       = {traj_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {log_interval}

continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")
        
        print(f"MDP files generated in {self.mdp_dir}")
        print(f"  - Energy minimization: {em_config['nsteps']} steps")
        print(f"  - NVT equilibration: {sim_config['equilibration_nvt_time_ps']} ps")
        print(f"  - NPT equilibration: {sim_config['equilibration_npt_time_ps']} ps")
        print(f"  - Production: {sim_config['production_time_ns']} ns")
        
        # Generate SMD MDP file if PMF is enabled
        if self.config.get('pmf', {}).get('enabled', False):
            pmf_config = self.config['pmf']
            print("\n=== Generating SMD MDP File ===")
            
            # Create PMF calculator to generate SMD MDP
            from .pmf import PMFCalculator
            pmf_calc = PMFCalculator(pmf_config, self.output_dir, "dummy_md_dir")
            smd_mdp_path = pmf_calc.generate_smd_mdp(self.mdp_dir)
            print(f"SMD MDP file generated: {smd_mdp_path}")
    
    def cleanup(self):
        """Clean up temporary files"""
        print("\n=== Cleaning up temporary files ===")
        
        # Cleanup directories for both GAFF and OpenFF
        cleanup_dirs = ["forcefield", "temp_openff"]
        
        for dir_name in cleanup_dirs:
            cleanup_dir = os.path.join(self.output_dir, dir_name)
            if os.path.exists(cleanup_dir):
                if self.ligand_forcefield == 'gaff':
                    temp_patterns = ["*.frcmod", "*.prep", "*.prmtop", "*.rst7", "*.log", "*.in",
                                   "ANTECHAMBER*", "ATOMTYPE*", "PREP*", "NEWPDB*", "sqm*", "leap*"]
                    
                    for pattern in temp_patterns:
                        for file_path in Path(cleanup_dir).glob(pattern):
                            try:
                                os.remove(file_path)
                            except:
                                pass
                else:
                    # For OpenFF, we might want to keep the directory clean
                    # or remove it entirely if it's temporary
                    try:
                        shutil.rmtree(cleanup_dir)
                    except:
                        pass
        
        print("Cleanup completed")
    
    def save_config(self):
        """Save the current configuration to a file"""
        config_file = os.path.join(self.output_dir, "prism_config.yaml")
        with open(config_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)
        print(f"Configuration saved to: {config_file}")
    
    def run(self):
        """Run the complete workflow"""
        print(f"\n{'='*60}")
        print("Starting PRISM Builder workflow")
        print(f"{'='*60}")
        
        try:
            # Save configuration for reference
            self.save_config()
            
            # Step 1: Generate ligand force field
            self.generate_ligand_forcefield()
            
            # Step 2: Clean protein
            cleaned_protein = self.clean_protein()
            
            # Step 3: Build model
            model_dir = self.build_model(cleaned_protein)
            if not model_dir:
                raise RuntimeError("Failed to build model")
            
            # Step 4: Generate MDP files
            self.generate_mdp_files()
            
            # Step 5: Cleanup
            self.cleanup()
            
            # Step 6: PMF Analysis (if enabled)
            pmf_results_dir = None
            if self.config.get('pmf', {}).get('enabled', False):
                print(f"\n{'='*60}")
                print("PMF Analysis Enabled - Starting PMF Calculation")
                print(f"{'='*60}")
                
                md_results_dir = os.path.join(self.output_dir, "GMX_PROLIG_MD")
                pmf_results_dir = run_pmf_analysis(
                    self.config['pmf'], 
                    self.output_dir, 
                    md_results_dir
                )
                
                if pmf_results_dir:
                    print(f"PMF analysis completed. Results in: {pmf_results_dir}")
                else:
                    print("PMF analysis failed or was disabled")
            
            print(f"\n{'='*60}")
            print("PRISM Builder workflow completed successfully!")
            print(f"{'='*60}")
            print(f"\nOutput files are in: {self.output_dir}")
            print(f"MD system files are in: {os.path.join(self.output_dir, 'GMX_PROLIG_MD')}")
            print(f"MDP files are in: {self.mdp_dir}")
            print(f"Configuration saved in: {os.path.join(self.output_dir, 'prism_config.yaml')}")
            if pmf_results_dir:
                print(f"PMF results are in: {pmf_results_dir}")
            print(f"\nProtein force field used: {self.forcefield['name']}")
            print(f"Ligand force field used: {self.ligand_forcefield.upper()}")
            print(f"Water model used: {self.water_model['name']}")
            print(f"\nYou can now run MD simulations using the generated files.")
            
            return self.output_dir
            
        except Exception as e:
            print(f"\nError during PRISM Builder workflow: {e}")
            import traceback
            traceback.print_exc()
            raise


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="PRISM Builder - Build protein-ligand systems for GROMACS with multiple force field support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Using GAFF force field (default)
  python prism_builder.py protein.pdb ligand.mol2 -o output_dir
  
  # Using OpenFF force field
  python prism_builder.py protein.pdb ligand.sdf -o output_dir --ligand-forcefield openff
  
  # With custom configuration
  python prism_builder.py protein.pdb ligand.mol2 -o output_dir --config config.yaml
  
  # Override specific parameters
  python prism_builder.py protein.pdb ligand.mol2 -o output_dir --forcefield 3 --water 1 --ligand-forcefield gaff
        """
    )
    
    parser.add_argument("protein", help="Path to protein PDB file")
    parser.add_argument("ligand", help="Path to ligand file (MOL2/SDF)")
    parser.add_argument("--output", "-o", default="prism_output",
                       help="Output directory (default: prism_output)")
    parser.add_argument("--ligand-forcefield", "-lff", choices=['gaff', 'openff'], default='gaff',
                       help="Force field for ligand (default: gaff)")
    parser.add_argument("--config", "-c", 
                       help="Path to configuration YAML file")
    parser.add_argument("--forcefield", "-ff", type=int,
                       help="Protein force field index (overrides config)")
    parser.add_argument("--water", "-w", type=int,
                       help="Water model index (overrides config)")
    parser.add_argument("--overwrite", "-f", action="store_true",
                       help="Overwrite existing files (overrides config)")
    
    args = parser.parse_args()
    
    # Create and run PRISM Builder
    builder = PRISMBuilder(
        args.protein,
        args.ligand,
        args.output,
        ligand_forcefield=args.ligand_forcefield,
        config_path=args.config,
        forcefield=args.forcefield,
        water_model=args.water,
        overwrite=args.overwrite
    )
    builder.run()


if __name__ == "__main__":
    main()
    