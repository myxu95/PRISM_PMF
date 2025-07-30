# # prism/__init__.py
# """
# PRISM - Protein Receptor Interaction Simulation Modeler

# A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
# """

# __version__ = "1.0.0"
# __author__ = "PRISM Development Team"

# from .builder import PRISMBuilder

# __all__ = ["PRISMBuilder"]

# # ===========================
# # prism/forcefield/__init__.py
# """
# Force field generators for PRISM
# """

# from forcefield.base import ForceFieldGeneratorBase
# from forcefield.gaff import GAFFForceFieldGenerator
# from forcefield.openff import OpenFFForceFieldGenerator

# __all__ = ["ForceFieldGeneratorBase", "GAFFForceFieldGenerator", "OpenFFForceFieldGenerator"]

# # ===========================
# # prism/utils/__init__.py
# """
# Utility functions for PRISM
# """

# # Future utilities can be added here


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder
from .core import PRISMSystem
from .pmf import PMFCalculator, run_pmf_analysis

# High-level API functions
def system(protein_path, ligand_path, config=None, **kwargs):
    """
    Create a protein-ligand system for MD simulation.
    
    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str  
        Path to ligand file (MOL2/SDF)
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)
    
    Returns:
    --------
    PRISMSystem
        A system object with methods to build and run simulations
        
    Examples:
    ---------
    >>> import prism as pm
    >>> system = pm.system("protein.pdb", "ligand.mol2")
    >>> system.build()
    >>> system.generate_mdp_files()
    
    >>> # With custom configuration
    >>> system = pm.system("protein.pdb", "ligand.sdf", 
    ...                    config="config.yaml",
    ...                    ligand_forcefield="openff")
    >>> system.build()
    """
    return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
    """
    Build a complete protein-ligand system (one-step function).
    
    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)  
    output_dir : str
        Output directory for generated files
    **kwargs : optional
        Additional parameters
        
    Returns:
    --------
    str
        Path to the output directory
        
    Examples:
    ---------
    >>> import prism as pm
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    return system_obj.build()

def run_pmf(protein_path, ligand_path, output_dir="prism_output", pmf_config=None, **kwargs):
    """
    Build a system and run PMF analysis (one-step function).
    
    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)
    output_dir : str
        Output directory for generated files
    pmf_config : dict, optional
        PMF configuration dictionary
    **kwargs : optional
        Additional parameters
        
    Returns:
    --------
    str
        Path to the PMF results directory
        
    Examples:
    ---------
    >>> import prism as pm
    >>> pmf_results = pm.run_pmf("protein.pdb", "ligand.mol2")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    system_obj.build()
    return system_obj.run_pmf_analysis(pmf_config)

def analyze_pmf(md_results_dir, output_dir="pmf_results", pmf_config=None):
    """
    Run PMF analysis on existing MD results.
    
    Parameters:
    -----------
    md_results_dir : str
        Directory containing MD simulation results
    output_dir : str
        Output directory for PMF results
    pmf_config : dict, optional
        PMF configuration dictionary
        
    Returns:
    --------
    str
        Path to the PMF results directory
        
    Examples:
    ---------
    >>> import prism as pm
    >>> pmf_results = pm.analyze_pmf("md_results", "pmf_output")
    """
    if pmf_config is None:
        # Default PMF configuration
        pmf_config = {
            'enabled': True,
            'method': 'umbrella',
            'reaction_coordinate': 'distance',
            'reference_group': 'Protein',
            'moving_group': 'LIG',
            'distance': {
                'start': 0.3,
                'end': 2.0,
                'step': 0.1,
                'force_constant': 1000.0
            },
            'umbrella': {
                'n_windows': 20,
                'equilibration_time_ps': 1000,
                'production_time_ps': 5000,
                'sampling_interval_ps': 10
            },
            'analysis': {
                'method': 'wham',
                'temperature': 310.0,
                'tolerance': 1e-6,
                'max_iterations': 50000
            },
            'output': {
                'plot_format': 'png',
                'data_format': 'txt',
                'include_error_bars': True
            }
        }
    
    return run_pmf_analysis(pmf_config, output_dir, md_results_dir)

# Export main classes and functions
__all__ = [
    "PRISMBuilder", 
    "PRISMSystem",
    "PMFCalculator",
    "system", 
    "build_system",
    "run_pmf",
    "analyze_pmf",
    "run_pmf_analysis",
    "__version__"
]