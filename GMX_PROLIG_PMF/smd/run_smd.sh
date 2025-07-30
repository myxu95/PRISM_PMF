#!/bin/bash
#SBATCH -J smd
#SBATCH -p quick
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

# Decide the software version
export LD_LIBRARY_PATH=/public/software/lib/:$LD_LIBRARY_PATH
source /public/software/profile.d/apps_gromacs_2023.2.sh

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

echo "Gromacs info:"
gmx --version

######################################################
# SMD SIMULATION SCRIPT
# Distance: 0.30 -> 2.00 nm
# Pull rate: 0.005 nm/ps
# Simulation time: 340.0 ps
######################################################

set -e
echo "=== SMD Simulation ==="

# Create necessary directories
mkdir -p results
mkdir -p analysis

echo "Step 1: Generating TPR file..."
gmx grompp -f smd.mdp -c md.gro -n index.ndx -p topol.top -o smd.tpr -maxwarn 10

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
