#!/bin/bash
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
