#!/bin/bash
######################################################
# SMD PREPARATION SCRIPT
######################################################

set -e
echo "=== SMD Preparation ==="

cd /public/home/xmy/PRISM/test/PRISM/prism/pmf_test/pmf_output/pmf

# Step 1: Create index file
if [ ! -f index.ndx ]; then
    echo "Creating index file..."
    echo "r LIG\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx
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
