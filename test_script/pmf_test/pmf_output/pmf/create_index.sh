#!/bin/bash
# Create index file for SMD groups
echo "Creating index file for groups..."

cd /public/home/xmy/PRISM/test/PRISM/prism/pmf_test/pmf_output/pmf

# Create index file
echo "r LIG\nq" | gmx make_ndx -f solv_ions.gro -o index.ndx

echo "✓ Index file created: index.ndx"
