#!/bin/bash

# Set variables
PREFIX=
TPR_FILE="md.tpr"
INDEX_FILE="index.ndx"
OUTPUT_XTC="${PREFIX}.rewrapped.pf1000ps.xtc"
OUTPUT_XVG="${PREFIX}.rmsd.xvg"
GMX="gmx"

# Create an index file with molecule groups
# $GMX make_ndx -s "$TPR_FILE" -n "$INDEX_FILE"

# Process the trajectory by adjusting PBC, centering it and saving every 1000th frame
$GMX trjconv -s "$TPR_FILE" -f "${PREFIX}.xtc" -n "$INDEX_FILE" -pbc mol -center -ur compact -dt 1000 -o "step1-$OUTPUT_XTC"
$GMX trjconv -s "$TPR_FILE" -f "step1-$OUTPUT_XTC" -n "$INDEX_FILE" -fit rot+trans -o "$OUTPUT_XTC"

# Calculate RMSD for the trajectory
$GMX rms -s "$TPR_FILE" -f "${OUTPUT_XTC}" -n "$INDEX_FILE" -o "$OUTPUT_XVG"

rm "step1-$OUTPUT_XTC"
