#!/usr/bin/env python3
"""
assign_chain_ids.py

Usage: python assign_chain_ids.py [filename]

This script reads a PDB file, identifies different molecules by finding breaks in residue IDs, 
and assigns chain IDs to them. It first assigns 'A', 'B', and 'C' to the first three separate 
molecules and 'X' to all the rest. A 'break' is identified when the residue ID is not the same 
or not one greater than the previous residue ID.

Prerequisites:
1. Python 3.x
2. A valid PDB file with properly sorted atoms. The script assumes that each ATOM line of the PDB 
file has the residue ID at columns 23-26 (1-based counting), and the chain ID at column 22. Please 
ensure that the ATOM lines are sorted by residue ID.

This script does not require any external Python libraries and uses only Python's built-in modules.

Written by: OpenAI's ChatGPT-4
"""

import sys

def assign_chain_ids(filename):
    # Open the PDB file
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Define variable to hold previous residue id, chain counter, and a list for the new file lines
    prev_res_id = None
    chain_counter = 0
    chains = ['A', 'B', 'C']
    new_lines = []

    for line in lines:
        # Only consider lines starting with 'ATOM'
        if line.startswith('ATOM'):
            # Extract residue id from the line
            res_id = int(line[22:26].strip())

            # If there is a break in the sequence of residue ids or the residue id decreases, a new molecule starts
            if prev_res_id is not None and (res_id != prev_res_id and res_id != prev_res_id + 1):
                chain_counter += 1

            # Assign the chain id
            chain_id = chains[chain_counter] if chain_counter < len(chains) else 'X'

            # Update the line with the new chain id
            new_line = line[:21] + chain_id + line[22:]
            new_lines.append(new_line)

            # Update the previous residue id
            prev_res_id = res_id

        # If the line does not start with 'ATOM', it is copied without changes
        else:
            new_lines.append(line)

    # Write the new file
    with open('new_' + filename, 'w') as f:
        f.writelines(new_lines)

# Get the filename from the command line arguments
if len(sys.argv) != 2:
    print("Usage: python script.py <input_PDB_file>")
    exit(1)
filename = sys.argv[1]

# Call the function with your filename
assign_chain_ids(filename)
