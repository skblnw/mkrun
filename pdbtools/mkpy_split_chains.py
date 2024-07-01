import os
import sys
from pathlib import Path

# Variables
PREFIX = "ionized"
INPUT_FILE = f"{PREFIX}.gro"
OUTPUT_FILE = "output_with_chains.pdb"
OUTPUT_DIR = "chains"

# Convert gro to pdb if needed
if INPUT_FILE.endswith(".gro"):
    os.system(f"gmx editconf -f {INPUT_FILE} -o {PREFIX}.pdb")
    INPUT_FILE = f"{PREFIX}.pdb"

def increment_chain(chain):
    return chr(ord(chain) + 1)

# Assign chain identifiers to repeating residue number ranges
with open(INPUT_FILE, "r") as inp, open(OUTPUT_FILE, "w") as out:
    current_chain = "A"
    last_residue = -1

    for line in inp:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue = int(line[22:26].strip())
            resname = line[17:20].strip()

            if resname == "SOL":
                break
            elif residue < last_residue:
                current_chain = increment_chain(current_chain)

            line = line[:21] + current_chain + line[22:]
            last_residue = residue
            out.write(line)

# Create output directory
Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Split the output pdb file into separate files based on chain names
with open(OUTPUT_FILE, "r") as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            output_file = f"{OUTPUT_DIR}/chain_{chain}.pdb"
            with open(output_file, "a") as out_f:
                out_f.write(line)

