#!/usr/bin/env python3
"""
This script processes PDB and PSF files to manipulate the chain information. It provides options to:
- Read one PDB file or a combination of PSF and PDB files.
- Include only specific chains in the output PDB file.

Usage:
- For a single PDB file: python script.py input.pdb
- For PSF and PDB files: python script.py input.pdb --psf input.psf
- To include specific chains: python script.py input.pdb --chains "A B C"
"""

import sys
import os
import subprocess
import argparse

def main(args):

    # Assign input arguments
    PDB = args.pdb
    PSF = args.psf
    chains = args.chains

    # Validation for file existence
    if not os.path.isfile(PDB):
        print(f"{PDB} \nStructure not found!")
        sys.exit(0)

    # Write initial TCL commands to load the structure
    with open("tcl", "w") as file:
        print(f"Loading {PDB}...")
        file.write(f"mol new {PDB} waitfor all\n")
        if PSF:
            print(f"Adding {PSF}...")
            file.write(f"mol addfile {PSF} waitfor all\n")

        # TCL commands to manipulate chain information
        tcl_commands = """
set alphabets {0 A B C D E F G H I J K L M N O P Q R S T U V W X Y Z}
set selall [atomselect top protein]
set total_num_of_chains [llength [lsort -unique [$selall get segname]]]
if { $total_num_of_chains == 1 } {
    set total_num_of_chains [llength [lsort -unique [$selall get chain]]]
    set chainlist [lsort -unique [$selall get chain]]
} else {
    set chainlist [lsort -unique [$selall get segname]]
}
foreach ii $chainlist cc [lrange $alphabets 1 $total_num_of_chains] {
    set sel [atomselect top "segname $ii"]
    $sel set chain $cc
}
set sel [atomselect top "resname HSD HSE HSP"]; $sel set resname HIS
set sel [atomselect top "resname ILE and name CD"]; $sel set name CD1
set sel [atomselect top "name OT1"]; $sel set name O
set sel [atomselect top "name OT2"]; $sel set name OXT
set sel [atomselect top "resname TIP"]; $sel set resname SOL
set sel [atomselect top "resname SOL and name H1"]; $sel set name HW1
set sel [atomselect top "resname SOL and name H2"]; $sel set name HW2
set sel [atomselect top "resname SOL and name OH2"]; $sel set name OW
set sel [atomselect top "all"]
$sel writepdb temp.pdb
quit
"""

        file.write(tcl_commands)

    print("Running VMD command...")
    subprocess.run(["vmd", "-dispdev", "text", "-e", "tcl"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Optionally filter the chains in the output PDB file
    if chains:
        print(f"Filtering chains: {chains}...")
        with open("temp.pdb", "r") as infile, open("fixed.pdb", "w") as outfile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21]
                    if chain_id in chains.split():
                        outfile.write(line)
        os.remove("temp.pdb")
    else:
        os.rename("temp.pdb", "fixed.pdb")

    os.remove("tcl")
    print("Script completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processing PDB and PSF files.")
    parser.add_argument("pdb", type=str, help="Path to the PDB file.")
    parser.add_argument("--psf", type=str, help="Path to the PSF file (optional).")
    parser.add_argument("--chains", type=str, help="Chains to include in the output, separated by spaces (e.g., 'A B C').")

    args = parser.parse_args()
    main(args)
