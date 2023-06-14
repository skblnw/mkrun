"""
# PDB Chain Identifier and Atom Name Fixer

This script takes an input PDB file and performs the following operations:

1. Assigns chain identifiers based on repeating residue number ranges.
2. Skips specified resnames: TIP3, SOD, POT, and CLA.
3. Fixes some atom names:
   - Changes HSD, HSE, and HSP to HIS.
   - Changes ILE CD to CD1.
   - Changes OT1 to O.
   - Changes OT2 to OXT.
4. Skips all atom names starting with "H".

## Usage

To run the script, use the following command:
```
python script.py <input_PDB_file>
```

The script will create an output file named "fixed.pdb" containing the updated chain identifiers, atom names, and filtered resnames without the specified hydrogen atoms.

## Requirements

This script is written in Python and requires Python 3 to run.

"""

import sys

def increment_chain(chain):
    return chr(ord(chain) + 1)

def assign_chain_identifiers(input_file, output_file):
    with open(input_file, "r") as inp, open(output_file, "w") as out:
        current_chain = "A"
        last_residue = -1
        skip_resnames = {"TIP", "SOD", "POT", "CLA"}

        for line in inp:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()

                # Skip specified resnames
                if resname in skip_resnames:
                    continue

                # Skip atom names starting with "H"
                if atom_name.startswith("H"):
                    continue

                # Fix some atom names
                if resname == "HSD" or resname == "HSE" or resname == "HSP":
                    line = line[:17] + "HIS" + line[20:]
                if resname == "ILE" and atom_name == "CD":
                    line = line[:12] + " CD1" + line[16:]
                if atom_name == "OT1":
                    line = line[:12] + " O  " + line[16:]
                if atom_name == "OT2":
                    line = line[:12] + " OXT" + line[16:]

                residue = int(line[22:26].strip())

                if residue < last_residue:
                    current_chain = increment_chain(current_chain)

                line = line[:21] + current_chain + line[22:]
                last_residue = residue
            out.write(line)

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_PDB_file>")
        exit(1)

    input_file = sys.argv[1]
    output_file = "fixed.pdb"
    assign_chain_identifiers(input_file, output_file)

if __name__ == "__main__":
    main()
