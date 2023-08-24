import argparse

def parse_pdb(pdb_filename):
    """
    Parse a PDB file and extract atom details.
    """
    atoms = []
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                atom = {
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:21].strip(),
                    "seg_name": line[72:76].strip(),
                    "element": line[76:78].strip(),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                    "residue": line[22:26].strip()
                }
                atoms.append(atom)
    return atoms

def detailed_lipid_counts(atoms, zcenter):
    """
    Count the number of unique residues for each lipid type in both the 'UP' and 'DOWN' regions.
    """
    lipid_residue_counts = {}
    lipid_residues = {atom["res_name"] for atom in atoms if atom["atom_name"] == "P"}
    for residue in lipid_residues:
        up_atoms = [atom for atom in atoms if atom["res_name"] == residue and atom["z"] > zcenter and atom["atom_name"] == "P"]
        down_atoms = [atom for atom in atoms if atom["res_name"] == residue and atom["z"] < zcenter and atom["atom_name"] == "P"]
        lipid_residue_counts[residue] = {
            "UP": len(up_atoms),
            "DOWN": len(down_atoms)
        }
    return lipid_residue_counts

def is_non_standard_amino_acid(res_name):
    """
    Check if a residue name does not belong to the 20 standard amino acids.
    """
    standard_amino_acids = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL"
    }
    return res_name not in standard_amino_acids

def count_residues_by_resname(atoms):
    """
    Count the number of residues for each res_name, with special handling for 'SOL'.
    """
    residues_by_res_name = {}
    # Define a list pairing res_name and atom_name
    atom_name_for_res = {
        "SOL": "OW",
        "NA": "NA",
        "CL": "CL"
    }
    
    # Count residues using specific atom names
    for res_name, atom_name in atom_name_for_res.items():
        residues_by_res_name[res_name] = sum(1 for atom in atoms if atom["res_name"] == res_name and atom["atom_name"] == atom_name)
    
    # Count unique residues for other residues
    for atom in atoms:
        res_name = atom["res_name"]
        atom_name = atom["atom_name"]
        
        # If res_name has already been counted, skip
        if res_name in atom_name_for_res:
            continue
        
        # For lipid and other residues, count unique residues
        if is_non_standard_amino_acid(res_name):
            residue_key = (res_name, atom["residue"])
            if res_name not in residues_by_res_name:
                residues_by_res_name[res_name] = set()
            residues_by_res_name[res_name].add(residue_key)
    
    # Convert sets to counts
    for res_name in residues_by_res_name.keys():
        if isinstance(residues_by_res_name[res_name], set):
            residues_by_res_name[res_name] = len(residues_by_res_name[res_name])
    return residues_by_res_name

def print_counts(atoms, residues_by_res_name, detailed_counts):
    """
    Modified print function to display counts of residues by res_name.
    """
    print("Counts of all unique residues:")
    for res_name, count in residues_by_res_name.items():
        print(f"{res_name}: {count} residues")
    
    # Print detailed lipid counts
    print("\nDetailed counts for each lipid type:")
    for lipid, counts in detailed_counts.items():
        print(f"{lipid}:")
        print(f"  UPPER: {counts['UP']} residues")
        print(f"  LOWER: {counts['DOWN']} residues")

def main(pdb_filename):
    """
    Modified main function to count residues by res_name and lipid residues in both the 'UP' and 'DOWN' regions.
    """
    atoms = parse_pdb(pdb_filename)
    z_values = [atom["z"] for atom in atoms]
    zcenter = sum(z_values) / len(z_values)
    residues_by_res_name = count_residues_by_resname(atoms)
    detailed_counts = detailed_lipid_counts(atoms, zcenter)
    print_counts(atoms, residues_by_res_name, detailed_counts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count lipid residues in both the 'UP' and 'DOWN' regions of a PDB file.")
    parser.add_argument("pdb_filename", type=str, help="Path to the PDB file.")
    args = parser.parse_args()
    main(args.pdb_filename)
