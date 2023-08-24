import argparse

def parse_pdb(pdb_filename):
    """
    Parse a PDB file and extract atom details.
    
    Args:
    - pdb_filename (str): Path to the PDB file.
    
    Returns:
    - List[Dict]: A list of dictionaries where each dictionary contains details of an atom.
    """
    
    atoms = []
    
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                atom = {
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "seg_name": line[72:76].strip(),
                    "element": line[76:78].strip(),
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                    "residue": int(line[22:26].strip())
                }
                atoms.append(atom)
                
    return atoms

def main(args):
    atoms = parse_pdb(args.pdb_file)
    print(atoms[:1])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PDB parser')
    parser.add_argument('pdb_file', type=str, help='PDB')
  
    args = parser.parse_args()
    main(args)
