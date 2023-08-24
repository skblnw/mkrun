import os
import subprocess
import argparse
import tempfile

# General-purpose functions

def BASH(command):
    """Run an external command."""
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        raise Exception(f"Command '{command}' failed with error: {result.stderr.decode('utf-8')}")
    return result.stdout.decode('utf-8')

def BASHz(cmd):
    subprocess.run(cmd, check=True, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def BASHv(cmd):
    subprocess.run(cmd, check=True, shell=True)

def check_files_existence(filenames):
    """Check if multiple files exist."""
    missing_files = [f for f in filenames if not os.path.exists(f)]
    if missing_files:
        raise FileNotFoundError(", ".join(missing_files) + " not found")

def combine_pdb_files(proteins, lipids):
    """Combine two pdb files into combined.pdb."""
    with open(proteins, 'r') as f1, open(lipids, 'r') as f2:
        proteins_lines = [line for line in f1 if line.startswith("ATOM")]
        lipids_lines = [line for line in f2 if line.startswith("ATOM")]

    with open("combined.pdb", 'w') as out_file:
        out_file.writelines(proteins_lines + lipids_lines)

def remove_overlapping_lipids():
    """Label and remove lipids that overlap in space with the proteins using a TCL script."""
    tcl_script_content = """
    mol new combined.pdb
    set seltext_for_lipids "resname POPC"
    set badlipids [atomselect top "$seltext_for_lipids and same resid as {$seltext_for_lipids and within 0.6 of protein}"]
    set sel [atomselect top "not index [$badlipids get index]"]
    $sel writepdb prot_memb_hole.pdb
    quit
    """

    # Write the TCL script to a temporary file
    with tempfile.NamedTemporaryFile(suffix=".tcl", delete=False) as temp_file:
        temp_file_name = temp_file.name
        temp_file.write(tcl_script_content.encode())

    # Execute the TCL script using vmd
    vmd_path = '"/Applications/VMD 1.9.4a51-arm64-Rev9.app/Contents/Resources/VMD.app/Contents/MacOS/VMD"'
    BASH(f"rlwrap {vmd_path} -dispdev text -e {temp_file_name}")

    # Clean up temporary file
    os.remove(temp_file_name)

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

def count_waters_and_lipids(atoms):
    """
    Count the number of waters and lipids in the atom list.
    """
    waters = sum(1 for atom in atoms if atom["res_name"] == "SOL" and atom["atom_name"] == "OW")
    lipids = sum(1 for atom in atoms if atom["res_name"] == "POPC" and atom["atom_name"] == "P")
    return waters, lipids

def extract_c1_z_coordinates(atoms):
    """
    Extract z-coordinates for "C1" atoms of lipids.
    """
    c1_z_coordinates = [atom["z"] for atom in atoms if atom["res_name"] == "POPC" and atom["atom_name"] == "C1"]
    return c1_z_coordinates

def determine_water_boundaries(atoms):
    """
    Determine the boundary z-coordinates for "bad" waters based on lipid "C1" atoms.
    """
    c1_z_coordinates = extract_c1_z_coordinates(atoms)
    average_z = sum(c1_z_coordinates) / len(c1_z_coordinates)
    
    upper_c1s = [z for z in c1_z_coordinates if z > average_z]
    lower_c1s = [z for z in c1_z_coordinates if z <= average_z]
    
    upper_boundary = min(upper_c1s)
    lower_boundary = max(lower_c1s)
    
    return lower_boundary, upper_boundary

def remove_water_between_lipids(lower_boundary, upper_boundary):
    """Label and remove lipids that overlap in space with the proteins using a TCL script."""
    tcl_script_content = f"""
    mol new solvate.pdb
    set badwater [atomselect top "resname SOL and same residue as z < {upper_boundary} and z > {lower_boundary}"]
    set sel [atomselect top "not index [$badwater get index]"]
    $sel writepdb solvate_delw.pdb
    quit
    """

    # Write the TCL script to a temporary file
    with tempfile.NamedTemporaryFile(suffix=".tcl", delete=False) as temp_file:
        temp_file_name = temp_file.name
        temp_file.write(tcl_script_content.encode())

    # Execute the TCL script using vmd
    vmd_path = '"/Applications/VMD 1.9.4a51-arm64-Rev9.app/Contents/Resources/VMD.app/Contents/MacOS/VMD"'
    BASH(f"rlwrap {vmd_path} -dispdev text -e {temp_file_name}")

    # Clean up temporary file
    os.remove(temp_file_name)

def modify_topol(lipids, waters):
    """Modify the topol.top file according to the specified criteria."""
    with open("topol.top", 'r') as file:
        lines = file.readlines()

    # Identify the line containing "[ molecules ]"
    index = lines.index("[ molecules ]\n")

    # Remove all lines following this that start with "SOL"
    lines = lines[:index+1] + [line for line in lines[index+1:] if not line.startswith("SOL")]

    # Append the lines "POPC {lipids}" and "SOL {waters}"
    lines.append(f"POPC {lipids}\n")
    lines.append(f"SOL {waters}\n")

    # Write the modified lines back to topol.top
    with open("topol.top", 'w') as file:
        file.writelines(lines)

def add_ions():
    mdp_content = f"""
    integrator  = md
    dt          = 0.001
    nsteps      = 50000
    nstlist         = 1
    cutoff-scheme   = Verlet
    ns_type         = grid
    coulombtype     = PME
    rcoulomb        = 1.0
    rvdw            = 1.0
    pbc             = xyz
    """

    # Write the TCL script to a temporary file
    with tempfile.NamedTemporaryFile(suffix=".mdp", delete=False) as temp_file:
        temp_file_name = temp_file.name
        temp_file.write(mdp_content.encode())

    BASH(f"gmx grompp -f {temp_file_name} -c solvate_delw.pdb -o ions.tpr -p topol.top -maxwarn 1")
    check_files_existence(["ions.tpr"])

    BASHv(f"gmx genion -s ions.tpr -o ionized.pdb -conc 0.15 -neutral -p topol.top")
    check_files_existence(["ionized.pdb"])

def main(args):

    proteins = args.proteins
    lipids = args.lipids

    # Check existence of initial files
    check_files_existence([proteins, lipids])

    # Combine proteins and lipids into combined.pdb
    combine_pdb_files(proteins, lipids)

    remove_overlapping_lipids()

    BASH(f"gmx editconf -f prot_memb_hole.pdb -o editconf.gro -box {args.box_size}")
    check_files_existence(["editconf.gro"])

    check_files_existence(["topol.top"])
    BASH(f"gmx solvate -cp editconf.gro -o solvate.pdb -p topol.top")
    check_files_existence(["solvate.pdb"])

    # Parse the PDB file after gmx editconf
    atoms = parse_pdb("solvate.pdb")

    # Determine the boundary z-coordinates for "bad" waters
    lower_boundary, upper_boundary = determine_water_boundaries(atoms)

    remove_water_between_lipids(lower_boundary, upper_boundary)

    # Count the number of waters and lipids in the new PDB file
    cleaned_atoms = parse_pdb("solvate_delw.pdb")
    waters, lipids = count_waters_and_lipids(cleaned_atoms)

    modify_topol(lipids, waters)

    add_ions()

    BASH(f"rm \#* mdout.mdp combined.pdb prot_memb_hole.pdb editconf.gro solvate.pdb solvate_delw.pdb ions.tpr")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Prepare a pdb file.")
    parser.add_argument("proteins", help="Path to proteins.pdb")
    parser.add_argument("lipids", help="Path to lipids.pdb")
    parser.add_argument("box_size", type=str, help="Box size for gmx editconf, e.g., '10.12176 10.12176 10'")
    args = parser.parse_args()

    main(args)
