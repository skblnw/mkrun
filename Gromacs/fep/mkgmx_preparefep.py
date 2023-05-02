import os
import subprocess
from pathlib import Path
import shutil
import glob

# Set variables
CHAIN_DIR = "./split/chains"
LOG_GROMPP = "LOG_grompp"
IONS_MDP = "ions.mdp"

# Remove the existing symlink if it exists
Path("chains").unlink(missing_ok=True)

# Create symlinks
Path("chains").symlink_to(CHAIN_DIR, target_is_directory=True)

def cleanup():
    files_to_remove = [
        "#*", "editconf.pdb", "solvate.pdb", "ions.mdp", "ions.tpr", "index.ndx", "mdout.mdp", "LOG_grompp",
        "conf.pdb", "free.pdb", "ionized.pdb", "mol.pdb", "peptide_mutate.pdb", "X.pdb"
    ]
    for file in files_to_remove:
        Path(file).unlink(missing_ok=True)

# Remove previous topol.top file
Path("topol.top").unlink(missing_ok=True)

# Run gmx to get conf.gro
subprocess.run(["gmx", "pdb2gmx", "-f", "chains/chain_C.pdb", "-ff", "charmm36-feb2021", "-water", "tip3p"])

# Mutate conf.gro using pmx
subprocess.run(["pmx", "mutate", "-f", "conf.gro", "-ff", "charmm36m-mut", "-o", "peptide_mutate.pdb", "--script", "mut"])

# Divide (and fix the names of) complex PDB into multiple fragments using segname
tcl_script = '''
mol new peptide_mutate.pdb
set sel [atomselect top "all"]
$sel set chain X
$sel writepdb X.pdb
quit
'''

with open("tcl", "w") as tcl_file:
    tcl_file.write(tcl_script)

subprocess.run(["vmd", "-dispdev", "text", "-e", "tcl"])
Path("tcl").unlink(missing_ok=True)

# Get free.pdb
with open("X.pdb", "r") as f, open("mol.pdb", "w") as mol_f:
    for line in f:
        if line.startswith("ATOM"):
            mol_f.write(line)

subprocess.run(["gmx", "pdb2gmx", "-f", "mol.pdb", "-ff", "charmm36m-mut", "-water", "tip3p", "-o", "free.pdb"])

# Get conf.pdb (complex) and topol.top (complex)
for chain in ["A", "B"]:
    with open(f"chains/chain_{chain}.pdb", "r") as chain_f, open("mol.pdb", "a") as mol_f:
        for line in chain_f:
            if line.startswith("ATOM"):
                mol_f.write(line)

subprocess.run(["gmx", "pdb2gmx", "-f", "mol.pdb", "-ff", "charmm36m-mut", "-water", "tip3p", "-o", "conf.pdb"])

# Generate hybrid topology using pmx
subprocess.run(["pmx", "gentop", "-p", "topol.top", "-ff", "charmm36m-mut"])

with open("pmxtop.top", "r") as f, open("pmxtop_free.top", "w") as free_f:
    for line in f:
        if not line.startswith("Protein_chain_"):
            free_f.write(line)

# Create ions.mdp file
ions_mdp_content = '''
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
'''

with open(IONS_MDP, "w") as ions_mdp_file:
    ions_mdp_file.write(ions_mdp_content)

# Remove and create directories
try:
    shutil.rmtree("free")
except FileNotFoundError:
    pass
try:
    shutil.rmtree("complex")
except FileNotFoundError:
    pass
Path("free/pdb2gmx").mkdir(parents=True, exist_ok=True)
Path("complex/pdb2gmx").mkdir(parents=True, exist_ok=True)


# Function for creating box, solvating, and adding ions
def prepare_system(system_name, input_pdb, output_top):
    # Create the box
    subprocess.run(f"echo 1 | gmx editconf -f {input_pdb} -o editconf.pdb -princ -d 1.0 -bt cubic", shell=True, check=True, text=True)

    # Solvate the box
    subprocess.run(["gmx", "solvate", "-cp", "editconf.pdb", "-o", "solvate.pdb", "-p", output_top])

    # Add ions to make it neutral and of 0.15 M NaCl
    # If you want KCl, add -pname K
    subprocess.run(["gmx", "grompp", "-f", IONS_MDP, "-c", "solvate.pdb", "-o", "ions.tpr", "-p", output_top, "-maxwarn", "1"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    subprocess.run(f"echo 13 | gmx genion -s ions.tpr -o ionized.pdb -conc 0.15 -neutral -p {output_top}", shell=True, check=True, text=True)

    # Create a local copy of ionized.pdb specific to the system
    shutil.copy("ionized.pdb", f"{system_name}_ionized.pdb")

    # Create index file
    with subprocess.Popen(f"gmx make_ndx -f {system_name}_ionized.pdb", shell=True, stdin=subprocess.PIPE, text=True) as proc:
        proc.communicate(f"q\n")

    # Copy files to the respective directory
    files_to_copy = [f"{system_name}_ionized.pdb", "index.ndx"] + glob.glob("pmx_*") + glob.glob("posre_*")
    for file in files_to_copy:
        shutil.copy(file, f"{system_name}/pdb2gmx/{file}")
        Path(file).unlink(missing_ok=True)

    Path(output_top).replace(f"{system_name}/pdb2gmx/topol.top")

# Prepare complex system
prepare_system("complex", "conf.pdb", "pmxtop.top")

# Prepare free system
# Get the box size first
subprocess.run(["gmx", "editconf", "-f", "ionized.pdb", "-o", "tmp.gro"])
with open("tmp.gro", "r") as f:
    bs = f.readlines()[-1]

Path("ionized.pdb").unlink(missing_ok=True)
Path("tmp.gro").unlink(missing_ok=True)

prepare_system("free", "free.pdb", "pmxtop_free.top")

# Cleanup temporary files
cleanup()
