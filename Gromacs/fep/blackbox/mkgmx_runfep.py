import os
import sys
import subprocess
import shutil
import re

def create_directories(rpath, positions, residues, num_dirs, prefix, free=True, complex=True):
    for pos in positions:
        for res in residues:
            create_pdb2gmx_directories(rpath, pos, res)
            create_windows(rpath, pos, res, num_dirs, prefix, free, complex)
            change_directory_and_run_script(rpath, pos, res)

def create_pdb2gmx_directories(rpath, pos, res):
    pdb2gmx_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")

    if os.path.exists(pdb2gmx_dir):
        shutil.rmtree(pdb2gmx_dir)
    os.makedirs(pdb2gmx_dir, exist_ok=True)

    with open(os.path.join(pdb2gmx_dir, "mut"), "w") as mut_file:
        mut_file.write(f"{pos} {res}")

    src_path = os.path.join(os.path.expanduser("~"), "charmm36-feb2021.ff")
    os.symlink(src_path, os.path.join(pdb2gmx_dir, "charmm36-feb2021.ff"))

    files_to_symlink = ['mkgmx_pdb2fep.py', 'ionized.gro']
    for file in files_to_symlink:
        src_path = os.path.join(rpath, file)
        os.symlink(src_path, os.path.join(pdb2gmx_dir, file))

    # free_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx/free/pdb2gmx")
    # complex_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx/complex/pdb2gmx")

    # os.makedirs(free_dir, exist_ok=True)
    # os.makedirs(complex_dir, exist_ok=True)

def create_windows(rpath, position, residue, num_dirs, prefix="win24.t1", free=True, complex=True):
    original_free_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/free")
    original_complex_dir = os.path.join(rpath, f"pos{position}/{residue}/pdb2gmx/complex")

    options = []
    if free:
        options.append(('free', original_free_dir))
    if complex:
        options.append(('complex', original_complex_dir))

    for option, original_dir in options:
        target_dir = os.path.join(rpath, f"pos{position}/{residue}/{prefix}/{option}")

        if os.path.exists(target_dir):
            shutil.rmtree(target_dir)
        os.makedirs(target_dir, exist_ok=True)

        os.symlink(os.path.join(original_dir, "pdb2gmx"), os.path.join(target_dir, "pdb2gmx"))
        shutil.copytree(f"{rpath}/mdp", f"{rpath}/pos{position}/{residue}/{prefix}/{option}", dirs_exist_ok=True)
        create_mdp_files(num_dirs, rpath, os.path.join(rpath, f"pos{position}/{residue}/{prefix}/{option}"))
        create_run_files(rpath, position, residue, option)

def create_pbs_files(rpath, pos, res, prefix):
    with open(os.path.join(rpath, "pbs"), "r") as pbs_file:
        pbs_content = pbs_file.read()
    pbs_content = re.sub(r'#PBS -N.*', f"#PBS -N p{pos}-{res}-{prefix}", pbs_content)
    with open(os.path.join(rpath, f"pos{pos}/{res}/win24.t1/{prefix}/pbs"), "w") as pbs_file:
        pbs_file.write(pbs_content)
    # subprocess.run(["qsub", os.path.join(rpath, f"pos{pos}/{res}/win24.t1/{prefix}/pbs")])

def create_run_files(rpath, pos, res, prefix):
    with open(os.path.join(rpath, "run"), "r") as run_file:
        run_content = run_file.read()
    with open(os.path.join(rpath, f"pos{pos}/{res}/win24.t1/{prefix}/run"), "w") as run_file:
        run_file.write(run_content)
    # subprocess.run(["bash", os.path.join(rpath, f"pos{pos}/{res}/win24.t1/{prefix}/run")])

def create_mdp_files(num_dirs, rpath, base_dir):
    for i in range(num_dirs):
        os.makedirs(os.path.join(base_dir, f"dir{i}"), exist_ok=True)
        for prefix in ['em', 'nvt', 'md']:
            with open(os.path.join(rpath, f"mdp/{prefix}.mdp"), 'r') as mod_mdp_file:
                mdp_content = mod_mdp_file.read()
            mdp_content = re.sub(r'; init_lambda_state.*',
                                 r'; init_lambda_state        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23',
                                 mdp_content)
            mdp_content = re.sub(r'^fep-lambdas.*',
                                 r'fep-lambdas              = 0 0.00001 0.0001 0.001 0.01 0.02 0.04 0.06 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.94 0.96 0.98 0.99 0.999 0.9999 0.99999 1.00',
                                 mdp_content, flags=re.MULTILINE)
            with open(os.path.join(base_dir, f"dir{i}/{prefix}.mdp"), 'w') as mdp_file:
                mdp_file.write(mdp_content.replace('init_lambda_state = 0', f'init_lambda_state = {i}'))

def change_directory_and_run_script(rpath, pos, res):
    original_dir = rpath
    target_dir = os.path.join(rpath, f"pos{pos}/{res}/pdb2gmx")

    script = 'mkgmx_pdb2fep.py'

    try:
        os.chdir(target_dir)
        subprocess.run(["python", script], text=True)
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":

    # Check if required files exist
    required_files = [
        "ionized.gro",
        "mkgmx_pdb2fep.py"
    ]
    for file in required_files:
        if not os.path.isfile(file):
            print(f"Error: Required file '{file}' not found.")
            sys.exit(1)

    # Check if required directories exist
    required_directories = [
        "mdp"
    ]
    for directory in required_directories:
        if not os.path.isdir(directory):
            print(f"Error: Required directory '{directory}' not found.")
            sys.exit(1)

    rpath = os.getcwd()
    positions = [2]
    residues = ["ALA"]
    create_directories(rpath, positions, residues, 24, 'win24.t1', free=True, complex=True)
