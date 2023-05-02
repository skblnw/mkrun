import os
import subprocess

NWINDOWS = 28
start = 0
end = NWINDOWS - 1

def run_command(cmd, logfile=None):
    if logfile:
        with open(logfile, "a") as log:
            subprocess.run(cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, shell=True, check=True)

for ii in range(start, end + 1):
    run_command(f"gmx_mpi grompp -f dir{ii}/em.mdp -c pdb2gmx/ionized.pdb -r pdb2gmx/ionized.pdb -p pdb2gmx/topol.top -o dir{ii}/em.tpr -maxwarn 3", "grompp_em.log")

run_command(f"mpirun -np {NWINDOWS} gmx_mpi mdrun -v -deffnm em -multidir {' '.join(['dir'+str(i) for i in range(NWINDOWS)])}")

for ii in range(start, end + 1):
    run_command(f"gmx_mpi grompp -f dir{ii}/nvt.mdp -c dir{ii}/em.gro -r dir{ii}/em.gro -p pdb2gmx/topol.top -o dir{ii}/nvt.tpr -maxwarn 2", "grompp_nvt.log")

run_command(f"mpirun -np {NWINDOWS} gmx_mpi mdrun -v -deffnm nvt -multidir {' '.join(['dir'+str(i) for i in range(NWINDOWS)])} -nb gpu")

for ii in range(start, end + 1):
    run_command(f"gmx_mpi grompp -f dir{ii}/md.mdp -c dir{ii}/nvt.gro -p pdb2gmx/topol.top -o dir{ii}/md.tpr -maxwarn 1", "grompp_md.log")

run_command(f"mpirun -np {NWINDOWS} gmx_mpi mdrun -v -deffnm md -multidir {' '.join(['dir'+str(i) for i in range(NWINDOWS)])} -replex 1000 -nb gpu")
