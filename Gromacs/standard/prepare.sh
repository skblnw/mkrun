# Create a cubic box of 2 nm on each side and place the system to the center
gmx editconf -f conf.gro -o editconf.gro -bt cubic -d 2 -c -princ
# Solvate the box
gmx solvate -cp editconf.gro -o solvate.gro -p topol.top
# Add ions to make it neutral and of 0.15M
gmx grompp -f ions.mdp -c solvate.gro -o ions.tpr -p topol.top
gmx genion -s ions.tpr -o ionized.gro -conc 0.15 -neutral -p topol.top