# Create a cubic box of 2 nm on each side and place the system to the center and align its principle axes to the reference axes
# -princ usually helps to reduce your system size
gmx editconf -f conf.gro -o editconf.pdb -bt cubic -d 2 -c -princ
# Solvate the box
gmx solvate -cp editconf.pdb -o solvate.gro -p topol.top
# Add ions to make it neutral and of 0.15 M NaCl
# If you want KCl, add -pname K
gmx grompp -f ions.mdp -c solvate.gro -o ions.tpr -p topol.top
gmx genion -s ions.tpr -o ionized.pdb -conc 0.15 -neutral -p topol.top
