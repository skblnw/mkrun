rm -f editconf.pdb solvate.gro ions.tpr ionized.pdb
# Create a cubic box of 2 nm on each side and place the system to the center and align its principle axes to the reference axes
# -princ usually helps to reduce your system size
gmx editconf -f conf.gro -o editconf.pdb -d 2 -princ
# Solvate the box
gmx solvate -cp editconf.pdb -o solvate.gro -p topol.top
# Add ions to make it neutral and of 0.15 M NaCl
# If you want KCl, add -pname K
rm -f ions.tpr
gmx grompp -f ions.mdp -c solvate.gro -o ions.tpr -p topol.top -maxwarn 1
gmx genion -s ions.tpr -o ionized.pdb -conc 0.15 -neutral -p topol.top
