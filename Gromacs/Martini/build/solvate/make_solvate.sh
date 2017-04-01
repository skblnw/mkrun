#########################################
## Description: A bash make for solvating martini protein (and lipids) with BMW waters
## Author: Kevin@CityUHK 31 Mar 2016
## Usage: 
## Input:
## Output:
## Units: 
## Other Notes: 
#########################################

# Some global variables
GMX_PATH=
itp_dir=/home/kevin/Dropbox/QWD/scripts/MD/martini/itp/BMW

output_top=system.top
output_gro=system.gro

pdb_dir=../insane_membrane
input_pdb=$pdb_dir/protein_membrane.gro
# Usually input_top specify only lipid compositions
top_dir=../insane_membrane
input_top=$top_dir/insane.top
# Usually TMP_TOP specify itps names
TMP_TOP=./templates/template.top

water_gro=/home/kevin/Dropbox/QWD/scripts/MD/martini/itp/BMW/bmw.gro # Where BMW.gro is
mol_name=W
radius=0.105
scale=1.14
maxsol=120000
d=5
conc=0.18

# Prepare output_top before adding water
sed '' $TMP_TOP > $output_top
sed '1,/\[ molecules \]/d' $input_top >> $output_top

# Solvate
# Step 1
gmx_mpi editconf -f $input_pdb -box 12 40 60 -center 6 15 30 -o step1.gro
wait
# Step 2
gmx_mpi solvate -cp step1.gro -cs $water_gro -scale $scale -maxsol $maxsol -p $output_top -o step2.gro >& step2.log
wait

# BMW-modifications to itp files according to the force-field
sed -i 's/AC1/C1/g' Protein_A.itp Protein_B.itp
sed -i 's/AC2/C2/g' Protein_A.itp Protein_B.itp

if ask ">>Do ionize? "; then
    cp step2.gro $output_gro
    rm -f \#* mdout.mdp *.tpr step*
    exit 1
else
    cp itps/* .
    cp mdps/* .
    # Ionize
    # Step 3
    gmx_mpi grompp -f mini-ionize.mdp -c step2.gro -p $output_top -o mini-ionize.tpr
    wait
    gmx_mpi genion -s mini-ionize.tpr -p $output_top -o step3.gro -pname NA+ -pq 1 -nname CL- -nq -1 -conc $conc
    wait
    
    sed -i 's/NA+/ION/g' step3.gro
    sed -i 's/CL-/ION/g' step3.gro
    sed -i 's/NA+/NA/g' $output_top
    sed -i 's/CL-/CL/g' $output_top
    #sed -i 's/ NA/NA+/g' $output_gro3
    #sed -i 's/ CL/CL-/g' $output_gro3
    
    cp step3.gro $output_gro
    rm -f \#* mdout.mdp *.tpr step*
fi
