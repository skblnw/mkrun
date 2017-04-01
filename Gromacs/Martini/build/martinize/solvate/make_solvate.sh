GMX_PATH=
itp_martini_name=martini_v2.2.itp
itp_dir=/home/kevin/Dropbox/QWD/scripts/MD/Martini/itp/particle_definitions/

input_dir=../
input_pdb=$input_dir/50_renumbered.pdb
input_top=$input_dir/3F_cg.top
d=2
output_gro1=tmp.gro
water_gro=/home/kevin/Dropbox/QWD/scripts/MD/Martini/gro/water-box-CG_303K-1bar.gro
radius=0.21
output_gro2=solvated.gro
output_top=system.top

ionize=false

output_gro3=system.gro
conc=0.18


# Solvate
gmx_mpi editconf -f $input_pdb -d $d -bt cubic -o $output_gro1
wait
gmx_mpi solvate -cp $output_gro1 -cs $water_gro -radius $radius -o $output_gro2 >& solvate.log
wait

NW=`grep "W (   1 atoms):" TMP | awk '{print $5}'`
echo -e "#include \"$itp_martini_name\"" > $output_top
echo -e "#include \"martini_v2.0_ions.itp\"" >> $output_top
sed '1d' $input_top >> $output_top
echo -e "\nW $NW" >> $output_top

# Ionize
if $ionize; then
    #cp $itp_dir/$itp_martini_name .
    cp $input_dir/*.itp .

    gmx_mpi grompp -f mini-ionize.mdp -c $output_gro2 -p $output_top -o mini-ionize.tpr
    wait
    gmx_mpi genion -s mini-ionize.tpr -p $output_top -o $output_gro3 -pname NA+ -pq 1 -nname CL- -nq -1 -conc $conc
    wait

    sed -i 's/NA+/ION/g' $output_gro3
    sed -i 's/CL-/ION/g' $output_gro3
    sed -i 's/ NA/NA+/g' $output_gro3
    sed -i 's/ CL/CL-/g' $output_gro3

    rm -f \#* mdout.mdp $output_gro1 $output_gro2
else
    cp solvated.gro system.gro
fi
