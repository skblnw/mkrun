#!/bin/bash

# Definitions of variables
VMD="vmd"
SEGNAME="PROC"
SELECT_TEXT="segname $SEGNAME"
MD_PDB_FILE="md.pdb"
FF_FILES=("readcharmmtop1.2/top_all36_prot.rtf" "readcharmmtop1.2/top_all36_hybrid.inp" "top_all36_propatch.rtf")
EQ_FILES=("fep.tcl" "fep.eq.namd" "fep.namd")


# /-------------------/
# /     Functions     /
# /-------------------/


# Function to check if a directory exists, and if not, make it
function check_and_make_dir() {
  if [ ! -d $1 ]; then
    mkdir -p $1
  fi
}

# Function to move files
function move_files() {
  mv ionized.fep $1
  mv ionized.p* cell_size.str $1/pdb2namd/vmd_solvate
}

# Function to rsync files
function sync_files() {
  local dest=$1
  for file in ${EQ_FILES[*]}; do
    rsync -aL ../../eq/$file $dest/eq/$file
  done
}

# Function to check if all EQ_FILES exist
function check_files_exist() {
  for file in ${EQ_FILES[*]}; do
    if [ ! -f "./eq/$file" ]; then
      echo "Error: File ./eq/$file does not exist!"
      exit 1
    fi
  done
  for file in ${FF_FILES[*]}; do
    if [ ! -f "./$file" ]; then
      echo "Error: File ./$file does not exist!"
      exit 1
    fi
  done
}

# Function to run namd3
function run_namd3() {
  if [[ "$1" == "-run" ]]; then 
    cd eq
    echo "mknamd> Running eq..."
    namd3 +p8 +devices 1 fep.eq.namd >& LOG_eq
    cd ..
  fi
}

# Function to create run.sh file
function create_run_file() {
  local filename=$1
  local dev=$2

  cat << EOF > $filename
#!/bin/bash

NAMD="namd3"

cd eq
\$NAMD +p8 +devices $dev fep.eq.namd >& LOG_eq
cd ..

for ii in \$(seq 1 1); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial\$ii
    \$NAMD +p1 +devices $dev trial\$ii/fep.namd >& trial\$ii/LOG_fep
done
EOF

  # Make the file executable
  chmod +x $filename
}

# Function to create slurm file
function create_slurm_file() {
  local filename=$1

  cat << EOF > $filename
#!/bin/bash
#SBATCH -J <job name>
#SBATCH -p single
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1

echo "Start time: \$(date)"
echo "SLURM_JOB_NODELIST: \$SLURM_JOB_NODELIST"
echo "hostname: \$(hostname)"
echo "CUDA_VISIBLE_DEVICES: \$CUDA_VISIBLE_DEVICES"
echo "Job directory: \$(pwd)"

# Decide the software version
source /public/software/profile.d/apps_namd-3.0alpha9.sh

NAMD="/public/software/apps/NAMD_3.0alpha9/namd3"

\$NAMD +p4 +devices 0 eq/fep.eq.namd >& eq/LOG_eq

for ii in \$(seq 1 1); do
    rsync -avh eq/fep.namd eq/fep.tcl eq/equilibrate.coor eq/equilibrate.vel eq/equilibrate.xsc trial\$ii
    \$NAMD +p1 +devices 0 trial\$ii/fep.namd >& trial\$ii/LOG_fep
done
EOF
}

psfgen () {
  cat > tcl <<EOF
package require psfgen

mol new $MD_PDB_FILE
foreach ii {A B} {
    set sel [atomselect top "segname PRO\$ii"]
    \$sel writepdb chains/PRO\$ii.pdb
}
set selsc [atomselect top "$SELECT_TEXT and resid $1 and not name C CA N O HN HA CB"]
foreach idx [\$selsc get index] {
  set sel [atomselect top "index \$idx"]
  set name [\$sel get name]
  \$sel set name \${name}A
}
set sel [atomselect top "$SELECT_TEXT"]
\$sel writepdb chains/mutant.pdb
EOF
  cat >> tcl <<'EOF'
resetpsf
topology readcharmmtop1.2/top_all36_prot.rtf
topology readcharmmtop1.2/top_all36_hybrid.inp
topology top_all36_propatch.rtf
# topology toppar_water_ions_namd.str

# Aliases borrowed from AutoPSF
  pdbalias atom ILE CD1 CD
  pdbalias atom SER HG HG1
  pdbalias residue HIS HSD

# Water aliases
  pdbalias residue HOH TIP3
  pdbalias atom TIP3 O OH2

# Ion aliases
  pdbalias residue K POT
  pdbalias atom K K POT
  pdbalias residue ICL CLA
  pdbalias atom ICL CL CLA
  pdbalias residue INA SOD
  pdbalias atom INA NA SOD
  pdbalias residue CA CAL
  pdbalias atom CA CA CAL

# Other aliases
  pdbalias atom LYS 1HZ HZ1
  pdbalias atom LYS 2HZ HZ2
  pdbalias atom LYS 3HZ HZ3

  pdbalias atom ARG 1HH1 HH11
  pdbalias atom ARG 2HH1 HH12
  pdbalias atom ARG 1HH2 HH21
  pdbalias atom ARG 2HH2 HH22

  pdbalias atom ASN 1HD2 HD21
  pdbalias atom ASN 2HD2 HD22

EOF

# /---------------------------------------------/
# /     Check if it's the terminal residue      /
# /       If yes, apply terminal patches        /
# /---------------------------------------------/

  if [[ $ii -eq 1 ]] && [[ $resname == "GLY" ]]; then
    cat >> tcl <<EOF
segment MUT {
  pdb chains/mutant.pdb
  mutate $ii ${mutation[$resname]}
  first none
  last CTER
}
EOF
  elif [[ $ii -eq $length_of_peptide ]] && [[ $resname == "GLY" ]]; then
    cat >> tcl <<EOF
segment MUT {
  pdb chains/mutant.pdb
  mutate $ii ${mutation[$resname]}
  first NTER
  last none
}
EOF
  else
    cat >> tcl <<EOF
segment MUT {
  pdb chains/mutant.pdb
  mutate $ii ${mutation[$resname]}
}
EOF
  fi

# /--------------------------------/
# /     Check if it's Proline      /
# /     If yes, apply patches      /
# /--------------------------------/

  if [[ $ii -ne 1 ]] && [[ $resname == "PRO" ]]; then
    jj=$((ii-1))
    cat >> tcl <<EOF
patch AABP MUT:$jj
patch AASP MUT:$ii
coordpdb chains/mutant.pdb MUT
EOF
  else
    cat >> tcl <<EOF
coordpdb chains/mutant.pdb MUT
EOF
  fi
}

psfgen_free () {
  cat >> tcl <<'EOF'
regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
  $VMD -dispdev text -e tcl >> LOG_vmd
}

psfgen_bound () {
  cat >> tcl <<'EOF'
segment PROA {
  pdb chains/PROA.pdb
}
patch DISU PROA:101 PROA:164
patch DISU PROA:203 PROA:259
coordpdb chains/PROA.pdb PROA

segment PROB {
  pdb chains/PROB.pdb
}
patch DISU PROB:25 PROB:80
coordpdb chains/PROB.pdb PROB

regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF
  $VMD -dispdev text -e tcl >> LOG_vmd
}

psfmerge () {
  cat > tcl2 <<EOF
package require topotools
mol new mutant.psf
mol addfile mutant.pdb
mol new $MD_PDB_FILE
set sel1 [atomselect 0 all]
set sel2 [atomselect 1 "segname PROA PROB"]
set mol [::TopoTools::selections2mol "\$sel1 \$sel2"]
animate write psf prot.psf \$mol
animate write pdb prot.pdb \$mol
quit
EOF
  $VMD -dispdev text -e tcl2 >& LOG_vmd
}

solvate () {
  cat > tcl3 <<'EOF'
package require Orient
mol new prot.pdb
set sel [atomselect top "all"]
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set I [Orient::calc_principalaxes $sel]
set A [Orient::orient $sel [lindex $I 1] {0 1 0}]
$sel move $A
$sel moveby [vecinvert [measure center $sel]]
$sel writepdb tmp.pdb

package require solvate
solvate prot.psf tmp.pdb -minmax {{-40 -40 -52} {40 40 52}} -o solvated
package require autoionize
autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -o ionized

mol new ionized.pdb type pdb waitfor all
set all [atomselect top "all"]
set fout [open "cell_size.str" w]
set all [atomselect top water] 
set minmax [measure minmax $all] 
set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
puts $fout "cellBasisVector1 [lindex $vec 0] 0 0" 
puts $fout "cellBasisVector2 0 [lindex $vec 1] 0" 
puts $fout "cellBasisVector3 0 0 [lindex $vec 2]" 
set center [measure center $all] 
puts $fout "cellOrigin $center" 
close $fout
quit
EOF
  $VMD -dispdev text -e tcl3 >> LOG_vmd
}

markfep () {
  cat > tcl4 <<EOF
mol new ionized.pdb
set sel [atomselect top all]
\$sel set beta 0
set sel [atomselect top "resname \".*2.*\" and not name CA CB HA HB and name \".*A\""]
\$sel set beta -1
set sel [atomselect top "resname \".*2.*\" and not name CA CB HA HB and name \".*B\""]
\$sel set beta 1
set sel [atomselect top all]
\$sel writepdb ionized.fep
quit
EOF
  $VMD -dispdev text -e tcl4 >> LOG_vmd
}


# /-------------------/
# /     Main body     /
# /-------------------/

[ $# -lt 1 ] && { echo "mknamd> Usage: $0 [all|RESID] [-run]"; echo "mknamd> Default peptide selection is: $SELECT_TEXT"; echo "mknamd> If apply multiple RESIDs, use e.g. \"1 2 3\""; exit 1; }

[ ! -f $MD_PDB_FILE ] && { echo "md.pdb does not exist!"; exit 1; }

check_files_exist

length_of_peptide=`grep $SEGNAME $MD_PDB_FILE | grep "CA" -c`
sequence=`grep $SEGNAME $MD_PDB_FILE | grep "CA" | awk '{print $4}'`
# echo -e "mknamd> the antigen sequence is: $sequence\nmknamd> Now working on..."

if [[ "$1" == "all" ]]; then
  list=$(seq 1 $length_of_peptide)
else
  list="$1"
fi

declare -A mutation=( ["ARG"]="R2A" \
                      ["ASN"]="N2A" \
                      ["ASP"]="D2A" \
                      ["CYS"]="C2A" \
                      ["GLN"]="Q2A" \
                      ["GLU"]="E2A" \
                      ["GLY"]="G2A" \
                      ["HIS"]="H2A" \
                      ["HSD"]="H2A" \
                      ["HSE"]="H2A" \
                      ["HSP"]="H2A" \
                      ["ILE"]="I2A" \
                      ["LEU"]="L2A" \
                      ["LYS"]="K2A" \
                      ["MET"]="M2A" \
                      ["PHE"]="F2A" \
                      ["PRO"]="P2A" \
                      ["SER"]="S2A" \
                      ["THR"]="T2A" \
                      ["TRP"]="W2A" \
                      ["TYR"]="Y2A" \
                      ["VAL"]="V2A" \
                      )

ii=0
for resname in $sequence
do
  ii=$((ii+1))
  echo "mknamd> $ii $resname"

  if [[ "$resname" == "ALA" ]]; then
    echo "mknamd> Do nothing"; continue
  fi

  if [[ ! $list =~ $ii ]]; then
    echo "mknamd> Do nothing"; continue
  else

    if [ -d pos$ii ]; then
      cp -r pos$ii pos$ii.BAK
      rm -r pos$ii
    fi

    check_and_make_dir "chains"
    check_and_make_dir "pos$ii"

    # /---------------------/
    # /     Bound State     /
    # /---------------------/

    check_and_make_dir "pos$ii/bound"
    check_and_make_dir "pos$ii/bound/eq"
    check_and_make_dir "pos$ii/bound/pdb2namd"
    check_and_make_dir "pos$ii/bound/pdb2namd/vmd_solvate"

    psfgen $ii
    psfgen_bound
    solvate
    markfep

    move_files "pos$ii/bound"

    cd pos$ii/bound
    sync_files "."
    #create_run_file "run.sh" 1
    #run_namd3 "$2"
    create_slurm_file "slurm"
    cd ../..

    # /--------------------/
    # /     Free State     /
    # /--------------------/

    check_and_make_dir "pos$ii/free"
    check_and_make_dir "pos$ii/free/eq"
    check_and_make_dir "pos$ii/free/pdb2namd"
    check_and_make_dir "pos$ii/free/pdb2namd/vmd_solvate"

    psfgen $ii
    psfgen_free
    solvate
    markfep

    move_files "pos$ii/free"

    cd pos$ii/free
    sync_files "."
    #create_run_file "run.sh" 1
    #run_namd3 "$2"
    create_slurm_file "slurm"
    cd ../..

    rm -rf chains tcl* tmp.pdb prot.p* solvated.*
  fi

  res_prev=$resname
done
