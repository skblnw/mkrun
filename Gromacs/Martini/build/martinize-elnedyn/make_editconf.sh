GMX_PATH=

input_pdb=cg_4nsw-elnedyn.pdb
output_pdb=cg_4nsw-elnedyn-oriented.pdb

gmx_mpi editconf -f $input_pdb -ndef -princ -rotate 90 0 0 -o $output_pdb
