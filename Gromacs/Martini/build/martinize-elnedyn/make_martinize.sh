input_pdb=final-modified.pdb
output_top=mdff_cg.top
output_pdb=mdff_cg.pdb
dssp_bin=/home/kevin/Dropbox/Source/DSSP/dssp-2.0.4-linux-amd64
# true or false
elastic=true
ef=500
el=0.5
eu=0.9
ea=0
ep=0
# martini22 or elnedyn22
ff=elnedyn22

if $elastic; then
    python martinize.py -f $input_pdb -o $output_top -x $output_pdb -dssp $dssp_bin -elastic -ef $ef -el $el -eu $eu -ea $ea -ep $ep -p backbone -ff $ff
else
    python martinize.py -f $input_pdb -o $output_top -x $output_pdb -dssp $dssp_bin -p backbone -ff $ff
fi
