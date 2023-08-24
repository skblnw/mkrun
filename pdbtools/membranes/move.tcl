set sel [atomselect 0 all]
set selprot [atomselect 1 all]
set selmemb [atomselect 2 all]
$selprot moveby [vecsub [measure center $sel] [measure center $selprot]]
$selmemb moveby [vecsub [measure center $sel] [measure center $selmemb]]
$selprot writepdb prot_move.pdb
$selmemb writepdb lipids_move.pdb
