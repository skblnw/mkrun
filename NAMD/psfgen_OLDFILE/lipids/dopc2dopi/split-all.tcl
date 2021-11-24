mol new chainA_fill_aligned.pdb
foreach chain {A B} {
	set sel [atomselect top "protein and chain $chain"]
	$sel writepdb chainA-$chain.pdb
	$sel delete
}