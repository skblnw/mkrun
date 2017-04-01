package require psfgen
resetpsf

topology toppar/top_all36_prot.rtf
topology toppar/toppar_water_ions.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
# Water aliases
pdbalias residue HOH TIP3
pdbalias atom TIP3 O OH2

foreach seg {A B} {
	segment $seg {
		pdb chainA-$seg.pdb
	}
	coordpdb chainA-$seg.pdb $seg
	regenerate angles dihedrals
}

segment WC {
	auto none
	pdb water_w5p.pdb
}
coordpdb water_w5p.pdb WC

guesscoord
writepdb chainA_wc_psfgen.pdb
writepsf chainA_wc_psfgen.psf

quit
