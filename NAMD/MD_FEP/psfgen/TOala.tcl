resetpsf

package require alchemify
package require psfgen
package require readcharmmtop
package require mutator 1.5
topology d:/VMD/readcharmmtop1.2/top_all36_prot.rtf
topology d:/VMD/readcharmmtop1.2/top_all36_lipid.rtf
topology d:/VMD/readcharmmtop1.2/top_all36_na.rtf
topology d:/VMD/readcharmmtop1.2/top_all36_carb.rtf
topology d:/VMD/readcharmmtop1.2/top_all36_cgenff.rtf
topology d:/VMD/readcharmmtop1.2/toppar_water_ions_namd.str
topology d:/VMD/readcharmmtop1.2/toppar_all36_carb_glycopeptide.str

topology d:/VMD/readcharmmtop1.2/top_all36_prot.rtf
topology d:/VMD/readcharmmtop1.2/toppar_water_ions_namd.str
set notprot [atomselect top "segid HLAA "]
$notprot writepdb chainA1.pdb
$notprot writepsf chainA1.psf
set notprot [atomselect top "segid HLAB "]
$notprot writepdb chainB1.pdb
$notprot writepsf chainB1.psf
set notprot [atomselect top "segid MYO"]
$notprot writepdb chainC1.pdb
$notprot writepsf chainC1.psf
set notprot [atomselect top "not protein and not segid MYO"]
$notprot writepdb notprotein.pdb
$notprot writepsf notprotein.psf


 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb chainC1.pdb

mutate -1 S2A}
coordpdb chainC1.pdb
regenerate angles dihedrals
guesscoord
writepdb S-1A.pdb
writepsf S-1A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb S-1A.pdb

mutate 1 G2A}
coordpdb S-1A.pdb
regenerate angles dihedrals
guesscoord
writepdb G1A.pdb
writepsf G1A.psf
 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb G1A.pdb

mutate 2 E2A}
coordpdb G1A.pdb
regenerate angles dihedrals
guesscoord
writepdb E2A.pdb
writepsf E2A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb E2A.pdb

mutate 3 I2A}
coordpdb E2A.pdb
regenerate angles dihedrals
guesscoord
writepdb I3A.pdb
writepsf I3A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb I3A.pdb

mutate 4 L2A}
coordpdb I3A.pdb
regenerate angles dihedrals
guesscoord
writepdb L4A.pdb
writepsf L4A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb L4A.pdb

mutate 5 L2A}
coordpdb L4A.pdb
regenerate angles dihedrals
guesscoord
writepdb L5A.pdb
writepsf L5A.psf
 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb L5A.pdb

mutate 6 Y2A}
coordpdb L5A.pdb
regenerate angles dihedrals
guesscoord
writepdb Y6A.pdb
writepsf Y6A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb Y6A.pdb

mutate 7 F2A}
coordpdb Y6A.pdb
regenerate angles dihedrals
guesscoord
writepdb F7A.pdb
writepsf F7A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb F7A.pdb

mutate 8 I2A}
coordpdb F7A.pdb
regenerate angles dihedrals
guesscoord
writepdb I8A.pdb
writepsf I8A.psf
 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb I8A.pdb

mutate 9 H2A}
coordpdb I8A.pdb
regenerate angles dihedrals
guesscoord
writepdb H9A.pdb
writepsf H9A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb H9A.pdb

mutate 10 F2A}
coordpdb H9A.pdb
regenerate angles dihedrals
guesscoord
writepdb F10A.pdb
writepsf F10A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb F10A.pdb

mutate 11 S2A}
coordpdb F10A.pdb
regenerate angles dihedrals
guesscoord
writepdb S11A.pdb
writepsf S11A.psf

 resetpsf  
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
        
    segment MYO {
                pdb S11A.pdb

mutate 12 R2A}
coordpdb S11A.pdb
regenerate angles dihedrals
guesscoord
writepdb R12A.pdb
writepsf R12A.psf

resetpsf
package require psfgen
topology  d:/VMD/readcharmmtop1.2/top_all36_prot.rtf
topology d:/VMD/readcharmmtop1.2/toppar_water_ions_namd.str
topology d:/VMD/readcharmmtop1.2/top_all36_hybrid.inp
proc psfalias {}  {

	#### CHARMM PROTEIN ####

	pdbalias residue HIS HSE
	pdbalias atom ILE CD1 CD


}

psfalias
segment HLAA {
    pdb chainA1.pdb
}
segment HLAB {
    pdb chainB1.pdb
}
segment MYO {
	pdb R12A.pdb
}
coordpdb chainA1.pdb
coordpdb chainB1.pdb
coordpdb R12A.pdb MYO
readpsf notprotein.psf
coordpdb notprotein.pdb
guesscoord 	 
writepdb MHC-iab-myo-ala.pdb 	 
writepsf MHC-iab-myo-ala.psf
