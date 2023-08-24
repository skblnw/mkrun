mol delete all
package require psfgen
resetpsf
readpsf popc.psf
coordpdb popc.pdb
readpsf ../5kxi.psf
coordpdb ../5kxi_1.pdb
writepsf 5kxi_popc_raw.psf
writepdb 5kxi_popc_raw.pdb


