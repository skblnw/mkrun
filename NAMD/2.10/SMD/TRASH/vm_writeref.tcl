#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
$all set occupancy 0
set sel [atomselect top "protein"]
$sel set beta 1
$all writepdb protein.ref

quit
