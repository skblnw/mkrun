#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

set PREFIX cons

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "backbone or element P"]
#set sel [atomselect top "backbone or name P1 P2 P3"]
$sel set beta 1
$all writepdb ${PREFIX}_bb_and_P.pdb

$all set beta 0
set sel [atomselect top "name CA or element P"]
#set sel [atomselect top "backbone or name P1 P2 P3"]
$sel set beta 1
$all writepdb ${PREFIX}_CA_and_P.pdb

quit
