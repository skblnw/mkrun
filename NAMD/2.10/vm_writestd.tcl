#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "backbone or element P"]
$sel set beta 1
$all writepdb cons_bb_and_P.pdb

$all set beta 0
set sel [atomselect top "name CA or element P"]
$sel set beta 1
$all writepdb cons_CA_and_P.pdb

$all set beta 0
set sel [atomselect top {not segname "WT.*" ION}]
$sel set beta 1
$all writepdb fix_solute.pdb

$all set beta 0
set sel [atomselect top "backbone or element P or segname W1 W2"]
$sel set beta 1
$all writepdb fix_heavy.pdb

#$all set beta 0
#set sel [atomselect top "water"]
#$sel set beta 1
#$all writepdb fix_water.pdb

quit
