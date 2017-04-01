#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

set PREFIX fix

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "not {water or ion} or segname W1 W2"]
$sel set beta 1
$all writepdb ${PREFIX}_solute.pdb

#$all set beta 0
#set sel [atomselect top "backbone or element P or segname W1 W2"]
#$sel set beta 1
#$all writepdb ${PREFIX}_heavy.pdb

#$all set beta 0
#set sel [atomselect top "water"]
#$sel set beta 1
#$all writepdb ${PREFIX}_water.pdb

quit