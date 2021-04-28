#!/usr/bin/tclsh
#################################
## TCL script for writepdb from selection
## ukevi Apr 2014
## USAGE: vmd -dispdev text -e writepdb.tcl
#################################

set EXTENSION ref

mol new ../ionized.pdb type pdb waitfor all
set all [atomselect top "all"]

$all set beta 0
set sel [atomselect top "resname DOPS and name N C12 C13 O13A O13B C11 P O13 O14 O12 O11 and z > 0"]
$sel set beta 1
$all writepdb restraints/dops_head_upper.$EXTENSION

$all set beta 0
set sel [atomselect top "resname DOPS and name N C12 C13 O13A O13B C11 P O13 O14 O12 O11 and z < 0"]
$sel set beta 1
$all writepdb restraints/dops_head_lower.$EXTENSION

$all set beta 0
set sel [atomselect top "resname DOPC and name N C12 C13 C14 C15 C11 P O13 O14 O12 O11 and z > 0"]
$sel set beta 1
$all writepdb restraints/dopc_head_upper.$EXTENSION

$all set beta 0
set sel [atomselect top "resname DOPC and name N C12 C13 C14 C15 C11 P O13 O14 O12 O11 and z < 0"]
$sel set beta 1
$all writepdb restraints/dopc_head_lower.$EXTENSION

quit
