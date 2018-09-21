#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

# Protocol:
# for numframes
#     for selections
#         open file for writing
#         [for segments]
#             call measure scripts
#         [end]
#     [end]
# [end]

source proc_calc-dist-btw-memb.tcl

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
#eval file delete [glob output/*.dat]

# Load packages for calculating principle axis (if needed)
#package require Orient 
#namespace import Orient::orient 

set OUTPUT_DIR com
exec mkdir -p $OUTPUT_DIR
set OUTPUT_NAME smd1

# Load your structure and frames
set molnum [mol new ../combine_dcd/initial-noW.psf waitfor all]
mol addfile ../combine_dcd/initial-noW.pdb waitfor all
mol addfile ../combine_dcd/smd5-pf10ps-1.dcd waitfor all

set total_frame [molinfo $molnum get numframes]
for {set nn 0} {$nn < $total_frame} {incr nn} {
    set nframe [expr $nn + 0]
    puts "Frame $nframe"

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment if you need PA
#    set selseg [atomselect top "protein and name CA" frame $nn]
#    set Ip [Orient::calc_principalaxes $selseg]

    # Reset index of selections
    set nsel 1
	
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
	# This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      set outf [open $OUTPUT_DIR/$OUTPUT_NAME.dat "a"]
      # Write TIME at the very first of a line
      set out_line [format "%d" $nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
	  # This determines <number of columns in one output file>
      foreach {seg1 seg2} {1 1} {
          # Call calc funtion you like
          lappend out_line [calc_dist $nn]
      }
      # Write to file
      puts $outf "$out_line"
      # Remember to close the file
      close $outf
	  
	  # Increase index of selections by 1
      incr nsel
    }
}

quit
