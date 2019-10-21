#
# nsimakov@psc.edu  2010
#
puts "This script removes water from stk trajectory"
puts "usage:"
puts "vmd -dispdev text -e nowater.tcl -args system.mae workdir/run2.stk nowater"
puts "  system.mae - system"
puts "  workdir/run2.stk - stk trajectory"
puts "  nowater - is the prefix for the  directory, waterless pdb, system.mae and stk trajectory that is created"



proc remove_water {SysMae TrjSTK Prefix Selection} {
    puts "###################################################"
    puts "Starting Script remove_water"
    puts "SysMae $SysMae"
    puts "TrjSTK $TrjSTK"
    puts "Prefix $Prefix"
    puts "Selection $Selection"
    
    mol new $SysMae
    
    set nowatsel [atomselect top "$Selection"]
    
    if { [file exists ${Prefix}_trj] } {
        puts "ERROR: ${Prefix}_trj exists, delete it or specify a diffrent prefix"
        exit
    }
    file mkdir ${Prefix}_trj
    
    
    puts "Writing mae file without waters"
    if { [file exists ${Prefix}_trj/${Prefix}.mae] } {
        puts "ERROR: ${Prefix}_trj/${Prefix}.mae exists, delete it or specify a diffrent prefix"
        exit
    }
    animate write mae ${Prefix}_trj/${Prefix}.mae beg 0 end 0 sel $nowatsel 0
    #puts [molinfo 0 get numframes]
    puts "Reading/writing stk file"
    
    
    set outstk [open "${Prefix}_trj/${Prefix}.stk" w]
    
    set instk [open $TrjSTK r]
    set idtr 0
    while { [gets $instk dtr] >= 0 } {
        puts "working on $dtr"
        
        set JobStep [string range [file tail [file dirname $dtr]] 0 5]
        puts "Jobstep: $JobStep"
        
        animate read dtr $dtr waitfor all 0
        set frames [expr [molinfo 0 get numframes] - 1]
        puts "read $frames frames"
        
        animate write dtr ${Prefix}_trj/$JobStep.dtr beg 1 end $frames sel $nowatsel 0
        
        puts $outstk "${Prefix}_trj/$JobStep.dtr"
        
        animate delete beg 1 end $frames
        #puts [molinfo 0 get numframes]
        incr idtr
    }
    close $instk
    close $outstk
    puts "Number of job steps read: $idtr"
    
    #
    # Also report it in an external file
    #
    #set outfile [open "report.out" w]
    #puts $outfile "Number of lines: $number"
    #close $outfile
    
    puts "Script remove_water is done"
    puts "###################################################"
}



if { [llength $argv] != 3} {
puts "There are $argc arguments to this script"
puts "The name of this script is $argv0"
if { [llength $argv] > 0} {puts "The other arguments are: $argv" }

exit
}

set SysMae [lindex $argv 0];
set TrjSTK [lindex $argv 1];
set Prefix [lindex $argv 2];

#This is where the water gets selected for VMD, if VMD does not
#remove the waters in your system you need to look at your system and this line:
remove_water $SysMae $TrjSTK $Prefix "not water"

exit