##############################################################
# FEP SCRIPT
# Jerome Henin <jerome.henin@uhp-nancy.fr>
#
# Example:
#
# source fep.tcl
#
# fep                  on
# fepFile              system.fep
# fepCol               B
# fepOutFreq           10
# fepOutFile           system.fepout
# fepEquilSteps        500
#
# set nSteps      5000
# set init {0 0.05 0.1}
# set end {0.9 0.95 1.0}
#
# runFEPlist $init $nSteps
# runFEP 0.1 0.9 0.1 $nSteps
# runFEPlist $end $nSteps
##############################################################

##############################################################
# proc runFEPlist { lambdaList nSteps }
#
# Run n FEP windows joining (n + 1) lambda-points
##############################################################

proc runFEPlist { lambdaList nSteps } {
    # Keep track of window number
    global win
    if {![info exists win]} {
	set win 1
    }

    set l1 [lindex $lambdaList 0]
    foreach l2 [lrange $lambdaList 1 end] {

	print [format "Running FEP window %3s: Lambda1 %-6s Lambda2 %-6s \[dLambda %-6s\]"\
		$win $l1 $l2 [expr $l2 - $l1]]
	firsttimestep 0
	alchlambda       $l1
	alchlambda2      $l2
	run $nSteps

	set l1 $l2
	incr win
    }
}


##############################################################
# proc runFEP { start stop dLambda nSteps }
#
# FEP windows of width dLambda between values start and stop
##############################################################

proc runFEP { start stop dLambda nSteps } {
    set ll " $start "
    set l2 [expr $start + $dLambda]

    # A small workaround for numerical rounding errors
    while { [expr {$l2 <= ($stop + 1e-15) } ] } {
	lappend ll $l2
	set l2 [expr {$l2 + $dLambda} ]
    } 

    runFEPlist $ll $nSteps
}
