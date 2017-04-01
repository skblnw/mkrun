#!/bin/bash

# It is an example for how to run gromacs (4.x.x, single core) 
# for DOPC membrane start from optimization to equilibration,
# as a regular process before the final production run.
# Each line without "#" can be used as individual comment line in Linux.

# For all the example mdp file below, the length of the optimization
# or simulation might be too short. You could certainly increase
# the total number of steps to make it longer.

# optimize the system with flexible water
grompp -c dopc.gro -f bmw_minimize.mdp -p flexible -n index -o em1

mdrun -s em1 -c opt1 -g outem1

# optimize the system again with rigid water
grompp -c opt1.gro -f bmw_minimize.mdp -p bmw_membrane -n index -o em2

mdrun -s em2 -c opt2 -g outem2

# If the results are not optimized good enough, you might want to
# use double precision of gromacs to optimize it again from flexible
# water to rigid water. Make sure the final energy & force is not
# too large

# start equilibration with NVT
grompp -c opt2.gro -f bmw_nvt.mdp -p bmw_membrane.top -n index.ndx -o equ1

mdrun -s equ1 -c endequ1 -g outequ1 -x equ1 -cpo equ1

# run another equilibration run with semiisotropic pressure coupling
# (NPxyPzT).
grompp -c endequ1.gro -f bmw_membrane.mdp -p bmw_membrane.top -n index.ndx -o equ2

mdrun -s equ2 -c endequ2 -g outequ2 -e eequ2 -x equ2 -cpo equ2

# Now you can start the production run with endequ2.gro and equ2.cpt

exit 0
