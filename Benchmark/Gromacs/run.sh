#!/bin/bash

bmk=true

if $bmk; then
    gmx mdrun -ntomp 8 -v -deffnm run_md -nsteps 10000
    mv run_md.log 150k.log
    rm -f run_md.*
fi

