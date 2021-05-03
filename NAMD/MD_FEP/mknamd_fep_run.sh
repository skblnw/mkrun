#!/bin/bash

/opt/namd/NAMD_3.0alpha8_Linux-x86_64-multicore-CUDA/namd3 +p4 equilibrate.namd >& LOG_eq
sleep 10
/opt/namd/NAMD_3.0alpha8_Linux-x86_64-multicore-CUDA/namd3 +p4 fep.soft.namd >& LOG_fep
