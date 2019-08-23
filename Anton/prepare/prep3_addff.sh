#!/bin/bash -i
prefix=md

garden av viparr-ff
garden av v2software
garden load v2software/1.40.0c7/bin
garden load msys/1.7.248c7/bin
garden load viparr/4.6.9-02-psc-c7/bin
garden load viparr/4.6.9-02-psc-c7/lib
garden load viparr/4.6.9-02-psc-c7/lib-python
garden load viparr/4.6.9-02-psc-c7/lib-python3
viparr -f amber99SBstar-ILDN -f tip3p ${prefix}_out.dms out.dms
