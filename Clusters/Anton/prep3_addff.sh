#!/bin/bash -i
prefix=md

#garden av viparr-ffpublic
#garden av v2software
#garden load v2software/1.40.0c7/bin
#garden load msys/1.7.248c7/bin
#garden load viparr/4.7.12c7/bin
#garden load viparr/4.7.12c7/lib
#garden load viparr/4.7.12c7/lib-python
#garden load viparr/4.7.12c7/lib-python3

garden with -m viparr/4.7.12c7/bin -m viparr-ffpublic/1.0.3c7/data -- viparr -f aa.amber.ff14SB -f na.amber.OL15 -f na.amber.OL3 -f water.tip3p -f ions.amber1234lm_anton.tip3p ${prefix}_out.dms out.dms
#garden with -m viparr/4.7.12c7/bin -m viparr-ffpublic/1.0.3c7/data -- viparr -f aa.amber.ff99SB-ILDN -f water.tip3p ${prefix}_out.dms out.dms
#garden with -m viparr/4.7.12c7/bin -m viparr-ffpublic/1.0.3c7/data -- viparr -f aa.charmm.c36m -f lipid.charmm.c36 -f water.tip3p -f ions.charmm36 ${prefix}_out.dms out.dms
