#/bin/bash -u -e -x -o pipefail

JOBNAME=
DMSFILE=out.dms
#EXTRACFG="--cfg anton.tune.last_time=50000"
WD=/anton2fs/raw/$USER/$JOBNAME/workdir.1
ASW=v2software/1.40.0c7/bin

anton2 create $JOBNAME --workdir $WD --sw $ASW --include psc.ark  --cfg boot.file=$DMSFILE  $EXTRACFG
