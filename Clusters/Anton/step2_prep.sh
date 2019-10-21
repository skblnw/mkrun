#/bin/bash -u -e -x -o pipefail

JOBNAME=6I1K
WD=/anton2fs/raw/$USER/$JOBNAME/workdir.1
ASW=v2software/1.40.0c7/bin

anton2 prep $WD
