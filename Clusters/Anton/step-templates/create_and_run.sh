#!/bin/bash

PREFIX=
NN=

cat >> step1_create.sh << EOF
#/bin/bash -u -e -x -o pipefail

JOBNAME=$PREFIX
DMSFILE=out.dms
WD=/anton2fs/raw/\$USER/\$JOBNAME/workdir.$NN
ASW=v2software/1.40.0c7/bin

anton2 create \$JOBNAME --workdir \$WD --sw \$ASW --include psc.ark  --cfg boot.file=\$DMSFILE  \$EXTRACFG
EOF

cat >> step2_prepar.sh << EOF
#/bin/bash -u -e -x -o pipefail

JOBNAME=$PREFIX
WD=/anton2fs/raw/\$USER/\$JOBNAME/workdir.$NN
ASW=v2software/1.40.0c7/bin

anton2 prep \$WD
EOF

cat >> step3_submit.sh << EOF
#/bin/bash -u -e -x -o pipefail

JOBNAME=$PREFIX

WD=/anton2fs/raw/\$USER/\$JOBNAME/workdir.$NN

ASW=v2software/1.40.0c7/bin
anton2 submit \$WD --send-email=always
EOF

chmod +x step*.sh
#./step1_create.sh
#./step2_prepar.sh
#./step3_submit.sh
