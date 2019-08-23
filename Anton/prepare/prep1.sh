#/bin/bash
prefix=md
module load gromacs/2016.1
source  /usr/local/packages/gromacs/2016.1/bin/GMXRC.bash
module load desmond_maestro/2016.3
module load python
python /usr/local/packages/InterMol/intermol/convert.py --gro_in $prefix.gro topol.top --desmond
