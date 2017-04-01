#########################################
## Description: A bash script for insane.py
## Author: Kevin@CityUHK 31 Mar 2016
## Usage: 
## Input:
## Output:
## Units: 
## Other Notes: 
#########################################

input_pdb=
output_gro=system.gro
output_top=insane.top
pbc=cubic
# e.g. 10,10,10
box=30,30,30
composition="-l DOPC:4 -l DOPE:3 -l POPS:2 -l PAP2:1"
dm=7
conc=0.18

# Membrane only
python insane_bmw.py -o $output_gro -p $output_top -pbc rectangular -box 50,12,30 $composition

# Membrane + Protein
python insane.py -f $input_pdb -o $output_gro -p $output_top -pbc $pbc -box $box $compostion -center -dm $dm -sol W -salt $conc -charge auto