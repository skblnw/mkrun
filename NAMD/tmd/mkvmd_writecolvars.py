#!/bin/python

import glob, os
import subprocess
import argparse
from argparse import RawDescriptionHelpFormatter

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
ap.add_argument('psf', type=str, help='PSF file')
ap.add_argument('pdb', type=str, help='Initial PDB file')
ap.add_argument('outfile', type=str, help='Output filename')
ap.add_argument('--tmdlist1', type=int, nargs='+', help='List of TMD residues 1')
ap.add_argument('--tmdlist2', type=int, nargs='+', help='List of TMD residues 2')
ap.add_argument('--reslist1', type=int, nargs='+', help='List of constrained residues 1')
ap.add_argument('--reslist2', type=int, nargs='+', help='List of constrained residues 2')
try:
  cmd = ap.parse_args()
except:
  ap.print_help()
  os._exit(0)

filenames = [cmd.psf, cmd.pdb]
for filename in filenames:
  if not os.path.isfile(filename):
    print(f"{filename} does not exist")
    os._exit(0)

os.system('mkdir -p restraints')

f=open(cmd.outfile,'w')
f.write("""colvarsTrajFrequency 1000
colvarsRestartFrequency 100

"""
  )
f.close()

for nn, (ii, jj) in enumerate(zip(cmd.tmdlist1,cmd.tmdlist2)):
  print(nn)
  f=open('tcl'+str(nn),'w')
  f.write("""
mol new step5_input.psf
mol addfile step7_3.restart.coor
set selall [atomselect 0 all]
set sel0 [atomselect 0 "segname PROP and resid """+str(ii)+""" and name N O"]
$selall set beta 0
$sel0 set beta 1
$selall writepdb restraints/"""+str(nn+1)+"""a.pdb
set sel0 [atomselect 0 "segname PROP and resid """+str(jj)+""" and name N O"]
$selall set beta 0
$sel0 set beta 1
$selall writepdb restraints/"""+str(nn+1)+"""b.pdb
quit"""
    )
  f.close()
  os.system('vmd -dispdev text -e tcl'+str(nn))

  f=open(cmd.outfile,'a')
  f.write("""
colvar {
name pair"""+str(nn+1)+"""
  coordnum {
    group1 {
      atomsFile restraints/"""+str(nn+1)+"""a.pdb
      atomsCol B
    }
    group2 {
      atomsFile restraints/"""+str(nn+1)+"""b.pdb
      atomsCol B
    }
    cutoff 3.3
  }
}

harmonic {
  colvars pair"""+str(nn+1)+"""
  centers 1.2
  targetCenters 2.7
  targetNumSteps """+str((nn+1)*2)+"""0000
  forceConstant 1000
  outputCenters on
}

  """
    )
  f.close()

for nn, (ii, jj) in enumerate(zip(cmd.reslist1,cmd.reslist2)):
  print(nn)
  f=open('tcl'+str(nn),'w')
  f.write("""
mol new step5_input.psf
mol addfile step7_3.restart.coor
set selall [atomselect 0 all]
set sel0 [atomselect 0 "segname PROP and resid """+str(ii)+""" and name N O"]
$selall set beta 0
$sel0 set beta 1
$selall writepdb restraints/res"""+str(nn+1)+"""a.pdb
set sel0 [atomselect 0 "segname PROP and resid """+str(jj)+""" and name N O"]
$selall set beta 0
$sel0 set beta 1
$selall writepdb restraints/res"""+str(nn+1)+"""b.pdb
quit"""
    )
  f.close()
  os.system('vmd -dispdev text -e tcl'+str(nn))

  f=open(cmd.outfile,'a')
  f.write("""
colvar {
name respair"""+str(nn+1)+"""
  coordnum {
    group1 {
      atomsFile restraints/res"""+str(nn+1)+"""a.pdb
      atomsCol B
    }
    group2 {
      atomsFile restraints/res"""+str(nn+1)+"""b.pdb
      atomsCol B
    }
    cutoff 3.3
  }
}

harmonic {
  colvars respair"""+str(nn+1)+"""
  centers 2.4
  forceConstant 50
  outputCenters on
}

  """
    )
  f.close()

for f in glob.glob("tcl*"):
    os.remove(f)
