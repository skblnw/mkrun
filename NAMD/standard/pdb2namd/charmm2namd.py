#!/usr/bin/python

# This script converts a CHARMm psf file to a NAMD (= XPLOR) psf file. 
# USAGE: charmm2namd psf-file rtf-file  
#
# EPF 04-04-2000
#

import sys,struct,os,time

if len(sys.argv) < 3:
  sys.stderr.write ("USAGE: charmm2namd psf-file rtf-file > psf-file\n")
  sys.stderr.write ("       The rtf-file must contain the correct MASS entries for all\n")
  sys.stderr.write ("       atom types present in the psf-file.\n")  
  os._exit(0)
 
psffile = sys.argv[1]
rtffile = sys.argv[2]

# Read the MASS entries from the rtf-file.

f=open(rtffile)

buffer = " "
potentials = []

while buffer != "":
  buffer = f.readline()
  if buffer != "":
    if buffer[0:4] == "MASS":
      potentials.append(buffer.split(None,4))  
     
f.close()

# Now open the psf file and do the conversion of the potential types. 

f=open(psffile)


buffer = f.readline()
# search for the line stating number of title lines
while buffer != "" and "!NTITLE" not in buffer:
  sys.stdout.write(buffer)
  buffer = f.readline()
 
titlelines = int(buffer[:buffer.find("!NTITLE")])
print '%8d !NTITLE' % (titlelines+1) 

# The title lines
for i in range(titlelines):
  buffer = f.readline()
  sys.stdout.write(buffer)

# add extra title line
print '* Converted to NAMD format using '+rtffile+' on '+time.ctime(time.time()) 

# search for the line stating number of atoms
while buffer != "" and "!NATOM" not in buffer:
  buffer = f.readline()
  sys.stdout.write(buffer)
      
numofatoms = float(buffer[:buffer.find("!NATOM")])


j = 0;
while buffer != "":
  buffer = f.readline()
  if buffer != "" and j < numofatoms:
    line = buffer.split(None)
   
    # now search through potential list 

    oldpot = line[5]
    for p in potentials:
      if p[1] == line[5]:
        line[5] = p[2]
    if oldpot == line[5]:
      sys.stderr.write('\nFatal error: no potential for atom type '+oldpot+'\n\n')
      os._exit(0)

    linetuple=(int(line[0]),line[1],line[2],line[3],line[4],line[5],float(line[6]),float(line[7]),int(line[8]))        
             
    psfline =  '   %5d %-4s %-4s %-4s %-4s %4s  %9.6f      %8.4f           %d' % linetuple
    print psfline
   
    j += 1
  else:
    sys.stdout.write(buffer)
             
f.close()  
sys.stderr.write('\nConverted '+str(j)+' atoms.\n\n')
