#!/home/PHARMACY/chan.773/anaconda3/bin/python3.6
import parmed as pmd
gmx = pmd.load_file('topol.top', xyz='ionized.pdb')
gmx.save('ionized.psf')
