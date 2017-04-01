#!/usr/bin/env python
import sqlite3
import shutil
import argparse
from argparse import RawTextHelpFormatter

def read_pdb_coor(pdbfile):
    fin = open(pdbfile, 'r')
    coor = []
    box = []
    for line in fin.readlines():
        line = line.strip('\n')
        linee = line.split()
        if line.startswith("ATOM") or line.startswith("HET"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coor.append([x, y, z])
        if line.startswith("CRYST1"):
            boxx = float(linee[1])
            boxy = float(linee[2])
            boxz = float(linee[3])
            box = [boxx, boxy, boxz]
    return coor, box
 

if __name__ == "__main__":
    des= """update the coordinates in a dms file using coordinates form a PDB file"""
    parser = argparse.ArgumentParser(description=des, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-dms', action="store", dest="dms", required=True,\
            help="input dms")
    parser.add_argument('-o', action="store", dest="outputdms", required=True,\
            help="output dms")
    parser.add_argument('-pdb', action="store", dest="pdb", required=True,\
            help="PDB file where the new coordinates are saved")

    inputarg = parser.parse_args()
    
    dms = inputarg.dms
    outputdms = inputarg.outputdms
    pdb = inputarg.pdb
    shutil.copy(dms, outputdms)

    pdbcoor, boxsize = read_pdb_coor(pdb)

    conn = sqlite3.connect(outputdms)
    c = conn.cursor()

    if len(boxsize) > 0:
        c.execute("""update global_cell set x= %f where id=%d"""%(boxsize[0], 1))
        c.execute("""update global_cell set y= %f where id=%d"""%(boxsize[1], 2))
        c.execute("""update global_cell set z= %f where id=%d"""%(boxsize[2], 3))

    c.execute("select * from particle")
    particle = c.fetchall()
    if len(particle) != len(pdbcoor):
        print "number of atom does not match. In dms %d, in PDB %d"%(len(particle), len(pdbcoor))
        sys.exit()

    for i in range(0, len(pdbcoor)):
        c.execute("""update particle set x= %f, y= %f, z= %f where id=%d"""%(pdbcoor[i][0], pdbcoor[i][1], pdbcoor[i][2], i))
    #update posre_harm
    #CREATE TABLE posre_harm_term (p0 integer, 'x0' float, 'y0' float, 'z0' float, param integer not null);
    c.execute("select p0 from posre_harm_term")
    posharm_id = c.fetchall()
    for ids in posharm_id:
        id = ids[0]
        c.execute("""update posre_harm_term set x0=%f, y0=%f, z0=%f where p0=%d"""%(pdbcoor[id][0], pdbcoor[id][1], pdbcoor[id][2], id ))



    conn.commit()

