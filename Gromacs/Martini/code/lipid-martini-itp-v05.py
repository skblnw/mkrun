#!/usr/bin/env python

"""lipid-martini-itp.py creates a customized Martini lipid topologies, use lipid-martini-itp.py -h for description"""

__author__  = "Helgi I. Ingolfsson, and Tsjerk A. Wassenaar"
__status__  = "Development"
__version__ = "0.5"
__email__   = "h.i.ingolfsson@rug.nl"

import sys,math

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]


# Description
desc = """
This scripts creates a customized Martini lipid topology based on the head, linker and
tail specification strings provided. The topology follows the standard Martini 2.0 lipid
definitions, lipids with these topologies have been explore e.g. in:
 -S.J. Marrink, A.H. de Vries, A.E. Mark.
  Coarse grained model for semi-quantitative lipid simulations.
  JPC-B, 108:750-760, 2004.
 -S.J. Marrink, H.J. Risselada, S. Yefimov, D.P. Tieleman, A.H. de Vries.
  The MARTINI force field: coarse grained model for biomolecular simulations.
  JPC-B, 111:7812-7824, 2007.
 -S.J. Marrink, A.H. de Vries, T.A. Harroun, J. Katsaras, S.R. Wassall.
  Cholesterol shows preference for the interior of polyunsaturated lipid membranes.
  JACS, 130:10-11, 2008.
 -H.J. Risselada, S.J. Marrink.
  The molecular face of lipid rafts in model membranes.
  PNAS, 105:17367-17372, 2008.
 -S. Baoukina, L. Monticelli, H.J. Risselada, S.J. Marrink, D.P. Tieleman.
  The molecular mechanism of lipid monolayer collapse.
  PNAS, 105:10803-10808, 2008.
 -H.I. Ingolfsson, M.N. Melo, F.J. van Eerden, C. Arnarez, C.A. Lopez, T.A. Wassenaar, 
  X. Periole, A.H. De Vries, D.P. Tieleman, S.J. Marrink. 
  Lipid organization of the plasma membrane. 
  JACS, 136:14554-14559, 2014. 
 -C.A. Lopez, Z. Sovova, F.J. van Eerden, A.H. de Vries, S.J. Marrink. 
  Martini force field parameters for glycolipids. 
  JCTC, 9:1694-1708, 2013.

Description of this script can be found in the following manuscript, please cite:
 -T.A. Wassenaar, H.I. Ingolfsson, R.A. Bockmann, D.P. Tieleman, S.J. Marrink. 
  Computational lipidomics with insane: a versatile tool for generating custom membranes for molecular simulations. 
  JCTC, 150410125128004, 2015.

WARNING:
  This script can generate topologies for numerous lipids many of which are unrealistic
  and untested, please use with discretion. E.g. none of the sphingomyelin lipids have been thoroughly tested.


The lipid descriptions supported are as follows:

Heads (-alhead): 
  Please provide a list of lipid head beads. The left most bead will be on top, they are
  connected in a sequence from left to right and the right most bead is connected to the
  first bead in the linker. Each bead is connected with a bond (R_b = 0.47 and K_b = 1250). 
  <Warning some charged Martini lipids have used a shorter bond, R_b = 0.37, but here all 
  lipids use the same R_b = 0.47 bond between headgroup beads.>
  There are no angles between the different head beads, but there is an angle
  [last head bead / first linker bead / first bead of first tail] Theta_a = 180 and
  K_a = 25 that helps orient the head with respect to the rest of the lipid. If an empty
  string is provided the lipid will have no head and starts with the linker beads. Spaces
  should separate the different beads; extra spaces are ignored.

  head bead types supported:
    C = NC3 = Choline      - bead Q0, charge +1
    E = NH3 = Ethanolamine - bead Qd, charge +1
    G = GL0 = Glycerol     - bead P4, charge  0
    S = CNO = Serine       - bead P5, charge  0
    P = PO4 = Phosphate    - bead Qa, charge -1
    O = PO4 = Phosphate    - bead Qa, charge -2

  Examples of lipid heads:
    "C P" -> 'NC3 PO4' - PC - PhosphatidylCholine
    "E P" -> 'NH3 PO4' - PE - PhosphatidylEthanolamine
    "G P" -> 'GL0 PO4' - PG - PhosphatidylGlycerol
    "S P" -> 'CNO PO4' - PS - PhosphatidylSerine 
    "P"   -> 'PO4 ---' - PA - Phosphatidic acid (charge -1, use "O" for charge -2) 
    "O"   -> 'PO4 ---' - PA - Phosphatidic acid, one bond, not protonated (charge -2) 
    ""    -> '--- ---' - DG - No head, Diacyl Glycerols if x2 Gly linkers are used (DAG)

  Current version also supports writing PI, PIP, PIP2 and PIP3 headgrous (called PI, P1, 
  P2 and P3). These headgroups are hardcoded based on Cesar Lopez parameters see Lopez 
  et al. 2013 JCTC.  

Linkers (-allink):
  Currently only lipids with Glycerol linker ir sogubgisube backboneare supported. Each 
  linker is connected with a bond R_b = 0.37 and K_b 1250. The number of linkers and tails 
  provided have to match and each linker is connected with its corresponding tail with a 
  bond R_b = 0.47 and K_b = 1250. Additionally if more than one linker is provided an angle
  [last head bead / first linker bead / second linker bead] Theta_a = 120 and K_a = 25 is
  added to support the head / linker / tail orientation. 
  The sphingosine linking (for ceramide and sphingomyeline lipids) is not tested and 
  should be used with care. Only works with x2 linker beads AM1 and AM2. They are defined 
  with the standard bond between them and linker angle to the head. Just remember that AM1 
  and AM2 contain parts of the lipid tail, and normally tail A (AM1 tail) should start
  with a T bead (containing the inital trans double bond) and be shorter than normal.

  linker beads types supported:
    G = GLY = Glycerols      - bead Na with tail but P1 without, charge  0
    A = Sphingosine backbone - always x2 - AM1 with P1 bead and AM2 with P5 bead, charge 0

  Examples of lipid linkers:
    "G G"   -> 'GLY GLY ---' - A glycerol linker
    "A A"   -> 'AM1 AM2'     - A sphingosine backbone

Tails (-altail):
  One lipid tail definition should be provided for each linker, separated with a space;
  extra spaces are ignored. Each tail can have an arbitrary number of tail beads. Tails
  are connected by bonds, first bead to the tail's linker bead and then from left to
  right (R_b = 0.47 and K_b 1250). To fix the tails orientation/dynamics an angle is
  added for each tail bead [tail bead - 1 or linker bead / tail bead / tail bead + 1]
  the Theta_a and K_a depend on the tail definition.  

  tail bead types supported:
    C = corresponding roughly to a linear combination of x4 CH2 groups (CH2-CH2-CH2-CH2).
        Represented with a C1 bead. Angle [X / C / X] Theta_a = 180 and K_a = 25.
    D = corresponding roughly to a linear combination of x2 CH2 and x2 CH groups
        (CH2-CH=CH-CH2), where the double bound is in a cis bond. Represented with a
        C3 bead, except if there is another D before or after then use a C4 bead. For
        the angles the standard [X / D / X] is a Theta_a = 120 and K_a = 45, except if
        the next bead is also D [X / D / D] then use Theta_a = 100 and K_a = 10.
    T = Same as D except with a trans double bond. The angle [X / T / X] is set with
        equilibrium bond angle Theta_a = 180 and K_a = 45. Represented with a C3 bead.
        
  Examples of tails:
    Lyso tails:
    "- CCCC       " - C16-18:0 Lyso
    "- CCCCC      " - C20-22:0 Lyso
    "- CCCCCC     " - C24-26:0 Lyso
    Saturated tails:
    "CC     CC    " - C08-10:0 - diOctanoyl or diDecanoyl
    "CCC    CCC   " - C12-14:0 - diLauric acid or diMyristoyl
    "CCCC   CCCC  " - C16-18:0 - diPalmitic acid or diStearoyl
    "CCCCC  CCCCC " - C20-22:0 - diArachidoyl or diBehenoyl
    "CCCCCC CCCCCC" - C24-26:0 - diLignoceroyl or diHexacosanoyl
    Unsaturated tails:
    "CDC    CDC   " - C14:1(9c) - diMyristoleoyl 
    "CDCC   CDCC  " - C16:1(9c) - diPalmitoleoyl / C18:1(9c) - diOleoyl 
    "CDDC   CDDC  " - C18:2(9c,12c) - diLinoleoyl 
    "CCDCC  CCDCC " - C20:1(11c) - diGondic acid / C22:1(11c) - diErucoyl
    "DDDDC  DDDDC " - C20:4(5c,8c,11c,14c) - diArachidonoyl 
    "DDDDDD DDDDDD" - C22:6(4c,7c,10c,13c,16c,19c) - diDocosahexaenoyl
    "CCCDCC CCCDCC" - C24:1(15c) - diNervonoyl
    Mixed tails:
    "CCCC   CCC   " - C14:0/C16:0 - MP / C14:0/C18:0 - MS
    "CDCC   CCCC  " - C16:0/C18:1(9c) - PO / C18:0/C18:1(9c) - SO
    "CDDC   CCCC  " - C16:0/C18:2(9c,12c) / C18:0/C18:2(9c,12c)
    "DDDDC  CCCC  " - C16:0/C20:4(5c,8c,11c,14c) - PA / C18:0/C20:4(5c,8c,11c,14c) - SA
    "DDDDDD CCCC  " - C16:0/C22:6(4c,7c,10c,13c,16c,19c) / C18:0/C22:6(4c,7c,10c,13c,16c,19c)
    Trans tails:
    "CTCC   CTCC  " - C18:1(9t) - dielaidoyl
    "TCC    CCCC  " - palmytoyl sphingomyeline tail (AM1 contains the start of the tail)

    NOTE: the first tail (tail A) is connected to linker 1 closer to head (this is sn-2 for GLY linker lipids), which is reverse order
    compared to how regular lipid names are written. The second tail is tail B (for GLY linker lipids this is sn-1)

Use:
  ./lipid-martini-itp-v05.py -alhead 'C P' -allink 'G G' -altail "CDCC CCCC" -alname POPC -o POPC-lipid.itp
"""

# Options
options = [
"""
Options:""",
("-o",       Option(str,    1,        "Martini-lipid.itp", "Output speciffic Martini lipid topology")),
("-alname",  Option(str,    1,        "POPC", "Four letter lipid name")),
("-alhead",  Option(str,    1,        "C P", "Lipid heads, see description")),
("-allink",  Option(str,    1,        "G G", "Lipid linkers, see description")),
("-altail",  Option(str,    1,        "CDCC CCCC", "Lipid tails, see description")),
("-name",    Option(str,    1,        "POPC", "A common name of the lipid, only use in comments")),
("-desc",    Option(str,    1,        "This is a ...", "A general description of what the FF is / represents, only use in comments")),
("-keyw",    Option(str,    1,        "", "List of keywords, only use in comments")),
("-parm",    Option(str,    1,        "Was modeled on ...", "Fow the FF was parameterized, only use in comments")),
("-refs",    Option(str,    1,        "", "List of references for the FF, only use in comments")),
("-crea",    Option(str,    1,        "", "FF created on, only use in comments")),
("-auth",    Option(str,    1,        "", "FF author, only use in comments")),
("-modi",    Option(str,    1,        "", "List of modifications to the FF, only use in comments")),
("-area",    Option(str,    1,        "", "Reference area per lipid, only use in comments")),
("-warn",    Option(str,    1,        "", "Warning(s)/Note(s) for the FF, only use in comments"))
          ]

# Define supported lipid head beads
# Lists all supported head bead types. One letter name mapped to type, atom name and charge
headMapp = {
    "C":  ['Q0', 'NC3', '1.0'],  # NC3 = Choline
    "E":  ['Qd', 'NH3', '1.0'],  # NH3 = Ethanolamine 
    "G":  ['P4', 'GL0', '0.0'],  # GL0 = Glycerol
    "S":  ['P5', 'CNO', '0.0'],  # CNO = Serine
    "P":  ['Qa', 'PO4', '-1.0'], # PO4 = Phosphate
    "O":  ['Qa', 'PO4', '-2.0']  # PO4 = Phosphate (one bond x2 charges can be used e.g. when making unprotonated PA lipids)
    }

# Define possible bond lengths and forces
defBlength = '0.47'
defShortBlength = '0.37'
defBforce = '1250'

# Define possible angles and forces
defAngle1 = '100.0'
defAngle2 = '120.0'
defAngle3 = '180.0'
defAforce1 = '10.0'
defAforce2 = '25.0'
defAforce3 = '45.0'


# Get arguments    
args = sys.argv[1:]

# Print help 
if '-h' in args or '--help' in args:
    print "\n",__file__
    print desc
    for thing in options:
        print type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing
    print
    sys.exit()

# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])

# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])
    
# Get ouput .itp file name
itpFileName  = options["-o"].value

# Get lipid description
lipidHead  = options["-alhead"].value
lipidLinker  = options["-allink"].value
lipidTail  = options["-altail"].value
if lipidLinker==None or lipidLinker==None or lipidTail==None:
    print >>sys.stderr, "You have to provide a header, linker and tail lipid description, if one should be missing provide an empty string"
    sys.exit()
lipidName = options["-alname"].value

lipidCommonName = options["-name"].value
lipidDesc = options["-desc"].value
lipidParm = options["-parm"].value
if lipidCommonName==None or lipidDesc==None or lipidParm==None:
    print >>sys.stderr, "You have to provide a common name, description and list how the FF was parameterized."
    sys.exit()
lCharge = 0  # Update when adding charged beads

progString = "The Martini lipid itp generator version " + __version__ + "  Args are: -o %s -alname %s -alhead '%s' -allink '%s' -altail '%s'" % (itpFileName, lipidName, lipidHead, lipidLinker, lipidTail)
print progString

headsArray = lipidHead.split()
linkersArray = lipidLinker.split()
linkersIndex = []
tailsArray = lipidTail.split()
if len(tailsArray)>len(linkersArray):
    print >>sys.stderr, "A linker definition has to be provided for each tail"
    sys.exit()


bondsArray = []
anglesArray = []
beadArray = []
dihedralsArray = []
constraintsArray = []
exclusionsArray = []

# If speciall head insert now all beads, bonds, angles, dihedrals, constreints etc 
index = 1
if len(headsArray)>0 and headsArray[0]=='PI': # Add PI head 
    # This is from the head of Cesars DPPI parameaters (in the glycolipids.itp)
    # Modified: - bead 4 CP name changed to PO4 - HII
    # Modified: - switch on / of constraints - Helgi, Xavier
    beadArray.append([1, 'P1', 1, lipidName, 'C1 ', 1, 0,    ''])
    beadArray.append([2, 'P4', 1, lipidName, 'C2 ', 2, 0,    ''])
    beadArray.append([3, 'P4', 1, lipidName, 'C3 ', 3, 0,    ''])
    beadArray.append([4, 'Qa', 1, lipidName, 'PO4', 4, -1.0, '; Name changed from CP to PO4'])
    index += 4
    lCharge += -1.0 # Keep track of overall lipid charge
    beadArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    bondsArray.append([-2, '#ifdef FLEXIBLE'])    
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, '0.40', '30000', ''])
    bondsArray.append([1, 3, '0.40', '30000', ''])
    bondsArray.append([2, 3, '0.40', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([1, 4, 0.35, defBforce, ''])
    bondsArray.append([4, 5, defBlength, defBforce, ''])  # This links the head to the linker
    bondsArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    anglesArray.append([3, 1, 4, '133.0', '100.0', ''])
    anglesArray.append([2, 1, 4, '100.0',  '70.0', ''])
    anglesArray.append([-1, 'Orient Head'])
    anglesArray.append([1, 4, 5, '140.0', '30.0', '; link to lipid'])
    anglesArray.append([-1, '4    5    6         2 120.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, '4    5    7         2 180.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])  

    dihedralsArray.append([-1, '3  1  4  5  1  -30.0  5.0  1 ; Removed as it was unstable - WARNING has not been tested'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])    
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.40', ''])
    constraintsArray.append([1, 3, '0.40', ''])
    constraintsArray.append([2, 3, '0.40', ''])
    constraintsArray.append([-2, '#endif'])
elif len(headsArray)>0 and headsArray[0]=='P1': # Add PIP_1 head 
    # This is from the head of Cesars PIP PI3 parameaters (in the glycolipids.itp)
    # Modified: - bead 4 CP name changed to PO4 - Helgi
    # Modified: - bead 2 type changed from Na to P1, oct 2013 - Helgi, Siewert
    # Modified: - switch on / of constraints - Helgi, Xavier
    beadArray.append([1, 'P1', 1, lipidName, 'C1 ', 1, 0,    ''])
    beadArray.append([2, 'P1', 1, lipidName, 'C2 ', 2, 0,    '; corrected particle type (P1 instead of Na), oct 2013']) 
    beadArray.append([3, 'P4', 1, lipidName, 'C3 ', 3, 0,    ''])
    beadArray.append([4, 'Qa', 1, lipidName, 'PO4', 4, -1.0, '; Name changed from CP to PO4'])
    beadArray.append([5, 'Qa', 1, lipidName, 'P1 ', 5, -2.0, ''])
    index += 5
    lCharge += -3.0 # Keep track of overall lipid charge
    beadArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    bondsArray.append([-2, '#ifdef FLEXIBLE'])    
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, '0.40', '30000', ''])
    bondsArray.append([1, 3, '0.40', '30000', ''])
    bondsArray.append([2, 3, '0.40', '30000', ''])
    bondsArray.append([1, 5, '0.40', '25000', ''])
    bondsArray.append([2, 5, '0.30', '30000', ''])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([1, 4, '0.35', defBforce, ''])
    bondsArray.append([4, 6, defBlength, defBforce, ''])  # This links the head to the linker
    bondsArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    anglesArray.append([-1, 'Here we have less angles than in PI, replaced by bonds/constraints'])
    anglesArray.append([-1, 'Orient Head'])
    anglesArray.append([1, 4, 6, '140.0', '25.0', '; link to lipid - PI has 30'])
    anglesArray.append([-1, '4    6    7         2 120.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, '4    6    8         2 180.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])  

    dihedralsArray.append([-1, '3  1  4  6  1  -30.0  5.0  1 ; Removed as it was unstable - WARNING has not been tested'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])    
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.40', ''])
    constraintsArray.append([1, 3, '0.40', ''])
    constraintsArray.append([2, 3, '0.40', ''])
    constraintsArray.append([1, 5, '0.40', ''])
    constraintsArray.append([2, 5, '0.30', ''])
    constraintsArray.append([-2, '#endif'])

elif len(headsArray)>0 and headsArray[0]=='P2': # Add PIP_2 head 
    # This is from the head of Cesars PIP2(3,4) parameaters (in the glycolipids.itp)
    # Modified: - bead 4 CP name changed to PO4 - Helgi
    # Modified: - switch on / of constraints - Helgi, Xavier
    beadArray.append([1, 'P1', 1, lipidName, 'C1 ', 1, 0,    ''])
    beadArray.append([2, 'Na', 1, lipidName, 'C2 ', 2, 0,    '']) 
    beadArray.append([3, 'P4', 1, lipidName, 'C3 ', 3, 0,    ''])
    beadArray.append([4, 'Qa', 1, lipidName, 'PO4', 4, -1.0, '; Name changed from CP to PO4'])
    beadArray.append([5, 'Qa', 1, lipidName, 'P1 ', 5, -2.0, ''])
    beadArray.append([6, 'Qa', 1, lipidName, 'P2 ', 6, -2.0, ''])
    index += 6
    lCharge += -5.0 # Keep track of overall lipid charge
    beadArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    bondsArray.append([-2, '#ifdef FLEXIBLE'])    
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, '0.40', '30000', ''])
    bondsArray.append([1, 3, '0.40', '30000', ''])
    bondsArray.append([2, 3, '0.40', '30000', ''])
    bondsArray.append([2, 5, '0.30', '25000', ''])
    bondsArray.append([2, 6, '0.35', '30000', ''])
    bondsArray.append([1, 5, '0.40', '25000', ''])
    bondsArray.append([3, 6, '0.31', '30000', ''])
    bondsArray.append([-1, '5  6  1  0.60  25000 ; Always keep as bond for stability'])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([5, 6, '0.60', '25000', '; Always keep as bond for stability'])
    bondsArray.append([1, 4, '0.35', defBforce, ''])
    bondsArray.append([4, 7, defBlength, defBforce, ''])  # This links the head to the linker
    bondsArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    anglesArray.append([-1, 'Here we have less angles than in PI, replaced by bonds/constraints'])
    anglesArray.append([-1, 'Orient Head'])
    anglesArray.append([1, 4, 7, '140.0', '25.0', '; link to lipid - PI has 30'])
    anglesArray.append([-1, '4  7  8  2  120.0 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, '4  7  9  2  180.0 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])  

    dihedralsArray.append([-1, '3  1  4  7  1  -30.0  5.0  1 ; Removed as it was unstable - WARNING has not been tested'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])    
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.40', ''])
    constraintsArray.append([1, 3, '0.40', ''])
    constraintsArray.append([2, 3, '0.40', ''])
    constraintsArray.append([2, 5, '0.30', ''])
    constraintsArray.append([2, 6, '0.35', ''])
    constraintsArray.append([1, 5, '0.40', ''])
    constraintsArray.append([3, 6, '0.31', ''])
    constraintsArray.append([-1, '5     6         1 0.60 ; Always keep as bond for stability'])
    constraintsArray.append([-2, '#endif'])

elif len(headsArray)>0 and headsArray[0]=='P3': # Add PIP_3 head 
    # This is derived from the heads of Cesars PIP2(3,4) and PIP1 parameaters (in the glycolipids.itp)
    # The third phosphate was placed by Helgi and Xavier  
    # Modified: - bead 4 CP name changed to PO4 - Helgi
    # Modified: - switch on / of constraints - Helgi, Xavier
    beadArray.append([1, 'P1', 1, lipidName, 'C1 ', 1, 0,    ''])
    beadArray.append([2, 'Na', 1, lipidName, 'C2 ', 2, 0,    '']) 
    beadArray.append([3, 'P4', 1, lipidName, 'C3 ', 3, 0,    ''])
    beadArray.append([4, 'Qa', 1, lipidName, 'PO4', 4, -1.0, '; Name changed from CP to PO4'])
    beadArray.append([5, 'Qa', 1, lipidName, 'P1 ', 5, -2.0, ''])
    beadArray.append([6, 'Qa', 1, lipidName, 'P2 ', 6, -2.0, ''])
    beadArray.append([7, 'Qa', 1, lipidName, 'P3 ', 7, -2.0, '; New P3 bead'])
    index += 7
    lCharge += -7.0 # Keep track of overall lipid charge
    beadArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    bondsArray.append([-2, '#ifdef FLEXIBLE'])    
    bondsArray.append([-1, 'Using bonds not constraints'])
    bondsArray.append([1, 2, '0.40', '30000', ''])
    bondsArray.append([1, 3, '0.40', '30000', ''])
    bondsArray.append([2, 3, '0.40', '30000', ''])
    bondsArray.append([2, 5, '0.30', '25000', ''])
    bondsArray.append([2, 6, '0.35', '30000', ''])
    bondsArray.append([1, 5, '0.40', '25000', ''])
    bondsArray.append([3, 6, '0.31', '30000', ''])
    bondsArray.append([-1, '5  6  1  0.60  25000 ; Always keep as bond for stability'])
    bondsArray.append([1, 7, '0.40', '25000', '; New to add P3 bead - not messured just placed'])
    bondsArray.append([3, 7, '0.30', '30000', '; New to add P3 bead - not messured just placed'])
    bondsArray.append([-1, '6  7  1  0.60  25000 ; New to add P3 bead - not messured just placed'])
    bondsArray.append([-2, '#endif'])
    bondsArray.append([5, 6, '0.60', '25000', '; Always keep as bond for stability'])
    bondsArray.append([6, 7, '0.60', '25000', '; New to add P3 bead - not messured just placed, added just as bonds to increase stability'])
    bondsArray.append([1, 4, '0.35', defBforce, ''])
    bondsArray.append([4, 8, defBlength, defBforce, ''])  # This links the head to the linker
    bondsArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])

    anglesArray.append([-1, 'Here we have less angles than in PI, replaced by bonds/constraints'])
    anglesArray.append([-1, 'Orient Head'])
    anglesArray.append([1, 4, 8, '140.0', '25.0', '; link to lipid - PI has 30'])
    anglesArray.append([-1, '4    8    9         2 120.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, '4    8   10         2 180.00 25.0 ; These are part of the default lipids rules but not used here'])  
    anglesArray.append([-1, 'Tail part (uses standar Martini v2.0 tail rules)'])  

    dihedralsArray.append([-1, '3  1  4  8  1  -30.0  5.0  1 ; Removed as it was unstable - WARNING has not been tested'])

    constraintsArray.append([-2, '#ifndef FLEXIBLE'])    
    constraintsArray.append([-1, 'Using constraints not bonds'])
    constraintsArray.append([1, 2, '0.40', ''])
    constraintsArray.append([1, 3, '0.40', ''])
    constraintsArray.append([2, 3, '0.40', ''])
    constraintsArray.append([2, 5, '0.30', ''])
    constraintsArray.append([2, 6, '0.35', ''])
    constraintsArray.append([1, 5, '0.40', ''])
    constraintsArray.append([3, 6, '0.31', ''])
    constraintsArray.append([-1, '5  6  1  0.60 ; Always keep as bond for stability'])
    constraintsArray.append([1, 7, '0.40', '; New to add P3 bad - not messured just placed'])
    constraintsArray.append([3, 7, '0.30', '; New to add P3 bad - not messured just placed'])
    constraintsArray.append([-1, '6  7  1  0.60 ; New to add P3 bad - not messured just placed, added just as bonds to increase stability'])
    constraintsArray.append([-2, '#endif'])

elif len(headsArray)>0: # Else build head (works for simple PC, PE, PA, PG, etc heads)
    # Start head beads
    head_charge = 0
    for cHead in headsArray:
        head_charge += int(float(headMapp[cHead][2]))

    for cHead in headsArray:
        beadArray.append([index, headMapp[cHead][0], 1, lipidName, headMapp[cHead][1], index, headMapp[cHead][2], ''])
        #print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%s' % (index, headMapp[cHead][0], 1, lipidName, headMapp[cHead][1], index, headMapp[cHead][2])
        lCharge += float(headMapp[cHead][2]) # Keep track of overall lipid charge
        if index > 1: # link head beads
            if head_charge==0:
                bondsArray.append([index - 1, index, defBlength, defBforce, ''])
            else:
                bondsArray.append([index - 1, index, defBlength, defBforce, ''])
                # Keep same 0.47 distance for now
                #bondsArray.append([index - 1, index, defShortBlength, defBforce, ''])
        index += 1
    
    # This links the head to the linker (that has to be some linker beads so the next bead is a linker)
    if head_charge == 0:
        bondsArray.append([index - 1, index, defBlength, defBforce, ''])
    else:
        bondsArray.append([index - 1, index, defBlength, defBforce, ''])
        # Keep same 0.47 distance for now
        #bondsArray.append([index - 1, index, defShortBlength, defBforce])

    # Orient lipid head, add angles between linkers and head + head-first linker-first tail
    # Add angle between the last head bead and first x2 linkers
    if len(linkersArray) > 1:
        if len(tailsArray[0]) == 1 and tailsArray[0][0]=='-':  # First tail is missing (like in LPC) keep <head - linker_1 - linker_2> straight
            anglesArray.append([index - 1, index, index + 1, defAngle3, defAforce2, ''])
        else:  # All normal cases <head - linker_1 - linker_2> have a defAngle2
            anglesArray.append([index - 1, index, index + 1, defAngle2, defAforce2, ''])
    # Add angle down to tail
    if len(linkersArray) > 1 and len(tailsArray[0]) == 1 and tailsArray[0][0]=='-':  
        # Special case first tail is missing (like in LPC) in stead of aligning <last head, linker_1, first tail> allign <linker_1, linker_2, first tail>
        anglesArray.append([index, index+1, index+len(linkersArray), defAngle3, defAforce2, ''])
    else: 
        # Normal case, add angle between the last head bead, the first linker and the first tail
        anglesArray.append([index - 1, index, index+len(linkersArray), defAngle3, defAforce2, ''])

else: #If -alhead was empty len(headsArray)==0 then no head to add (DAG, CER etc)
    pass

# End Add heads 
headIndex = index - 1 # To know what was the last headbead


for linkerBeadIndex in range(0, len(linkersArray)):
    cLinker = linkersArray[linkerBeadIndex]
    if (cLinker != "G") and (cLinker != "A"): # GLY or AM*
        print >>sys.stderr, "This script currently only supports GLY or AM linkers"
        sys.exit()
   
    if (cLinker == "G"): # GLY
        if len(tailsArray[linkerBeadIndex]) == 1 and tailsArray[linkerBeadIndex][0]=='-':  # No tail given for current linker bead so free OH group Na changed to P1
            beadArray.append([index, 'P1', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0, ''])  # Example of this is in LPC (it as not tail A so GL1 changes to P1 type
        else:  # For all other lipids
            #print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%i' % (index, 'Na', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0)
            beadArray.append([index, 'Na', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0, ''])
    else: # A linker
        if linkerBeadIndex == 0:
            #print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%i' % (index, 'P1', 1, lipidName, "AM1", index, 0)
            beadArray.append([index, 'P1', 1, lipidName, "AM1", index, 0, ''])
        elif linkerBeadIndex == 1:
            #print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%i' % (index, 'P5', 1, lipidName, "AM2", index, 0)
            beadArray.append([index, 'P5', 1, lipidName, "AM2", index, 0, ''])
        else:
            print >>sys.stderr, "This script only supports x2 A linkers (making a sphingosine backbone)"
            sys.exit()
    #lCharge += 0 # Keep track of overall lipid charge, not needed as current linkers are all uncharged

    linkersIndex.append(index)
    if index > 1: # There have been beads before (don't do anything if this was the first linker bead and no head)
        if (index -1) > headIndex:  # Then this a linker - linker bond
            bondsArray.append([index - 1, index, defShortBlength, defBforce, ''])
    index += 1
# End linkersArray loop


tailIndex = 0
indexToLetter = "A B C D E F G H I J K L M N".split()
for cTail in tailsArray:   

    # If tail is empty (-) skip this section (like in LPC with no GL1 tail A)
    if len(cTail) == 1 and cTail[0]=='-':
        tailIndex += 1
        continue

    # Add bond from current tail to assosiated linker
    bondsArray.append([linkersIndex[tailIndex], index, defBlength, defBforce, ''])

    for cTailBeadIndex in range(0, len(cTail)):
        cTailBead = cTail[cTailBeadIndex]
        if cTailBead=='C':
            bType = 'C1'
        elif cTailBead=='T':
            bType = 'C3'
        elif cTailBead=='D':
            if (cTailBeadIndex > 0 and cTail[cTailBeadIndex - 1] == 'D') or (len(cTail) > (cTailBeadIndex + 1) and cTail[cTailBeadIndex + 1] == 'D'):
                # Another D before or after
                bType = 'C4'
            else:
                bType = 'C3'
        else:
            print >>sys.stderr, "Tail definition \"%s\" not recognized" % (cTailBead)
            sys.exit()
        
        atomName = cTailBead + str(cTailBeadIndex+1) + indexToLetter[tailIndex]       
        #print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%i' % (index, bType, 1, lipidName, atomName, index, 0)
        beadArray.append([index, bType, 1, lipidName, atomName, index, 0, ''])
        #lCharge += 0 # Keep track of overall lipid charge, not needed as current tails are all uncharged

        # Add bond between tail beads
        if cTailBeadIndex > 0:
            bondsArray.append([index-1, index, defBlength, defBforce, ''])

        # Add angles to support the tails (only add angles when considdering the middle bead in the angle)
        if (cTailBeadIndex + 1) < len(cTail): # Else we are looking at the last bead (which can't have an angle)
            # Set angle constreins (default regular 'C' angle)
            cdefAngle = defAngle3
            cdefAforce = defAforce2
            if cTail[cTailBeadIndex]=='D': 
                if cTail[cTailBeadIndex + 1]=='D':  # has another D right after
                    cdefAngle = defAngle1
                    cdefAforce = defAforce1                    
                else: # does not have another D after
                    cdefAngle = defAngle2
                    cdefAforce = defAforce3
            elif cTail[cTailBeadIndex]=='T':  # Trans double bond
                cdefAngle = defAngle3
                cdefAforce = defAforce3

            if cTailBeadIndex == 0: # top angle connecting to the tail linker [linker / head tail / head+1 tail]
                anglesArray.append([linkersIndex[tailIndex], index, index + 1, cdefAngle, cdefAforce, ''])
            else: # regular angle connecting to the tail linker [current-1 tail / current tail / current+1 tail bead]
                anglesArray.append([index - 1, index, index + 1, cdefAngle, cdefAforce, ''])             
        # end angle stuff
         
        index += 1
    tailIndex += 1
# End tailsArray loop    

# Make .itp headder
itpFile = open(itpFileName,"w")
print >>itpFile, ';;;;;; Martini lipid topology for ' + lipidCommonName + ', generated using:'
print >>itpFile, '; ' + progString
print >>itpFile, '; WARNING: Lipids topology was generated following the Martini 2.0 guidelines but this specific lipid type might'
print >>itpFile, ';          not have been tested and should therefore be used with care. \n;'
print >>itpFile, '; Description:'
print >>itpFile, ';   ' + lipidDesc.replace('\\n','\n;  ')
current = options["-keyw"].value
if current!=None and len(current) > 0:
    print >>itpFile, ';@Keywords: '+current
print >>itpFile, '; Parameterization:'
print >>itpFile, ';   ' + lipidParm.replace('\\n','\n;  ')
current = options["-refs"].value
if current!=None and len(current) > 0:
    print >>itpFile, '; Reference(s): '
    print >>itpFile, ';   ' + current.replace('\\n','\n;  ')
current = options["-crea"].value
if current!=None and len(current) > 0:
    print >>itpFile, '; Created: ' + current
current = options["-auth"].value
if current!=None and len(current) > 0:
    print >>itpFile, '; Author(s): ' + current
current = options["-modi"].value
if current!=None and len(current) > 0:
    print >>itpFile, '; Modified:'
    print >>itpFile, ';   ' + current.replace('\\n','\n;  ')
lArea = options["-area"].value
if lArea!=None and len(lArea) > 0:
    print >>itpFile, '; Reference area per lipid: ' + lArea + ' nm^2'
current = options["-warn"].value
if current!=None and len(current) > 0:
    print >>itpFile, '; Warning(s)/Note(s):'
    print >>itpFile, ';   ' + current.replace('\\n','\n;  ')
print >>itpFile, ';'

# Add INSANE input string
current = '@INSANE alhead='+lipidHead+', allink='+lipidLinker+', altail='+lipidTail+', alname='+lipidName
if lCharge!=None:
    current += ', charge='+str(lCharge)
if lArea!=None and len(lArea) > 0:
    current += ', area='+lArea
print >>itpFile, ';' + current

# Add @RESNTEST, test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
resntest = ""
cutoflen = 3
if len(lipidName) < cutoflen: # fix test for short lipid names
    cutoflen = len(lipidName)
if len(headsArray)>0 and headsArray[0]=='PI': # Add PI head 
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[4]==GL1'
elif len(headsArray)>0 and headsArray[0]=='P1': # Add PIP_1 head 
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[4]==P1'
elif len(headsArray)>0 and headsArray[0]=='P2': # Add PIP_2 head 
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[5]==P2'
elif len(headsArray)>0 and headsArray[0]=='P3': # Add PIP_3 head 
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4] +' and atoms[6]==P3'
elif len(headsArray)>0: # Else build head (works for simple PC, PE, PA, PG, etc heads)
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4]
else: #If -alhead was empty len(headsArray)==0 then no head to add (DAG, CER etc)
    resntest=lipidName[0:cutoflen]+'=='+lipidName+' if: atoms[0]=='+ beadArray[0][4]
# Add test if using x3 bead resnames how to fine the last letter (e.g. POP is it POPC, POPE, POPS etc)
print >>itpFile, ';@RESNTEST '+resntest 

# Add all beads
sBeads = ""
beadNameDict = {}
for cBead in beadArray:
    if cBead[0] > 0:
        sBeads += cBead[4].strip()+" "
        beadNameDict[cBead[0]] = cBead[4]
print >>itpFile, ';@BEADS '+sBeads

# Add all bonds
sBonds = ""
for cBond in bondsArray:
    if cBond[0] > 0:
        sBonds += beadNameDict[cBond[0]].strip()+"-"+beadNameDict[cBond[1]].strip()+" "
print >>itpFile, ';@BONDS '+sBonds

print >>itpFile, ';'
print >>itpFile, ''
print >>itpFile, '[moleculetype]'
print >>itpFile, '; molname      nrexcl'
print >>itpFile, '  ' + lipidName + '          1'

# Write beads
print >>itpFile, '\n[atoms]'
print >>itpFile, '; id 	type 	resnr 	residu 	atom 	cgnr 	charge'
for cBead in beadArray:
    if cBead[0] > 0:
        print >>itpFile, '  %2i \t%s \t%2i \t%s \t%s \t%2i \t%s \t%s' % (cBead[0], cBead[1], cBead[2], cBead[3], cBead[4], cBead[5], cBead[6], cBead[7])
    elif cBead[0] == -1:  # Regular comment 
        print >>itpFile, '; ' + cBead[1]
    elif cBead[0] == -2:  # gromacs system line
        print >>itpFile, cBead[1]

# Write lipid bonds
print >>itpFile, '\n[bonds]'
print >>itpFile, ';  i  j 	funct 	length 	force.c.'
for cBond in bondsArray:
    if cBond[0] > 0:
        print >>itpFile, '  %2i %2i \t1 \t%s \t%s \t%s' % (cBond[0], cBond[1], cBond[2], cBond[3], cBond[4])
    elif cBond[0] == -1:  # Regular comment 
        print >>itpFile, '; ' + cBond[1]
    elif cBond[0] == -2:  # gromacs system line
        print >>itpFile, cBond[1]
    
# Write lipid angles
print >>itpFile, '\n[angles]'
print >>itpFile, ';  i  j  k 	funct 	angle 	force.c.'
for cAngle in anglesArray:
    if cAngle[0] > 0:
         print >>itpFile, '  %2i %2i %2i \t2 \t%s \t%s \t%s' % (cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4], cAngle[5])    
    elif cAngle[0] == -1:  # Regular comment 
        print >>itpFile, '; ' + cAngle[1]
    elif cAngle[0] == -2:  # gromacs system line
        print >>itpFile, cAngle[1]
    
# Write lipid dihedrals
if len(dihedralsArray) > 0: 
    print >>itpFile, '\n[dihedrals]'
    print >>itpFile, ';  i  j  k  l 	funct 	angle 	force.c.'
    for cDihedral in dihedralsArray:
        if cDihedral[0] > 0:
            print >>itpFile, '  %2i %2i %2i %2i \t%i \t%s \t%s \t%s' % (cDihedral[0], cDihedral[1], cDihedral[2], cDihedral[3], cDihedral[4], cDihedral[5], cDihedral[6], cDihedral[7])    
        elif cDihedral[0] == -1:  # Regular comment 
            print >>itpFile, '; ' + cDihedral[1]
        elif cDihedral[0] == -2:  # gromacs system line
            print >>itpFile, cDihedral[1]

# Write lipid constraints
if len(constraintsArray) > 0: 
    print >>itpFile, '\n[constraints]'
    print >>itpFile, ';  i  j  k 	funct 	length'
    for cConstraint in constraintsArray:
        if cConstraint[0] > 0:
            print >>itpFile, '  %2i %2i \t1 \t%s \t%s' % (cConstraint[0], cConstraint[1], cConstraint[2], cConstraint[3])    
        elif cConstraint[0] == -1:  # Regular comment 
            print >>itpFile, '; ' + cConstraint[1]
        elif cConstraint[0] == -2:  # gromacs system line
            print >>itpFile, cConstraint[1]        

# Write lipid exclusions
if len(exclusionsArray) > 0: 
    print >>itpFile, '\n[exclusions]'
    print >>itpFile, ';  i  j'
    for cExclusion in exclusionsArray:
        if cExclusion[0] > 0:
            print >>itpFile, '  %2i %2i \t%s' % (cExclusion[0], cExclusion[1], cExclusion[2])    
        elif cExclusion[0] == -1:  # Regular comment 
            print >>itpFile, '; ' + cExclusion[1]
        elif cExclusion[0] == -2:  # gromacs system line
            print >>itpFile, cExclusion[1]

print >>itpFile, ''
itpFile.close()
# End lipid-martini-itp



