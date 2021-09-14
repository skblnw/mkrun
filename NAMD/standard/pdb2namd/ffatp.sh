#!/bin/zsh
source ~/.zshrc

REF_RTF='drawing_3D/drawing_3D.rtf'
LIGAND_NAME='LIG'
LIGAND_CHARGE='-4.00'
LIST="C8 C9 H8"

cat > hybrid.rtf <<'EOF'
* ---- 
* Built RTF for drawing_3D.mol2 
*    by user vzoete     Wed Sep  8 13:04:31 UTC 2021 
* ---- 
*
   22    0

EOF
grep "^MASS" $REF_RTF >> hybrid.rtf
cat >> hybrid.rtf <<'EOF'
!hydrogens
MASS   256 HGA1     1.00800  ! alphatic proton, CH
MASS   257 HGA2     1.00800  ! alphatic proton, CH2
MASS   258 HGA3     1.00800  ! alphatic proton, CH3
MASS   259 HGA4     1.00800  ! alkene proton; RHC=
MASS   260 HGA5     1.00800  ! alkene proton; H2C=CR
MASS   261 HGA6     1.00800  ! aliphatic H on fluorinated C, monofluoro
MASS   262 HGA7     1.00800  ! aliphatic H on fluorinated C, difluoro
MASS   263 HGAAM0   1.00800  ! aliphatic H, NEUTRAL trimethylamine (#)
MASS   264 HGAAM1   1.00800  ! aliphatic H, NEUTRAL dimethylamine (#)
MASS   265 HGAAM2   1.00800  ! aliphatic H, NEUTRAL methylamine (#)
!(#) EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
!on NEUTRAL METHYLAMINE groups, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens
MASS   266 HGP1     1.00800  ! polar H
MASS   267 HGP2     1.00800  ! polar H, +ve charge
MASS   268 HGP3     1.00800  ! polar H, thiol
MASS   269 HGP4     1.00800  ! polar H, neutral conjugated -NH2 group (NA bases)
MASS   270 HGP5     1.00800  ! polar H on quarternary ammonium salt (choline)
MASS   271 HGPAM1   1.00800  ! polar H, NEUTRAL dimethylamine (#), terminal alkyne H
MASS   272 HGPAM2   1.00800  ! polar H, NEUTRAL methylamine (#)
MASS   273 HGPAM3   1.00800  ! polar H, NEUTRAL ammonia (#)
!(#) EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
!on NEUTRAL METHYLAMINE groups, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens
MASS   274 HGR51    1.00800  ! nonpolar H, neutral 5-mem planar ring C, LJ based on benzene
MASS   275 HGR52    1.00800  ! Aldehyde H, formamide H (RCOH); nonpolar H, neutral 5-mem planar ring C adjacent to heteroatom or + charge
MASS   276 HGR53    1.00800  ! nonpolar H, +ve charge HIS he1(+1)
MASS   277 HGR61    1.00800  ! aromatic H
MASS   278 HGR62    1.00800  ! nonpolar H, neutral 6-mem planar ring C adjacent to heteroatom
MASS   279 HGR63    1.00800  ! nonpolar H, NAD+ nicotineamide all ring CH hydrogens
MASS   280 HGR71    1.00800  ! nonpolar H, neutral 7-mem arom ring, AZUL, azulene, kevo
!carbons
MASS   281 CG1T1   12.01100  ! internal alkyne R-C#C
MASS   282 CG1T2   12.01100  ! terminal alkyne H-C#C
MASS   283 CG1N1   12.01100  ! C for cyano group
MASS   284 CG2D1   12.01100  ! alkene; RHC= ; imine C
MASS   285 CG2D2   12.01100  ! alkene; H2C=
MASS   286 CG2D1O  12.01100  ! double bond carbon adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC1.
MASS   287 CG2D2O  12.01100  ! double bond carbon adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC2.
MASS   288 CG2DC1  12.01100  ! conjugated alkenes, R2C=CR2
MASS   289 CG2DC2  12.01100  ! conjugated alkenes, R2C=CR2
MASS   290 CG2DC3  12.01100  ! conjugated alkenes, H2C=
MASS   291 CG2N1   12.01100  ! conjugated C in guanidine/guanidinium
MASS   292 CG2N2   12.01100  ! conjugated C in amidinium cation
MASS   293 CG2O1   12.01100  ! carbonyl C: amides
MASS   294 CG2O2   12.01100  ! carbonyl C: esters, [neutral] carboxylic acids
MASS   295 CG2O3   12.01100  ! carbonyl C: [negative] carboxylates
MASS   296 CG2O4   12.01100  ! carbonyl C: aldehydes
MASS   297 CG2O5   12.01100  ! carbonyl C: ketones
MASS   298 CG2O6   12.01100  ! carbonyl C: urea, carbonate
MASS   299 CG2O7   12.01100  ! CO2 carbon
MASS   300 CG2R51  12.01100  ! 5-mem ring, his CG, CD2(0), trp
MASS   301 CG2R52  12.01100  ! 5-mem ring, double bound to N, PYRZ, pyrazole
MASS   302 CG2R53  12.01100  ! 5-mem ring, double bound to N and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo
MASS   303 CG2R57  12.01100  ! 5-mem ring, bipyrroles
MASS   304 CG2R61  12.01100  ! 6-mem aromatic C
MASS   305 CG2R62  12.01100  ! 6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)
MASS   306 CG2R63  12.01100  ! 6-mem aromatic amide carbon (NA) (and other 6-mem aromatic carbonyls?)
MASS   307 CG2R64  12.01100  ! 6-mem aromatic amidine and guanidine carbon (between 2 or 3 Ns and double-bound to one of them), NA, PYRM
MASS   308 CG2R66  12.01100  ! 6-mem aromatic carbon bound to F
MASS   309 CG2R67  12.01100  ! 6-mem aromatic carbon of biphenyl
MASS   310 CG2RC0  12.01100  ! 6/5-mem ring bridging C, guanine C4,C5, trp
MASS   311 CG2R71  12.01100  ! 7-mem ring arom C, AZUL, azulene, kevo
MASS   312 CG2RC7  12.01100  ! sp2 ring connection with single bond(!), AZUL, azulene, kevo
MASS   313 CG301   12.01100  ! aliphatic C, no hydrogens, neopentane
MASS   314 CG302   12.01100  ! aliphatic C, no hydrogens, trifluoromethyl
MASS   315 CG311   12.01100  ! aliphatic C with 1 H, CH
MASS   316 CG312   12.01100  ! aliphatic C with 1 H, difluoromethyl
MASS   317 CG314   12.01100  ! aliphatic C with 1 H, adjacent to positive N (PROT NTER) (+)
MASS   318 CG321   12.01100  ! aliphatic C for CH2
MASS   319 CG322   12.01100  ! aliphatic C for CH2, monofluoromethyl
MASS   320 CG323   12.01100  ! aliphatic C for CH2, thiolate carbon
MASS   321 CG324   12.01100  ! aliphatic C for CH2, adjacent to positive N (piperidine) (+)
MASS   322 CG331   12.01100  ! aliphatic C for methyl group (-CH3)
MASS   323 CG334   12.01100  ! aliphatic C for methyl group (-CH3), adjacent to positive N (PROT NTER) (+)
MASS   324 CG3AM0  12.01100  ! aliphatic C for CH3, NEUTRAL trimethylamine methyl carbon (#)
MASS   325 CG3AM1  12.01100  ! aliphatic C for CH3, NEUTRAL dimethylamine methyl carbon (#)
MASS   326 CG3AM2  12.01100  ! aliphatic C for CH3, NEUTRAL methylamine methyl carbon (#)
!(#) EXTREME care is required when doing atom typing on compounds that look like this. Use ONLY
!on NEUTRAL METHYLAMINE groups, NOT ETHYL, NOT Schiff Bases, but DO use on 2 out of 3 guanidine nitrogens
MASS   327 CG3C31  12.01100  ! cyclopropyl carbon
MASS   328 CG3C41  12.01100  ! cyclobutyl carbon
MASS   329 CG3C50  12.01100  ! 5-mem ring aliphatic quaternary C (cholesterol, bile acids)
MASS   330 CG3C51  12.01100  ! 5-mem ring aliphatic CH  (proline CA, furanoses)
MASS   331 CG3C52  12.01100  ! 5-mem ring aliphatic CH2 (proline CB/CG/CD, THF, deoxyribose)
MASS   332 CG3C53  12.01100  ! 5-mem ring aliphatic CH  adjacent to positive N (proline.H+ CA) (+)
MASS   333 CG3C54  12.01100  ! 5-mem ring aliphatic CH2 adjacent to positive N (proline.H+ CD) (+)
MASS   334 CG3RC1  12.01100  ! bridgehead in bicyclic systems containing at least one 5-membered or smaller ring
!(+) Includes protonated Shiff base (NG3D5, NG2R52 in 2HPP) but NOT amidinium (NG2R52 in IMIM), guanidinium
!nitrogens
MASS   335 NG1T1   14.00700  ! N for cyano group
!MASS   336 NG1D1   14.00700  ! terminal N in azides, lsk
MASS   337 NG2D1   14.00700  ! N for neutral imine/Schiff's base (C=N-R, acyclic amidine, gunaidine)
MASS   338 NG2S0   14.00700  ! N,N-disubstituted amide, proline N (CO=NRR')
MASS   339 NG2S1   14.00700  ! peptide nitrogen (CO=NHR)
MASS   340 NG2S2   14.00700  ! terminal amide nitrogen (CO=NH2)
MASS   341 NG2S3   14.00700  ! external amine ring nitrogen (planar/aniline), phosphoramidate
!MASS   342 NG2S4   14.00700  ! neutral hydroxamic acid
MASS   343 NG2O1   14.00700  ! NITB, nitrobenzene
MASS   344 NG2P1   14.00700  ! N for protonated imine/Schiff's base (C=N(+)H-R, acyclic amidinium, guanidinium)
MASS   345 NG2R43  14.00700  ! amide in 4-memebered ring (planar), AZDO, lsk
MASS   346 NG2R50  14.00700  ! double bound neutral 5-mem planar ring, purine N7
MASS   347 NG2R51  14.00700  ! single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)
MASS   348 NG2R52  14.00700  ! protonated schiff base, amidinium, guanidinium in 5-membered ring, HIS, 2HPP, kevo
MASS   349 NG2R53  14.00700  ! amide in 5-memebered NON-SP2 ring (slightly pyramidized), 2PDO, kevo
MASS   350 NG2R57  14.00700  ! 5-mem ring, bipyrroles
MASS   351 NG2R60  14.00700  ! double bound neutral 6-mem planar ring, pyr1, pyzn
MASS   352 NG2R61  14.00700  ! single bound neutral 6-mem planar ring imino nitrogen; glycosyl linkage
MASS   353 NG2R62  14.00700  ! double bound 6-mem planar ring with heteroatoms in o or m, pyrd, pyrm
MASS   354 NG2R67  14.00700  ! 6-mem planar ring substituted with 6-mem planar ring (N-phenyl pyridinones etc.)
MASS   355 NG2RC0  14.00700  ! 6/5-mem ring bridging N, indolizine, INDZ, kevo
MASS   356 NG301   14.00700  ! neutral trimethylamine nitrogen
MASS   357 NG311   14.00700  ! neutral dimethylamine nitrogen
MASS   358 NG321   14.00700  ! neutral methylamine nitrogen
MASS   359 NG331   14.00700  ! neutral ammonia nitrogen
MASS   360 NG3C51  14.00700  ! secondary sp3 amine in 5-membered ring
MASS   361 NG3N1   14.00700  ! N in hydrazine, HDZN
MASS   362 NG3P0   14.00700  ! quarternary N+, choline
MASS   363 NG3P1   14.00700  ! tertiary NH+ (PIP)
MASS   364 NG3P2   14.00700  ! secondary NH2+ (proline)
MASS   365 NG3P3   14.00700  ! primary NH3+, phosphatidylethanolamine
!oxygens
MASS   366 OG2D1   15.99940  ! carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, uera
MASS   367 OG2D2   15.99940  ! carbonyl O: negative groups: carboxylates, carbonate
MASS   368 OG2D3   15.99940  ! carbonyl O: ketones
MASS   369 OG2D4   15.99940  ! 6-mem aromatic carbonyl oxygen (nucleic bases)
MASS   370 OG2D5   15.99940  ! CO2 oxygen
MASS   371 OG2N1   15.99940  ! NITB, nitrobenzene
MASS   372 OG2P1   15.99940  ! =O in phosphate or sulfate
MASS   373 OG2R50  15.99940  ! FURA, furan
MASS   374 OG3R60  15.99940  ! O in 6-mem cyclic enol ether (PY01, PY02) or ester
MASS   375 OG301   15.99940  ! ether -O- !SHOULD WE HAVE A SEPARATE ENOL ETHER??? IF YES, SHOULD WE MERGE IT WITH OG3R60???
MASS   376 OG302   15.99940  ! ester -O-
MASS   377 OG303   15.99940  ! phosphate/sulfate ester oxygen
MASS   378 OG304   15.99940  ! linkage oxygen in pyrophosphate/pyrosulphate
MASS   379 OG311   15.99940  ! hydroxyl oxygen
MASS   380 OG312   15.99940  ! ionized alcohol oxygen
MASS   381 OG3C31  15.99940  ! epoxide oxygen, 1EOX, 1BOX, sc
MASS   382 OG3C51  15.99940  ! 5-mem furanose ring oxygen (ether)
MASS   383 OG3C61  15.99940  ! DIOX, dioxane, ether in 6-membered ring !SHOULD WE MERGE THIS WITH OG3R60???
!sulphurs
MASS   384 SG2D1   32.06000  ! thiocarbonyl S
MASS   385 SG2R50  32.06000  ! THIP, thiophene
MASS   386 SG311   32.06000  ! sulphur, SH, -S-
MASS   387 SG301   32.06000  ! sulfur C-S-S-C type
MASS   388 SG302   32.06000  ! thiolate sulfur (-1)
MASS   389 SG3O1   32.06000  ! sulfate -1 sulfur
MASS   390 SG3O2   32.06000  ! neutral sulfone/sulfonamide sulfur
MASS   391 SG3O3   32.06000  ! neutral sulfoxide sulfur
!halogens
MASS   392 CLGA1   35.45300  ! CLET, DCLE, chloroethane, 1,1-dichloroethane
MASS   393 CLGA3   35.45300  ! TCLE, 1,1,1-trichloroethane
MASS   394 CLGR1   35.45300  ! CHLB, chlorobenzene
MASS   395 BRGA1   79.90400  ! BRET, bromoethane
MASS   396 BRGA2   79.90400  ! DBRE, 1,1-dibromoethane
MASS   397 BRGA3   79.90400  ! TBRE, 1,1,1-dibromoethane
MASS   398 BRGR1   79.90400  ! BROB, bromobenzene
MASS   399 IGR1   126.90447  ! IODB, iodobenzene
MASS   400 FGA1    18.99800  ! aliphatic fluorine, monofluoro
MASS   401 FGA2    18.99800  ! aliphatic fluorine, difluoro
MASS   402 FGA3    18.99800  ! aliphatic fluorine, trifluoro
MASS   403 FGP1    18.99800  ! anionic F, for ALF4 AlF4-
MASS   404 FGR1    18.99800  ! aromatic flourine
!miscellaneous
MASS   405 PG0     30.97380  ! neutral phosphate
MASS   406 PG1     30.97380  ! phosphate -1
MASS   407 PG2     30.97380  ! phosphate -2
MASS   408 ALG1    26.98154  ! Aluminum, for ALF4, AlF4-

MASS   409 CG25C1  12.01100  ! same as CG2DC1 but in 5-membered ring with exocyclic double bond
MASS   410 CG25C2  12.01100  ! same as CG2DC2 but in 5-membered ring with exocyclic double bond
MASS   411 CG251O  12.01100  ! same as CG2D1O but in 5-membered ring with exocyclic double bond
MASS   412 CG252O  12.01100  ! same as CG2D2O but in 5-membered ring with exocyclic double bond

!MASS   413 HGTIP3   1.00800  ! polar H, TIPS3P WATER HYDROGEN
!MASS   414 OGTIP3  15.99940  ! TIPS3P WATER OXYGEN
!MASS   415 DUM      0.00000  ! dummy atom
!MASS   416 HE       4.00260  ! helium
!MASS   417 NE      20.17970  ! neon
MASS   429 DUMM      0.001    ! Dummy with mass
EOF
cat >> hybrid.rtf <<EOF

AUTOGENERATE ANGLES DIHE
DEFA FIRS NONE LAST NONE

RESI $LIGAND_NAME  $LIGAND_CHARGE ! LIG-ATP linked by a soft bond between <backbone C of an unnatural AA> and <alpha P of an ATP>
GROUP
EOF

#################################################

grep "^ATOM" drawing_3D/drawing_3D.rtf | awk '{$2 = ($2 != "C8" && 
                                                     $2 != "C9" && 
                                                     $2 != "H8" \
                                                     ? $2"X" : $2"A")} 1' | column -t >> hybrid.rtf

cat >> hybrid.rtf <<'EOF'
ATOM  H8B   HCMM  -0.0730
EOF

cat >> hybrid.rtf <<'EOF'
GROUP
ATOM C4'  CG3C51    0.11  !                       H61  H62
ATOM H4'  HGA1      0.09  !                         \  /
ATOM O4'  OG3C51   -0.40  !                          N6
ATOM C1'  CG3C51    0.11  !                          |
ATOM H1'  HGA1      0.09  !                          C6
GROUP                     !                        //  \
ATOM C5   CG2RC0    0.28  !                        N1   C5--N7\\
ATOM N7   NG2R50   -0.71  !                        |    ||     C8-H8
ATOM C8   CG2R53    0.34  !                        C2   C4--N9/
ATOM H8   HGR52     0.12  !                       / \\ /     \
ATOM N9   NG2R51   -0.05  !                     H2   N3       \
                          !                                    \
ATOM N1   NG2R62   -0.74  !                                     \
ATOM C2   CG2R64    0.50  !                                      \
ATOM H2   HGR62     0.13  ! (-)O3G   O2B     O1A    H5' H4'  O4'  \
ATOM N3   NG2R62   -0.75  !    |      |       |      |    \ /   \  \
ATOM C4   CG2RC0    0.43  !O1G=PG-O3B-PB-O3A-PA-O5'-C5'---C4'    C1'
ATOM C6   CG2R64    0.46  !    |      |       |      |     \     / \
                          ! (-)O2G (-)O1B (-)O2A    H5''   C3'--C2' H1'
ATOM N6   NG2S3    -0.77  !                               / \   / \
ATOM H61  HGP4      0.38  !                            O3' H3' O2' H2''
ATOM H62  HGP4      0.38  !                             |       |
GROUP                     !                            H3T     H2'
ATOM C2'  CG3C51    0.14
ATOM H2'' HGA1      0.09
ATOM O2'  OG311    -0.65
ATOM H2'  HGP1      0.42
GROUP
ATOM C3'  CG3C51    0.14
ATOM H3'  HGA1      0.09
ATOM O3'  OG311    -0.65
ATOM H3T  HGP1      0.42
GROUP
ATOM C5'  CG321    -0.08
ATOM H5'  HGA2      0.09
ATOM H5'' HGA2      0.09
ATOM O5'  OG303    -0.62
ATOM PA   PG1       1.50
ATOM O1A  OG2P1    -0.82
ATOM O2A  OG2P1    -0.82
ATOM O3A  OG304    -0.74
ATOM PB   PG1       1.50
ATOM O1B  OG2P1    -0.82
ATOM O2B  OG2P1    -0.82
ATOM O3B  OG304    -0.86 ! charge adjusted to yield total triP of -4.0
ATOM PG   PG2       1.10
ATOM O1G  OG2P1    -0.90
ATOM O2G  OG2P1    -0.90
ATOM O3G  OG2P1    -0.90
EOF

list=`grep "^ATOM" drawing_3D/drawing_3D.rtf | awk '$2 !~ /^C8$|^C9$|^H8$/ {print $2}'`
#################################################

grep "^BOND" drawing_3D/drawing_3D.rtf > tmp
for ii in $LIST; do
    sed -i '' 's/ '$ii' / '$ii'A/g' tmp
done
for ii in $list; do
    sed -i '' 's/ '$ii' / '$ii'X/g' tmp
done
cat >> tmp <<'EOF'
BOND  C4X  H8B
EOF

cat >> tmp <<'EOF'
BOND O5'  C5'       O5'  PA        PA   O1A       PA   O2A       PA   O3A
BOND O3A  PB        PB   O1B       PB   O2B       PB   O3B       O3B  PG
BOND PG   O1G       PG   O2G       PG   O3G
BOND C5'  C4'       C4'  O4'       C4'  C3'       O4'  C1'
BOND C1'  N9        C1'  C2'       N9   C4        N9   C8        C4   N3
BOND C2   N1        C6   N6
BOND N6   H61       N6   H62       C6   C5        C5   N7
BOND C2'  C3'       C2'  O2'       O2'  H2'       C3'  O3'       O3'  H3T
BOND C1'  H1'       C2'  H2''      C3'  H3'       C4'  H4'       C5'  H5'
BOND C5'  H5''      C8   H8        C2   H2
EOF

cat >> tmp <<'EOF'
BOND  CX   PA
EOF

cat >> tmp <<'EOF'
DOUBLE N1   C6    N3   C2    C4   C5        N7   C8
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IMPH" drawing_3D/drawing_3D.rtf > tmp
for ii in $LIST; do
    sed -i '' 's/ '$ii' / '$ii'A/g' tmp
done
for ii in $list; do
    sed -i '' 's/ '$ii' / '$ii'X/g' tmp
done
cat >> tmp <<'EOF'
IMPH  C4X  C3X  C5X  H8B   
EOF

cat >> tmp <<'EOF'
IMPR C6 C5  N1  N6
IMPR N6 H62 H61 C6
DONO H61  N6
DONO H62  N6
DONO H2'  O2'
ACCE N3
ACCE N7
ACCE N1
ACCE O1A  PA
ACCE O2A  PA
ACCE O2'
ACCE O3'
ACCE O4'
ACCE O5'
ACCE O3A
ACCE O2B
ACCE O1B
ACCE O3B
ACCE O3G
ACCE O1G
ACCE O2G
EOF

column -t tmp >> hybrid.rtf 
#################################################

grep "^IC" drawing_3D/drawing_3D.rtf > tmp
for ii in $LIST; do
    sed -i '' 's/ '$ii' / '$ii'A/g' tmp
done
for ii in $list; do
    sed -i '' 's/ '$ii' / '$ii'X/g' tmp
done
cat >> tmp <<'EOF'
IC  C2X  C3X   C4X  H8B   1.39   119.97  169.96  120.00  1.02  
IC  C6X  C5X   C4X  H8B   1.39   120.00 -179.96  119.97  1.02  
IC  O1X  C5X   C4X  H8B   1.35   120.04    0.01  119.97  1.02  
IC  H8B  C4X   C3X  H5X   1.02   120.00   -0.09  120.01  1.02  
IC  C3X  C5X  *C4X  H8B   0.00     0.00  180.00    0.00  0.00  
EOF

cat >> tmp <<'EOF'
IC PA   O3A  PB   O3B    0.0000  000.00  180.0   000.00   0.0000
IC PB   O3A  PA   O5'    0.0000  000.00  120.0   000.00   0.0000
IC O5'  O3A  *PA  O1A    0.0000  000.00  120.0   000.00   0.0000
IC O5'  O3A  *PA  O2A    0.0000  000.00 -120.0   000.00   0.0000
IC O3A  O3B  *PB  O1B    0.0000  000.00  120.0   000.00   0.0000
IC O3A  O3B  *PB  O2B    0.0000  000.00 -120.0   000.00   0.0000
IC O3A  PB   O3B  PG     0.0000  000.00  150.0   000.00   0.0000
IC PB   O3B  PG   O1G    0.0000  000.00  150.0   000.00   0.0000
IC O3B  O1G  *PG  O2G    0.0000  000.00  120.0   000.00   0.0000
IC O3B  O1G  *PG  O3G    0.0000  000.00 -120.0   000.00   0.0000
IC O3A  PA   O5'  C5'    0.0000  000.00  -30.0   000.00   0.0000
!IC PA   O5'  C5'  C4'    1.5996  119.00 -151.39  110.04   1.5160 !Not-so-stable minimum
IC PA   O5'  C5'  C4'    0.0000  000.00 -110.0   000.00   0.0000
!IC O5'  C5'  C4'  C3'    1.4401  108.83 -179.85  116.10   1.5284 !Not-so-stable minimum
IC O5'  C5'  C4'  C3'    0.0000  000.00   70.0   000.00   0.0000
IC C5'  C4'  C3'  O3'    1.5160  116.10   76.70  115.12   1.4212
IC H3T  O3'  C3'  C4'    0.9650  105.47   38.18  111.98   1.5386
IC O4'  C3'  *C4' C5'    1.4572  104.06 -120.04  116.10   1.5160
IC C2'  C4'  *C3' O3'    1.5284  100.16 -124.08  115.12   1.4212
IC C4'  C3'  C2'  C1'    1.5284  100.16   39.58  102.04   1.5251
!IC C3'  C2'  C1'  N9     1.5284  101.97  144.39  113.71   1.4896 !144.39 can't be right
IC O4'  C2'  *C1' N9     0.0       0.0   120.0     0.0    0.0
!IC O4'  C1'  N9   C4     1.5251  113.71  -96.00  125.97   1.3703 !clash after fixing N9
IC O4'  C1'  N9   C4     0.0       0.0  -150.0     0.0    0.0
IC C1'  C4   *N9  C8     1.4896  125.97 -179.94  106.0    1.367
IC C4   N9   C8   N7     1.376   106.0     0.0   113.6    1.312
IC C8   N9   C4   C5     1.367   106.0     0.0   105.6    1.382
IC C8   N7   C5   C6     0.0       0.0   180.0     0.0    0.0
IC N7   C5   C6   N1     0.0       0.0   180.0     0.0    0.0
IC C5   C6   N1   C2     0.0       0.0     0.0     0.0    0.0
IC N9   C5   *C4  N3     1.376   105.6  -180.0   126.9    1.342
IC C5   N1   *C6  N6     1.409   117.6  -180.0   121.2    1.337
IC N1   C6   N6   H61    1.337   121.2     0.0   119.0    1.01
IC H61  C6   *N6  H62    1.01    119.0   180.0   119.00   1.01
IC C5   N1   *C6  N6     1.409   117.6  -180.0   119.0    1.337
IC N1   C6   N6   H61    1.337   119.0     0.0   119.0    1.01
IC H61  C6   *N6  H62    1.01    119.0   180.0   121.00   1.01
IC N9   N7   *C8   H8    0.0       0.0   180.0     0.0    0.0
IC N1   N3   *C2   H2    0.0       0.0   180.0     0.0    0.0
IC C1'  C3'  *C2' O2'    1.5284  102.04 -114.67  110.81   1.4212
IC H2'  O2'  C2'  C3'    0.9600  114.97  148.63  111.92   1.5284
IC O4'  C2'  *C1' H1'    0.0       0.0  -115.0     0.0    0.0
IC C1'  C3'  *C2' H2''   0.0       0.0   115.0     0.0    0.0
IC C2'  C4'  *C3' H3'    0.0       0.0   115.0     0.0    0.0
IC C3'  O4'  *C4' H4'    0.0       0.0  -115.0     0.0    0.0
IC C4'  O5'  *C5' H5'    0.0       0.0  -115.0     0.0    0.0
IC C4'  O5'  *C5' H5''   0.0       0.0   115.0     0.0    0.0
EOF

column -t tmp >> hybrid.rtf 
#################################################

cat > tcl <<'EOF'
mol new drawing_3D/drawing_3D.pdb
set sel [atomselect top "name C CA O N OXT"]
$sel set resname LIG
$sel set resid 1
$sel writepdb 1.pdb
mol new anp.pdb
set sel [atomselect top all]
$sel set resname LIG
$sel set resid 1
$sel writepdb 2.pdb

package require psfgen
resetpsf
topology readcharmmtop1.2/top_all36_prot.rtf
topology readcharmmtop1.2/top_all36_hybrid.inp
topology toppar_water_ions_namd.str
topology hybrid.rtf
segment LIG {
    pdb 1.pdb
    foreach name {C CA O N OXT H HA} {
        pdbalias atom LIG ${name} ${name}X
    }
}
coordpdb 1.pdb LIG
coordpdb 2.pdb LIG
regenerate angles dihedrals
guesscoord
writepsf prot.psf
writepdb prot.pdb
quit
EOF

vmd -dispdev text -e tcl
rm tcl tmp