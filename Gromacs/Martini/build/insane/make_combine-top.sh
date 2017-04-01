martinize_top=../../dimer/martinize-elnedyn/dimer.top
insane_top=insane.top

itp_martini_name=martini_v2.2.itp
itp_lipid_name=martini_v2.0_lipids_all_201506.itp
itp_ion_name=martini_v2.0_ions.itp
itp_selfdef_name="martini_v2.0_PAP2.itp"

echo -e "#include \"$itp_martini_name\"" > TMP_itp
sed '/\[/,$d' $martinize_top | sed '1d' >> TMP_itp
echo -e "#include \"$itp_lipid_name\"" >> TMP_itp
echo -e "#include \"$itp_ion_name\"" >> TMP_itp
for ii in $itp_selfdef_name; do
    echo -e "#include \"$ii\"" >> TMP_itp
done
echo -e "\n\n\n" >> TMP_itp

sed -n '/\]/,/\[/p' $martinize_top | sed '$d' > TMP_system
sed -n '/\]/,/\[/p' $insane_top | sed '1,2d' | sed '$d' >> TMP_system

sed -n '/\[ molecules \]/,$p' $martinize_top > TMP_molecules
echo "" >> TMP_molecules
sed -n '/\[ molecules \]/,$p' $insane_top | sed '1,3d' >> TMP_molecules

cat TMP_itp TMP_system TMP_molecules > system.top
rm -f TMP_itp TMP_system TMP_molecules
