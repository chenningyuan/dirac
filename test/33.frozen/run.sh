#
pam -rmwrk -utcmo -run std O O
mv DFCOEF cf.O
ln -f -s cf.O DFPROJ
pam -rmwrk -copy DFPROJ -nice -run f1s2s_geo H2O H2O
rm DFPROJ cf.O
exit 0
