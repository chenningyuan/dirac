#!/bin/sh
pam -rmwrk -utcmo  H H
mv DFCOEF cf.H
ln -s cf.H DFHAT1
ln -s cf.H DFHAT2
ln -s cf.H DFHAT3
ln -s cf.H DFHAT4
pam -utcmo -rmwrk C C
mv DFCOEF cf.C
ln -s cf.C DFCAT
pam -rmwrk -copy "DFCAT DFHAT1 DFHAT2 DFHAT3 DFHAT4" CH4 CH4
exit 0
