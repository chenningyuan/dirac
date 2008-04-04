#!/bin/bash

for basis in 4s 5s 6s 7s 4p 5p 6p 7p 3d 4d 5d 6d 4f 5f; do
   file=$basis"_archive.tar.gz"
   wget http://dirac.chem.sdu.dk/basisarchives/dyall/$file
   tar zxvf $file
   rm $file
done
