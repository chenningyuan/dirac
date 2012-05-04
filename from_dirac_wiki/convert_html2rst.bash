#!/bin/bash

while read html_line
do a=$(($a+1));
#echo ${html_line##*/} | sed 's/^[\. _]//'
file_name=`echo ${html_line##*/} | sed 's/^[\. _]//'`
echo $a  $html_line   $file_name
pandoc --from=html  $html_line -o $file_name.rst
done < "all_wiki_pages"

echo "Final line count is: $a";
