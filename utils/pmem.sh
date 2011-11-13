#!/bin/bash
# Small script to get the static memory of an executable. Thanks to Simen Reine and Pawel Salek !
objdump -hw $1| perl -p -e '@a=split;if($a[1] eq ".bss") {$s="0x".$a[2];} END{print hex($s)/(1024*1024), "MB\n"}'
exit 0
