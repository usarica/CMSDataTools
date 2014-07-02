#!/bin/sh

if [ $# == 0 ]
	then
	echo "Usage: ./pushjson.sh <Directory of LHC_*TeV templates> <TemplateBuilder/run location>"
	exit
elif [ $# == 1 ] 
    then 
    direct=$1
    tbloc=../TemplateBuilder/run
elif [ $# == 2 ] 
    then
    direct=$1
    tbloc=$2
fi

ls HiggsWidth/*/*tpl > jsonlist

while read line; do
	echo $line
	filebase=`echo "${line:0:${#line}-4}"`
	sed -e 's|<DIR>|'$direct'|g' < $line > $filebase.json
done < jsonlist

cp -rp HiggsWidth $tbloc/.
rm $tbloc/HiggsWidth/*/*tpl

rm jsonlist