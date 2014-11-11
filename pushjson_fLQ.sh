#!/bin/sh

if [ $# == 0 ]
	then
	echo "Usage: ./pushjson_fLQ.sh <Location of LHC_*TeV directories that contain templates with trees> <Use Djet(0/1)> <TemplateBuilder/run location>"
	exit
elif [ $# == 1 ] 
    then 
    direct=$1
    Djet=1
    tbloc=../TemplateBuilder/run
elif [ $# == 2 ] 
    then
    direct=$1
    Djet=$2
    tbloc=../TemplateBuilder/run
elif [ $# == 3 ]
	then
	direct=$1
	Djet=$2
	tbloc=$3
fi

if [ [ $Djet != 0 ] && [ $Djet != 1 ] ]
	then
	echo "Djet option only accepts 0 or 1"
	exit
fi

ls HiggsWidth/*/*_fLQ.tpl > jsonlist

systematics=( Nominal SysUp_ggQCD SysDown_ggQCD SysUp_ggPDF SysDown_ggPDF )

while read line; do
	echo $line
	filebase=`echo "${line:0:${#line}-4}"`
	for i in `seq 0 4`;
	do
		filebase_new=`echo "${filebase}_${systematics[i]}"`
		echo $filebase_new
		sed -e 's|<DIR>|'$direct'|g' < $line > $filebase_new.json
		sed -i 's|<SYST>|'${systematics[i]}'|g' $filebase_new.json
		sed -i 's|<JET>||g' $filebase_new.json
		filebase_new_djet=`echo "${filebase}_${systematics[i]}_Djet"`
		echo $filebase_new_djet
		sed -e 's|<DIR>|'$direct'|g' < $line > $filebase_new_djet.json
		sed -i 's|<SYST>|'${systematics[i]}'|g' $filebase_new_djet.json
		sed -i 's|<JET>|_Djet|g' $filebase_new_djet.json
		filebase_new_nondjet=`echo "${filebase}_${systematics[i]}_nonDjet"`
		echo $filebase_new_nondjet
		sed -e 's|<DIR>|'$direct'|g' < $line > $filebase_new_nondjet.json
		sed -i 's|<SYST>|'${systematics[i]}'|g' $filebase_new_nondjet.json
		sed -i 's|<JET>|_nonDjet|g' $filebase_new_nondjet.json
	done
done < jsonlist

cp -rp HiggsWidth $tbloc/.
rm $tbloc/HiggsWidth/*/*tpl

rm jsonlist