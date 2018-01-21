#!/bin/sh

version=$1
stage=$2
tpldir="output/LHC_13TeV/Templates/"$version"/"

targetdir=$CMSSW_BASE"/src/HZZ4l/CreateWidthDatacards/templates2D/"$version"/"
mkdir -p $targetdir"SM/13TeV/../../a2/13TeV/../../a3/13TeV/../../L1/13TeV/"

for f in $(ls $tpldir | grep -e "_Stage"$stage".root" | grep -e "Check"); do
	newname=$f
	GEN=""
	for g in POWHEG MCFM;do
		if [[ "$newname" == *"$g"* ]];then
			GEN=$g
			newname=${newname/"_"$g/""}
		fi
	done
	if [[ "$GEN" == "" ]];then
		continue
	fi
	ACdir=""
	for AC in SM a2 a3 L1; do
		if [[ "$newname" == *"Check"$AC"Discriminants"* ]];then
			ACdir=$AC
			newname=${newname/"_Check"$AC"Discriminants_Stage"$stage/""}
		fi
	done
	if [[ "$ACdir" == "" ]];then
		continue
	fi

	let doCP=0
	if [[ "$GEN" == "MCFM" ]];then
	if [[ "$newname" == *"ggZZ"* ]]; then
	if [[ "$newname" == *"Inclusive"* ]]; then
		let doCP=1
	fi;fi;fi
	if [[ "$GEN" != "MCFM" ]];then
	if [[ "$newname" == *"ggZZ"* ]]; then
	if [[ "$newname" != *"Inclusive"* ]]; then
		let doCP=1
	fi;fi;fi
	if [[ "$GEN" != "MCFM" ]];then
	if [[ "$newname" != *"ggZZ"* ]]; then
		let doCP=1
	fi;fi
	if [ $doCP -eq 1 ];then
		newname=${newname/"Dn.root"/"Down.root"}
		echo $newname
		ln -sf $(pwd)"/"$tpldir$f $targetdir$ACdir"/13TeV/"$newname
	fi
	
done
