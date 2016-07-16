#!/bin/bash/

dir="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/CMSSW_8_0_12/src/Analysis/ICHEP2016/"
todaysdate=$1
energy=$2
channel=$3
xbinning=$4
syst=$5
djet=$6
#if [[ "$djet" != "" ]];then
#  djet="_"$djet
#endif

tplbldr="/work-zfs/lhc/usarica/CMS-related/Analysis/MassWidth_2016/CMSSW_8_0_12/src/Analysis/TemplateBuilder/buildTemplate.exe"
tplfile="templates_ggTo4l_V1.tpl"
jsonfile="templates_ggTo4l_"$energy"TeV_"$channel"_"$syst$djet".json"

cp $tplfile $jsonfile

SEDCOMMAND="s.<DIR>."$dir".g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<TODAYSDATE>/"$todaysdate"/g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<ENERGY>/"$energy"/g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<CHANNEL>/"$channel"/g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<XBINNING>/"$xbinning"/g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<SYST>/"$syst"/g"
sed -i $SEDCOMMAND $jsonfile
SEDCOMMAND="s/<DJET>/"$djet"/g"
sed -i $SEDCOMMAND $jsonfile

$tplbldr $jsonfile


