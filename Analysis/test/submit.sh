#!/bin/sh

FILE=$1
FCN=$2
FCNARGS=$3

mkdir -p ./output/Logs

extLog=$FCN
if [[ "$FCNARGS" != "" ]];then
  fcnargname=${FCNARGS//\"}
  fcnargname=${fcnargname//"("}
  fcnargname=${fcnargname//")"}
  fcnargname=${fcnargname//","/"_"}
  extLog=$FCN"_"$fcnargname
fi

if [[ -f $FCN".c" ]]; then
  echo "File "$FCN".c"" already exists."
else
  cp $FILE $FCN".c"
fi

root -l -b -q -e "gROOT->ProcessLine(\".x loadLib.C\"); gROOT->ProcessLine(\".L "$FCN".c+\");"

bsub -q 2nd -C 0 -o "./output/Logs/lsflog_"$extLog".txt" -e "./output/Logs/lsferr_"$extLog".err" submit.lsf.sh $FCN $FCNARGS

