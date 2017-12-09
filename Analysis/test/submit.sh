#!/bin/sh

FILE=$1
FCN=$2
FCNARGS=$3

mkdir -p ./output/Logs

extLog=""
if [[ "$FCNARGS" != "" ]];then
  fcnargname=${FCNARGS//\"}
  fcnargname=${FCNARGS//(}
  fcnargname=${FCNARGS//)}
  fcnargname=${FCNARGS//,/_}
  extLog=$FCN"_"$fcnargname
else
  extLog=$FCN
fi

bsub -q 2nd -C 0 -o ./output/Logs/lsflog_"$extLog".txt" -e ./output/Logs/lsferr_"$extLog".err submit.lsf.sh $FCN $FCNARGS

