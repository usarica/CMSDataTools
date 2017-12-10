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

hname=$(hostname)
echo $hname
if [[ "$hname" == *"lxplus"* ]];then
  echo "Host is on LXPLUS, so need to use LXBATCH"
  bsub -q 2nd -C 0 -o "./output/Logs/lsflog_"$extLog".txt" -e "./output/Logs/lsferr_"$extLog".err" submit.lsf.sh $FCN $FCNARGS
elif [[ "$hname" == *"login-node"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch"
  sbatch --output="./output/Logs/lsflog_"$extLog".txt" --error="./output/Logs/lsferr_"$extLog".err" submit.slurm.sh $FCN $FCNARGS
fi


