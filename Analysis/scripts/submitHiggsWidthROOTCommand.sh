#!/bin/sh

SCRIPT=$1
FCN=$2
FCNARGS=$3
QUEUE=$4

CMSENVDIR=$CMSSW_BASE
if [[ "$CMSENVDIR" == "" ]];then
  echo "Set up CMSSW first!"
  exit 1
fi

if [[ "$QUEUE" == "" ]];then
  QUEUE="default"
fi

echo "Calling "$SCRIPT"::"$FCN"("$FCNARGS")"


mkdir -p ./output/Logs

extLog=$SCRIPT"_"$FCN
extLog="${extLog/./_}"
if [[ "$FCNARGS" != "" ]];then
  fcnargname=${FCNARGS//\"}
  fcnargname=${fcnargname//"("}
  fcnargname=${fcnargname//")"}
  fcnargname=${fcnargname//"\\"}
  fcnargname=${fcnargname//","/"_"}
  fcnargname="${fcnargname/./_}"
  extLog=$extLog"_"$fcnargname
fi
echo "Log files will be appended "$extLog

if [[ -f $SCRIPT ]]; then
  echo "File "$SCRIPT" exists."
  SOFILE=$SCRIPT
  SOFILE="${SOFILE/./_}"
  SOFILE=$SOFILE".so"

  if [[ ! -f $SOFILE ]]; then
    echo "Compiling "$SCRIPT"..."
    root -l -b -q -e "gROOT->ProcessLine(\".x loadLib.C\"); gROOT->ProcessLine(\".L "$SCRIPT"+\");"
  fi

  hname=$(hostname)
  echo $hname
  if [[ "$hname" == *"lxplus"* ]];then
    echo "Host is on LXPLUS, so need to use LXBATCH"
    THEQUEUE="2nd"
    if [[ "$QUEUE" != "default" ]];then
      THEQUEUE=$QUEUE
    fi
    bsub -q $THEQUEUE -C 0 -o "./output/Logs/lsflog_"$extLog".txt" -e "./output/Logs/lsferr_"$extLog".err" submitHiggsWidthROOTCommand.lsf.sh $CMSENVDIR $SCRIPT $FCN $FCNARGS
  elif [[ "$hname" == *"login-node"* ]]; then
    echo "Host is on MARCC, so need to use SLURM batch"
    THEQUEUE="lrgmem"
    if [[ "$QUEUE" != "default" ]];then
      THEQUEUE=$QUEUE
    fi
    sbatch --output="./output/Logs/lsflog_"$extLog".txt" --error="./output/Logs/lsferr_"$extLog".err" --partition=$THEQUEUE submitHiggsWidthROOTCommand.slurm.sh $CMSENVDIR $SCRIPT $FCN $FCNARGS
  fi

fi

