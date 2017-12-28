#!/bin/sh

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

runfile=$1
extcmd="()"
if [[ "$2" != "" ]];then
  extcmd=$2
fi

cmd=$runfile".c+"$extcmd

root -b -l -q loadLib.C $cmd


