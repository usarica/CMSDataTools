#!/bin/sh

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

eval `scram runtime -sh`

echo $CMSSW_VERSION

infile=$1
runfile=$2
extcmd="()"
if [[ "$3" != "" ]];then
  extcmd=$3
fi

#rm -rf $runfile".c"
if [[ -f $runfile".c" ]]; then
  echo "File "$runfile".c"" already exists."
else
  cp $infile $runfile".c"
fi
cmd=$runfile".c+"$extcmd

root -b -l -q loadLib.C $cmd


