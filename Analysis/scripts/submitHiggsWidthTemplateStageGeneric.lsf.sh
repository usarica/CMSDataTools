#!/bin/sh


RUNDIR=${LS_SUBCWD}
cd $RUNDIR
echo "LSF job running in: " `pwd`

CMSENVDIR=$1
cd $CMSENVDIR
eval `scram runtime -sh`
cd $RUNDIR

echo $CMSSW_VERSION

runfile=$2
extcmd="()"
if [[ "$3" != "" ]];then
  extcmd=$3
fi

cmd=$runfile".c+"$extcmd

root -b -l -q loadLib.C $cmd


