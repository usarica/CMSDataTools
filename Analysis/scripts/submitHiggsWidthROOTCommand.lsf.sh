#!/bin/sh

RUNDIR=${LS_SUBCWD}
cd $RUNDIR
echo "LSF job running in: " `pwd`

CMSENVDIR=$1
cd $CMSENVDIR
eval `scram runtime -sh`
cd $RUNDIR

echo $CMSSW_VERSION

echo "Host name: "$(hostname)

runHiggsWidthROOTCommand.py --loadlib loadLib.C --script $2 --function $3 --command $4

