#!/bin/sh

cd ${LS_SUBCWD}

echo "LSF job running in: " `pwd`

CMSDIR=$(readlink cmsdir)
if [[ "$CMSDIR" != "" ]];
  pushd $CMSDIR
fi
eval `scram runtime -sh`
if [[ "$CMSDIR" != "" ]];
  popd
fi

echo $CMSSW_VERSION

echo "Host name: "$(hostname)

runHiggsWidthROOTCommand.py --loadlib loadLib.C --script $1 --function $2 --command $3

