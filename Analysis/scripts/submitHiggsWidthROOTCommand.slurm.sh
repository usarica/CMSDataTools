#!/bin/sh

#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=lrgmem
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=usarica1@jhu.edu

cd ${SLURM_SUBMIT_DIR}
echo "SLURM job running in: " `pwd`

# ROOT
# xpm (needed for ROOT)
export CPATH="/work-zfs/lhc/usarica/libXpm-3.5.11/include:$CPATH"
source /work-zfs/lhc/usarica/ROOT/bin/thisroot.sh

# CMSSW
source /work-zfs/lhc/cms9/cmsset_default.sh
module load boost/1.60.0
export LIBRARY_PATH=$LIBRARY_PATH:/cm/shared/apps/boost/1.60.0/lib
export CPATH=$CPATH:/cm/shared/apps/boost/1.60.0/include

eval `scram runtime -sh`

echo $CMSSW_VERSION

echo "Host name: "$(hostname)

runHiggsWidthROOTCommand.py --loadlib loadLib.C --script $1 --function $2 --command $3

