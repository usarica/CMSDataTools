#!/bin/bash


INPUTDIR=$1
FILENAME=$2
OUTPUTSITE=$3 # e.g. 't2.ucsd.edu'
OUTPUTDIR=$4 # Must be absolute path
RENAMEFILE=$FILENAME
if [[ "$5" != "" ]];then
  RENAMEFILE=$5
fi


echo "Copy from Condor is called with"
echo "INPUTDIR: ${INPUTDIR}"
echo "FILENAME: ${FILENAME}"
echo "OUTPUTSITE: ${OUTPUTSITE}"
echo "OUTPUTDIR: ${OUTPUTDIR}"
echo "RENAMEFILE: ${RENAMEFILE}"


if [[ "$INPUTDIR" == "" ]];then #Input directory is empty, so assign pwd
  INPUTDIR=$(pwd)
elif [[ "$INPUTDIR" != "/"* ]];then # Input directory is a relative path
  INPUTDIR=$(pwd)/${INPUTDIR}
fi

if [[ "$OUTPUTDIR" != "/"* ]];then # Output directory must be an absolute path!
  echo "Output directory must be an absolute path! Cannot transfer the file..."
  exit 1
fi


if [[ ! -z ${FILENAME} ]];then
  echo -e "\n--- begin copying output ---\n"

  echo "Sending output file ${FILENAME}"

  if [[ ! -e ${INPUTDIR}/${FILENAME} ]]; then
    echo "ERROR! Output ${FILENAME} doesn't exist"
    exit 1
  fi

  echo "Time before copy: $(date +%s)"

  COPY_SRC="file://${INPUTDIR}/${FILENAME}"
  COPY_DEST=""
  echo "OUTPUTFILE: ${OUTPUTDIR}/${RENAMEFILE}"
  if [[ "${OUTPUTSITE}" == *"t2.ucsd.edu"* ]]; then
    if [[ ! -s ${INPUTDIR}/${FILENAME} ]]; then
      echo "Running the gsiftp endpoint..."
      COPY_DEST="gsiftp://gftp.${OUTPUTSITE}${OUTPUTDIR}/${RENAMEFILE}"
    else
      echo "Running the xrootd endpoint..."
      COPY_DEST="davs://redirector.t2.ucsd.edu:1094${OUTPUTDIR}/${RENAMEFILE}"
      COPY_DEST=${COPY_DEST/'/hadoop/cms'/''}
    fi
  elif [[ "${OUTPUTSITE}" == *"eoscms.cern.ch"* ]]; then
    COPY_DEST="root://eoscms.cern.ch${OUTPUTDIR}/${RENAMEFILE}"
    COPY_DEST=${COPY_DEST/'/eos/cms'/''}
  elif [[ "${OUTPUTSITE}" == *"iihe.ac.be"* ]]; then
    COPY_DEST="srm://maite.iihe.ac.be:8443${OUTPUTDIR}/${RENAMEFILE}"
    #COPY_DEST=${COPY_DEST/'/pnfs/iihe/cms'/''}
  elif [[ "${OUTPUTSITE}" == *"ihep.ac.cn"* ]]; then
    COPY_DEST="srm://srm.ihep.ac.cn:8443${OUTPUTDIR}/${RENAMEFILE}"
    #COPY_DEST=${COPY_DEST/'/data/cms'/''}
  elif [[ "${OUTPUTSITE}" == *"m45.ihep.su"* ]]; then
    COPY_DEST="srm://dp0015.m45.ihep.su:8443${OUTPUTDIR}/${RENAMEFILE}"
    #COPY_DEST=${COPY_DEST/'/pnfs/m45.ihep.su/data/cms'/''}
  fi

  runCmd="env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 14400 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}"
  echo "Running: ${runCmd}"
  declare -i itry=0
  declare -i COPY_STATUS=-1
  declare -i REMOVE_STATUS=-1
  while [[ $itry -lt 5 ]]; do
    echo " - Trial ${itry}:"
    ${runCmd}
    COPY_STATUS=$?
    if [[ $COPY_STATUS -eq 0 ]]; then
      break
    fi
    (( itry += 1 ))
  done
  if [[ $COPY_STATUS -ne 0 ]]; then
    echo "Removing output file because gfal-copy crashed with code $COPY_STATUS"
    env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-rm -t 14400 --verbose ${COPY_DEST}
    REMOVE_STATUS=$?
    if [[ $REMOVE_STATUS -ne 0 ]]; then
        echo "gfal-copy crashed and then the gfal-rm also crashed with code $REMOVE_STATUS"
        echo "You probably have a corrupt file sitting on ${OUTPUTDIR} now."
        exit $REMOVE_STATUS
    fi
    exit $COPY_STATUS
  else
    echo "Time after copy: $(date +%s)"
    echo "Copied successfully!"
  fi

  echo -e "\n--- end copying output ---\n"
else
  echo "File name is not specified!"
  exit 1
fi
