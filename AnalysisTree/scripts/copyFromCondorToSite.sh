#!/bin/bash

INPUTDIR=$1
if [[ "$INPUTDIR" == "" ]];then #Input directory is empty, so assign pwd
  INPUTDIR=$(pwd)
elif [[ "$INPUTDIR" != "/"* ]];then # Input directory is a relative path
  INPUTDIR=$(pwd)/${INPUTDIR}
fi

FILENAME=$2

OUTPUTSITE=$3 # e.g. 't2.ucsd.edu'

OUTPUTDIR=$4 # Must be absolute path
if [[ "$OUTPUTDIR" != "/"* ]];then # Output directory must be an absolute path!
  echo "output directory must be an absolute path! Cannot transfer the files..."
  exit 1
fi

RENAMEFILE=$FILENAME
if [[ "$5" != "" ]];then
  RENAMEFILE=$5
fi


if [ ! -z ${FILENAME} ];then
  echo -e "\n--- begin copying output ---\n"

  echo "Sending output file ${FILENAME}"

  if [ ! -e "$FILENAME" ]; then
    echo "ERROR! Output $FILENAME doesn't exist"
    exit 1
  fi

  echo "time before copy: $(date +%s)"

  COPY_SRC="file://${INPUTDIR}/${FILENAME}"
  COPY_DEST="gsiftp://gftp.${OUTPUTSITE}${OUTPUTDIR}/${RENAMEFILE}"
  echo "Running: env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}"
  env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}
  COPY_STATUS=$?
  if [[ $COPY_STATUS != 0 ]]; then
    echo "Removing output file because gfal-copy crashed with code $COPY_STATUS"
    env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-rm -t 7200 --verbose ${COPY_DEST}
    REMOVE_STATUS=$?
    if [[ $REMOVE_STATUS != 0 ]]; then
        echo "gfal-copy crashed and then the gfal-rm also crashed with code $REMOVE_STATUS"
        echo "You probably have a corrupt file sitting on ${OUTPUTDIR} now."
        exit $REMOVE_STATUS
    fi
    exit $COPY_STATUS
  else
    echo "Copied successfully!"
  fi

  echo -e "\n--- end copying output ---\n"
fi
