#!/bin/sh

fixedDate=""
if [[ "$1" != "" ]];then
  fixedDate=$1
fi

script="makeQQBkgTemplatesFromPOWHEG"

. submit.sh $script".cc" $script"_one" "(k4e,Inclusive,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k4mu,Inclusive,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k2e2mu,Inclusive,\"Nominal\",\""$fixedDate"\")"
