#!/bin/sh

fixedDate=""
if [[ "$1" != "" ]];then
  fixedDate=$1
fi

script="makeQQBBKGTemplatesFromPOWHEG"

. submit.sh $script".cc" $script"_checkstage" "(k4e,Inclusive,kSM,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k4mu,Inclusive,kSM,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k2e2mu,Inclusive,kSM,\"Nominal\",2,\""$fixedDate"\")"

. submit.sh $script".cc" $script"_checkstage" "(k4e,Inclusive,kL1,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k4mu,Inclusive,kL1,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k2e2mu,Inclusive,kL1,\"Nominal\",2,\""$fixedDate"\")"

. submit.sh $script".cc" $script"_checkstage" "(k4e,Inclusive,kA2,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k4mu,Inclusive,kA2,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k2e2mu,Inclusive,kA2,\"Nominal\",2,\""$fixedDate"\")"

. submit.sh $script".cc" $script"_checkstage" "(k4e,Inclusive,kA3,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k4mu,Inclusive,kA3,\"Nominal\",2,\""$fixedDate"\")"
. submit.sh $script".cc" $script"_checkstage" "(k2e2mu,Inclusive,kA3,\"Nominal\",2,\""$fixedDate"\")"
