#!/bin/sh

fixedDate=""
if [[ "$1" != "" ]];then
  fixedDate=$1
fi

script="makeGGTemplatesFromPOWHEG"

rm -f $script".cc"

. submit.sh $script".cc" $script"_one" "(k4e,Inclusive,kSM,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k4mu,Inclusive,kSM,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k2e2mu,Inclusive,kSM,\"Nominal\",\""$fixedDate"\")"

. submit.sh $script".cc" $script"_one" "(k4e,Inclusive,kL1,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k4mu,Inclusive,kL1,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k2e2mu,Inclusive,kL1,\"Nominal\",\""$fixedDate"\")"

. submit.sh $script".cc" $script"_one" "(k4e,Inclusive,kA2,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k4mu,Inclusive,kA2,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k2e2mu,Inclusive,kA2,\"Nominal\",\""$fixedDate"\")"

. submit.sh $script".cc" $script"_one" "(k4e,Inclusive,kA3,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k4mu,Inclusive,kA3,\"Nominal\",\""$fixedDate"\")"
. submit.sh $script".cc" $script"_one" "(k2e2mu,Inclusive,kA3,\"Nominal\",\""$fixedDate"\")"
