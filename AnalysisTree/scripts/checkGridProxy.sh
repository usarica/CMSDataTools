#!/bin/bash

USER_ID=$(id -u)
proxy_valid=$(voms-proxy-info --timeleft)
if [[ "$proxy_valid" == "" ]];then
  voms-proxy-init --rfc --voms cms --hours 72 --valid 72:00
  proxy_valid=$(voms-proxy-info --timeleft)
fi
proxy_file="x509up_u$USER_ID"
if [ $proxy_valid > 10 ];then
   echo "GRID proxy found, validity: $proxy_valid s"
   if [ ! -z $X509_USER_PROXY ];then 
     if [[ "$X509_USER_PROXY" != "~/$proxy_file" ]];then
       echo "Copying proxy file $X509_USER_PROXY -> ~/$proxy_file"
       cp $X509_USER_PROXY "~/$proxy_file"
     fi
   elif [ -e "/tmp/$proxy_file" ];then
     echo "Copying proxy file /tmp/$proxy_file -> ~/$proxy_file"
     cp "/tmp/$proxy_file" ~/
   fi
else # Last attempt: Try to see if a valid proxy already exists in ~
   if [ -e "~/$proxy_file" ];then
      export X509_USER_PROXY="~/$proxy_file"
      proxy_valid=${voms-proxy-info --timeleft}
      if [ $proxy_valid > 10 ];then
         echo "GRID proxy found in ~, validity: $proxy_valid s"
      fi
   fi 

   if [ $proxy_valid < 10 ];then
      echo "Error: no valid GRID proxy found."
      exit
   fi
fi
