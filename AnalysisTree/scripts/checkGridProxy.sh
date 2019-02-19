#!/bin/bash

USER_ID=$(id -u)
proxy_valid=$(voms-proxy-info --timeleft)
if [ -z $proxy_valid ] || [ $proxy_valid -le 10 ];then
  voms-proxy-init --rfc --voms cms --hours 72 --valid 72:00
  proxy_valid=$(voms-proxy-info --timeleft)
fi
proxy_file="x509up_u${USER_ID}"
if [ $proxy_valid -gt 10 ];then
   echo "GRID proxy found, validity: $proxy_valid s"
   if [ ! -z $X509_USER_PROXY ];then 
     if [[ "$X509_USER_PROXY" != "~/$proxy_file" ]];then
       let doCopy=1
       if [ -f ~/$proxy_file ];then
         cmp $X509_USER_PROXY ~/$proxy_file &> /dev/null
         if [ $? -eq 0 ];then
           let doCopy=0
           echo "Proxy files $X509_USER_PROXY and ~/$proxy_file are the same."
         fi
       fi
       if [ $doCopy -eq 1 ];then
         echo "Copying proxy file $X509_USER_PROXY -> ~/$proxy_file"
         cp $X509_USER_PROXY ~/
       fi
     fi
   elif [ -f "/tmp/$proxy_file" ];then
     let doCopy=1
     if [ -f ~/$proxy_file ];then
       cmp /tmp/$proxy_file ~/$proxy_file &> /dev/null
       if [ $? -eq 0 ];then
         let doCopy=0
         echo "Proxy files /tmp/$proxy_file and ~/$proxy_file are the same."
       fi
     else
       echo "Proxy file ~/$proxy_file does not exist."
     fi
     if [ $doCopy -eq 1 ];then
       echo "Copying proxy file /tmp/$proxy_file -> ~/$proxy_file"
       cp /tmp/$proxy_file ~/
     fi
   fi
else # Last attempt: Try to see if a valid proxy already exists in ~
   if [ -f ~/$proxy_file ];then
      export X509_USER_PROXY="~/$proxy_file"
      proxy_valid=${voms-proxy-info --timeleft}
      if [ $proxy_valid -gt 10 ];then
         echo "GRID proxy found in ~, validity: $proxy_valid s"
      fi
   fi 

   if [ $proxy_valid -le 10 ];then
      echo "Error: no valid GRID proxy found."
      exit
   fi
fi
