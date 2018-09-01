#include "HostHelpers.h"
#include <iostream>


HostHelpers::Hosts HostHelpers::GetHostLocation(){
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  TString strhost = hostname;
  if (strhost.Contains("lxplus") || strhost.Contains("cern")) return kLXPLUS;
  else if (strhost.Contains("login-node") || strhost.Contains("bc-login") || strhost.Contains("gateway") || strhost.Contains("compute") || strhost.Contains("bigmem")) return kMARCC;
  else return kUNKNOWN;
}
bool HostHelpers::DirectoryExists(const char* dirname){
  struct stat sb;
  return (stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode));
}

time_t HostHelpers::GetTimestamp(const char* fname){
  struct stat st;
  int ierr=stat(fname, &st);
  if (ierr!=0) return 0;
  else return st.st_mtime;
}
TString HostHelpers::GetTimestampConverted(const char* fname){
  time_t t = GetTimestamp(fname);
  struct tm* tt=localtime(&t);
  char buf[200];
  strftime(buf, sizeof(buf), "%d.%m.%Y %H:%M:%S", tt);
  TString res=buf;
  return res;
}

TString HostHelpers::GetCJLSTSamplesDirectory(const TString proddir){
  const Hosts host = GetHostLocation();

  TString theDir;
  if (host==kLXPLUS) theDir = "/eos/user/u/usarica/CJLST/4l";
  else if (host==kMARCC) theDir = "/work-zfs/lhc/CJLSTtrees";

  if (theDir!=""){
    TString testdir = theDir + "/" + proddir;
    if (DirectoryExists(testdir.Data())) return theDir;
    else if (host==kLXPLUS){
      theDir = "root://lxcms03//data3/Higgs";
      return theDir; // There is no easy way to check if proddir exists over xrootd
    }
  }

  std::cout << "CJLST samples directory could not be found!" << std::endl;
  assert(0);
  return "";
}
