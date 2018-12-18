#include "HostHelpersCore.h"


HostHelpers::Hosts HostHelpers::GetHostLocation(){
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  TString strhost = hostname;
  if (strhost.Contains("lxplus") || strhost.Contains("cern")) return kLXPLUS;
  else if (strhost.Contains("login-node") || strhost.Contains("bc-login") || strhost.Contains("gateway") || strhost.Contains("compute") || strhost.Contains("bigmem")) return kMARCC;
  else if (strhost.Contains("t2.ucsd.edu")) return kUCSDT2;
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
