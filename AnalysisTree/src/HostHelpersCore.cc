#include "HostHelpersCore.h"
#include "TFile.h"
#include "TSystem.h"
#if defined(R__UNIX)
#include "TUnixSystem.h"
#define ROOT_GSYSTEM_TYPE TUnixSystem
#elif defined(R__WIN32)
#include "TWinNTSystem.h"
#define ROOT_GSYSTEM_TYPE TWinNTSystem
#endif
#include "MELAStreamHelpers.hh"


HostHelpers::Hosts HostHelpers::GetHostLocation(){
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  TString strhost = hostname;
  if (strhost.Contains("lxplus") || strhost.Contains("cern")) return kLXPLUS;
  else if (strhost.Contains("login-node") || strhost.Contains("bc-login") || strhost.Contains("gateway") || strhost.Contains("compute") || strhost.Contains("bigmem")) return kMARCC;
  else if (strhost.Contains("t2.ucsd.edu")) return kUCSDT2;
  else return kUNKNOWN;
}
int HostHelpers::GetCurrentDirectory(TString& dirname){
  char cpath[FILENAME_MAX];
  if (!getcwd(cpath, FILENAME_MAX)) return errno;
  dirname = cpath;
  return 0;
}

bool HostHelpers::DirectoryExists(const char* dirname){
  if (!dirname) return false;
  struct stat sb;
  return (stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode));
}
bool HostHelpers::FileExists(const char* fname){
  if (!fname) return false;
  struct stat sb;
  return (stat(fname, &sb) == 0 && S_ISREG(sb.st_mode));
}
bool HostHelpers::FileReadable(const char* fname){
  if (!fname) return false;
  if (FileExists(fname)) return (access(fname, R_OK)==0);
  else if (strstr(fname, "root://")){
    MELAStreamHelpers::MELAout
      << "HostHelpers::FileReadable: "
      << fname << " is remote, so will check if it is readable by explicitly opening it first..."
      << std::endl;

    bool res = false;
    const char strext_root[]=".root";
    const char* strext_isroot = strstr(fname, strext_root);
    if (strext_isroot && strcmp(strext_isroot, strext_root)==0){
      TFile* ftmp = TFile::Open(fname, "read");
      if (ftmp){
        if (ftmp->IsOpen()){
          if (!ftmp->IsZombie()) res = true;
          ftmp->Close();
        }
        else delete ftmp;
      }
    }
    return res;
  }
  else return false;
}
int HostHelpers::ExecuteCommand(const char* strCmd){
  if (system(nullptr)) return system(strCmd);
  else{
    MELAStreamHelpers::MELAerr << "HostHelpers::ExecuteCommand failed because the processor is not available!" << std::endl;
    return -1;
  }
}
int HostHelpers::ExecuteCommand(TString strCmd){ return ExecuteCommand(strCmd.Data()); }

time_t HostHelpers::GetTimestamp(const char* fname){
  struct stat sb;
  int ierr = stat(fname, &sb);
  if (ierr != 0) return 0;
  else return sb.st_mtime;
}
TString HostHelpers::GetTimestampConverted(const char* fname){
  time_t t = GetTimestamp(fname);
  struct tm* tt = localtime(&t);
  char buf[200];
  strftime(buf, sizeof(buf), "%d.%m.%Y %H:%M:%S", tt);
  TString res = buf;
  return res;
}
