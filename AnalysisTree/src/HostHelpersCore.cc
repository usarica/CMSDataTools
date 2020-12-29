#include <cassert>
#include <regex>
#include <signal.h>
#include "HostHelpersCore.h"
#include "HelperFunctionsCore.h"
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
  return GetHostLocation(strhost);
}
HostHelpers::Hosts HostHelpers::GetHostLocation(TString const& strhost){
  if (strhost.Contains("marcc")) return kMARCC;
  else if (strhost.Contains("t2.ucsd.edu")) return kUCSDT2;
  else if (strhost.Contains("iihe.ac.be")) return kIIHET2;
  else if (strhost.Contains("eos.cms")) return kEOSCMS; // Wouldn't really happen, but anyway...
  else if (strhost.Contains("cern")) return kLXPLUS; // Let this come after eos.
  else return kUNKNOWN;
}

TString HostHelpers::GetHostLocalRedirector(Hosts const& host, bool isForFileOps){
  switch (host){
  case kUCSDT2:
    return "davs://redirector.t2.ucsd.edu:1094/"; // "root://redirector.t2.ucsd.edu/" also works but is slower.
  case kIIHET2:
    return "srm://maite.iihe.ac.be:8443"; // No extra '/' at the end in srm protocols!
  case kEOSCMS:
    return "root://eoscms.cern.ch/";
  default:
    MELAStreamHelpers::MELAerr << "HostHelpers::GetHostRemoteConnectionSpecifier: Host " << host << " might not support remote connection to read files. Returning the widest but slowest possible option." << std::endl;
    return (isForFileOps ? "" : "root://cms-xrd-global.cern.ch/");
  }
}
bool HostHelpers::CheckContainsURL(const char* name){
  return name && strstr(name, "://");
}
TString HostHelpers::GetHostPathToStore(Hosts const& host){
  switch (host){
  case kUCSDT2:
    return "/hadoop/cms";
  case kIIHET2:
    return "/pnfs/iihe/cms"; // But we don't take this out in srm protocols.
  case kEOSCMS:
    return "/eos/cms";
  default:
    MELAStreamHelpers::MELAerr << "HostHelpers::GetHostPathToStore: Host " << host << " is not supported. Returning empty string..." << std::endl;
    return "";
  }
}
TString HostHelpers::GetStandardHostPathToStore(const char* name, Hosts const& host){
  if (!name) return "";
  TString res = name;
  // Only convert paths that refer to /store and if they don't already have a redirector.
  if (res.Contains("/store") && !CheckContainsURL(name)){
    TString strstorepath = GetHostPathToStore(host);
    if (host != GetHostLocation()){
      TString strlocalredir = GetHostLocalRedirector(host, true);
      // The path should not include the portion up to '/store' unless the protocol is srm.
      if (!strlocalredir.Contains("srm")) HelperFunctions::replaceString<TString, TString const>(res, strstorepath, "");
      else if (!res.Contains(strstorepath)) res = strstorepath + "/" + res;
      res = strlocalredir + res;
    }
    else if (!res.Contains(strstorepath) && res.BeginsWith("/store")) res = strstorepath + res;
  }
  return res;
}

TString HostHelpers::GetX509Proxy(){
  TString res;
  struct stat sb;
  const char* proxy_by_env = getenv("X509_USER_PROXY");
  if (proxy_by_env) res = proxy_by_env;
  if (stat(res.Data(), &sb) == 0) return res;
  unsigned int uid = geteuid(); // Actually the effective user id, this is what 'id -u' does.
  const char* homedir = getenv("HOME");
  if (homedir && uid) res = Form("%s/x509up_u%i", homedir, uid);
  if (stat(res.Data(), &sb) == 0) return res;
  if (uid) res = Form("/tmp/x509up_u%u", uid);
  if (stat(res.Data(), &sb) == 0) return res;
  MELAStreamHelpers::MELAerr << "HostHelpers::GetX509Proxy: No grid proxy found Please put one in ~, /tmp, or define through the environment variable X509_USER_PROXY." << std::endl;
  assert(0);
  return "";
}

int HostHelpers::GetCurrentDirectory(TString& dirname){
  char cpath[FILENAME_MAX];
  if (!getcwd(cpath, FILENAME_MAX)) return errno;
  dirname = cpath;
  return 0;
}

bool HostHelpers::DirectoryExists(const char* dirname){
  if (!dirname) return false;
  if (CheckContainsURL(dirname)){
    TString strCmd = Form("[[ ! -z $(env -i X509_USER_PROXY=%s gfal-stat %s | grep -i -e \"directory\") ]]", GetX509Proxy().Data(), dirname);
    return ExecuteCommand(strCmd.Data())==0;
  }
  struct stat sb;
  return (stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode));
}
bool HostHelpers::FileExists(const char* fname){
  if (!fname) return false;
  if (CheckContainsURL(fname)){
    TString strCmd = Form("[[ ! -z $(env -i X509_USER_PROXY=%s gfal-stat %s | grep -i -e \"regular file\") ]]", GetX509Proxy().Data(), fname);
    ExpandEnvironmentVariables(strCmd);
    return ExecuteCommand(strCmd.Data())==0;
  }
  struct stat sb;
  return (stat(fname, &sb) == 0 && S_ISREG(sb.st_mode));
}
bool HostHelpers::FileReadable(const char* fname){
  if (!fname) return false;
  if (CheckContainsURL(fname)){
    bool res = false;
    const char strext_root[]=".root";
    const char* strext_isroot = strstr(fname, strext_root);
    if (strext_isroot && strcmp(strext_isroot, strext_root)==0){
      MELAStreamHelpers::MELAout
        << "HostHelpers::FileReadable: "
        << fname << " is remote, so will check if it is readable by explicitly opening it first..."
        << std::endl;
      TFile* ftmp = TFile::Open(fname, "read");
      if (ftmp){
        if (ftmp->IsOpen()){
          if (!ftmp->IsZombie()) res = true;
          ftmp->Close();
        }
        else delete ftmp;
      }
    }
    else{
      TString strCmd = Form("[[ ! -z $(env -i X509_USER_PROXY=%s gfal-stat %s | grep -i -e \"access\" | grep -i -e \"/-r\") ]]", GetX509Proxy().Data(), fname);
      ExpandEnvironmentVariables(strCmd);
      res = (ExecuteCommand(strCmd.Data())==0);
    }
    return res;
  }
  else if (FileExists(fname)) return (access(fname, R_OK)==0);
  else return false;
}
int HostHelpers::ExecuteCommand(const char* strCmd){
  int sys_ret = system(nullptr);
  struct sigaction childact;
  sigaction(SIGCHLD, NULL, &childact);
  if (sys_ret==1 || childact.sa_handler==SIG_IGN) return system(strCmd);
  else{
    MELAStreamHelpers::MELAerr << "HostHelpers::ExecuteCommand failed because the processor is not available (exit code " << sys_ret << ")!" << std::endl;
    return -1;
  }
}
int HostHelpers::ExecuteCommand(TString strCmd){ return ExecuteCommand(strCmd.Data()); }
void HostHelpers::ExpandEnvironmentVariables(TString& str){
  std::string ss = str.Data();
  ExpandEnvironmentVariables(ss);
  str = ss.data();
}
void HostHelpers::ExpandEnvironmentVariables(std::string& str){
  static std::regex env("\\$\\{([^}]+)\\}");
  std::smatch match;
  while (std::regex_search(str, match, env)){
    const char* s = getenv(match[1].str().c_str());
    const std::string var(s == NULL ? "" : s);
    str.replace(match[0].first, match[0].second, var);
  }
}


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
