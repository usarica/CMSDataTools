#ifndef HOSTHELPERSCORE_H
#define HOSTHELPERSCORE_H

#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <errno.h>
#include <cstdio>
#include <ctime>
#include <cassert>
#include <string>
#include "TString.h"

namespace HostHelpers{

  enum Hosts{
    kLXPLUS,
    kMARCC,
    kUCSDT2,
    kEOSCMS,
    kUNKNOWN
  };

  Hosts GetHostLocation();
  TString GetHostLocalRedirector(Hosts const& host, bool isForFileOps);
  TString GetHostPathToStore(Hosts const& host);

  TString GetX509Proxy();

  int GetCurrentDirectory(TString& dirname);

  bool CheckContainsURL(const char* name);
  bool DirectoryExists(const char* dirname);
  bool FileExists(const char* fname);
  bool FileReadable(const char* fname);
  int ExecuteCommand(const char* strCmd);
  int ExecuteCommand(TString strCmd);
  void ExpandEnvironmentVariables(TString& str);
  void ExpandEnvironmentVariables(std::string& str);

  time_t GetTimestamp(const char* fname);
  TString GetTimestampConverted(const char* fname);

}


#endif
