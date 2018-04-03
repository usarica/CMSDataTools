#ifndef HOSTHELPERS_H
#define HOSTHELPERS_H

#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <cstdio>
#include <ctime>
#include "TString.h"
#include <cassert>

namespace HostHelpers{

  enum Hosts{
    kLXPLUS,
    kMARCC,
    kUNKNOWN
  };

  Hosts GetHostLocation();
  bool DirectoryExists(const char* dirname);

  time_t GetTimestamp(const char* fname);
  TString GetTimestampConverted(const char* fname);

  TString GetCJLSTSamplesDirectory(const TString proddir);

}


#endif
