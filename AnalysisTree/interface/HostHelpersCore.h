#ifndef HOSTHELPERSCORE_H
#define HOSTHELPERSCORE_H

#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <cstdio>
#include <ctime>
#include <cassert>
#include "TString.h"

namespace HostHelpers{

  enum Hosts{
    kLXPLUS,
    kMARCC,
    kUCSDT2,
    kUNKNOWN
  };

  Hosts GetHostLocation();
  bool DirectoryExists(const char* dirname);

  time_t GetTimestamp(const char* fname);
  TString GetTimestampConverted(const char* fname);

}


#endif
