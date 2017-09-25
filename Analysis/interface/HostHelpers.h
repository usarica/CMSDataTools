#ifndef HOSTHELPERS_H
#define HOSTHELPERS_H

#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <cstdio>
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

  TString GetCJLSTSamplesDirectory(const TString proddir);

}


#endif
