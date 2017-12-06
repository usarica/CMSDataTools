#ifndef STDEXTENSIONS_H
#define STDEXTENSIONS_H

#include <string>
#include <functional>
#include "TString.h"

namespace std{

  template<> struct hash<TString>{
    typedef TString argument_type;
    typedef size_t result_type;
    result_type operator()(argument_type const& arg) const{ return hash<string>{}(arg.Data()); }
  };

}

#endif
