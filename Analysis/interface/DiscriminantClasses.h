#ifndef DISCRIMINANTCLASSES_H
#define DISCRIMINANTCLASSES_H

#include "SimpleDiscriminant.h"
#include "PA2PB1Discriminant.h"
#include "PA2PB2Discriminant.h"
#include "PA1PB1PBp1Discriminant.h"
#include "PA1PB2PBp2Discriminant.h"
#include "JJEWQCDBkgDiscriminant.h"
#include "SimpleInterferenceDiscriminant.h"

namespace DiscriminantClasses{
  enum Type{
    kDbkgkin,
    kDbkgdec,

    kDggbkgkin,
    kDggint,

    kDjVBF,
    kDjjVBF,

    kDjjZH,
    kDjjWH,

    kDbkgjjEWQCD,

    kNTypes
  };

  Type getKDType(const TString name);
  Discriminant* constructKDFromType(
    const Type type,
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="", const float gscale=1
  );
  Discriminant* constructKDFromType(
    const TString name,
    const TString cfilename="", const TString splinename="",
    const TString gfilename="", const TString gsplinename="", const float gscale=1
  );

  std::vector<TString> getKDVars(const Type type);
  std::vector<TString> getKDVars(const TString name);

};


#endif
