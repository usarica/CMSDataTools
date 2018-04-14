#ifndef DISCRIMINANTCLASSES_H
#define DISCRIMINANTCLASSES_H

#include "HelperFunctionsCore.h"
#include "SimpleDiscriminant.h"
#include "PA2PB1Discriminant.h"
#include "PA2PB2Discriminant.h"
#include "PA1PB1PBp1Discriminant.h"
#include "PA1PB2PBp2Discriminant.h"
#include "VHProdDecACDiscriminant.h"
#include "JJEWQCDBkgDiscriminant.h"
#include "JJEWQCDBkgM4LDiscriminant.h"
#include "SimpleInterferenceDiscriminant.h"
#include "VHProdIntACDiscriminant.h"
#include "JJEWQCDInterferenceDiscriminant.h"
#include "SimpleInterferenceTrigPhase.h"
#include "SimpleAverageInterferenceTrigPhase.h"
#include "VHProdIntTrigPhaseACDiscriminant.h"


namespace DiscriminantClasses{
  enum Type{
    kDbkgkin,
    kDbkgdec,
    kDbkgm4l,

    kDggbkgkin,
    kDggint,
    kCggint,

    kDjVBF,
    kDjjVBF,
    kDjjZH,
    kDjjWH,

    kDjjVBFL1,
    kDjjZHL1,
    kDjjWHL1,

    kDjjVBFa2,
    kDjjZHa2,
    kDjjWHa2,

    kDjjVBFa3,
    kDjjZHa3,
    kDjjWHa3,

    kDbkgjjEWQCD,
    kDbkgm4ljjEWQCD,
    kDintjjEWQCD,
    kCjjVBFint,
    kCjjVHint,

    kDL1dec,
    kDL1decint,
    kCL1decint,
    kDa2dec,
    kDa2decint,
    kCa2decint,
    kDa3dec,
    kDa3decint,
    kCa3decint,

    kDL1jjVBFdec,
    kDL1jjVBFint,
    kCL1jjVBFint,
    kDa2jjVBFdec,
    kDa2jjVBFint,
    kCa2jjVBFint,
    kDa3jjVBFdec,
    kDa3jjVBFint,
    kCa3jjVBFint,

    kDL1jjVHdec,
    kDL1jjVHint,
    kCL1jjVHint,
    kDa2jjVHdec,
    kDa2jjVHint,
    kCa2jjVHint,
    kDa3jjVHdec,
    kDa3jjVHint,
    kCa3jjVHint,

    kNTypes
  };

  struct KDspecs{
    TString KDname;
    Type KDtype;
    std::vector<TString> KDvars;

    Discriminant* KD;

    KDspecs();
    KDspecs(TString strname);
    KDspecs(DiscriminantClasses::Type type);
    bool isValid() const;
  };

  extern const std::unordered_map<TString, DiscriminantClasses::Type> mapKDNameType;
  std::unordered_map<TString, DiscriminantClasses::Type> getKDNameTypeMap();

  DiscriminantClasses::Type getKDType(const TString name);
  TString getKDName(DiscriminantClasses::Type type);
  TString getKDLabel(DiscriminantClasses::Type type);
  TString getKDLabel(const TString name);
  float getKDWP(DiscriminantClasses::Type type);
  float getKDWP(const TString name);

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

  bool isCPSensitive(const Type type);
  bool isCPSensitive(const TString name);

};


#endif
