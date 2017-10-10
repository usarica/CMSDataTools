#ifndef DISCRIMINANTCLASSES_H
#define DISCRIMINANTCLASSES_H

#include "SimpleDiscriminant.h"
#include "PA2PB1Discriminant.h"
#include "PA2PB2Discriminant.h"
#include "PA1PB1PBp1Discriminant.h"
#include "JJEWQCDBkgDiscriminant.h"

namespace DiscriminantClasses{
  enum Type{
    kDbkgkin,
    kDbkgdec,
    kDjVBF,
    kDjjVBF,
    kDjjVH,
    kDbkgjjEWQCD,
    kNTypes
  };

  Type getKDType(const TString name);
  Discriminant* constructKDFromType(const Type type, const TString filename="", const TString splinename="");
  Discriminant* constructKDFromType(const TString name, const TString filename="", const TString splinename="");

};


#endif
