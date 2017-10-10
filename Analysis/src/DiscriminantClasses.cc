#include "DiscriminantClasses.h"


DiscriminantClasses::Type DiscriminantClasses::getKDType(const TString name){
  if (name=="Dbkgkin") return kDbkgkin;
  else if (name=="Dbkgdec") return kDbkgdec;
  else if (name=="DjVBF") return kDjVBF;
  else if (name=="DjjVBF") return kDjjVBF;
  else if (name=="DjjZH" || name=="DjjWH") return kDjjVH;
  else if (name=="DbkgjjEWQCD") return kDbkgjjEWQCD;
  else return kNTypes;
}
Discriminant* DiscriminantClasses::constructKDFromType(const DiscriminantClasses::Type type, const TString filename, const TString splinename){
  Discriminant* res=nullptr;
  switch (type){
  case kDbkgkin:
    return new Dbkgkin_t(filename, splinename);
  case kDbkgdec:
    return new Dbkgdec_t(filename, splinename);
  case kDjVBF:
    return new DjVBF_t(filename, splinename);
  case kDjjVBF:
    return new DjjVBF_t(filename, splinename);
  case kDjjVH:
    return new DjjVH_t(filename, splinename);
  case kDbkgjjEWQCD:
    return new DbkgjjEWQCD_t(filename, splinename);
  default:
    return res;
  };
}
Discriminant* DiscriminantClasses::constructKDFromType(const TString name, const TString filename, const TString splinename){ return constructKDFromType(getKDType(name), filename, splinename); }

