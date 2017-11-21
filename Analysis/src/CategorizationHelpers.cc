#include <cmath>
#include "CategorizationHelpers.h"


using namespace CategorizationHelpers;


TString CategorizationHelpers::getCategoryName(CategorizationHelpers::Category category){
  switch(category){
  case Inclusive:
    return "Inclusive";
  case Untagged:
    return "Untagged";
  case JJVBFTagged:
    return "JJVBFTagged";
  case HadVHTagged:
    return "HadVHTagged";
  default:
    return "";
  }
}

CategorizationHelpers::Category CategorizationHelpers::getCategory(const float& DjjVBF, const float& DjjZH, const float& DjjWH, const bool forceUntagged){
  if (
    forceUntagged
    ||
    (DjjVBF<0. && DjjZH<0. && DjjWH<0.)
    ) return Untagged;
  else if (
    DjjVBF>=0.5
    ) return JJVBFTagged;
  else if (
    std::max(DjjZH, DjjWH)>=0.5
    ) return HadVHTagged;
  else return Untagged;
}

