#include <cmath>
#include "CategorizationHelpers.h"



namespace CategorizationHelpers{
  CategorizationScheme globalCategorizationScheme=UntaggedOrJJVBFOrHadVH;
}

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

void CategorizationHelpers::setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme){ globalCategorizationScheme=scheme; }

CategorizationHelpers::Category CategorizationHelpers::getCategory(const float& DjjVBF, const float& DjjZH, const float& DjjWH, const bool forceUntagged){
  if (
    forceUntagged
    ||
    (DjjVBF<0. && DjjZH<0. && DjjWH<0.)
    ) return Untagged;
  if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH){
    if (
      DjjVBF>=0.5
      ) return JJVBFTagged;
    else if (
      std::max(DjjZH, DjjWH)>=0.5
      ) return HadVHTagged;
    else return Untagged;
  }
  else if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH_Arbitrated){
    float DjjVH=std::max(DjjZH, DjjWH);
    if (
      DjjVBF>=0.5 && DjjVBF>=DjjVH
      ) return JJVBFTagged;
    else if (
      DjjVH>=0.5 && DjjVBF<DjjVH
      ) return HadVHTagged;
    else return Untagged;
  }
  else/* if (globalCategorizationScheme==UntaggedOrJJVBF)*/{
    if (
      DjjVBF>=0.5
      ) return JJVBFTagged;
    else return Untagged;
  }
}

