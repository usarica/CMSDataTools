#include <cmath>
#include <cassert>
#include "CategorizationHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


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

TString CategorizationHelpers::getCategoryLabel(CategorizationHelpers::Category category){
  switch (category){
  case Inclusive:
    return "Inclusive";
  case Untagged:
    return "Untagged";
  case JJVBFTagged:
    return "VBF 2 jets";
  case HadVHTagged:
    return "Had. VH";
  default:
    return "";
  }
}

TString CategorizationHelpers::getCategoryLabelForKDs(CategorizationHelpers::Category category){
  switch (category){
  case JJVBFTagged:
    return "VBF";
  case HadVHTagged:
    return "VH";
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
std::vector<CategorizationHelpers::Category> CategorizationHelpers::getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme){
  std::vector<CategorizationHelpers::Category> res;
  res.push_back(Inclusive);
  switch (globalCategorizationScheme){
  case UntaggedOrJJVBFOrHadVH:
  case UntaggedOrJJVBFOrHadVH_Arbitrated:
    res.push_back(HadVHTagged);
  case UntaggedOrJJVBF:
    res.push_back(JJVBFTagged);
    res.push_back(Untagged);
  default:
    break;
  }
  return res;
}
bool CategorizationHelpers::testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory){
  std::vector<Category> cats=getAllowedCategories(globalCategorizationScheme);
  for (Category const& cat:cats){ if (cat==theCategory) return true; }
  return false;
}

TString CategorizationHelpers::getMassRegionName(CategorizationHelpers::MassRegion massregion){
  switch (massregion){
  case kOnshell:
    return "Onshell";
  case kOffshell:
    return "Offshell";
  default:
    MELAerr << "CategorizationHelpers::getMassRegionName: Mass region " << massregion << " has undefined name!" << endl;
    assert(0);
    return "";
  }
}

