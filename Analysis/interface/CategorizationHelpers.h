#ifndef CATEGORIZATIONHELPERS_H
#define CATEGORIZATIONHELPERS_H

#include "TString.h"
#include "DiscriminantClasses.h"


namespace CategorizationHelpers{
  enum Category{
    Inclusive,
    Untagged,
    JJVBFTagged,
    HadVHTagged,
    nCategories
  };
  enum MassRegion{
    kOnshell,
    kOffshell,
    NMassRegions
  };
  enum CategorizationScheme{
    UntaggedOrJJVBF,
    UntaggedOrJJVBFOrHadVH,
    UntaggedOrJJVBFOrHadVH_Arbitrated
  };

  extern CategorizationScheme globalCategorizationScheme;
  void setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme);
  std::vector<CategorizationHelpers::Category> getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme);
  bool testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory);

  TString getCategoryName(CategorizationHelpers::Category category);
  TString getCategoryLabel(CategorizationHelpers::Category category);
  TString getCategoryLabelForKDs(CategorizationHelpers::Category category);

  CategorizationHelpers::Category getCategory(const float& DjjVBF, const float& DjjZH, const float& DjjWH, const bool forceUntagged);

  TString getMassRegionName(CategorizationHelpers::MassRegion massregion);

}


#endif
