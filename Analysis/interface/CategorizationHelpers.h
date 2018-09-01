#ifndef CATEGORIZATIONHELPERS_H
#define CATEGORIZATIONHELPERS_H

#include <vector>
#include "TString.h"
#include "SimpleEntry.h"
#include "DiscriminantClasses.h"
#include "SystematicVariations.h"


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
    UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag,
    UntaggedOrJJVBFOrHadVH_Arbitrated
  };

  extern CategorizationScheme globalCategorizationScheme;
  void setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme);
  std::vector<CategorizationHelpers::Category> getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme);
  bool testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory);

  TString getCategoryName(CategorizationHelpers::Category category);
  TString getCategoryLabel(CategorizationHelpers::Category category);
  TString getCategoryLabelForKDs(CategorizationHelpers::Category category);

  CategorizationHelpers::Category getCategory(SimpleEntry const& sel_vars, const bool forceUntagged);

  TString getMassRegionName(CategorizationHelpers::MassRegion massregion);

  template<typename var_t> void getExtraCategorizationVariables(CategorizationHelpers::CategorizationScheme scheme, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& res);
  template<> void getExtraCategorizationVariables<short>(CategorizationHelpers::CategorizationScheme scheme, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& res);

}


#endif
