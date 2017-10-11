#ifndef CATEGORIZATIONHELPERS_H
#define CATEGORIZATIONHELPERS_H

#include "TString.h"


namespace CategorizationHelpers{
  enum WidthCategory{
    Inclusive,
    Untagged,
    JJVBFTagged,
    HadVHTagged,
    nCategoriesMor17
  };

  TString getWidthCategoryName(CategorizationHelpers::WidthCategory category);

}


#endif
