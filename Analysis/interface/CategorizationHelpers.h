#ifndef CATEGORIZATIONHELPERS_H
#define CATEGORIZATIONHELPERS_H

#include "TString.h"


namespace CategorizationHelpers{
  enum Category{
    Inclusive,
    Untagged,
    JJVBFTagged,
    HadVHTagged,
    nCategories
  };

  TString getCategoryName(CategorizationHelpers::Category category);

}


#endif
