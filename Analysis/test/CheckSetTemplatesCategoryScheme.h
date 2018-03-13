#ifndef CHECKSETTEMPLATESCATEGORYSCHEME_H
#define CHECKSETTEMPLATESCATEGORYSCHEME_H

#include "CategorizationHelpers.h"

using namespace std;
using namespace CategorizationHelpers;

// Set categorization scheme for templates and check the current category validity
bool CheckSetTemplatesCategoryScheme(CategorizationHelpers::Category category){
  CategorizationHelpers::setGlobalCategorizationScheme(UntaggedOrJJVBFOrHadVH);
  return testCategoryAgainstGlobalScheme(category);
}


#endif
