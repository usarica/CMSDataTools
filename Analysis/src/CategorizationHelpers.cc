#include "CategorizationHelpers.h"


using namespace CategorizationHelpers;


enum WidthCategory{
  Inclusive,
  Untagged,
  JJVBFTagged,
  HadVHTagged,
  nCategoriesMor17
};

TString CategorizationHelpers::getWidthCategoryName(CategorizationHelpers::WidthCategory category){
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

