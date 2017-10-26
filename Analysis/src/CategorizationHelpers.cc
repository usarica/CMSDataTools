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

