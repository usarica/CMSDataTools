#ifndef ACHYPOTHESISHELPERS_H
#define ACHYPOTHESISHELPERS_H

#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"


namespace ACHypothesisHelpers{
  enum ACHypothesis{
    kSM=0,
    kL1,
    kA2,
    kA3,
    nACHypotheses
  };

  TString getACHypothesisName(ACHypothesisHelpers::ACHypothesis hypo);

  std::vector<DiscriminantClasses::Type> getACHypothesisKDSet(ACHypothesisHelpers::ACHypothesis hypo, CategorizationHelpers::Category category, CategorizationHelpers::MassRegion massregion);
  std::vector<TString> getACHypothesisKDNameSet(ACHypothesisHelpers::ACHypothesis hypo, CategorizationHelpers::Category category, CategorizationHelpers::MassRegion massregion);

  float getACHypothesisMEHZZGVal(ACHypothesisHelpers::ACHypothesis hypo);
  float getACHypothesisHZZGVal(ACHypothesisHelpers::ACHypothesis hypo);

}

#endif
