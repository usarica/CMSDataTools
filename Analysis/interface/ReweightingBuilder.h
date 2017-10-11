#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "CJLSTTree.h"


class ReweightingBuilder{
public:
  enum BWScheme{
    kDoNotConsider,
    kDivideSampleBW,
    kMultiplySampleBW
  };

  ReweightingBuilder(TString inStrWeight, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);
  ReweightingBuilder(std::vector<TString> inStrWeights, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);

protected:
  const BWScheme theScheme;
  std::vector<TString> strWeights;

  virtual float eval(CJLSTTree* theTree);
};


#endif
