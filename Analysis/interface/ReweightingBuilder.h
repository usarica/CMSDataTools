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

  ReweightingBuilder(CJLSTTree* intree, TString inStrWeight, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);
  ReweightingBuilder(CJLSTTree* intree, std::vector<TString> inStrWeights, ReweightingBuilder::BWScheme inscheme = ReweightingBuilder::kDoNotConsider);

protected:
  CJLSTTree* const theTree;
  const BWScheme theScheme;
  std::vector<TString> strWeights;

  virtual float eval();
};


#endif
