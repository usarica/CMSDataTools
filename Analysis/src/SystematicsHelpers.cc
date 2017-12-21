#include "SystematicsHelpers.h"


SystematicsHelpers::ProcessSystematic::ProcessSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&)) : rule(infcn), evaluators(inEvaluators) {}
float SystematicsHelpers::ProcessSystematic::eval(CJLSTTree* theTree) const{
  float weight=0;
  if (rule) weight=rule(theTree, evaluators);
  return weight;
}

float SystematicsHelpers::getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  float weight=1;
  for (auto const& builder:builders) weight *= builder->getPostThresholdWeight(theTree);
  return weight;
}
float SystematicsHelpers::getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  float weight=0;
  assert(builders.size()==2 && builders.at(0) && builders.at(1));
  float sumWeights = builders.at(0)->getSumPostThresholdWeights(theTree);
  if (sumWeights!=0.) weight = builders.at(0)->getPostThresholdWeight(theTree)*builders.at(1)->getSumPostThresholdWeights(theTree)/sumWeights;
  return weight;
}
