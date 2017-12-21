#ifndef SYSTEMATICSHELPERS_H
#define SYSTEMATICSHELPERS_H

#include "ReweightingBuilder.h"


namespace SystematicsHelpers{

  class SystematicsClass{
  public:
    SystematicsClass(){}
    virtual float eval(CJLSTTree* theTree) const=0;
  };

  class YieldSystematic : public SystematicsClass{
  protected:
    float(*rule)(CJLSTTree*, const std::vector<ReweightingBuilder*>&);
    std::vector<ReweightingBuilder*> evaluators;

  public:
    YieldSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&));
    virtual float eval(CJLSTTree* theTree) const;
  };

  float getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);
  float getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);


  class PerLeptonSystematic : public SystematicsClass{
  protected:
    float(*rule)(std::vector<short>* const&, std::vector<std::vector<float>*> const&);
    TString strLepId;
    std::vector<TString> strVars;
    std::unordered_map<CJLSTTree*,
      std::pair<
      std::vector<short>*, std::vector<std::vector<float>*>
      >
    > componentRefs;
  public:
    PerLeptonSystematic(const TString inStrLepId, const std::vector<TString>& inStrVars, float(*infcn)(std::vector<short>* const&, std::vector<std::vector<float>*> const&));
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
  };


}

#endif
