#ifndef SYSTEMATICSHELPERS_H
#define SYSTEMATICSHELPERS_H

#include "ReweightingBuilder.h"


namespace SystematicsHelpers{

  class SystematicsClass{
  public:
    SystematicsClass(){}
    virtual ~SystematicsClass(){}
    virtual float eval(CJLSTTree* theTree) const=0;
  };

  class YieldSystematic : public SystematicsClass{
  protected:
    float(*rule)(CJLSTTree*, const std::vector<ReweightingBuilder*>&);
    std::vector<ReweightingBuilder*> evaluators;

  public:
    YieldSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&));
    virtual ~YieldSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
  };

  float getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);
  float getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders);


  class PerLeptonSystematic : public SystematicsClass{
  protected:
    std::pair<float, float> (*rule)(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const&);
    std::vector<TString> strVars;
    bool doUp;
    std::unordered_map<CJLSTTree*,
      std::pair<
      std::vector<short*>, std::vector<std::vector<float>*>
      >
    > componentRefs;
  public:
    PerLeptonSystematic(const std::vector<TString>& inStrVars, std::pair<float, float> (*infcn)(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const&), bool doUp_);
    virtual ~PerLeptonSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
  };

  std::pair<float, float> getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const& LepVars);


}

#endif
