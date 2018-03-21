#ifndef SYSTEMATICSHELPERS_H
#define SYSTEMATICSHELPERS_H

#include "SystematicVariations.h"
#include "SampleHelpers.h"
#include "ReweightingBuilder.h"
#include "DiscriminantClasses.h"
#include "ProcessHandler.h"
#include "CategorizationHelpers.h"


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
    std::pair<float, float> (*rule)(short const&, short const&, std::vector<std::vector<float>*> const&, unsigned int const);
    std::vector<TString> strVars;
    unsigned int const id_requested;
    bool doUp;
    std::unordered_map<CJLSTTree*,
      std::pair<
      std::vector<short*>, std::vector<std::vector<float>*>
      >
    > componentRefs;
  public:
    PerLeptonSystematic(
      const std::vector<TString>& inStrVars,
      std::pair<float, float>(*infcn)(short const&, short const&, std::vector<std::vector<float>*> const&, unsigned int const),
      unsigned int const id_requested_,
      bool doUp_
    );
    virtual ~PerLeptonSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
  };

  std::pair<float, float> getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const& LepVars, unsigned int const idreq);


  int convertSystematicVariationTypeToInt(SystematicsHelpers::SystematicVariationTypes type);

  std::vector<SystematicsHelpers::SystematicVariationTypes> getProcessSystematicVariations(
    CategorizationHelpers::Category const category,
    SampleHelpers::Channel const channel,
    ProcessHandler::ProcessType const proc,
    TString strGenerator
  );

  bool systematicAllowed(
    CategorizationHelpers::Category const category,
    SampleHelpers::Channel const channel,
    ProcessHandler::ProcessType const proc,
    SystematicsHelpers::SystematicVariationTypes const syst,
    TString strGenerator
  );

  SystematicsClass* constructSystematic(
    CategorizationHelpers::Category const category,
    SampleHelpers::Channel const channel,
    ProcessHandler::ProcessType const proc,
    SystematicsHelpers::SystematicVariationTypes const syst,
    std::vector<CJLSTTree*> trees,
    std::vector<ReweightingBuilder*>& extraEvaluators,
    TString strGenerator
  );

  void adjustDiscriminantJECVariables(SystematicsHelpers::SystematicVariationTypes const syst, std::vector<DiscriminantClasses::KDspecs>& KDlist);

  TString getSystematicsName(SystematicsHelpers::SystematicVariationTypes const syst);

  TString getSystematicsCombineName(
    CategorizationHelpers::Category const category,
    SampleHelpers::Channel const channel,
    ProcessHandler::ProcessType const proc,
    SystematicsHelpers::SystematicVariationTypes const syst
  );

}

#endif
