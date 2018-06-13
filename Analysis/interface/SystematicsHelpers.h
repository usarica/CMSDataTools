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
  public:
    typedef std::pair<float, float> (*PerLeptonSystematicFunction_t)(short const&, short const&, std::vector<std::vector<float>* const*> const&, unsigned int const);

  protected:
    PerLeptonSystematicFunction_t rule;
    std::vector<TString> strVars;
    unsigned int const id_requested;
    bool doUp;
    std::unordered_map<CJLSTTree*,
      std::pair<
      std::vector<short*>, std::vector<std::vector<float>* const*>
      >
    > componentRefs;
  public:
    PerLeptonSystematic(
      const std::vector<TString>& inStrVars,
      PerLeptonSystematicFunction_t infcn,
      unsigned int const id_requested_,
      bool doUp_
    );
    virtual ~PerLeptonSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
  };

  std::pair<float, float> getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq);


  class PerLeptonScaleSystematic : public SystematicsClass{
  public:
    typedef std::pair<float, float>(*PerLeptonScaleSystematicFunction_t)(short const&, short const&, float const&, std::vector<std::vector<float>* const*> const&, unsigned int const);

  protected:
    struct componentData{
      std::vector<short*> ref_shorts;
      std::vector<float*> ref_floats;
      std::vector<std::vector<float>* const*> ref_vfloats;
      componentData(){}
      componentData(std::vector<short*> const& rs, std::vector<float*> const& rf, std::vector<std::vector<float>* const*> const& rvf) : ref_shorts(rs), ref_floats(rf), ref_vfloats(rvf){}
      componentData(componentData const& other) : ref_shorts(other.ref_shorts), ref_floats(other.ref_floats), ref_vfloats(other.ref_vfloats){}
    };

    PerLeptonScaleSystematicFunction_t rule;
    std::vector<TString> strVars;
    unsigned int const id_requested;
    bool doUp;
    std::unordered_map<CJLSTTree*, componentData> componentRefs;
  public:
    PerLeptonScaleSystematic(
      const std::vector<TString>& inStrVars,
      PerLeptonScaleSystematicFunction_t infcn,
      unsigned int const id_requested_,
      bool doUp_
    );
    virtual ~PerLeptonScaleSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
  };
  class PerLeptonResSystematic : public SystematicsClass{
  public:
    typedef std::pair<float, float>(*PerLeptonResSystematicFunction_t)(short const&, short const&, float const&, float const&, std::vector<std::vector<float>* const*> const&, unsigned int const);

  protected:
    struct componentData{
      std::vector<short*> ref_shorts;
      std::vector<float*> ref_floats;
      std::vector<std::vector<float>* const*> ref_vfloats;
      componentData(){}
      componentData(std::vector<short*> const& rs, std::vector<float*> const& rf, std::vector<std::vector<float>* const*> const& rvf) : ref_shorts(rs), ref_floats(rf), ref_vfloats(rvf){}
      componentData(componentData const& other) : ref_shorts(other.ref_shorts), ref_floats(other.ref_floats), ref_vfloats(other.ref_vfloats){}
    };

    PerLeptonResSystematicFunction_t rule;
    std::vector<TString> strVars;
    unsigned int const id_requested;
    bool doUp;
    float centralValue;
    std::unordered_map<CJLSTTree*, componentData> componentRefs;
  public:
    PerLeptonResSystematic(
      const std::vector<TString>& inStrVars,
      PerLeptonResSystematicFunction_t infcn,
      unsigned int const id_requested_,
      bool doUp_
    );
    virtual ~PerLeptonResSystematic(){}
    virtual float eval(CJLSTTree* theTree) const;
    virtual void setup(CJLSTTree* theTree);
    virtual void setCentralValue(float centralValue_){ centralValue=centralValue_; }
  };

  std::pair<float, float> getLeptonScaleSystematic(short const& Z1Flav, short const& Z2Flav, float const& ZZMass, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq);
  std::pair<float, float> getLeptonResSystematic(short const& Z1Flav, short const& Z2Flav, float const& ZZMass, float const& centralValue, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq);


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
  TString getSystematicsLabel(SystematicsHelpers::SystematicVariationTypes const syst);

  TString getSystematicsCombineName(
    CategorizationHelpers::Category const category,
    SampleHelpers::Channel const channel,
    ProcessHandler::ProcessType const proc,
    SystematicsHelpers::SystematicVariationTypes const syst
  );

}

#endif
