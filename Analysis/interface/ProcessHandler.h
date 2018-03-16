#ifndef PROCESSHANDLER_H
#define PROCESSHANDLER_H

#include "HelperFunctions.h"
#include "ACHypothesisHelpers.h"
#include "ZXFakeRateHandler.h"


class ProcessHandler{
public:
  enum ProcessType{
    kGG,
    kVV,
    kVBF,
    kZH,
    kWH,
    kQQBkg,
    kZX,
    nProcessTypes
  };

  ProcessHandler(ProcessType proctype_, CategorizationHelpers::MassRegion massregion_);

  const ProcessHandler::ProcessType& getProcessType() const;
  const CategorizationHelpers::MassRegion& getProcessMassRegion() const;
  const TString& getProcessName() const;

  virtual void imposeTplPhysicality(std::vector<float>& vals) const;

protected:
  ProcessType const proctype;
  CategorizationHelpers::MassRegion const massregion;
  TString procname;

  void assignProcessName();

};


class GGProcessHandler : public ProcessHandler{
public:
  enum HypothesisType{
    GGBkg=0,
    GGSig=1, // fai=0
    GGBSI=2, // fai=0
    nGGSMTypes=3,

    GGSigBSM=3, // fai=1 sig.
    GGSigBSMSMInt=4, // fai=0.5 sig.
    GGBBI=5, // fai=1 BSI
    nGGTypes=6
  };
  enum TemplateType{
    GGTplBkg=0,
    GGTplSig=1, // fai=0
    GGTplInt_Re=2, // fai=0
    nGGTplSMTypes=3, // fai=0 int.

    GGTplSigBSM=3, // fai=1 sig.
    GGTplSigBSMSMInt_Re=4, // fai=0.5 sig.
    GGTplIntBSM_Re=5, // fai=1 int.
    nGGTplTypes=6
  };

  struct TemplateContributionList{
    GGProcessHandler::TemplateType type;
    float coefficient;
    std::vector<std::pair<GGProcessHandler::TemplateType, float>> TypePowerPair;
    TemplateContributionList(GGProcessHandler::TemplateType type_);
  };

  GGProcessHandler(CategorizationHelpers::MassRegion massregion_);

  TString getOutputTreeName(GGProcessHandler::HypothesisType type) const;
  TString getTemplateName(GGProcessHandler::TemplateType type) const;
  TString getMELAHypothesisWeight(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
  std::vector<GGProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
  TString getProcessLabel(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
  TString getProcessLabel(GGProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

  static int castHypothesisTypeToInt(GGProcessHandler::HypothesisType type);
  static int castTemplateTypeToInt(GGProcessHandler::TemplateType type);
  static GGProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
  static GGProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
  static bool isInterferenceContribution(GGProcessHandler::TemplateType const type);

  void imposeTplPhysicality(std::vector<float>& vals) const;
  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<typename T> void recombineTemplatesWithPhaseRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

};
template<> void GGProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

template<typename T> void GGProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
  if (vals.empty()) return;

  if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nGGSMTypes);
  else assert(vals.size()==nGGTypes);

  vector<vector<unsigned int>> divideByTpl;
  divideByTpl.assign(vals.size(), vector<unsigned int>());
  for (unsigned int t=0; t<vals.size(); t++){
    if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplInt_Re)){
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplBkg));
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSig));
    }
    else if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplSigBSMSMInt_Re)){
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSig));
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSigBSM));
    }
    else if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplIntBSM_Re)){
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplBkg));
      divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSigBSM));
    }
    else divideByTpl.at(t).push_back(t);
  }
  for (unsigned int t=0; t<vals.size(); t++){
    if (divideByTpl.at(t).size()==1) continue;
    vector<std::pair<T, float>> ctpls; ctpls.reserve(divideByTpl.at(t).size());
    for (unsigned int& ht:divideByTpl.at(t)) ctpls.push_back(std::pair<T, float>(vals.at(ht), 0.5));
    HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &ctpls);
  }
  for (unsigned int t=0; t<vals.size(); t++){
    if (divideByTpl.at(t).size()!=1) continue;
    HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
  }
}
template void GGProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
template void GGProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


class VVProcessHandler : public ProcessHandler{
public:
  enum HypothesisType{
    VVBkg=0,
    VVSig=1, // fai=0
    VVBSI=2, // fai=0
    nVVSMTypes=3,

    VVSigBSM=3, // fai=1 sig.
    VVSigBSMSMInt0p25=4, // fai=0.25 sig.
    VVSigBSMSMInt0p5=5, // fai=0.5 sig.
    VVSigBSMSMInt0p75=6, // fai=0.75 sig.
    VVBBI=7, // fai=1 BSI
    VVBMI=8, // fai=0.5 BSI
    nVVTypes=9
  };
  enum TemplateType{
    VVTplBkg=0, // ai**0 a1**0 B**2
    VVTplSig=1, // ai**0 a1**4 B**0
    VVTplInt_Re=2, // ai**0 a1**2 B**1
    nVVTplSMTypes=3,

    VVTplSigBSM=3, // ai**4 a1**0 B**0
    VVTplSigBSMSMInt_ai1_1_Re=4, // ai**1 a1**3 B**0
    VVTplSigBSMSMInt_ai1_2_PosDef=5, // ai**2 a1**2 B**0
    VVTplSigBSMSMInt_ai1_3_Re=6, // ai**3 a1**1 B**0
    VVTplIntBSM_ai1_1_Re=7, // ai**1 a1**1 B**1
    VVTplIntBSM_ai1_2_Re=8, // ai**2 a1**0 B**1
    nVVTplTypes=9
  };

  struct TemplateContributionList{
    VVProcessHandler::TemplateType type;
    float coefficient;
    std::vector<std::pair<VVProcessHandler::TemplateType, float>> TypePowerPair;
    TemplateContributionList(VVProcessHandler::TemplateType type_);
  };

  VVProcessHandler(CategorizationHelpers::MassRegion massregion_, ProcessHandler::ProcessType proctype_=ProcessHandler::kVV);

  TString getOutputTreeName(VVProcessHandler::HypothesisType type) const;
  TString getTemplateName(VVProcessHandler::TemplateType type) const;
  TString getMELAHypothesisWeight(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
  std::vector<VVProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
  TString getProcessLabel(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
  TString getProcessLabel(VVProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

  static int castHypothesisTypeToInt(VVProcessHandler::HypothesisType type);
  static int castTemplateTypeToInt(VVProcessHandler::TemplateType type);
  static VVProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
  static VVProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
  static bool isInterferenceContribution(VVProcessHandler::TemplateType const type);

  void imposeTplPhysicality(std::vector<float>& vals) const;
  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<typename T> void recombineTemplatesWithPhaseRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

};
template<> void VVProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineTemplatesWithPhaseRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

template<typename T> void VVProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
  if (vals.empty()) return;

  if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nVVSMTypes);
  else assert(vals.size()==nVVTypes);

  vector<vector<std::pair<T, float>>> divideByTpl;
  divideByTpl.assign(vals.size(), vector<std::pair<T, float>>());
  for (unsigned int t=0; t<vals.size(); t++){
    if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplInt_Re)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.5));
    }
    else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_1_Re)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.75));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.25));
    }
    else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_2_PosDef)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.5));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.5));
    }
    else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_3_Re)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.25));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.75));
    }
    else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplIntBSM_ai1_1_Re)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.25));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.25));
    }
    else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplIntBSM_ai1_2_Re)){
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
      divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.5));
    }
    else divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(t), 1));
  }
  for (unsigned int t=0; t<vals.size(); t++){
    if (divideByTpl.at(t).size()==1) continue;
    HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &(divideByTpl.at(t)));
  }
  for (unsigned int t=0; t<vals.size(); t++){
    if (divideByTpl.at(t).size()!=1) continue;
    HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
  }
}
template void VVProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
template void VVProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


class QQBkgProcessHandler : public ProcessHandler{
public:
  enum HypothesisType{
    QQBkg=0,
    nQQBkgTypes
  };
  enum TemplateType{
    QQBkgTpl=0,
    nQQBkgTplTypes
  };

  QQBkgProcessHandler(CategorizationHelpers::MassRegion massregion_);

  TString getOutputTreeName() const;
  TString getTemplateName() const;
  TString getMELAHypothesisWeight(unsigned int njets) const;
  TString getProcessLabel() const;
  std::vector<QQBkgProcessHandler::HypothesisType> getHypotheses() const;

  static int castHypothesisTypeToInt(QQBkgProcessHandler::HypothesisType type);
  static int castTemplateTypeToInt(QQBkgProcessHandler::TemplateType type);
  static QQBkgProcessHandler::HypothesisType castIntToHypothesisType(int type);
  static QQBkgProcessHandler::TemplateType castIntToTemplateType(int type);

  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals) const;

  template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

};
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals) const;
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals) const;
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals) const;

template<typename T> void QQBkgProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
  if (vals.empty()) return;
  for (T& hh:vals) HelperFunctions::conditionalizeHistogram(hh, iaxis);
}
template void QQBkgProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
template void QQBkgProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

class ZXProcessHandler : public ProcessHandler{
public:
  enum HypothesisType{
    ZX=0,
    nZXTypes
  };
  enum TemplateType{
    ZXTpl=0,
    nZXTplTypes
  };

  ZXProcessHandler(CategorizationHelpers::MassRegion massregion_);

  TString getOutputTreeName() const;
  TString getTemplateName() const;
  TString getProcessLabel() const;

  static int castHypothesisTypeToInt(ZXProcessHandler::HypothesisType type);
  static int castTemplateTypeToInt(ZXProcessHandler::TemplateType type);
  static ZXProcessHandler::HypothesisType castIntToHypothesisType(int type);
  static ZXProcessHandler::TemplateType castIntToTemplateType(int type);

  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals) const;

  template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

};
template<> void ZXProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals) const;
template<> void ZXProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals) const;
template<> void ZXProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals) const;

template<typename T> void ZXProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
  if (vals.empty()) return;
  for (T& hh:vals) HelperFunctions::conditionalizeHistogram(hh, iaxis);
}
template void ZXProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
template void ZXProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


#endif
