#ifndef PROCESSHANDLER_H
#define PROCESSHANDLER_H

#include "HelperFunctions.h"
#include "ACHypothesisHelpers.h"


class ProcessHandler{
public:
  enum ProcessType{
    kGG,
    kVV,
    kQQBkg,
    nProcessTypes
  };

  ProcessHandler(ProcessType proctype_, bool useOffshell_);

  const TString& getProcessName() const;
  const ProcessHandler::ProcessType& getProcessType() const;

protected:
  ProcessType const proctype;
  bool const useOffshell;
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

  GGProcessHandler(bool useOffshell_);

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

  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

};
template<> void GGProcessHandler::recombineHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;


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

  VVProcessHandler(bool useOffshell_);

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

  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

};
template<> void VVProcessHandler::recombineHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;


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

  QQBkgProcessHandler(bool useOffshell_);

  TString getOutputTreeName() const;
  TString getTemplateName() const;
  TString getMELAHypothesisWeight(unsigned int njets) const;
  TString getProcessLabel() const;

  static int castHypothesisTypeToInt(QQBkgProcessHandler::HypothesisType type);
  static int castTemplateTypeToInt(QQBkgProcessHandler::TemplateType type);
  static QQBkgProcessHandler::HypothesisType castIntToHypothesisType(int type);
  static QQBkgProcessHandler::TemplateType castIntToTemplateType(int type);

  template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals) const;

};
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals) const;
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals) const;
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals) const;

#endif
