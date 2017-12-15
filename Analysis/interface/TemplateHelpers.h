#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "DiscriminantClasses.h"
#include "SampleHelpers.h"
#include "CategorizationHelpers.h"
#include "ACHypothesisHelpers.h"


namespace TemplateHelpers{

  /*********************/
  /* General functions */
  /*********************/
  template<typename T> void setTemplateAxisLabels(T* histo);
  void getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void getCategorizationDiscriminants(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);
  void adjustDiscriminantJECVariables(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist);

  /****************/
  /* Gluon fusion */
  /****************/
  enum GGHypothesisType{
    GGBkg=0,
    GGSig=1, // fai=0
    GGBSI=2, // fai=0
    nGGSMTypes=3,

    GGSigBSM=3, // fai=1 sig.
    GGSigBSMSMInt=4, // fai=0.5 sig.
    GGBBI=5, // fai=1 BSI
    nGGTypes=6
  };
  enum GGTemplateType{
    GGTplBkg=0,
    GGTplSig=1, // fai=0
    GGTplInt_Re=2, // fai=0
    nGGTplSMTypes=3, // fai=0 int.

    GGTplSigBSM=3, // fai=1 sig.
    GGTplSigBSMSMInt_Re=4, // fai=0.5 sig.
    GGTplIntBSM_Re=5, // fai=1 int.
    nGGTplTypes=6
  };

  TString getGGProcessName(bool useOffshell);
  TString getGGOutputTreeName(TemplateHelpers::GGHypothesisType type, bool useOffshell);
  TString getGGTemplateName(TemplateHelpers::GGTemplateType type, bool useOffshell);
  TString getMELAGGHypothesisWeight(TemplateHelpers::GGHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo);
  std::vector<TemplateHelpers::GGHypothesisType> getGGHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo);
  TString getGGProcessLabel(TemplateHelpers::GGHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo);
  TString getGGProcessLabel(TemplateHelpers::GGTemplateType type, ACHypothesisHelpers::ACHypothesis hypo);

  int castGGHypothesisTypeToInt(TemplateHelpers::GGHypothesisType type);
  int castGGTemplateTypeToInt(TemplateHelpers::GGTemplateType type);
  GGHypothesisType castIntToGGHypothesisType(int type, bool useN=false);
  GGTemplateType castIntToGGTemplateType(int type, bool useN=false);

  template<typename T> void recombineGGHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineGGHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineGGHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineGGHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineGGHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);


  /****************/
  /* Gluon fusion */
  /****************/
  enum VVHypothesisType{
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
  enum VVTemplateType{
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

  TString getVVProcessName(bool useOffshell);
  TString getVVOutputTreeName(TemplateHelpers::VVHypothesisType type, bool useOffshell);
  TString getVVTemplateName(TemplateHelpers::VVTemplateType type, bool useOffshell);
  TString getMELAVVHypothesisWeight(TemplateHelpers::VVHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo);
  std::vector<TemplateHelpers::VVHypothesisType> getVVHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo);
  TString getVVProcessLabel(TemplateHelpers::VVHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo);
  TString getVVProcessLabel(TemplateHelpers::VVTemplateType type, ACHypothesisHelpers::ACHypothesis hypo);

  int castVVHypothesisTypeToInt(TemplateHelpers::VVHypothesisType type);
  int castVVTemplateTypeToInt(TemplateHelpers::VVTemplateType type);
  VVHypothesisType castIntToVVHypothesisType(int type, bool useN=false);
  VVTemplateType castIntToVVTemplateType(int type, bool useN=false);

  template<typename T> void recombineVVHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineVVHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineVVHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineVVHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);
  template<> void recombineVVHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo);


  /*****************/
  /* QQ background */
  /*****************/
  enum QQBkgHypothesisType{
    QQBkg=0,
    nQQBkgTypes
  };
  enum QQBkgTemplateType{
    QQBkgTpl=0,
    nQQBkgTplTypes
  };

  TString getQQBkgProcessName(bool useOffshell);
  TString getQQBkgOutputTreeName(bool useOffshell);
  TString getQQBkgTemplateName(bool useOffshell);
  TString getQQBkgProcessLabel();

  int castQQBkgHypothesisTypeToInt(TemplateHelpers::QQBkgHypothesisType type);
  int castQQBkgTemplateTypeToInt(TemplateHelpers::QQBkgTemplateType type);
  QQBkgHypothesisType castIntToQQBkgHypothesisType(int type);
  QQBkgTemplateType castIntToQQBkgTemplateType(int type);

  template<typename T> void recombineQQBkgHistogramsToTemplates(std::vector<T*>& vals);

}

#endif
