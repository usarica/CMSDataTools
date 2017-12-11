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
  void getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, std::vector<DiscriminantClasses::KDspecs>& KDlist);

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
