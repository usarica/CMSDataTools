#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ACHypothesisHelpers.h"


namespace TemplateHelpers{
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
    GGTplSigBSMSMInt=4, // fai=0.5 sig.
    GGTplIntBSM_Re=5, // fai=1 int.
    nGGTplTypes=6
  };

  TString getGGProcessName(bool useOffshell);
  TString getGGOutputTreeName(TemplateHelpers::GGHypothesisType type, bool useOffshell);
  TString getGGTemplateName(TemplateHelpers::GGTemplateType type, bool useOffshell);
  TString getMELAGGHypothesisWeight(TemplateHelpers::GGHypothesisType type);

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

}

#endif
