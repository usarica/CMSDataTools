#include "TemplateHelpers.h"


TString TemplateHelpers::getGGProcessName(bool useOffshell){ return (useOffshell ? "ggZZ_offshell" : "ggH"); }
TString TemplateHelpers::getGGOutputTreeName(TemplateHelpers::GGHypothesisType type, bool useOffshell){
  TString res;
  switch (type){
  case GGBkg:
    res="Bkg";
    break;
  case GGSig:
    res="Sig";
    break;
  case GGBSI:
    res="BSI";
    break;
  case GGSigBSM:
    res="SigBSM";
    break;
  case GGSigBSMSMInt:
    res="SigBSMSMInt";
    break;
  case GGBBI:
    res="BBI";
    break;
  default:
    break;
  };
  if (res!="") res = Form("T_%s_%s_Tree", TemplateHelpers::getGGProcessName(useOffshell).Data(), res.Data());
  return res;
}
TString TemplateHelpers::getGGTemplateName(TemplateHelpers::GGTemplateType type, bool useOffshell){
  TString res;
  switch (type){
  case GGTplBkg:
    res="Bkg";
    break;
  case GGTplSig:
    res="Sig";
    break;
  case GGTplInt_Re:
    res="Int_Re";
    break;
  case GGTplSigBSM:
    res="Sig_ai1_2";
  case GGTplSigBSMSMInt:
    res="Sig_ai1_1_Re";
  case GGTplIntBSM_Re:
    res="Int_ai1_1_Re";
  default:
    break;
  };
  if (res!="") res = Form("T_%s_%s", TemplateHelpers::getGGProcessName(useOffshell).Data(), res.Data());
  return res;
}
TString TemplateHelpers::getMELAGGHypothesisWeight(TemplateHelpers::GGHypothesisType type){
  TString strWeight;
  switch (type){
  case GGBkg:
    strWeight = "p_Gen_GG_BKG_MCFM";
    break;
  case GGSig:
    strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";
    break;
  case GGBSI:
    strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM";
    break;
  default:
    break;
  };
  return strWeight;
}



TString TemplateHelpers::getQQBkgProcessName(bool useOffshell){ return (useOffshell ? "qqZZ" : "bkg_qqzz"); }
TString TemplateHelpers::getQQBkgOutputTreeName(bool useOffshell){
  TString res = Form("T_%s_Tree", TemplateHelpers::getQQBkgProcessName(useOffshell).Data());
  return res;
}

TString TemplateHelpers::getQQBkgTemplateName(bool useOffshell){
  TString res = Form("T_%s", TemplateHelpers::getQQBkgProcessName(useOffshell).Data());
  return res;
}
