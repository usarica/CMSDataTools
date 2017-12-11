#include "TemplateHelpers.h"
#include "TemplateHelpers.hpp"


/*********************/
/* General functions */
/*********************/
void TemplateHelpers::getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  using namespace SampleHelpers;
  using namespace DiscriminantClasses;
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
    KDspecs KDbkg("Dbkgkin");
    KDbkg.KD = constructKDFromType(kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
    KDbkg.KDvars = getKDVars(kDbkgkin);
    KDlist.push_back(KDbkg);
    KDspecs KDggint("Dggint");
    //KDggint.KD = constructKDFromType(kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dggbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varReco_Constant_Smooth");
    //KDggint.KDvars = getKDVars(kDggint);
    KDlist.push_back(KDggint);
    KDspecs KDL1("DL1dec");
    KDL1.KD = constructKDFromType(kDL1dec, "", "", "../data/gConstant_HZZ2e2mu_L1.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1", 1e-4);
    KDL1.KDvars = getKDVars(kDL1dec);
    KDlist.push_back(KDL1);
    KDspecs KDa2("Da2dec");
    KDa2.KD = constructKDFromType(kDa2dec, "", "", "../data/gConstant_HZZ2e2mu_g2.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");
    KDa2.KDvars = getKDVars(kDa2dec);
    KDlist.push_back(KDa2);
    KDspecs KDa3("Da3dec");
    KDa3.KD = constructKDFromType(kDa3dec, "", "", "../data/gConstant_HZZ2e2mu_g4.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");
    KDa3.KDvars = getKDVars(kDa3dec);
    KDlist.push_back(KDa3);
  }
}

/****************/
/* Gluon fusion */
/****************/
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
    break;
  case GGTplSigBSMSMInt_Re:
    res="Sig_ai1_1_Re";
    break;
  case GGTplIntBSM_Re:
    res="Int_ai1_1_Re";
    break;
  default:
    break;
  };
  if (res!="") res = Form("T_%s_%s", TemplateHelpers::getGGProcessName(useOffshell).Data(), res.Data());
  return res;
}
TString TemplateHelpers::getMELAGGHypothesisWeight(TemplateHelpers::GGHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo){
  TString strWeight;
  if (type==GGBkg) strWeight = "p_Gen_GG_BKG_MCFM";
  else if (type==GGSig) strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";
  else if (type==GGBSI) strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM";
  else if (type==GGSigBSM){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==GGSigBSMSMInt){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz4_1_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==GGBBI){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz4_1_MCFM";
      break;
    default:
      break;
    };
  }
  return strWeight;
}
std::vector<TemplateHelpers::GGHypothesisType> TemplateHelpers::getGGHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo){
  std::vector<TemplateHelpers::GGHypothesisType> res;
  // Order matters!
  if (hypo==ACHypothesisHelpers::kSM){
    res.push_back(TemplateHelpers::GGBkg);
    res.push_back(TemplateHelpers::GGSig);
    res.push_back(TemplateHelpers::GGBSI);
  }
  else{
    res.push_back(TemplateHelpers::GGSigBSM);
    res.push_back(TemplateHelpers::GGSigBSMSMInt);
    res.push_back(TemplateHelpers::GGBBI);
  }
  return res;
}
TString TemplateHelpers::getGGProcessLabel(TemplateHelpers::GGHypothesisType type, ACHypothesisHelpers::ACHypothesis hypo){
  TString acname;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    acname="f_{#Lambda1}";
    break;
  case ACHypothesisHelpers::kA2:
    acname="f_{a2}";
    break;
  case ACHypothesisHelpers::kA3:
    acname="f_{a3}";
    break;
  default:
    break;
  };
  switch (type){
  case GGBkg:
    return "gg #rightarrow 4l bkg.";
  case GGSig:
    return "gg #rightarrow 4l SM sig.";
  case GGBSI:
    return "gg #rightarrow 4l SM sig.+bkg.";
  case GGSigBSM:
    return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=1");
  case GGSigBSMSMInt:
    return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=0.5");
  case GGBBI:
    return Form("gg #rightarrow 4l %s%s sig.+bkg.", acname.Data(), "=1");
  default:
    return "";
  };
}
TString TemplateHelpers::getGGProcessLabel(TemplateHelpers::GGTemplateType type, ACHypothesisHelpers::ACHypothesis hypo){
  TString acname;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    acname="f_{#Lambda1}";
    break;
  case ACHypothesisHelpers::kA2:
    acname="f_{a2}";
    break;
  case ACHypothesisHelpers::kA3:
    acname="f_{a3}";
    break;
  default:
    break;
  };
  switch (type){
  case GGTplBkg:
    return "gg #rightarrow 4l bkg.";
  case GGTplSig:
    return "gg #rightarrow 4l SM sig.";
  case GGTplInt_Re:
    return "gg #rightarrow 4l SM sig.-bkg. interference";
  case GGTplSigBSM:
    return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=1");
  case GGTplSigBSMSMInt_Re:
    return Form("gg #rightarrow 4l %s%s interference", acname.Data(), "=0.5");
  case GGTplIntBSM_Re:
    return Form("gg #rightarrow 4l %s%s sig.-bkg. interference", acname.Data(), "=1");
  default:
    return "";
  };
}

int TemplateHelpers::castGGHypothesisTypeToInt(TemplateHelpers::GGHypothesisType type){ return (int)type; }
int TemplateHelpers::castGGTemplateTypeToInt(TemplateHelpers::GGTemplateType type){ return (int) type; }
TemplateHelpers::GGHypothesisType TemplateHelpers::castIntToGGHypothesisType(int type, bool useN){
  switch (type){
  case 0:
    return GGBkg;
  case 1:
    return GGSig;
  case 2:
    return GGBSI;
  case 3:
    return (!useN ? GGSigBSM : nGGSMTypes);
  case 4:
    return GGSigBSMSMInt;
  case 5:
    return GGBBI;
  default:
    return nGGTypes;
  };
}
TemplateHelpers::GGTemplateType TemplateHelpers::castIntToGGTemplateType(int type, bool useN){
  switch (type){
  case 0:
    return GGTplBkg;
  case 1:
    return GGTplSig;
  case 2:
    return GGTplInt_Re;
  case 3:
    return (!useN ? GGTplSigBSM : nGGTplSMTypes);
  case 4:
    return GGTplSigBSMSMInt_Re;
  case 5:
    return GGTplIntBSM_Re;
  default:
    return nGGTplTypes;
  };
}

template<> void TemplateHelpers::recombineGGHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo){
  if (vals.empty()) return;
  std::vector<float> res; 
  res.assign(vals.size(), 0);
  if (hypo==ACHypothesisHelpers::kSM){
    assert(vals.size()==nGGSMTypes);
    const float invA[nGGSMTypes][nGGSMTypes]={
      { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, -1, 1 }
    };
    for (int ix=0; ix<(int)nGGSMTypes; ix++){ for (int iy=0; iy<(int) nGGSMTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
  }
  else{
    assert(vals.size()==nGGTypes);
    const float couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
    const float couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
    const float cscale = couplA/couplM;
    const float cscalesq = pow(cscale, float(2));
    const float invA[nGGTypes][nGGTypes]={
      { 1, 0, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0, 0 },
      { -1, -1, 1, 0, 0, 0 },
      { 0, 0, 0, cscalesq, 0, 0 },
      { 0, -cscale, 0, -cscale, cscale, 0 },
      { -cscale, 0, 0, -cscale, 0, cscale }
    };
    for (int ix=0; ix<(int)nGGTypes; ix++){ for (int iy=0; iy<(int)nGGTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
  }
  std::swap(vals, res);
}
template<> void TemplateHelpers::recombineGGHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo){
  if (vals.empty()) return;
  typedef TH1F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  for (int ix=1; ix<=nx; ix++){
    std::vector<float> binvals; binvals.assign(vals.size(), 0);
    std::vector<float>::iterator ih=binvals.begin();
    for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix); ih++; }
    TemplateHelpers::recombineGGHistogramsToTemplates<float>(binvals, hypo);
    ih=binvals.begin();
    for (htype_t*& hh:vals){ hh->SetBinContent(ix, *ih); ih++; }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    GGTemplateType type = castIntToGGTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplInt_Re || type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
        if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
      }
      else{
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
      }
    }
    for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
    for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

    HelperFunctions::wipeOverUnderFlows(hh);
    HelperFunctions::divideBinWidth(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getGGProcessLabel(type, hypo));
    setTemplateAxisLabels(hh);
  }
}
template<> void TemplateHelpers::recombineGGHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo){
  if (vals.empty()) return;
  typedef TH2F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  int const ny = vals.at(0)->GetNbinsY();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      std::vector<float> binvals; binvals.assign(vals.size(), 0);
      std::vector<float>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix, iy); ih++; }
      TemplateHelpers::recombineGGHistogramsToTemplates<float>(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, *ih); ih++; }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    GGTemplateType type = castIntToGGTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplInt_Re || type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
        if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
        if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
      }
      else{
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
      }
    }
    for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
    for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

    HelperFunctions::wipeOverUnderFlows(hh);
    HelperFunctions::divideBinWidth(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getGGProcessLabel(type, hypo));
    setTemplateAxisLabels(hh);
  }
}
template<> void TemplateHelpers::recombineGGHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo){
  if (vals.empty()) return;
  typedef TH3F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  int const ny = vals.at(0)->GetNbinsY();
  int const nz = vals.at(0)->GetNbinsZ();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      for (int iz=1; iz<=nz; iz++){
        std::vector<float> binvals; binvals.assign(vals.size(), 0);
        std::vector<float>::iterator ih=binvals.begin();
        for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix, iy, iz); ih++; }
        TemplateHelpers::recombineGGHistogramsToTemplates<float>(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, iz, *ih); ih++; }
      }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    GGTemplateType type = castIntToGGTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplInt_Re || type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
        if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
        if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) asymAxes.push_back(2);
        if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
      }
      else{
        if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
        if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) symAxes.push_back(2);
      }
    }
    for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
    for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

    HelperFunctions::wipeOverUnderFlows(hh);
    HelperFunctions::divideBinWidth(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getGGProcessLabel(type, hypo));
    setTemplateAxisLabels(hh);
  }
}


/*****************/
/* QQ background */
/*****************/
TString TemplateHelpers::getQQBkgProcessName(bool useOffshell){ return (useOffshell ? "qqZZ" : "bkg_qqzz"); }
TString TemplateHelpers::getQQBkgOutputTreeName(bool useOffshell){
  TString res = Form("T_%s_Tree", TemplateHelpers::getQQBkgProcessName(useOffshell).Data());
  return res;
}
TString TemplateHelpers::getQQBkgTemplateName(bool useOffshell){
  TString res = Form("T_%s", TemplateHelpers::getQQBkgProcessName(useOffshell).Data());
  return res;
}
TString TemplateHelpers::getQQBkgProcessLabel(){ return "q#bar{q} #rightarrow 4l"; }

int TemplateHelpers::castQQBkgHypothesisTypeToInt(TemplateHelpers::QQBkgHypothesisType type){ return (int) type; }
int TemplateHelpers::castQQBkgTemplateTypeToInt(TemplateHelpers::QQBkgTemplateType type){ return (int) type; }
TemplateHelpers::QQBkgHypothesisType TemplateHelpers::castIntToQQBkgHypothesisType(int type){ return (QQBkgHypothesisType) type; }
TemplateHelpers::QQBkgTemplateType TemplateHelpers::castIntToQQBkgTemplateType(int type){ return (QQBkgTemplateType) type; }

