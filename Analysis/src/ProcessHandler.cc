#include <cassert>
#include "MELAStreamHelpers.hh"
#include "ProcessHandler.h"
#include "TemplateHelpers.h"
#include "TemplateHelpers.hpp"


using namespace std;
using namespace MELAStreamHelpers;


ProcessHandler::ProcessHandler(ProcessType proctype_, bool useOffshell_) : proctype(proctype_), useOffshell(useOffshell_) {
  assignProcessName();
}
void ProcessHandler::assignProcessName(){
  switch (proctype){
  case kGG:
    procname=(useOffshell ? "ggZZ_offshell" : "ggH");
    break;
  case kVV:
    procname= (useOffshell ? "VVZZ_offshell" : "VBF");
    break;
  case kQQBkg:
    procname= (useOffshell ? "qqZZ" : "bkg_qqzz");
    break;
  default:
    procname="";
    assert(0);
    break;
  };
}
const TString& ProcessHandler::getProcessName() const{ return procname; }
const ProcessHandler::ProcessType& ProcessHandler::getProcessType() const{ return proctype; }


/****************/
/* Gluon fusion */
/****************/
GGProcessHandler::GGProcessHandler(bool useOffshell_) : ProcessHandler(ProcessHandler::kGG, useOffshell_)
{}

TString GGProcessHandler::getOutputTreeName(GGProcessHandler::HypothesisType type) const{
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
  if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
  return res;
}
TString GGProcessHandler::getTemplateName(GGProcessHandler::TemplateType type) const{
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
  if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
  return res;
}
TString GGProcessHandler::getMELAHypothesisWeight(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
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
std::vector<GGProcessHandler::HypothesisType> GGProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
  std::vector<GGProcessHandler::HypothesisType> res;
  // Order matters!
  if (hypo==ACHypothesisHelpers::kSM){
    res.push_back(GGBkg);
    res.push_back(GGSig);
    res.push_back(GGBSI);
  }
  else{
    res.push_back(GGSigBSM);
    res.push_back(GGSigBSMSMInt);
    res.push_back(GGBBI);
  }
  return res;
}
TString GGProcessHandler::getProcessLabel(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
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
TString GGProcessHandler::getProcessLabel(GGProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
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

int GGProcessHandler::castHypothesisTypeToInt(GGProcessHandler::HypothesisType type){ return (int) type; }
int GGProcessHandler::castTemplateTypeToInt(GGProcessHandler::TemplateType type){ return (int) type; }
GGProcessHandler::HypothesisType GGProcessHandler::castIntToHypothesisType(int type, bool useN){
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
GGProcessHandler::TemplateType GGProcessHandler::castIntToTemplateType(int type, bool useN){
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

template<> void GGProcessHandler::recombineHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
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
    for (int ix=0; ix<(int) nGGSMTypes; ix++){ for (int iy=0; iy<(int) nGGSMTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
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
    for (int ix=0; ix<(int) nGGTypes; ix++){ for (int iy=0; iy<(int) nGGTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
  }
  std::swap(vals, res);
}
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
  if (vals.empty()) return;
  typedef TH1F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  for (int ix=1; ix<=nx; ix++){
    std::vector<float> binvals; binvals.assign(vals.size(), 0);
    std::vector<float>::iterator ih=binvals.begin();
    for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix); ih++; }
    GGProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
    ih=binvals.begin();
    for (htype_t*& hh:vals){ hh->SetBinContent(ix, *ih); ih++; }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
  if (vals.empty()) return;
  typedef TH2F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  int const ny = vals.at(0)->GetNbinsY();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      std::vector<float> binvals; binvals.assign(vals.size(), 0);
      std::vector<float>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix, iy); ih++; }
      GGProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, *ih); ih++; }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}
template<> void GGProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
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
        GGProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, iz, *ih); ih++; }
      }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}


/****************/
/* EW VV fusion */
/****************/
VVProcessHandler::VVProcessHandler(bool useOffshell_) : ProcessHandler(ProcessHandler::kVV, useOffshell_)
{}

TString VVProcessHandler::getOutputTreeName(VVProcessHandler::HypothesisType type) const{
  TString res;
  switch (type){
  case VVBkg:
    res="Bkg";
    break;
  case VVSig:
    res="Sig";
    break;
  case VVBSI:
    res="BSI";
    break;
  case VVSigBSM:
    res="SigBSM";
    break;
  case VVSigBSMSMInt0p25:
    res="SigBSMSMInt0p25";
    break;
  case VVSigBSMSMInt0p5:
    res="SigBSMSMInt0p5";
    break;
  case VVSigBSMSMInt0p75:
    res="SigBSMSMInt0p75";
    break;
  case VVBBI:
    res="BBI";
    break;
  case VVBMI:
    res="BMI";
    break;
  default:
    break;
  };
  if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
  return res;
}
TString VVProcessHandler::getTemplateName(VVProcessHandler::TemplateType type) const{
  TString res;
  switch (type){
  case VVTplBkg:
    res="Bkg";
    break;
  case VVTplSig:
    res="Sig";
    break;
  case VVTplInt_Re:
    res="Int_Re";
    break;
  case VVTplSigBSM:
    res="Sig_ai1_4";
    break;
  case VVTplSigBSMSMInt_ai1_1_Re:
    res="Sig_ai1_1_Re";
    break;
  case VVTplSigBSMSMInt_ai1_2_PosDef:
    res="Sig_ai1_2_PosDef";
    break;
  case VVTplSigBSMSMInt_ai1_3_Re:
    res="Sig_ai1_3_Re";
    break;
  case VVTplIntBSM_ai1_1_Re:
    res="Int_ai1_1_Re";
    break;
  case VVTplIntBSM_ai1_2_Re:
    res="Int_ai1_2_Re";
    break;
  default:
    break;
  };
  if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
  return res;
}
TString VVProcessHandler::getMELAHypothesisWeight(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
  TString strWeight;
  if (type==VVBkg) strWeight = "p_Gen_JJEW_BKG_MCFM";
  else if (type==VVSig) strWeight = "p_Gen_JJEW_SIG_ghv1_1_MCFM";
  else if (type==VVBSI) strWeight = "p_Gen_JJEW_BSI_ghv1_1_MCFM";
  else if (type==VVSigBSM){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_SIG_ghv1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_JJEW_SIG_ghv2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_SIG_ghv4_1_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==VVSigBSMSMInt0p25){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_25E2_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      //strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p25_MCFM";
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p251_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p25_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==VVSigBSMSMInt0p5){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_50E2_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p5_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p5_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==VVSigBSMSMInt0p75){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_75E2_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p75_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p75_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==VVBBI){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_BSI_ghv1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_JJEW_BSI_ghv2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_BSI_ghv4_1_MCFM";
      break;
    default:
      break;
    };
  }
  else if (type==VVBMI){
    switch (hypo){
    case ACHypothesisHelpers::kL1:
      strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv1prime2_1E4_MCFM";
      break;
    case ACHypothesisHelpers::kA2:
      strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv2_1_MCFM";
      break;
    case ACHypothesisHelpers::kA3:
      strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv4_1_MCFM";
      break;
    default:
      break;
    };
  }
  return strWeight;
}
std::vector<VVProcessHandler::HypothesisType> VVProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
  std::vector<HypothesisType> res;
  // Order matters!
  if (hypo==ACHypothesisHelpers::kSM){
    for (int i=0; i<castHypothesisTypeToInt(nVVSMTypes); i++){
      res.push_back(castIntToHypothesisType(i, false));
    }
  }
  else{
    for (int i=castHypothesisTypeToInt(nVVSMTypes); i<castHypothesisTypeToInt(nVVTypes); i++){
      res.push_back(castIntToHypothesisType(i, false));
    }
  }
  return res;
}
TString VVProcessHandler::getProcessLabel(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
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
  case VVBkg:
    return "VV #rightarrow 4l bkg.";
  case VVSig:
    return "VV #rightarrow 4l SM sig.";
  case VVBSI:
    return "VV #rightarrow 4l SM sig.+bkg.";
  case VVSigBSM:
    return Form("VV #rightarrow 4l %s%s sig.", acname.Data(), "=1");
  case VVSigBSMSMInt0p25:
    return Form("VV #rightarrow 4l %s%s sig.", acname.Data(), "=0.059");
  case VVSigBSMSMInt0p5:
    return Form("VV #rightarrow 4l %s%s sig.", acname.Data(), "=0.2");
  case VVSigBSMSMInt0p75:
    return Form("VV #rightarrow 4l %s%s sig.", acname.Data(), "=0.36");
  case VVBBI:
    return Form("VV #rightarrow 4l %s%s sig.+bkg.", acname.Data(), "=1");
  case VVBMI:
    return Form("VV #rightarrow 4l %s%s sig.+bkg.", acname.Data(), "=0.5");
  default:
    return "";
  };
}
TString VVProcessHandler::getProcessLabel(VVProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
  TString acname;
  switch (hypo){
  case ACHypothesisHelpers::kL1:
    acname="#Lambda_{1}";
    break;
  case ACHypothesisHelpers::kA2:
    acname="a_{2}";
    break;
  case ACHypothesisHelpers::kA3:
    acname="a_{3}";
    break;
  default:
    break;
  };
  switch (type){
  case VVTplBkg:
    return "VV #rightarrow 4l bkg.";
  case VVTplSig:
    return "VV #rightarrow 4l SM sig.";
  case VVTplInt_Re:
    return "VV #rightarrow 4l SM sig.-bkg. interference";
  case VVTplSigBSM:
    return Form("VV #rightarrow 4l %s sig.", acname.Data());
  case VVTplSigBSMSMInt_ai1_1_Re:
    return Form("VV #rightarrow 4l %s%s interference", acname.Data(), "^{1}");
  case VVTplSigBSMSMInt_ai1_2_PosDef:
    return Form("VV #rightarrow 4l %s%s interference", acname.Data(), "^{2}");
  case VVTplSigBSMSMInt_ai1_3_Re:
    return Form("VV #rightarrow 4l %s%s interference", acname.Data(), "^{3}");
  case VVTplIntBSM_ai1_1_Re:
    return Form("VV #rightarrow 4l %s%s sig.-bkg. interference", acname.Data(), "^{1}");
  case VVTplIntBSM_ai1_2_Re:
    return Form("VV #rightarrow 4l %s%s sig.-bkg. interference", acname.Data(), "^{2}");
  default:
    return "";
  };
}

int VVProcessHandler::castHypothesisTypeToInt(VVProcessHandler::HypothesisType type){ return (int) type; }
int VVProcessHandler::castTemplateTypeToInt(VVProcessHandler::TemplateType type){ return (int) type; }
VVProcessHandler::HypothesisType VVProcessHandler::castIntToHypothesisType(int type, bool useN){
  switch (type){
  case 0:
    return VVBkg;
  case 1:
    return VVSig;
  case 2:
    return VVBSI;
  case 3:
    return (!useN ? VVSigBSM : nVVSMTypes);
  case 4:
    return VVSigBSMSMInt0p25;
  case 5:
    return VVSigBSMSMInt0p5;
  case 6:
    return VVSigBSMSMInt0p75;
  case 7:
    return VVBBI;
  case 8:
    return VVBMI;
  default:
    return nVVTypes;
  };
}
VVProcessHandler::TemplateType VVProcessHandler::castIntToTemplateType(int type, bool useN){
  switch (type){
  case 0:
    return VVTplBkg;
  case 1:
    return VVTplSig;
  case 2:
    return VVTplInt_Re;
  case 3:
    return (!useN ? VVTplSigBSM : nVVTplSMTypes);
  case 4:
    return VVTplSigBSMSMInt_ai1_1_Re;
  case 5:
    return VVTplSigBSMSMInt_ai1_2_PosDef;
  case 6:
    return VVTplSigBSMSMInt_ai1_3_Re;
  case 7:
    return VVTplIntBSM_ai1_1_Re;
  case 8:
    return VVTplIntBSM_ai1_2_Re;
  default:
    return nVVTplTypes;
  };
}

template<> void VVProcessHandler::recombineHistogramsToTemplates<float>(std::vector<float>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
  if (vals.empty()) return;
  std::vector<float> res;
  res.assign(vals.size(), 0);
  if (hypo==ACHypothesisHelpers::kSM){
    assert(vals.size()==nVVSMTypes);
    const float invA[nVVSMTypes][nVVSMTypes]={
      { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, -1, 1 }
    };
    for (int ix=0; ix<(int) nVVSMTypes; ix++){ for (int iy=0; iy<(int) nVVSMTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
  }
  else{
    assert(vals.size()==nVVTypes);
    const float couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
    const float couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
    const float c = couplA/couplM;
    const float c2 = pow(c, 2);
    const float c3 = pow(c, 3);
    const float c4 = pow(c, 4);
    const float invA[nVVTypes][nVVTypes]={
      { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      { -1, -1, 1, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, c4, 0, 0, 0, 0, 0 },
      { 0, float(-22./3.)*c, 0, float(-3./32.)*c, float(12.)*c, float(-6.)*c, float(4./3.)*c, 0, 0 },
      { 0, float(16.)*c2, 0, float(11./16.)*c2, float(-40.)*c2, float(32.)*c2, float(-8.)*c2, 0, 0 },
      { 0, float(-32./3.)*c3, 0, float(-3./2.)*c3, float(32.)*c3, float(-32.)*c3, float(32./3.)*c3, 0, 0 },
      { c, float(2.)*c, -c, float(29./32.)*c, float(-4.)*c, float(6.)*c, float(-4.)*c, -c, c },
      { -c2, 0, 0, -c2, 0, 0, 0, c2, 0 }
    };
    for (int ix=0; ix<(int) nVVTypes; ix++){ for (int iy=0; iy<(int) nVVTypes; iy++) res.at(ix) += invA[ix][iy]*vals.at(iy); }
  }
  std::swap(vals, res);
}
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
  if (vals.empty()) return;
  typedef TH1F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  for (int ix=1; ix<=nx; ix++){
    std::vector<float> binvals; binvals.assign(vals.size(), 0);
    std::vector<float>::iterator ih=binvals.begin();
    for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix); ih++; }
    VVProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
    ih=binvals.begin();
    for (htype_t*& hh:vals){ hh->SetBinContent(ix, *ih); ih++; }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
  if (vals.empty()) return;
  typedef TH2F htype_t;
  int const nx = vals.at(0)->GetNbinsX();
  int const ny = vals.at(0)->GetNbinsY();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      std::vector<float> binvals; binvals.assign(vals.size(), 0);
      std::vector<float>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){ *ih=hh->GetBinContent(ix, iy); ih++; }
      VVProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, *ih); ih++; }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}
template<> void VVProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
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
        VVProcessHandler::recombineHistogramsToTemplates<float>(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){ hh->SetBinContent(ix, iy, iz, *ih); ih++; }
      }
    }
  }
  for (int t=0; t<(int) vals.size(); t++){
    htype_t*& hh=vals.at(t);
    TemplateType type = castIntToTemplateType(t);
    std::vector<unsigned int> symAxes;
    std::vector<unsigned int> asymAxes;
    if (hypo==ACHypothesisHelpers::kA3){
      if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
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
    hh->SetTitle(getProcessLabel(type, hypo));
    TemplateHelpers::setTemplateAxisLabels(hh);
  }
}


/*****************/
/* QQ background */
/*****************/
QQBkgProcessHandler::QQBkgProcessHandler(bool useOffshell_) : ProcessHandler(ProcessHandler::kQQBkg, useOffshell_)
{}

TString QQBkgProcessHandler::getOutputTreeName() const{
  TString res = Form("T_%s_Tree", getProcessName().Data());
  return res;
}
TString QQBkgProcessHandler::getTemplateName() const{
  TString res = Form("T_%s", getProcessName().Data());
  return res;
}
TString QQBkgProcessHandler::getProcessLabel() const{ return "q#bar{q} #rightarrow 4l"; }
TString QQBkgProcessHandler::getMELAHypothesisWeight(unsigned int njets) const{
  switch (njets){
  case 2:
    return "p_Gen_JJQCD_BKG_MCFM";
  default:
    return "p_Gen_QQB_BKG_MCFM";
  }
}
int QQBkgProcessHandler::castHypothesisTypeToInt(QQBkgProcessHandler::HypothesisType type){ return (int) type; }
int QQBkgProcessHandler::castTemplateTypeToInt(QQBkgProcessHandler::TemplateType type){ return (int) type; }
QQBkgProcessHandler::HypothesisType QQBkgProcessHandler::castIntToHypothesisType(int type){ return (HypothesisType) type; }
QQBkgProcessHandler::TemplateType QQBkgProcessHandler::castIntToTemplateType(int type){ return (TemplateType) type; }

template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals) const{
  typedef TH1F T;
  if ((int) vals.size()!=castHypothesisTypeToInt(nQQBkgTypes)) return;
  for (T*& hh:vals){
    HelperFunctions::wipeOverUnderFlows<T>(hh);
    HelperFunctions::divideBinWidth<T>(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getProcessLabel());
    TemplateHelpers::setTemplateAxisLabels<T>(hh);
  }
}
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals) const{
  typedef TH2F T;
  if ((int) vals.size()!=castHypothesisTypeToInt(nQQBkgTypes)) return;
  for (T*& hh:vals){
    HelperFunctions::wipeOverUnderFlows<T>(hh);
    HelperFunctions::divideBinWidth<T>(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getProcessLabel());
    TemplateHelpers::setTemplateAxisLabels<T>(hh);
  }
}
template<> void QQBkgProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals) const{
  typedef TH3F T;
  if ((int) vals.size()!=castHypothesisTypeToInt(nQQBkgTypes)) return;
  for (T*& hh:vals){
    HelperFunctions::wipeOverUnderFlows<T>(hh);
    HelperFunctions::divideBinWidth<T>(hh);
    hh->Scale(xsecScale);
    hh->SetTitle(getProcessLabel());
    TemplateHelpers::setTemplateAxisLabels<T>(hh);
  }
}