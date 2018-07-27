#include <iterator>
#include "SystematicsHelpers.h"
#include "TemplateHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace SystematicsHelpers;


SystematicsHelpers::YieldSystematic::YieldSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&)) : rule(infcn), evaluators(inEvaluators) {}
float SystematicsHelpers::YieldSystematic::eval(CJLSTTree* theTree) const{
  float weight=0;
  if (rule) weight=rule(theTree, evaluators);
  if (!HelperFunctions::checkVarNanInf(weight)){
    MELAerr << "SystematicsHelpers::YieldSystematic::eval(" << theTree->sampleIdentifier << ") weight (" << weight << ") is nan/inf!" << endl;
    weight=0;
  }
  return weight;
}

float SystematicsHelpers::getRawSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  float weight=1;
  for (auto const& builder:builders) weight *= builder->getPostThresholdWeight(theTree);
  return weight;
}
float SystematicsHelpers::getNormalizedSystematic(CJLSTTree* theTree, const std::vector<ReweightingBuilder*>& builders){
  if (builders.size()==1) return SystematicsHelpers::getRawSystematic(theTree, builders);

  float weight=0;
  assert(builders.size()==2 && builders.at(0) && builders.at(1));
  float sumWeights = builders.at(0)->getSumAllPostThresholdWeights(theTree);
  if (sumWeights!=0.) weight = builders.at(0)->getPostThresholdWeight(theTree)*builders.at(1)->getSumAllPostThresholdWeights(theTree)/sumWeights;
  return weight;
}


SystematicsHelpers::PerLeptonSystematic::PerLeptonSystematic(
  const std::vector<TString>& inStrVars,
  SystematicsHelpers::PerLeptonSystematic::PerLeptonSystematicFunction_t infcn,
  unsigned int const id_requested_,
  bool doUp_
) : rule(infcn), strVars(inStrVars), id_requested(id_requested_), doUp(doUp_)
{}
float SystematicsHelpers::PerLeptonSystematic::eval(CJLSTTree* theTree) const{
  auto it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "PerLeptonSystematic::eval: Could not find the variables for tree " << theTree->sampleIdentifier << ". Call PerLeptonSystematic::setup first!" << endl;
    return 0;
  }
  std::pair<float, float> res = rule(*(it->second.first.at(0)), *(it->second.first.at(1)), it->second.second, id_requested);
  return (doUp ? res.second : res.first);
}
void SystematicsHelpers::PerLeptonSystematic::setup(CJLSTTree* theTree){
  if (!theTree || strVars.empty()) return;
  std::vector<short*> refId; refId.assign(2, nullptr);
  theTree->getValRef("Z1Flav", refId.at(0));
  theTree->getValRef("Z2Flav", refId.at(1));
  std::vector<std::vector<float>* const*> refVar;
  for (TString const& s:strVars){
    vector<float>* const* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (!refVar.empty()) componentRefs[theTree]=std::pair<std::vector<short*>, std::vector<std::vector<float>* const*>>(refId, refVar);
}

std::pair<float, float> SystematicsHelpers::getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq){
  std::pair<float, float> res(1, 1);
  if ((Z1Flav*Z2Flav)%(short)idreq != 0) return res;
  assert(LepVars.size()==5);

  std::vector<float>* const& LepPt = *(LepVars.at(0)); assert(LepPt);
  std::vector<float>* const& LepRecoSF = *(LepVars.at(1)); assert(LepRecoSF);
  std::vector<float>* const& LepRecoSF_Unc = *(LepVars.at(2)); assert(LepRecoSF_Unc);
  std::vector<float>* const& LepSelSF = *(LepVars.at(3)); assert(LepSelSF);
  std::vector<float>* const& LepSelSF_Unc = *(LepVars.at(4)); assert(LepSelSF_Unc);

  unsigned int const nleps = LepPt->size(); assert(nleps>1);
  float const& lep4pt = LepPt->back();
  for (unsigned int ilep=0; ilep<nleps; ilep++){
    short absZflav;
    if (ilep<2) absZflav=std::abs(Z1Flav);
    else absZflav=std::abs(Z2Flav);
    if (absZflav%(short) idreq==0){
      res.first *= (1.-LepRecoSF_Unc->at(ilep)/LepRecoSF->at(ilep))*(1.-LepSelSF_Unc->at(ilep)/LepSelSF->at(ilep));
      res.second *= (1.+LepRecoSF_Unc->at(ilep)/LepRecoSF->at(ilep))*(1.+LepSelSF_Unc->at(ilep)/LepSelSF->at(ilep));
    }
  }
  if (std::abs(Z1Flav)==169 && std::abs(Z2Flav)==169){
    if (lep4pt<7.){
      res.first *= 1.-0.032;
      res.second *= 1.+0.001;
    }
    else{
      res.first *= 1.-0.015;
      res.second *= 1.+0.001;
    }
  }
  else if (std::abs(Z1Flav)==121 && std::abs(Z2Flav)==121){
    if (lep4pt<12.){
      res.first *= 1.-0.011;
      res.second *= 1.+0.001;
    }
    else{
      res.first *= 1.-0.010;
      res.second *= 1.+0.005;
    }
  }
  else{
    if (lep4pt<7.){
      res.first *= 1.-0.080;
      res.second *= 1.+0.005;
    }
    else if (lep4pt<12.){
      res.first *= 1.-0.032;
      res.second *= 1.+0.005;
    }
    else{
      res.first *= 1.-0.010;
      res.second *= 1.+0.001;
    }
  }
  return res;
}


SystematicsHelpers::PerLeptonScaleSystematic::PerLeptonScaleSystematic(
  const std::vector<TString>& inStrVars,
  SystematicsHelpers::PerLeptonScaleSystematic::PerLeptonScaleSystematicFunction_t infcn,
  unsigned int const id_requested_,
  bool doUp_
) : rule(infcn), strVars(inStrVars), id_requested(id_requested_), doUp(doUp_)
{}
float SystematicsHelpers::PerLeptonScaleSystematic::eval(CJLSTTree* theTree) const{
  auto it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "PerLeptonScaleSystematic::eval: Could not find the variables for tree " << theTree->sampleIdentifier << ". Call PerLeptonScaleSystematic::setup first!" << endl;
    return 0;
  }
  std::pair<float, float> res = rule(*(it->second.ref_shorts.at(0)), *(it->second.ref_shorts.at(1)), *(it->second.ref_floats.at(0)), it->second.ref_vfloats, id_requested);
  return (doUp ? res.second : res.first);
}
void SystematicsHelpers::PerLeptonScaleSystematic::setup(CJLSTTree* theTree){
  if (!theTree || strVars.empty()) return;
  std::vector<short*> refId; refId.assign(2, nullptr);
  std::vector<float*> refMass; refMass.assign(1, nullptr);
  theTree->getValRef("Z1Flav", refId.at(0));
  theTree->getValRef("Z2Flav", refId.at(1));
  theTree->getValRef("ZZMass", refMass.at(0));
  std::vector<std::vector<float>* const*> refVar;
  for (TString const& s:strVars){
    vector<float>* const* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonScaleSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (!refVar.empty()) componentRefs[theTree]=componentData(refId, refMass, refVar);
}

std::pair<float, float> SystematicsHelpers::getLeptonScaleSystematic(short const& Z1Flav, short const& Z2Flav, float const& ZZMass, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq){
  std::pair<float, float> res(ZZMass, ZZMass);
  if ((Z1Flav*Z2Flav)%(short) idreq != 0) return res;
  // LepVars: LepPt, LepEta, LepPhi, LepScaleUnc, LepResUnc
  assert(LepVars.size()==5);
  std::vector<float>* const& LepPt = *(LepVars.at(0)); assert(LepPt);
  std::vector<float>* const& LepEta = *(LepVars.at(1)); assert(LepEta);
  std::vector<float>* const& LepPhi = *(LepVars.at(2)); assert(LepPhi);
  std::vector<float>* const& LepScaleUnc = *(LepVars.at(3)); assert(LepScaleUnc);
  std::vector<float>* const& LepResUnc = *(LepVars.at(4)); assert(LepResUnc);
  unsigned int const nleps = LepPt->size(); assert(nleps>1);
  unsigned int const nperms = pow(3, nleps);
  std::vector<std::vector<TLorentzVector>> LepP; LepP.assign(nleps, std::vector<TLorentzVector>()); for (auto& ls:LepP) ls.assign(3, TLorentzVector(0, 0, 0, 0));
  for (unsigned int ilep=0; ilep<nleps; ilep++){
    for (int isyst=0; isyst<3; isyst++){
      TLorentzVector& pL=LepP.at(ilep).at(isyst);
      pL.SetPtEtaPhiM(LepPt->at(ilep), LepEta->at(ilep), LepPhi->at(ilep), 0);
      short absZflav;
      if (ilep<2) absZflav=std::abs(Z1Flav);
      else absZflav=std::abs(Z2Flav);
      if (absZflav%(short) idreq==0 && isyst>0) pL = pL*pow(LepScaleUnc->at(ilep)*LepResUnc->at(ilep), 2*isyst-3);
    }
  }
  float scale=1;
  std::vector<float> ZZMassVars;
  std::vector<std::pair<float, int>> ZZMassDiffVars_neg; ZZMassDiffVars_neg.reserve(nperms);
  std::vector<std::pair<float, int>> ZZMassDiffVars_pos; ZZMassDiffVars_pos.reserve(nperms);
  for (unsigned int iperm=0; iperm<nperms; iperm++){
    std::vector<unsigned int> isyst; isyst.reserve(nleps);
    unsigned int jperm=iperm;
    for (unsigned int ilep=0; ilep<nleps; ilep++){ isyst.push_back(jperm%3); jperm/=3; }
    TLorentzVector sumP(0, 0, 0, 0);
    for (unsigned int ilep=0; ilep<nleps; ilep++) sumP += LepP.at(ilep).at(isyst.at(ilep));
    float mass=sumP.M();
    if (iperm==0) scale = ZZMass/mass; // Account for FSR and missing lepton masses
    mass *= scale;

    float massdiff=mass-ZZMass;
    if (massdiff<0.) HelperFunctions::addByLowest(ZZMassDiffVars_neg, std::pair<float, int>(fabs(massdiff), iperm), true);
    else HelperFunctions::addByLowest(ZZMassDiffVars_pos, std::pair<float, int>(fabs(massdiff), iperm), true);
    ZZMassVars.push_back(mass);
  }

  if (!ZZMassDiffVars_neg.empty()){
    unsigned int index=ZZMassDiffVars_neg.size()*0.684;
    res.first = ZZMassVars.at(ZZMassDiffVars_neg.at(index).second);
  }
  if (!ZZMassDiffVars_pos.empty()){
    unsigned int index=ZZMassDiffVars_pos.size()*0.684;
    res.second = ZZMassVars.at(ZZMassDiffVars_pos.at(index).second);
  }
  return res;
}

SystematicsHelpers::PerLeptonResSystematic::PerLeptonResSystematic(
  const std::vector<TString>& inStrVars,
  SystematicsHelpers::PerLeptonResSystematic::PerLeptonResSystematicFunction_t infcn,
  unsigned int const id_requested_,
  bool doUp_
) : rule(infcn), strVars(inStrVars), id_requested(id_requested_), doUp(doUp_), centralValue(-999.)
{}
float SystematicsHelpers::PerLeptonResSystematic::eval(CJLSTTree* theTree) const{
  auto it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "PerLeptonResSystematic::eval: Could not find the variables for tree " << theTree->sampleIdentifier << ". Call PerLeptonResSystematic::setup first!" << endl;
    return 0;
  }
  if (centralValue==-999.) MELAerr << "PerLeptonResSystematic::eval: Reference value is " << centralValue << endl;
  std::pair<float, float> res = rule(*(it->second.ref_shorts.at(0)), *(it->second.ref_shorts.at(1)), *(it->second.ref_floats.at(0)), centralValue, it->second.ref_vfloats, id_requested);
  return (doUp ? res.second : res.first);
}
void SystematicsHelpers::PerLeptonResSystematic::setup(CJLSTTree* theTree){
  if (!theTree || strVars.empty()) return;
  std::vector<short*> refId; refId.assign(2, nullptr);
  std::vector<float*> refMass; refMass.assign(1, nullptr);
  theTree->getValRef("Z1Flav", refId.at(0));
  theTree->getValRef("Z2Flav", refId.at(1));
  theTree->getValRef("ZZMass", refMass.at(0));
  std::vector<std::vector<float>* const*> refVar;
  for (TString const& s:strVars){
    vector<float>* const* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonResSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (!refVar.empty()) componentRefs[theTree]=componentData(refId, refMass, refVar);
}

std::pair<float, float> SystematicsHelpers::getLeptonResSystematic(short const& Z1Flav, short const& Z2Flav, float const& ZZMass, float const& centralValue, std::vector<std::vector<float>* const*> const& LepVars, unsigned int const idreq){
  std::pair<float, float> res(ZZMass, ZZMass);
  if ((Z1Flav*Z2Flav)%(short) idreq != 0) return res;
  // LepVars: LepPt, LepEta, LepPhi, LepScaleUnc, LepResUnc
  assert(LepVars.size()==5);
  std::vector<float>* const& LepPt = *(LepVars.at(0)); assert(LepPt);
  std::vector<float>* const& LepEta = *(LepVars.at(1)); assert(LepEta);
  std::vector<float>* const& LepPhi = *(LepVars.at(2)); assert(LepPhi);
  std::vector<float>* const& LepScaleUnc = *(LepVars.at(3)); assert(LepScaleUnc);
  std::vector<float>* const& LepResUnc = *(LepVars.at(4)); assert(LepResUnc);
  unsigned int const nleps = LepPt->size(); assert(nleps>1);
  unsigned int const nperms = pow(3, nleps);
  std::vector<std::vector<TLorentzVector>> LepP; LepP.assign(nleps, std::vector<TLorentzVector>()); for (auto& ls:LepP) ls.assign(3, TLorentzVector(0, 0, 0, 0));
  for (unsigned int ilep=0; ilep<nleps; ilep++){
    for (int isyst=0; isyst<3; isyst++){
      TLorentzVector& pL=LepP.at(ilep).at(isyst);
      pL.SetPtEtaPhiM(LepPt->at(ilep), LepEta->at(ilep), LepPhi->at(ilep), 0);
      short absZflav;
      if (ilep<2) absZflav=std::abs(Z1Flav);
      else absZflav=std::abs(Z2Flav);
      if (absZflav%(short) idreq==0 && isyst>0) pL = pL*pow(LepScaleUnc->at(ilep)*LepResUnc->at(ilep), 2*isyst-3);
    }
  }
  float scale=1;
  float nominaldiff=fabs(ZZMass-centralValue);
  std::vector<float> ZZMassVars; ZZMassVars.reserve(nperms);
  std::vector<std::pair<float, int>> ZZMassDiffVars_neg; ZZMassDiffVars_neg.reserve(nperms);
  std::vector<std::pair<float, int>> ZZMassDiffVars_pos; ZZMassDiffVars_pos.reserve(nperms);
  for (unsigned int iperm=0; iperm<nperms; iperm++){
    std::vector<unsigned int> isyst; isyst.reserve(nleps);
    unsigned int jperm=iperm;
    for (unsigned int ilep=0; ilep<nleps; ilep++){ isyst.push_back(jperm%3); jperm/=3; }
    TLorentzVector sumP(0, 0, 0, 0);
    for (unsigned int ilep=0; ilep<nleps; ilep++) sumP += LepP.at(ilep).at(isyst.at(ilep));
    float mass=sumP.M();
    if (iperm==0) scale = ZZMass/mass; // Account for FSR and missing lepton masses
    mass *= scale;
    float massdiff=fabs(mass-centralValue)-nominaldiff;
    if (massdiff<0.) HelperFunctions::addByLowest(ZZMassDiffVars_neg, std::pair<float, int>(fabs(massdiff), iperm), true);
    else HelperFunctions::addByLowest(ZZMassDiffVars_pos, std::pair<float, int>(fabs(massdiff), iperm), true);
    ZZMassVars.push_back(mass);
  }

  if (!ZZMassDiffVars_neg.empty()){
    unsigned int index=ZZMassDiffVars_neg.size()*0.684;
    res.first = ZZMassVars.at(ZZMassDiffVars_neg.at(index).second);
  }
  if (!ZZMassDiffVars_pos.empty()){
    unsigned int index=ZZMassDiffVars_pos.size()*0.684;
    res.second = ZZMassVars.at(ZZMassDiffVars_pos.at(index).second);
  }
  return res;
}


int SystematicsHelpers::convertSystematicVariationTypeToInt(SystematicsHelpers::SystematicVariationTypes type){ return (int) type; }
std::vector<SystematicsHelpers::SystematicVariationTypes> SystematicsHelpers::getProcessSystematicVariations(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  TString strGenerator
){
  std::vector<SystematicVariationTypes> res;

  res.push_back(sNominal);

  if (proc==ProcessHandler::kZX){
    res.push_back(eZXStatsDn);
    res.push_back(eZXStatsUp);
  }
  else{
    if (channel==SampleHelpers::k4e || channel==SampleHelpers::k2e2mu){
      res.push_back(eLepSFEleDn);
      res.push_back(eLepSFEleUp);
    }
    if (channel==SampleHelpers::k4mu || channel==SampleHelpers::k2e2mu){
      res.push_back(eLepSFMuDn);
      res.push_back(eLepSFMuUp);
    }
    if (
      (proc==ProcessHandler::kGG || proc==ProcessHandler::kVV || proc==ProcessHandler::kVBF || proc==ProcessHandler::kZH || proc==ProcessHandler::kWH || proc==ProcessHandler::kTT || proc==ProcessHandler::kBB)
      &&
      (strGenerator=="POWHEG" || strGenerator=="")
      ){
      if (channel==SampleHelpers::k4e || channel==SampleHelpers::k2e2mu){
        res.push_back(eLepScaleEleDn);
        res.push_back(eLepScaleEleUp);
        res.push_back(eLepResEleDn);
        res.push_back(eLepResEleUp);
      }
      if (channel==SampleHelpers::k4mu || channel==SampleHelpers::k2e2mu){
        res.push_back(eLepScaleMuDn);
        res.push_back(eLepScaleMuUp);
        res.push_back(eLepResMuDn);
        res.push_back(eLepResMuUp);
      }
    }

    res.push_back(tPDFScaleDn);
    res.push_back(tPDFScaleUp);
    res.push_back(tQCDScaleDn);
    res.push_back(tQCDScaleUp);
    res.push_back(tAsMZDn);
    res.push_back(tAsMZUp);
    res.push_back(tPDFReplicaDn);
    res.push_back(tPDFReplicaUp);

    if (proc==ProcessHandler::kQQBkg){
      res.push_back(tQQBkgEWCorrDn);
      res.push_back(tQQBkgEWCorrUp);
    }

    if (category!=CategorizationHelpers::Inclusive){
      if ((proc==ProcessHandler::kGG || proc==ProcessHandler::kVV || proc==ProcessHandler::kVBF || proc==ProcessHandler::kZH || proc==ProcessHandler::kWH || proc==ProcessHandler::kTT/* || proc==ProcessHandler::kBB*/) && (strGenerator=="POWHEG" || strGenerator=="")){
        res.push_back(tPythiaScaleDn);
        res.push_back(tPythiaScaleUp);
        res.push_back(tPythiaTuneDn);
        res.push_back(tPythiaTuneUp);
      }
      if (proc==ProcessHandler::kGG && (strGenerator=="POWHEG" || strGenerator=="")){
        res.push_back(tMINLODn);
        res.push_back(tMINLOUp);
      }
      if (
        CategorizationHelpers::globalCategorizationScheme==CategorizationHelpers::UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag
        &&
        (proc==ProcessHandler::kTT || proc==ProcessHandler::kBB) && (strGenerator=="POWHEG" || strGenerator=="")
        ){
        res.push_back(eBTagSFDn);
        res.push_back(eBTagSFUp);
      }
      res.push_back(eJECDn);
      res.push_back(eJECUp);
    }
  }
  return res;
}
bool SystematicsHelpers::systematicAllowed(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst,
  TString strGenerator
){
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedTypes = SystematicsHelpers::getProcessSystematicVariations(category, channel, proc, strGenerator);
  for (SystematicVariationTypes& st:allowedTypes){ if (st==syst) return true; }
  return false;
}
bool SystematicsHelpers::systematicHasMassRatio(SystematicsHelpers::SystematicVariationTypes const syst){
  return !(
    syst==sNominal
    || syst==eLepScaleEleDn || syst==eLepScaleEleUp
    || syst==eLepResEleDn || syst==eLepResEleUp
    || syst==eLepScaleMuDn || syst==eLepScaleMuUp
    || syst==eLepResMuDn || syst==eLepResMuUp
    );
}

SystematicsHelpers::SystematicsClass* SystematicsHelpers::constructSystematic(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst,
  std::vector<CJLSTTree*> trees,
  std::vector<ReweightingBuilder*>& extraEvaluators,
  TString strGenerator
){
  SystematicsClass* res=nullptr;
  if (!systematicAllowed(category, channel, proc, syst, strGenerator)) return res;

  ExtendedBinning binning("GenHMass");
  {
    ExtendedBinning binning_offshell = TemplateHelpers::getDiscriminantFineBinning(channel, CategorizationHelpers::Inclusive, ACHypothesisHelpers::kSM, "ZZMass", CategorizationHelpers::kOffshell);
    ExtendedBinning binning_onshell = TemplateHelpers::getDiscriminantCoarseBinning(channel, CategorizationHelpers::Inclusive, ACHypothesisHelpers::kSM, "ZZMass", CategorizationHelpers::kOnshell);
    for (double const& bb:binning_offshell.getBinningVector()) binning.addBinBoundary(bb);
    binning.addBinBoundary((binning_onshell.getMax() + binning_offshell.getMin())/2.);
    binning.addBinBoundary(70.);
    binning.addBinBoundary(85.);
    binning.addBinBoundary(100.);
    binning.addBinBoundary(110.);
    binning.addBinBoundary(120.);
    binning.addBinBoundary(140.);
  }

  ReweightingBuilder* rewgtbuilder=nullptr;
  ReweightingBuilder* normbuilder=nullptr;
  std::vector<TString> strVars, strVarsNorm;
  std::vector<ReweightingBuilder*> evaluators;
  ReweightingFunctions::ReweightingFunction_t computeFcn=nullptr;
  if (syst==eLepSFEleDn || syst==eLepSFEleUp || syst==eLepSFMuDn || syst==eLepSFMuUp){
    strVars.reserve(5);
    strVars.push_back("LepPt");
    strVars.push_back("LepRecoSF");
    strVars.push_back("LepRecoSF_Unc");
    strVars.push_back("LepSelSF");
    strVars.push_back("LepSelSF_Unc");
    for (CJLSTTree*& tree:trees){
      tree->bookBranch<short>("Z1Flav", 0);
      tree->bookBranch<short>("Z2Flav", 0);
      for (TString const& s:strVars) tree->bookBranch<vector<float>*>(s, nullptr);
    }
    res = new PerLeptonSystematic(
      strVars,
      SystematicsHelpers::getLeptonSFSystematic,
      (syst==eLepSFEleDn || syst==eLepSFEleUp ? 11 : 13),
      (syst==eLepSFEleUp || syst==eLepSFMuUp)
    );
    for (CJLSTTree*& tree:trees) ((PerLeptonSystematic*) res)->setup(tree);
  }
  else if (syst==eLepScaleEleDn || syst==eLepScaleEleUp || syst==eLepScaleMuDn || syst==eLepScaleMuUp){
    strVars.reserve(4);
    strVars.push_back("LepPt");
    strVars.push_back("LepEta");
    strVars.push_back("LepPhi");
    strVars.push_back("LepScale_Unc");
    strVars.push_back("LepSmear_Unc");
    for (CJLSTTree*& tree:trees){
      tree->bookBranch<short>("Z1Flav", 0);
      tree->bookBranch<short>("Z2Flav", 0);
      tree->bookBranch<float>("ZZMass", 0);
      for (TString const& s:strVars) tree->bookBranch<vector<float>*>(s, nullptr);
    }
    res = new PerLeptonScaleSystematic(
      strVars,
      SystematicsHelpers::getLeptonScaleSystematic,
      (syst==eLepScaleEleDn || syst==eLepScaleEleUp ? 11 : 13),
      (syst==eLepScaleEleUp || syst==eLepScaleMuUp)
    );
    for (CJLSTTree*& tree:trees) ((PerLeptonScaleSystematic*) res)->setup(tree);
  }
  else if (syst==eLepResEleDn || syst==eLepResEleUp || syst==eLepResMuDn || syst==eLepResMuUp){
    strVars.reserve(4);
    strVars.push_back("LepPt");
    strVars.push_back("LepEta");
    strVars.push_back("LepPhi");
    strVars.push_back("LepScale_Unc");
    strVars.push_back("LepSmear_Unc");
    for (CJLSTTree*& tree:trees){
      tree->bookBranch<short>("Z1Flav", 0);
      tree->bookBranch<short>("Z2Flav", 0);
      tree->bookBranch<float>("ZZMass", 0);
      for (TString const& s:strVars) tree->bookBranch<vector<float>*>(s, nullptr);
    }
    res = new PerLeptonResSystematic(
      strVars,
      SystematicsHelpers::getLeptonResSystematic,
      (syst==eLepResEleDn || syst==eLepResEleUp ? 11 : 13),
      (syst==eLepResEleUp || syst==eLepResMuUp)
    );
    for (CJLSTTree*& tree:trees) ((PerLeptonResSystematic*) res)->setup(tree);
  }
  else if (syst==tQQBkgEWCorrDn || syst==tQQBkgEWCorrUp){
    strVars.push_back("KFactor_EW_qqZZ");
    strVars.push_back("KFactor_EW_qqZZ_unc");
    computeFcn = (syst==tQQBkgEWCorrUp ? ReweightingFunctions::getOnePlusB1OverA1Weight : ReweightingFunctions::getOneMinusB1OverA1Weight);
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn);
    evaluators.push_back(rewgtbuilder);
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, SystematicsHelpers::getRawSystematic);
  }
  else if (syst==tPDFScaleDn || syst==tPDFScaleUp){
    computeFcn = ReweightingFunctions::getA1OverB1Weight;
    strVars.push_back((syst==tPDFScaleDn ? "LHEweight_QCDscale_muR1_muF0p5" : "LHEweight_QCDscale_muR1_muF2"));
    strVars.push_back("LHEweight_QCDscale_muR1_muF1");
    if (proc==ProcessHandler::kGG){
      strVarsNorm.push_back((syst==tPDFScaleDn ? "KFactor_QCD_ggZZ_PDFScaleDn" : "KFactor_QCD_ggZZ_PDFScaleUp"));
      strVarsNorm.push_back("KFactor_QCD_ggZZ_Nominal");
    }
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn); evaluators.push_back(rewgtbuilder);
    if (!strVarsNorm.empty()){ normbuilder = new ReweightingBuilder(strVarsNorm, computeFcn); evaluators.push_back(normbuilder); }
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, (normbuilder ? SystematicsHelpers::getNormalizedSystematic : SystematicsHelpers::getRawSystematic));
  }
  else if (syst==tQCDScaleDn || syst==tQCDScaleUp){
    computeFcn = ReweightingFunctions::getA1OverB1Weight;
    strVars.push_back((syst==tQCDScaleDn ? "LHEweight_QCDscale_muR0p5_muF1" : "LHEweight_QCDscale_muR2_muF1"));
    strVars.push_back("LHEweight_QCDscale_muR1_muF1");
    if (proc==ProcessHandler::kGG){
      strVarsNorm.push_back((syst==tQCDScaleDn ? "KFactor_QCD_ggZZ_QCDScaleDn" : "KFactor_QCD_ggZZ_QCDScaleUp"));
      strVarsNorm.push_back("KFactor_QCD_ggZZ_Nominal");
    }
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn); evaluators.push_back(rewgtbuilder);
    if (!strVarsNorm.empty()){ normbuilder = new ReweightingBuilder(strVarsNorm, computeFcn); evaluators.push_back(normbuilder); }
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, (normbuilder ? SystematicsHelpers::getNormalizedSystematic : SystematicsHelpers::getRawSystematic));
  }
  else if (syst==tAsMZDn || syst==tAsMZUp){
    computeFcn = ReweightingFunctions::getA1OverB1Weight;
    strVars.push_back((syst==tAsMZDn ? "LHEweight_AsMZ_Dn" : "LHEweight_AsMZ_Up"));
    strVars.push_back("LHEweight_QCDscale_muR1_muF1");
    if (proc==ProcessHandler::kGG){
      strVarsNorm.push_back((syst==tAsMZDn ? "KFactor_QCD_ggZZ_AsDn" : "KFactor_QCD_ggZZ_AsUp"));
      strVarsNorm.push_back("KFactor_QCD_ggZZ_Nominal");
    }
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn); evaluators.push_back(rewgtbuilder);
    if (!strVarsNorm.empty()){ normbuilder = new ReweightingBuilder(strVarsNorm, computeFcn); evaluators.push_back(normbuilder); }
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, (normbuilder ? SystematicsHelpers::getNormalizedSystematic : SystematicsHelpers::getRawSystematic));
  }
  else if (syst==tPDFReplicaDn || syst==tPDFReplicaUp){
    computeFcn = ReweightingFunctions::getA1OverB1Weight;
    strVars.push_back((syst==tPDFReplicaDn ? "LHEweight_PDFVariation_Dn" : "LHEweight_PDFVariation_Up"));
    strVars.push_back("LHEweight_QCDscale_muR1_muF1");
    if (proc==ProcessHandler::kGG){
      strVarsNorm.push_back((syst==tPDFReplicaDn ? "KFactor_QCD_ggZZ_PDFReplicaDn" : "KFactor_QCD_ggZZ_PDFReplicaUp"));
      strVarsNorm.push_back("KFactor_QCD_ggZZ_Nominal");
    }
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn); evaluators.push_back(rewgtbuilder);
    if (!strVarsNorm.empty()){ normbuilder = new ReweightingBuilder(strVarsNorm, computeFcn); evaluators.push_back(normbuilder); }
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, (normbuilder ? SystematicsHelpers::getNormalizedSystematic : SystematicsHelpers::getRawSystematic));
  }
  else if (theDataPeriod=="2017" && (syst==SystematicsHelpers::tPythiaScaleDn || syst==SystematicsHelpers::tPythiaScaleUp)){
    computeFcn = ReweightingFunctions::getSimpleWeight;
    strVars.push_back((syst==tPythiaScaleDn ? "PythiaWeight_isr_muR0p25" : "PythiaWeight_isr_muR4"));
    strVars.push_back((syst==tPythiaScaleDn ? "PythiaWeight_fsr_muR0p25" : "PythiaWeight_fsr_muR4"));
    for (CJLSTTree*& tree:trees){ for (TString const& s:strVars) tree->bookBranch<float>(s, 1); for (TString const& s:strVarsNorm) tree->bookBranch<float>(s, 1); }
    rewgtbuilder = new ReweightingBuilder(strVars, computeFcn); evaluators.push_back(rewgtbuilder);
    if (!strVarsNorm.empty()){ normbuilder = new ReweightingBuilder(strVarsNorm, computeFcn); evaluators.push_back(normbuilder); }
    for (ReweightingBuilder*& rb:evaluators){
      rb->setWeightBinning(binning);
      for (CJLSTTree*& tree:trees) rb->setupWeightVariables(tree, (normbuilder ? 1 : -1), 0);
    }
    res = new YieldSystematic(evaluators, (normbuilder ? SystematicsHelpers::getNormalizedSystematic : SystematicsHelpers::getRawSystematic));
  }

  MELAout << "SystematicsHelpers::constructSystematic: Systematics constructed with:\n"
    << "\t- Vars: " << strVars
    << "\t- Norm: " << strVarsNorm
    << endl;
  std::copy(evaluators.begin(), evaluators.end(), std::back_inserter(extraEvaluators));
  return res;
}
void SystematicsHelpers::adjustDiscriminantJECVariables(SystematicsHelpers::SystematicVariationTypes const syst, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  if (syst==eJECDn || syst==eJECUp){
    for (DiscriminantClasses::KDspecs& KD:KDlist){
      for (TString& var:KD.KDvars) HelperFunctions::replaceString<TString, const TString>(var, TString("JECNominal"), TString((syst==eJECDn ? "JECDn" : "JECUp")));
    }
  }
}
TString SystematicsHelpers::getSystematicsName(SystematicsHelpers::SystematicVariationTypes const syst){
  switch (syst){
  case sNominal:
    return "Nominal";
  case tPDFScaleDn:
    return "PDFScaleDn";
  case tPDFScaleUp:
    return "PDFScaleUp";
  case tQCDScaleDn:
    return "QCDScaleDn";
  case tQCDScaleUp:
    return "QCDScaleUp";
  case tAsMZDn:
    return "AsMZDn";
  case tAsMZUp:
    return "AsMZUp";
  case tPDFReplicaDn:
    return "PDFReplicaDn";
  case tPDFReplicaUp:
    return "PDFReplicaUp";
  case tPythiaScaleDn:
    return "PythiaScaleDn";
  case tPythiaScaleUp:
    return "PythiaScaleUp";
  case tPythiaTuneDn:
    return "PythiaTuneDn";
  case tPythiaTuneUp:
    return "PythiaTuneUp";
  case tMINLODn:
    return "MINLODn";
  case tMINLOUp:
    return "MINLOUp";
  case tQQBkgEWCorrDn:
    return "EWCorrDn";
  case tQQBkgEWCorrUp:
    return "EWCorrUp";
  case eLepSFEleDn:
    return "LepEffEleDn";
  case eLepSFEleUp:
    return "LepEffEleUp";
  case eLepSFMuDn:
    return "LepEffMuDn";
  case eLepSFMuUp:
    return "LepEffMuUp";
  case eLepScaleEleDn:
    return "LepScaleEleDn";
  case eLepScaleEleUp:
    return "LepScaleEleUp";
  case eLepScaleMuDn:
    return "LepScaleMuDn";
  case eLepScaleMuUp:
    return "LepScaleMuUp";
  case eLepResEleDn:
    return "LepResEleDn";
  case eLepResEleUp:
    return "LepResEleUp";
  case eLepResMuDn:
    return "LepResMuDn";
  case eLepResMuUp:
    return "LepResMuUp";
  case eJECDn:
    return "JECDn";
  case eJECUp:
    return "JECUp";
  case eBTagSFDn:
    return "BTagSFDn";
  case eBTagSFUp:
    return "BTagSFUp";
  case eZXStatsDn:
    return "ZJetsStatsDn";
  case eZXStatsUp:
    return "ZJetsStatsUp";
  default:
    return "";
  }
}
TString SystematicsHelpers::getSystematicsLabel(SystematicsHelpers::SystematicVariationTypes const syst){
  switch (syst){
  case sNominal:
    return "";
  case tPDFScaleDn:
    return "PDF scale down";
  case tPDFScaleUp:
    return "PDF scale up";
  case tQCDScaleDn:
    return "QCD scale dn";
  case tQCDScaleUp:
    return "QCD scale up";
  case tAsMZDn:
    return "#alpha_{s}(m_{Z}) down";
  case tAsMZUp:
    return "#alpha_{s}(m_{Z}) up";
  case tPDFReplicaDn:
    return "PDF replica down";
  case tPDFReplicaUp:
    return "PDF replica up";
  case tPythiaScaleDn:
    return "Pythia scale down";
  case tPythiaScaleUp:
    return "Pythia scale up";
  case tPythiaTuneDn:
    return "Pythia tune down";
  case tPythiaTuneUp:
    return "Pythia tune up";
  case tMINLODn:
    return "MINLO-HJ down";
  case tMINLOUp:
    return "MINLO-HJ up";
  case tQQBkgEWCorrDn:
    return "NLO EW down";
  case tQQBkgEWCorrUp:
    return "NLO EW up";
  case eLepSFEleDn:
    return "#epsilon_{e} down";
  case eLepSFEleUp:
    return "#epsilon_{e} up";
  case eLepSFMuDn:
    return "#epsilon_{#mu} down";
  case eLepSFMuUp:
    return "#epsilon_{#mu} up";
  case eJECDn:
    return "JEC down";
  case eJECUp:
    return "JEC up";
  case eBTagSFDn:
    return "b-tag SF down";
  case eBTagSFUp:
    return "b-tag SF up";
  case eZXStatsDn:
    return "Fake rate down";
  case eZXStatsUp:
    return "Fake rate up";
  default:
    assert(0);
    return "";
  }
}
TString SystematicsHelpers::getSystematicsCombineName_NoDownUp(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst
){
  TString systname;
  switch (syst){
  case sNominal:
    systname="Nominal";
    break;
  case tPDFScaleDn:
  case tPDFScaleUp:
    systname="QCDscale_fac_[process]";
    break;
  case tQCDScaleDn:
  case tQCDScaleUp:
    systname="QCDscale_ren_[process]";
    break;
  case tAsMZDn:
  case tAsMZUp:
    systname="pdf_asmz_[process]";
    break;
  case tPDFReplicaDn:
  case tPDFReplicaUp:
    systname="pdf_variation_[process]";
    break;
  case tPythiaScaleDn:
  case tPythiaScaleUp:
    systname="CMS_scale_pythia";
    break;
  case tPythiaTuneDn:
  case tPythiaTuneUp:
    systname="CMS_tune_pythia";
    break;
  case tMINLODn:
  case tMINLOUp:
  {
    // Not actually QCD scale but rather missing ME components
    // Also not the BNL prescription http://cms.cern.ch/iCMS/jsp/openfile.jsp?type=NOTE&year=2011&files=NOTE2011_005.pdf
    // Still use QCDscale_[process][Njets]in to correlate with older datacards
    if (proc==ProcessHandler::kGG || proc==ProcessHandler::kQQBkg) systname="QCDscale_[process]2in"; // ggH is NLO, so only has 1 jet from the ME
    else systname="QCDscale_[process]4in"; // VH and VBF already have 3 jets from the ME
    break;
  }
  case tQQBkgEWCorrDn:
  case tQQBkgEWCorrUp:
    systname="EWcorr_[process]";
    break;
  case eLepSFEleDn:
  case eLepSFEleUp:
    systname="CMS_eff_e";
    break;
  case eLepSFMuDn:
  case eLepSFMuUp:
    systname="CMS_eff_m";
    break;
  case eLepScaleEleDn:
  case eLepScaleEleUp:
    systname="CMS_scale_e";
    break;
  case eLepScaleMuDn:
  case eLepScaleMuUp:
    systname="CMS_scale_m";
    break;
  case eLepResEleDn:
  case eLepResEleUp:
    systname="CMS_res_e";
    break;
  case eLepResMuDn:
  case eLepResMuUp:
    systname="CMS_res_m";
    break;
  case eJECDn:
  case eJECUp:
    systname="CMS_scale_j_[sqrts]TeV_[year]";
    break;
  case eBTagSFDn:
  case eBTagSFUp:
    systname="CMS_btag_comb_[sqrts]TeV_[year]";
    break;
  case eZXStatsDn:
  case eZXStatsUp:
    systname="CMS_fake_[channel]";
    break;
  default:
    MELAerr << "SystematicsHelpers::getSystematicsCombineName_NoDownUp: Combine name for systematic " << getSystematicsName(syst) << " is not found! Aborting..." << endl;
    assert(0);
  }

  TString strProcess;
  switch (proc){
  case ProcessHandler::kGG:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ) strProcess="ggH";
    else strProcess="Higgs_gg";
    break;
  }
  case ProcessHandler::kTT:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ) strProcess="ttH";
    else strProcess="Higgs_gg";
    break;
  }
  case ProcessHandler::kBB:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ) strProcess="bbH";
    else strProcess="Higgs_gg";
    break;
  }
  case ProcessHandler::kVV:
  case ProcessHandler::kVBF:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ) strProcess="qqH";
    else strProcess="Higgs_qqbar";
    break;
  }
  case ProcessHandler::kZH:
  case ProcessHandler::kWH:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ) strProcess="VH";
    else strProcess="Higgs_qqbar";
    break;
  }
  case ProcessHandler::kQQBkg:
  {
    if (
      syst==tQCDScaleDn || syst==tQCDScaleUp
      ||
      syst==tPDFScaleDn || syst==tPDFScaleUp
      ||
      syst==tMINLODn || syst==tMINLOUp
      ||
      syst==tQQBkgEWCorrDn || syst==tQQBkgEWCorrUp
      ) strProcess="VV";
    else strProcess="qqbar";
    break;
  }
  case ProcessHandler::kZX:
    strProcess="zjets";
    break;
  default:
    MELAerr << "SystematicsHelpers::getSystematicsCombineName_NoDownUp: " << proc << " process name of systematic " << getSystematicsName(syst) << " is not found! Aborting..." << endl;
    assert(0);
  }
  assert(!systname.Contains("[process]") || strProcess!=""); HelperFunctions::replaceString(systname, "[process]", strProcess.Data());

  TString strSqrts = Form("%i", theSqrts);
  assert(!systname.Contains("[sqrts]") || strSqrts!=""); HelperFunctions::replaceString(systname, "[sqrts]", strSqrts.Data());

  TString strYear = theDataPeriod;
  assert(!systname.Contains("[year]") || strYear!=""); HelperFunctions::replaceString(systname, "[year]", strYear.Data());

  TString strChannel = SampleHelpers::getChannelName(channel);
  assert(!systname.Contains("[channel]") || strChannel!=""); HelperFunctions::replaceString(systname, "[channel]", strChannel.Data());

  TString strCategory = CategorizationHelpers::getCategoryName(category);
  assert(!systname.Contains("[category]") || strCategory!=""); HelperFunctions::replaceString(systname, "[category]", strCategory.Data());

  return systname;
}
TString SystematicsHelpers::getSystematicsCombineName(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst
){
  TString systname=SystematicsHelpers::getSystematicsCombineName_NoDownUp(
    category,
    channel,
    proc,
    syst
  );
  if (((int) syst)%2==1) systname += "Down";
  else if (syst!=sNominal) systname += "Up";
  return systname;
}
