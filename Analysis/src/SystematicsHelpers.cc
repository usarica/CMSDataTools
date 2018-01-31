#include <iterator>
#include "SystematicsHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace SystematicsHelpers;


SystematicsHelpers::YieldSystematic::YieldSystematic(const std::vector<ReweightingBuilder*>& inEvaluators, float(*infcn)(CJLSTTree*, const std::vector<ReweightingBuilder*>&)) : rule(infcn), evaluators(inEvaluators) {}
float SystematicsHelpers::YieldSystematic::eval(CJLSTTree* theTree) const{
  float weight=0;
  if (rule) weight=rule(theTree, evaluators);
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
  std::pair<float, float>(*infcn)(short const&, short const&, std::vector<std::vector<float>*> const&, unsigned int const),
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
  std::vector<std::vector<float>*> refVar;
  for (TString const& s:strVars){
    vector<float>* v = nullptr;
    theTree->getValRef(s, v);
    if (!v) MELAerr << "PerLeptonSystematic::setup: Could not get the reference for " << s << endl;
    refVar.push_back(v);
  }
  if (!refVar.empty()) componentRefs[theTree]=std::pair<std::vector<short*>, std::vector<std::vector<float>*>>(refId, refVar);
}

std::pair<float, float> SystematicsHelpers::getLeptonSFSystematic(short const& Z1Flav, short const& Z2Flav, std::vector<std::vector<float>*> const& LepVars, unsigned int const idreq){
  std::pair<float, float> res(1, 1);
  if ((Z1Flav*Z2Flav)%(short)idreq != 0) return res;
  assert(LepVars.size()==5);
  float const& lep4pt = LepVars.at(0)->back();
  for (unsigned int ilep=0; ilep<LepVars.at(0)->size(); ilep++){
    short absZflav;
    if (ilep<2) absZflav=std::abs(Z1Flav);
    else absZflav=std::abs(Z2Flav);
    if (absZflav%(short) idreq==0){
      res.first *= (1.-LepVars.at(2)->at(ilep)/LepVars.at(1)->at(ilep))*(1.-LepVars.at(4)->at(ilep)/LepVars.at(3)->at(ilep));
      res.second *= (1.+LepVars.at(2)->at(ilep)/LepVars.at(1)->at(ilep))*(1.+LepVars.at(4)->at(ilep)/LepVars.at(3)->at(ilep));
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


int SystematicsHelpers::convertSystematicVariationTypeToInt(SystematicsHelpers::SystematicVariationTypes type){ return (int) type; }
std::vector<SystematicsHelpers::SystematicVariationTypes> SystematicsHelpers::getProcessSystematicVariations(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const type
){
  std::vector<SystematicVariationTypes> res;

  res.push_back(sNominal);

  //res.push_back(eZJetsStatsDn);
  //res.push_back(eZJetsStatsUp);

  // FIXME: Z+X DOES NOT HAVE THESE, BUT WE NEED ZJETS PROCESS FIRST TO TEST
  if (channel==SampleHelpers::k4e || channel==SampleHelpers::k2e2mu){
    res.push_back(eLepSFEleDn);
    res.push_back(eLepSFEleUp);
  }
  if (channel==SampleHelpers::k4mu || channel==SampleHelpers::k2e2mu){
    res.push_back(eLepSFMuDn);
    res.push_back(eLepSFMuUp);
  }

  res.push_back(tPDFScaleDn);
  res.push_back(tPDFScaleUp);
  res.push_back(tQCDScaleDn);
  res.push_back(tQCDScaleUp);
  res.push_back(tAsMZDn);
  res.push_back(tAsMZUp);
  res.push_back(tPDFReplicaDn);
  res.push_back(tPDFReplicaUp);

  if (type==ProcessHandler::kQQBkg){
    res.push_back(tQQBkgEWCorrDn);
    res.push_back(tQQBkgEWCorrUp);
  }

  if (category==CategorizationHelpers::JJVBFTagged || category==CategorizationHelpers::HadVHTagged){
    res.push_back(eJECDn);
    res.push_back(eJECUp);
  }

  return res;
}
bool SystematicsHelpers::systematicAllowed(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst
){
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedTypes = SystematicsHelpers::getProcessSystematicVariations(category, channel, proc);
  for (SystematicVariationTypes& st:allowedTypes){ if (st==syst) return true; }
  return false;
}
SystematicsHelpers::SystematicsClass* SystematicsHelpers::constructSystematic(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proc,
  SystematicsHelpers::SystematicVariationTypes const syst,
  std::vector<CJLSTTree*> trees,
  std::vector<ReweightingBuilder*>& extraEvaluators
){
  SystematicsClass* res=nullptr;
  if (!systematicAllowed(category, channel, proc, syst)) return res;

  ExtendedBinning binning((theSqrts*1000.-70.)/10., 70., theSqrts*1000., "GenHMass");
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
  case eJECDn:
    return "JECDn";
  case eJECUp:
    return "JECUp";
  case eZJetsStatsDn:
    return "ZJetsStatsDn";
  case eZJetsStatsUp:
    return "ZJetsStatsUp";
  default:
    return "";
  }
}
