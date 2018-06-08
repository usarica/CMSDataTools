#include <cmath>
#include <cassert>
#include "CategorizationHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


namespace CategorizationHelpers{
  CategorizationScheme globalCategorizationScheme=UntaggedOrJJVBFOrHadVH;
}

using namespace CategorizationHelpers;


TString CategorizationHelpers::getCategoryName(CategorizationHelpers::Category category){
  switch(category){
  case Inclusive:
    return "Inclusive";
  case Untagged:
    return "Untagged";
  case JJVBFTagged:
    return "JJVBFTagged";
  case HadVHTagged:
    return "HadVHTagged";
  default:
    return "";
  }
}

TString CategorizationHelpers::getCategoryLabel(CategorizationHelpers::Category category){
  switch (category){
  case Inclusive:
    return "Inclusive";
  case Untagged:
    return "Untagged";
  case JJVBFTagged:
    return "VBF 2 jets";
  case HadVHTagged:
    return "Had. VH";
  default:
    return "";
  }
}

TString CategorizationHelpers::getCategoryLabelForKDs(CategorizationHelpers::Category category){
  switch (category){
  case JJVBFTagged:
    return "VBF";
  case HadVHTagged:
    return "VH";
  default:
    return "";
  }
}

void CategorizationHelpers::setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme){ globalCategorizationScheme=scheme; }

CategorizationHelpers::Category CategorizationHelpers::getCategory(SimpleEntry const& sel_vars, const bool forceUntagged){
  if (forceUntagged) return Untagged;

  if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag){
    float DjjVBF=-1; sel_vars.getNamedVal("DjjVBF", DjjVBF);
    float DjjZH=-1; sel_vars.getNamedVal("DjjZH", DjjZH);
    float DjjWH=-1; sel_vars.getNamedVal("DjjWH", DjjWH);
    short nExtraLep=0; sel_vars.getNamedVal("nExtraLep", nExtraLep);
    short nJets=0; sel_vars.getNamedVal("nJets", nJets);
    short nJets_BTagged=0; sel_vars.getNamedVal("nJets_BTagged", nJets_BTagged);
    float const DjjVH=std::max(DjjZH, DjjWH);
    if (
      DjjVBF>=0.5
      &&
      nExtraLep==0
      &&
      (((nJets==2 || nJets==3) && nJets_BTagged<=1) || (nJets>=4 && nJets_BTagged==0))
      ) return JJVBFTagged;
    else if (
      DjjVH>=0.5
      &&
      nExtraLep==0
      &&
      (nJets==2 || nJets==3 || (nJets>=4 && nJets_BTagged==0))
      ) return HadVHTagged;
    else return Untagged;
  }
  else if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH){
    float DjjVBF=-1; sel_vars.getNamedVal("DjjVBF", DjjVBF);
    float DjjZH=-1; sel_vars.getNamedVal("DjjZH", DjjZH);
    float DjjWH=-1; sel_vars.getNamedVal("DjjWH", DjjWH);
    float const DjjVH=std::max(DjjZH, DjjWH);
    if (
      DjjVBF>=0.5
      ) return JJVBFTagged;
    else if (
      DjjVH>=0.5
      ) return HadVHTagged;
    else return Untagged;
  }
  else if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH_Arbitrated){
    float DjjVBF=-1; sel_vars.getNamedVal("DjjVBF", DjjVBF);
    float DjjZH=-1; sel_vars.getNamedVal("DjjZH", DjjZH);
    float DjjWH=-1; sel_vars.getNamedVal("DjjWH", DjjWH);
    float const DjjVH=std::max(DjjZH, DjjWH);
    if (
      DjjVBF>=0.5 && DjjVBF>=DjjVH
      ) return JJVBFTagged;
    else if (
      DjjVH>=0.5 && DjjVBF<DjjVH
      ) return HadVHTagged;
    else return Untagged;
  }
  else/* if (globalCategorizationScheme==UntaggedOrJJVBF)*/{
    float DjjVBF=-1; sel_vars.getNamedVal("DjjVBF", DjjVBF);
    if (
      DjjVBF>=0.5
      ) return JJVBFTagged;
    else return Untagged;
  }
}
std::vector<CategorizationHelpers::Category> CategorizationHelpers::getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme){
  std::vector<CategorizationHelpers::Category> res;
  res.push_back(Inclusive);
  switch (scheme){
  case UntaggedOrJJVBFOrHadVH:
  case UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag:
  case UntaggedOrJJVBFOrHadVH_Arbitrated:
    res.push_back(HadVHTagged);
  case UntaggedOrJJVBF:
    res.push_back(JJVBFTagged);
    res.push_back(Untagged);
  default:
    break;
  }
  return res;
}
bool CategorizationHelpers::testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory){
  std::vector<Category> cats=getAllowedCategories(globalCategorizationScheme);
  for (Category const& cat:cats){ if (cat==theCategory) return true; }
  return false;
}

TString CategorizationHelpers::getMassRegionName(CategorizationHelpers::MassRegion massregion){
  switch (massregion){
  case kOnshell:
    return "Onshell";
  case kOffshell:
    return "Offshell";
  default:
    MELAerr << "CategorizationHelpers::getMassRegionName: Mass region " << massregion << " has undefined name!" << endl;
    assert(0);
    return "";
  }
}

template<> void CategorizationHelpers::getExtraCategorizationVariables<short>(CategorizationHelpers::CategorizationScheme scheme, SystematicsHelpers::SystematicVariationTypes syst, std::vector<TString>& res){
  if (scheme==CategorizationHelpers::UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag){
    res.push_back("nExtraLep");
    switch (syst){
    case SystematicsHelpers::eJECDn:
      res.push_back("nCleanedJetsPt30_jecDn");
      res.push_back("nCleanedJetsPt30BTagged_bTagSF_jecDn");
      break;
    case SystematicsHelpers::eJECUp:
      res.push_back("nCleanedJetsPt30_jecUp");
      res.push_back("nCleanedJetsPt30BTagged_bTagSF_jecUp");
      break;
    case SystematicsHelpers::eBTagSFDn:
      res.push_back("nCleanedJetsPt30");
      res.push_back("nCleanedJetsPt30BTagged_bTagSFDn");
      break;
    case SystematicsHelpers::eBTagSFUp:
      res.push_back("nCleanedJetsPt30");
      res.push_back("nCleanedJetsPt30BTagged_bTagSFUp");
      break;
    default:
      res.push_back("nCleanedJetsPt30");
      res.push_back("nCleanedJetsPt30BTagged_bTagSF");
      break;
    }
  }
}
