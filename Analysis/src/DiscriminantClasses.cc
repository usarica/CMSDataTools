#include "DiscriminantClasses.h"


namespace DiscriminantClasses{
  const std::unordered_map<TString, DiscriminantClasses::Type> mapKDNameType = DiscriminantClasses::getKDNameTypeMap();
}

DiscriminantClasses::KDspecs::KDspecs() : KD(nullptr) {}
DiscriminantClasses::KDspecs::KDspecs(TString strname) : KDname(strname), KD(nullptr) {}
bool DiscriminantClasses::KDspecs::isValid() const{ return (KD!=nullptr); }

std::unordered_map<TString, DiscriminantClasses::Type> DiscriminantClasses::getKDNameTypeMap(){
  std::unordered_map<TString, DiscriminantClasses::Type> res;

  res["Dbkgkin"] = kDbkgkin;
  res["Dbkgdec"] = kDbkgdec;

  res["Dggbkgkin"] = kDggbkgkin;
  res["Dggint"] = kDggint;

  res["DjVBF"] = kDjVBF;
  res["DjjVBF"] = kDjjVBF;

  res["DjjZH"] = kDjjZH;
  res["DjjWH"] = kDjjWH;

  res["DbkgjjEWQCD"] = kDbkgjjEWQCD;

  res["DL1dec"] = kDL1dec;
  res["DL1decint"] = kDL1decint;
  res["Da2dec"] = kDa2dec;
  res["Da2decint"] = kDa2decint;
  res["Da3dec"] = kDa3dec;
  res["Da3decint"] = kDa3decint;

  return res;
}

DiscriminantClasses::Type DiscriminantClasses::getKDType(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return it->second;
  else return kNTypes;
}
TString DiscriminantClasses::getKDName(DiscriminantClasses::Type type){
  for (auto it=mapKDNameType.cbegin(); it!=mapKDNameType.cend(); it++){ if (it->second==type) return it->first; }
  return "";
}

TString DiscriminantClasses::getKDLabel(DiscriminantClasses::Type type){
  switch (type){
  case kDbkgkin:
    return "D^{kin}_{bkg}";
  case kDbkgdec:
    return "D^{dec}_{bkg}";

  case kDggbkgkin:
    return "D^{gg}_{bkg}";
  case kDggint:
    return "D^{gg}_{int}";

  case kDjVBF:
    return "D^{VBF}_{j}";
  case kDjjVBF:
    return "D^{VBF}_{jj}";

  case kDjjZH:
    return "D^{ZH}_{jj}";
  case kDjjWH:
    return "D^{WH}_{jj}";

  case kDbkgjjEWQCD:
    return "D^{jj+dec}_{bkg}";

  case kDL1dec:
    return "D^{dec}_{#Lambda1}";
  case kDa2dec:
    return "D^{dec}_{0h+}";
  case kDa3dec:
    return "D^{dec}_{0-}";
  case kDL1decint:
    return "D^{dec}_{#Lambda1, int}";
  case kDa2decint:
    return "D^{dec}_{int}";
  case kDa3decint:
    return "D^{dec}_{CP}";

  default:
    return "";
  };
}
TString DiscriminantClasses::getKDLabel(TString name){ return getKDLabel(getKDType(name)); }

Discriminant* DiscriminantClasses::constructKDFromType(
  const DiscriminantClasses::Type type,
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename, const float gscale
){
  Discriminant* res=nullptr;
  switch (type){
  case kDbkgkin:
    return new Dbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDbkgdec:
    return new Dbkgdec_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDggbkgkin:
    return new Dggbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDggint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDjVBF:
    return new DjVBF_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDjjVBF:
    return new DjjVBF_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDjjZH:
  case kDjjWH:
    return new DjjVH_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDbkgjjEWQCD:
    return new DbkgjjEWQCD_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDL1dec:
  case kDa2dec:
  case kDa3dec:
    return new Dbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDL1decint:
  case kDa2decint:
  case kDa3decint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  default:
    return res;
  };
}
Discriminant* DiscriminantClasses::constructKDFromType(
  const TString name,
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename, const float gscale
){ return constructKDFromType(getKDType(name), cfilename, splinename, gfilename, gsplinename, gscale); }

std::vector<TString> DiscriminantClasses::getKDVars(const Type type){
  std::vector<TString> res;
  // In the following statements, JECNominal is to be replaced later
  switch (type){
  case kDbkgkin:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_QQB_BKG_MCFM");
    break;
  case kDbkgdec:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    res.push_back("p_QQB_BKG_MCFM");
    res.push_back("p_Const_GG_BKG_MCFM");
    res.push_back("p_Const_QQB_BKG_MCFM");
    break;

  case kDggbkgkin:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    break;
  case kDggint:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    res.push_back("p_GG_BSI_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_Const_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_Const_GG_BKG_MCFM");
    break;

  case kDjVBF:
    res.push_back("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBF:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;

  case kDjjZH:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;
  case kDjjWH:
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;

  case kDbkgjjEWQCD:
    /* Not yet implemented */
    break;

  case kDL1dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    break;
  case kDa2dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    break;
  case kDa3dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    break;

  case kDL1decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen");
    break;
  case kDa2decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen");
    break;
  case kDa3decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen");
    break;

  default:
    break;
  };
  return res;
}
std::vector<TString> DiscriminantClasses::getKDVars(const TString name){ return getKDVars(getKDType(name)); }

bool DiscriminantClasses::isCPSensitive(const Type type){
  bool res = (type==kDa3decint);
  return res;
}
bool DiscriminantClasses::isCPSensitive(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return isCPSensitive(it->second);
  else return false;
}
