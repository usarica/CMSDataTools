#include "DiscriminantClasses.h"


DiscriminantClasses::Type DiscriminantClasses::getKDType(const TString name){
  if (name=="Dbkgkin") return kDbkgkin;
  else if (name=="Dbkgdec") return kDbkgdec;
  else if (name=="Dggint") return kDggint;

  else if (name=="DjVBF") return kDjVBF;
  else if (name=="DjjVBF") return kDjjVBF;

  else if (name=="DjjZH") return kDjjZH;
  else if (name=="DjjWH") return kDjjWH;

  else if (name=="DbkgjjEWQCD") return kDbkgjjEWQCD;

  else return kNTypes;
}
Discriminant* DiscriminantClasses::constructKDFromType(const DiscriminantClasses::Type type, const TString filename, const TString splinename){
  Discriminant* res=nullptr;
  switch (type){
  case kDbkgkin:
    return new Dbkgkin_t(filename, splinename);
  case kDbkgdec:
    return new Dbkgdec_t(filename, splinename);
  case kDggint:
    return new Dintkin_t(filename, splinename);

  case kDjVBF:
    return new DjVBF_t(filename, splinename);
  case kDjjVBF:
    return new DjjVBF_t(filename, splinename);

  case kDjjZH:
  case kDjjWH:
    return new DjjVH_t(filename, splinename);

  case kDbkgjjEWQCD:
    return new DbkgjjEWQCD_t(filename, splinename);

  default:
    return res;
  };
}
Discriminant* DiscriminantClasses::constructKDFromType(const TString name, const TString filename, const TString splinename){ return constructKDFromType(getKDType(name), filename, splinename); }

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

  default:
    break;
  };
  return res;
}
std::vector<TString> DiscriminantClasses::getKDVars(const TString name){ return getKDVars(getKDType(name)); }
