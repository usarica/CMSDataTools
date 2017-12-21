#include "TemplateHelpers.h"
#include "TemplateHelpers.hpp"


/*********************/
/* General functions */
/*********************/
void TemplateHelpers::getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  using namespace SampleHelpers;
  using namespace DiscriminantClasses;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
    KDspecs KDbkg("Dbkgkin");
    KDbkg.KD = constructKDFromType(kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
    KDbkg.KDvars = getKDVars(kDbkgkin);
    KDlist.push_back(KDbkg);
    KDspecs KDbkgsigint("Dggint");
    //KDbkgsigint.KD = constructKDFromType(kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dggbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
    //KDbkgsigint.KDvars = getKDVars(kDggint);
    KDlist.push_back(KDbkgsigint);
    KDspecs KDL1("DL1dec");
    KDL1.KD = constructKDFromType(kDL1dec, "", "", "../data/gConstant_HZZ2e2mu_L1.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1));
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
  else if (category==CategorizationHelpers::JJVBFTagged){
    KDspecs KDbkg("Dbkgkin"); // REPLACE ME
    KDbkg.KD = constructKDFromType(kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
    KDbkg.KDvars = getKDVars(kDbkgkin);
    KDlist.push_back(KDbkg);
    KDspecs KDbkgsigint("Dggint");
    //KDbkgsigint.KD = constructKDFromType(kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dggbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
    //KDbkgsigint.KDvars = getKDVars(kDggint);
    KDlist.push_back(KDbkgsigint);
    KDspecs KDL1("DL1jjVBFdec");
    KDL1.KD = constructKDFromType(kDL1jjVBFdec, "", "", "", "", pow(1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1), 2));
    KDL1.KD->addAdditionalG("../data/gConstant_VBF_L1.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_L1");
    KDL1.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_L1.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");
    KDL1.KDvars = getKDVars(kDL1jjVBFdec);
    KDlist.push_back(KDL1);
    KDspecs KDa2("Da2jjVBFdec");
    KDa2.KD = constructKDFromType(kDa2jjVBFdec, "", "", "", "");
    KDa2.KD->addAdditionalG("../data/gConstant_VBF_g2.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g2");
    KDa2.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_g2.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");
    KDa2.KDvars = getKDVars(kDa2jjVBFdec);
    KDlist.push_back(KDa2);
    KDspecs KDa3("Da3jjVBFdec");
    KDa3.KD = constructKDFromType(kDa3jjVBFdec, "", "", "", "");
    KDa3.KD->addAdditionalG("../data/gConstant_VBF_g4.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g4");
    KDa3.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_g4.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");
    KDa3.KDvars = getKDVars(kDa3jjVBFdec);
    KDlist.push_back(KDa3);
  }
  else if (category==CategorizationHelpers::HadVHTagged){
    KDspecs KDbkg("Dbkgkin"); // REPLACE ME
    KDbkg.KD = constructKDFromType(kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
    KDbkg.KDvars = getKDVars(kDbkgkin);
    KDlist.push_back(KDbkg);
    KDspecs KDbkgsigint("Dggint");
    //KDbkgsigint.KD = constructKDFromType(kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dggbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
    //KDbkgsigint.KDvars = getKDVars(kDggint);
    KDlist.push_back(KDbkgsigint);
    KDspecs KDL1("DL1jjVHdec");
    KDL1.KD = constructKDFromType(kDL1jjVHdec, "", "", "", "", pow(1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1), 2));
    KDL1.KD->addAdditionalG("../data/gConstant_VH_L1.root", "sp_tgfinal_ZH_SM_plus_tgfinal_WH_SM_over_tgfinal_ZH_L1_plus_tgfinal_WH_L1");
    KDL1.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_L1.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");
    KDL1.KDvars = getKDVars(kDL1jjVHdec);
    KDlist.push_back(KDL1);
    KDspecs KDa2("Da2jjVHdec");
    KDa2.KD = constructKDFromType(kDa2jjVHdec, "", "", "", "");
    KDa2.KD->addAdditionalG("../data/gConstant_VH_g2.root", "sp_tgfinal_ZH_SM_plus_tgfinal_WH_SM_over_tgfinal_ZH_g2_plus_tgfinal_WH_g2");
    KDa2.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_g2.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");
    KDa2.KDvars = getKDVars(kDa2jjVHdec);
    KDlist.push_back(KDa2);
    KDspecs KDa3("Da3jjVHdec");
    KDa3.KD = constructKDFromType(kDa3jjVHdec, "", "", "", "");
    KDa3.KD->addAdditionalG("../data/gConstant_VH_g4.root", "sp_tgfinal_ZH_SM_plus_tgfinal_WH_SM_over_tgfinal_ZH_g4_plus_tgfinal_WH_g4");
    KDa3.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_g4.root", "sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");
    KDa3.KDvars = getKDVars(kDa3jjVHdec);
    KDlist.push_back(KDa3);
  }

  adjustDiscriminantJECVariables(strSystematics, KDlist);
}
void TemplateHelpers::getCategorizationDiscriminants(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  using namespace DiscriminantClasses;

  KDspecs KDjjVBF("DjjVBF");
  KDjjVBF.KD = constructKDFromType(kDjjVBF, "../data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth");
  KDjjVBF.KDvars = getKDVars(kDjjVBF);
  KDlist.push_back(KDjjVBF);
  KDspecs KDjjZH("DjjZH");
  KDjjZH.KD = constructKDFromType(kDjjZH, "../data/SmoothKDConstant_m4l_DjjZH_13TeV.root", "sp_gr_varReco_Constant_Smooth");
  KDjjZH.KDvars = getKDVars(kDjjZH);
  KDlist.push_back(KDjjZH);
  KDspecs KDjjWH("DjjWH");
  KDjjWH.KD = constructKDFromType(kDjjWH, "../data/SmoothKDConstant_m4l_DjjWH_13TeV.root", "sp_gr_varReco_Constant_Smooth");
  KDjjWH.KDvars = getKDVars(kDjjWH);
  KDlist.push_back(KDjjWH);

  KDspecs KDjjVBFL1("DjjVBFL1");
  KDjjVBFL1.KD = constructKDFromType(kDjjVBFL1, "../data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_VBF_L1.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_L1", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1));
  KDjjVBFL1.KD->setInvertG(true);
  KDjjVBFL1.KDvars = getKDVars(kDjjVBFL1);
  KDlist.push_back(KDjjVBFL1);
  KDspecs KDjjZHL1("DjjZHL1");
  KDjjZHL1.KD = constructKDFromType(kDjjZHL1, "../data/SmoothKDConstant_m4l_DjjZH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_ZH_L1.root", "sp_tgfinal_ZH_SM_over_tgfinal_ZH_L1", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1));
  KDjjZHL1.KD->setInvertG(true);
  KDjjZHL1.KDvars = getKDVars(kDjjZHL1);
  KDlist.push_back(KDjjZHL1);
  KDspecs KDjjWHL1("DjjWHL1");
  KDjjWHL1.KD = constructKDFromType(kDjjWHL1, "../data/SmoothKDConstant_m4l_DjjWH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_WH_L1.root", "sp_tgfinal_WH_SM_over_tgfinal_WH_L1", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1));
  KDjjWHL1.KD->setInvertG(true);
  KDjjWHL1.KDvars = getKDVars(kDjjWHL1);
  KDlist.push_back(KDjjWHL1);

  KDspecs KDjjVBFa2("DjjVBFa2");
  KDjjVBFa2.KD = constructKDFromType(kDjjVBFa2, "../data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_VBF_g2.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g2");
  KDjjVBFa2.KD->setInvertG(true);
  KDjjVBFa2.KDvars = getKDVars(kDjjVBFa2);
  KDlist.push_back(KDjjVBFa2);
  KDspecs KDjjZHa2("DjjZHa2");
  KDjjZHa2.KD = constructKDFromType(kDjjZHa2, "../data/SmoothKDConstant_m4l_DjjZH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_ZH_g2.root", "sp_tgfinal_ZH_SM_over_tgfinal_ZH_g2");
  KDjjZHa2.KD->setInvertG(true);
  KDjjZHa2.KDvars = getKDVars(kDjjZHa2);
  KDlist.push_back(KDjjZHa2);
  KDspecs KDjjWHa2("DjjWHa2");
  KDjjWHa2.KD = constructKDFromType(kDjjWHa2, "../data/SmoothKDConstant_m4l_DjjWH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_WH_g2.root", "sp_tgfinal_WH_SM_over_tgfinal_WH_g2");
  KDjjWHa2.KD->setInvertG(true);
  KDjjWHa2.KDvars = getKDVars(kDjjWHa2);
  KDlist.push_back(KDjjWHa2);

  KDspecs KDjjVBFa3("DjjVBFa3");
  KDjjVBFa3.KD = constructKDFromType(kDjjVBFa3, "../data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_VBF_g4.root", "sp_tgfinal_VBF_SM_over_tgfinal_VBF_g4");
  KDjjVBFa3.KD->setInvertG(true);
  KDjjVBFa3.KDvars = getKDVars(kDjjVBFa3);
  KDlist.push_back(KDjjVBFa3);
  KDspecs KDjjZHa3("DjjZHa3");
  KDjjZHa3.KD = constructKDFromType(kDjjZHa3, "../data/SmoothKDConstant_m4l_DjjZH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_ZH_g4.root", "sp_tgfinal_ZH_SM_over_tgfinal_ZH_g4");
  KDjjZHa3.KD->setInvertG(true);
  KDjjZHa3.KDvars = getKDVars(kDjjZHa3);
  KDlist.push_back(KDjjZHa3);
  KDspecs KDjjWHa3("DjjWHa3");
  KDjjWHa3.KD = constructKDFromType(kDjjWHa3, "../data/SmoothKDConstant_m4l_DjjWH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_WH_g4.root", "sp_tgfinal_WH_SM_over_tgfinal_WH_g4");
  KDjjWHa3.KD->setInvertG(true);
  KDjjWHa3.KDvars = getKDVars(kDjjWHa3);
  KDlist.push_back(KDjjWHa3);

  adjustDiscriminantJECVariables(strSystematics, KDlist);
}
void TemplateHelpers::adjustDiscriminantJECVariables(const TString strSystematics, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  if (strSystematics=="JECUp" || strSystematics=="JECDn"){
    for (DiscriminantClasses::KDspecs& KD:KDlist){
      for (TString& var:KD.KDvars) HelperFunctions::replaceString<TString, const TString>(var, TString("JECNominal"), strSystematics);
    }
  }
}
