#include "TemplateHelpers.h"
#include "TemplateHelpers.hpp"
#include "MELAStreamHelpers.hh"


using namespace MELAStreamHelpers;


/*********************/
/* General functions */
/*********************/
void TemplateHelpers::getLikelihoodDiscriminants(const SampleHelpers::Channel channel, const CategorizationHelpers::Category category, const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist){
  using namespace SampleHelpers;
  using namespace DiscriminantClasses;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
    if (channel==k4e || channel==k4mu || channel==k2e2mu){
      KDspecs KDbkgkin("Dbkgkin");
      KDbkgkin.KD = constructKDFromType(kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
      KDbkgkin.KDvars = getKDVars(kDbkgkin);
      KDlist.push_back(KDbkgkin);
      KDspecs KDbkgm4l("Dbkgm4l");
      KDbkgm4l.KD = constructKDFromType(kDbkgm4l, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
      KDbkgm4l.KDvars = getKDVars(kDbkgm4l);
      KDlist.push_back(KDbkgm4l);
      /*
      KDspecs KDbkgsigint("Dggint");
      KDbkgsigint.KD = constructKDFromType(kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dggbkgkin_", strChannel.Data(), "_13TeV.root"), "sp_gr_varReco_Constant_Smooth");
      KDbkgsigint.KDvars = getKDVars(kDggint);
      KDlist.push_back(KDbkgsigint);
      */
      KDspecs KCggint("Cggint");
      KCggint.KD = constructKDFromType(kCggint, "", "");
      KCggint.KDvars = getKDVars(kCggint);
      KDlist.push_back(KCggint);
    }

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
    KDspecs KDL1ZGs("DL1ZGsdec");
    KDL1ZGs.KD = constructKDFromType(kDL1ZGsdec, "", "", "../data/gConstant_HZZ2e2mu_L1Zgs.root", "sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1ZGs));
    KDL1ZGs.KDvars = getKDVars(kDL1ZGsdec);
    KDlist.push_back(KDL1ZGs);

    //KDspecs KDL1int("CL1decint");
    //KDL1int.KD = constructKDFromType(kCL1decint, "", "");
    //KDL1int.KDvars = getKDVars(kCL1decint);
    //KDlist.push_back(KDL1int);
    KDspecs KDa2int("Ca2decint");
    KDa2int.KD = constructKDFromType(kCa2decint, "", "");
    KDa2int.KDvars = getKDVars(kCa2decint);
    KDlist.push_back(KDa2int);
    KDspecs KDa3int("Ca3decint");
    KDa3int.KD = constructKDFromType(kCa3decint, "", "");
    KDa3int.KDvars = getKDVars(kCa3decint);
    KDlist.push_back(KDa3int);
    KDspecs KDL1ZGsint("CL1ZGsdecint");
    KDL1ZGsint.KD = constructKDFromType(kCL1ZGsdecint, "", "");
    KDL1ZGsint.KDvars = getKDVars(kCL1ZGsdecint);
    KDlist.push_back(KDL1ZGsint);
  }
  else if (category==CategorizationHelpers::JJVBFTagged){
    //getLikelihoodDiscriminants(channel, CategorizationHelpers::Inclusive, syst, KDlist);
    if (channel==k4e || channel==k4mu || channel==k2e2mu){
      KDspecs KDbkgjjEWQCD("DbkgjjEWQCD");
      KDbkgjjEWQCD.KD = constructKDFromType(
        kDbkgjjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDbkgjjEWQCD.KDvars = getKDVars(kDbkgjjEWQCD);
      KDlist.push_back(KDbkgjjEWQCD);
      KDspecs KDbkgm4ljjEWQCD("Dbkgm4ljjEWQCD");
      KDbkgm4ljjEWQCD.KD = constructKDFromType(
        kDbkgm4ljjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDbkgm4ljjEWQCD.KDvars = getKDVars(kDbkgm4ljjEWQCD);
      KDlist.push_back(KDbkgm4ljjEWQCD);
      /*
      KDspecs KDintjjEWQCD("DintjjEWQCD");
      KDintjjEWQCD.KD = constructKDFromType(
        kDintjjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDintjjEWQCD.KDvars = getKDVars(kDintjjEWQCD);
      KDlist.push_back(KDintjjEWQCD);
      */
      KDspecs KCjjVBFint("CjjVBFint");
      KCjjVBFint.KD = constructKDFromType(kCjjVBFint, "", "");
      KCjjVBFint.KDvars = getKDVars(kCjjVBFint);
      KDlist.push_back(KCjjVBFint);
    }

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
    KDspecs KDL1ZGs("DL1ZGsjjVBFdec");
    KDL1ZGs.KD = constructKDFromType(kDL1ZGsjjVBFdec, "", "", "", "", pow(1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1ZGs), 2));
    KDL1ZGs.KD->addAdditionalG("../data/gConstant_VBF_L1Zgs.root", "sp_tgfinal_VBF_SM_photoncut_over_tgfinal_VBF_L1Zgs");
    KDL1ZGs.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_L1Zgs.root", "sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");
    KDL1ZGs.KDvars = getKDVars(kDL1ZGsjjVBFdec);
    KDlist.push_back(KDL1ZGs);

    //KDspecs KDL1int("CL1jjVBFint");
    //KDL1int.KD = constructKDFromType(kCL1jjVBFint, "", "");
    //KDL1int.KDvars = getKDVars(kCL1jjVBFint);
    //KDlist.push_back(KDL1int);
    KDspecs KDa2int("Ca2jjVBFint");
    KDa2int.KD = constructKDFromType(kCa2jjVBFint, "", "");
    KDa2int.KDvars = getKDVars(kCa2jjVBFint);
    KDlist.push_back(KDa2int);
    KDspecs KDa3int("Ca3jjVBFint");
    KDa3int.KD = constructKDFromType(kCa3jjVBFint, "", "");
    KDa3int.KDvars = getKDVars(kCa3jjVBFint);
    KDlist.push_back(KDa3int);
    KDspecs KDL1ZGsint("CL1ZGsjjVBFint");
    KDL1ZGsint.KD = constructKDFromType(kCL1ZGsjjVBFint, "", "");
    KDL1ZGsint.KDvars = getKDVars(kCL1ZGsjjVBFint);
    KDlist.push_back(KDL1ZGsint);
  }
  else if (category==CategorizationHelpers::HadVHTagged){
    if (channel==k4e || channel==k4mu || channel==k2e2mu){
      KDspecs KDbkgjjEWQCD("DbkgjjEWQCD");
      KDbkgjjEWQCD.KD = constructKDFromType(
        kDbkgjjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDbkgjjEWQCD.KDvars = getKDVars(kDbkgjjEWQCD);
      KDlist.push_back(KDbkgjjEWQCD);
      KDspecs KDbkgm4ljjEWQCD("Dbkgm4ljjEWQCD");
      KDbkgm4ljjEWQCD.KD = constructKDFromType(
        kDbkgm4ljjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDbkgm4ljjEWQCD.KDvars = getKDVars(kDbkgm4ljjEWQCD);
      KDlist.push_back(KDbkgm4ljjEWQCD);
      /*
      KDspecs KDintjjEWQCD("DintjjEWQCD");
      KDintjjEWQCD.KD = constructKDFromType(
        kDintjjEWQCD,
        Form(
          "%s%s%s",
          (channel==k2e2mu ? "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_" : "../data/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_"),
          strCategory.Data(), "_13TeV.root"
        ),
        "sp_gr_varReco_Constant_Smooth"
      );
      KDintjjEWQCD.KDvars = getKDVars(kDintjjEWQCD);
      KDlist.push_back(KDintjjEWQCD);
      */
      KDspecs KCjjVHint("CjjVHint");
      KCjjVHint.KD = constructKDFromType(kCjjVHint, "", "");
      KCjjVHint.KDvars = getKDVars(kCjjVHint);
      KDlist.push_back(KCjjVHint);
    }

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
    KDspecs KDL1ZGs("DL1ZGsjjVHdec");
    KDL1ZGs.KD = constructKDFromType(kDL1ZGsjjVHdec, "", "", "", "", pow(1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1ZGs), 2));
    KDL1ZGs.KD->addAdditionalG("../data/gConstant_VH_L1Zgs.root", "sp_tgfinal_ZH_SM_photoncut_plus_tgfinal_WH_SM_over_tgfinal_ZH_L1Zgs");
    KDL1ZGs.KD->addAdditionalG("../data/gConstant_HZZ2e2mu_L1Zgs.root", "sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");
    KDL1ZGs.KDvars = getKDVars(kDL1ZGsjjVHdec);
    KDlist.push_back(KDL1ZGs);

    //KDspecs KDL1int("CL1jjVHint");
    //KDL1int.KD = constructKDFromType(kCL1jjVHint, "", "");
    //KDL1int.KDvars = getKDVars(kCL1jjVHint);
    //KDlist.push_back(KDL1int);
    KDspecs KDa2int("Ca2jjVHint");
    KDa2int.KD = constructKDFromType(kCa2jjVHint, "", "");
    KDa2int.KDvars = getKDVars(kCa2jjVHint);
    KDlist.push_back(KDa2int);
    KDspecs KDa3int("Ca3jjVHint");
    KDa3int.KD = constructKDFromType(kCa3jjVHint, "", "");
    KDa3int.KDvars = getKDVars(kCa3jjVHint);
    KDlist.push_back(KDa3int);
    KDspecs KDL1ZGsint("CL1ZGsjjVHint");
    KDL1ZGsint.KD = constructKDFromType(kCL1ZGsjjVHint, "", "");
    KDL1ZGsint.KDvars = getKDVars(kCL1ZGsjjVHint);
    KDlist.push_back(KDL1ZGsint);
  }

  SystematicsHelpers::adjustDiscriminantJECVariables(syst, KDlist);
}
void TemplateHelpers::getCategorizationDiscriminants(const SystematicsHelpers::SystematicVariationTypes syst, std::vector<DiscriminantClasses::KDspecs>& KDlist){
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

  KDspecs KDjjVBFL1ZGs("DjjVBFL1ZGs");
  KDjjVBFL1ZGs.KD = constructKDFromType(kDjjVBFL1ZGs, "../data/SmoothKDConstant_m4l_DjjVBF_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_VBF_L1Zgs.root", "sp_tgfinal_VBF_SM_photoncut_over_tgfinal_VBF_L1Zgs", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1ZGs));
  KDjjVBFL1ZGs.KD->setInvertG(true);
  KDjjVBFL1ZGs.KDvars = getKDVars(kDjjVBFL1ZGs);
  KDlist.push_back(KDjjVBFL1ZGs);
  KDspecs KDjjZHL1ZGs("DjjZHL1ZGs");
  KDjjZHL1ZGs.KD = constructKDFromType(kDjjZHL1ZGs, "../data/SmoothKDConstant_m4l_DjjZH_13TeV.root", "sp_gr_varReco_Constant_Smooth", "../data/gConstant_ZH_L1Zgs.root", "sp_tgfinal_ZH_SM_photoncut_over_tgfinal_ZH_L1Zgs", 1./ACHypothesisHelpers::getACHypothesisMEHZZGVal(ACHypothesisHelpers::kL1ZGs));
  KDjjZHL1ZGs.KD->setInvertG(true);
  KDjjZHL1ZGs.KDvars = getKDVars(kDjjZHL1ZGs);
  KDlist.push_back(KDjjZHL1ZGs);

  // Set WPs
  for (auto& KDspec:KDlist) KDspec.KD->setWP(DiscriminantClasses::getKDWP(KDspec.KDtype));

  SystematicsHelpers::adjustDiscriminantJECVariables(syst, KDlist);
}
ExtendedBinning TemplateHelpers::getDiscriminantFineBinning(const SampleHelpers::Channel /*channel*/, const CategorizationHelpers::Category category, ACHypothesisHelpers::ACHypothesis hypo, TString const strKD, CategorizationHelpers::MassRegion const massregion){
  MELAout << "TemplateHelpers::getDiscriminantFineBinning: Inquiring binning for variable " << strKD << endl;
  ExtendedBinning res(strKD);
  if (strKD=="ZZMass"){
    switch (massregion){
    case CategorizationHelpers::kOffshell:
    {
      res.addBinBoundary(theSqrts*1000.);
      res.addBinBoundary(220);
      res.addBinBoundary(230);
      res.addBinBoundary(240);
      res.addBinBoundary(250);
      res.addBinBoundary(260);
      res.addBinBoundary(280);
      res.addBinBoundary(310);
      res.addBinBoundary(340);
      res.addBinBoundary(370);
      res.addBinBoundary(400);
      res.addBinBoundary(475);
      res.addBinBoundary(550);
      res.addBinBoundary(625);
      res.addBinBoundary(700);
      res.addBinBoundary(800);
      res.addBinBoundary(900);
      res.addBinBoundary(1000);
      res.addBinBoundary(1200);
      res.addBinBoundary(1600);
      res.addBinBoundary(2000);
      res.addBinBoundary(3000);
      break;
    }
    case CategorizationHelpers::kOnshell:
    {
      double vlow=105, vhigh=140;
      if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
        unsigned int nbins = 35;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
      else if (category==CategorizationHelpers::JJVBFTagged || category==CategorizationHelpers::HadVHTagged){
        unsigned int const nbins = 7;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
      break;
    }
    default:
      assert(0);
    }
  }
  else if (strKD.Contains("int")){
    unsigned int nbins=20;
    if (
      (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged)
      &&
      (massregion==CategorizationHelpers::kOffshell || hypo==ACHypothesisHelpers::kSM)
      ) nbins += 10;
    double boundary=1; if (category==CategorizationHelpers::HadVHTagged && strKD.Contains(DiscriminantClasses::getKDName(DiscriminantClasses::kDintjjEWQCD))) boundary=0.4;
    double stepsize=2.*boundary/double(nbins);
    for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(-boundary+double(i)*stepsize);
  }
  else{
    unsigned int nbins=20;
    if (
      (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged)
      &&
      (massregion==CategorizationHelpers::kOffshell || hypo==ACHypothesisHelpers::kSM)
      ) nbins += 10;
    double stepsize=1./double(nbins);
    for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(double(i)*stepsize);
  }
  return res;
}
ExtendedBinning TemplateHelpers::getDiscriminantCoarseBinning(const SampleHelpers::Channel /*channel*/, const CategorizationHelpers::Category category, ACHypothesisHelpers::ACHypothesis hypo, TString const strKD, CategorizationHelpers::MassRegion const massregion){
  ExtendedBinning res(strKD);
  if (strKD=="ZZMass"){
    switch (massregion){
    case CategorizationHelpers::kOffshell:
    {
      res.addBinBoundary(220);
      res.addBinBoundary(theSqrts*1000.);
      if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
        res.addBinBoundary(230);
        res.addBinBoundary(240);
        res.addBinBoundary(250);
        res.addBinBoundary(260);
        res.addBinBoundary(280);
        res.addBinBoundary(310);
        res.addBinBoundary(350);
        res.addBinBoundary(430);
        res.addBinBoundary(590);
        res.addBinBoundary(710);
        res.addBinBoundary(1000);
      }
      else if (category==CategorizationHelpers::JJVBFTagged || category==CategorizationHelpers::HadVHTagged){
        res.addBinBoundary(260);
        res.addBinBoundary(310);
        res.addBinBoundary(400);
        res.addBinBoundary(750);
      }
      else{
        MELAerr << "TemplateHelpers::getDiscriminantCoarseBinning: Category " << CategorizationHelpers::getCategoryName(category)  << " not yet implemented!" << endl;
        assert(0);
      }
      break;
    }
    case CategorizationHelpers::kOnshell:
    {
      double vlow=105, vhigh=140;
      if (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged){
        unsigned int nbins = 7;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
      else if (category==CategorizationHelpers::JJVBFTagged || category==CategorizationHelpers::HadVHTagged){
        unsigned int const nbins = 5;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
      break;
    }
    default:
      assert(0);
    }
  }
  else if (strKD.Contains("int")){
    unsigned int nbins=10;
    if (
      (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged)
      &&
      (massregion==CategorizationHelpers::kOffshell || hypo==ACHypothesisHelpers::kSM)
      ) nbins += 10;
    double boundary=1; if (category==CategorizationHelpers::HadVHTagged && strKD.Contains(DiscriminantClasses::getKDName(DiscriminantClasses::kDintjjEWQCD))) boundary=0.4;
    double stepsize=2.*boundary/double(nbins);
    for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(-boundary+double(i)*stepsize);
  }
  else{
    unsigned int nbins=10;
    if (
      (category==CategorizationHelpers::Inclusive || category==CategorizationHelpers::Untagged)
      &&
      (massregion==CategorizationHelpers::kOffshell || hypo==ACHypothesisHelpers::kSM)
      ) nbins += 10;
    double stepsize=1./double(nbins);
    for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(double(i)*stepsize);
  }
  return res;
}
float TemplateHelpers::getDiscriminantSmearingStrengthCoefficient(CategorizationHelpers::Category category, ACHypothesisHelpers::ACHypothesis hypo, TString KDname, ProcessHandler::ProcessType proctype, CategorizationHelpers::MassRegion massregion){
  float res=0;
  if (massregion==CategorizationHelpers::kOnshell){
    if (KDname=="ZZMass") res=2;
    else{
      DiscriminantClasses::Type KDtype=DiscriminantClasses::getKDType(KDname);
      switch (KDtype){
      case DiscriminantClasses::kDL1dec:
      case DiscriminantClasses::kDL1decint:
      case DiscriminantClasses::kCL1decint:
      case DiscriminantClasses::kDa2dec:
      case DiscriminantClasses::kDa2decint:
      case DiscriminantClasses::kCa2decint:
      case DiscriminantClasses::kDa3dec:
      case DiscriminantClasses::kDa3decint:
      case DiscriminantClasses::kCa3decint:
      case DiscriminantClasses::kDL1ZGsdec:
      case DiscriminantClasses::kDL1ZGsdecint:
      case DiscriminantClasses::kCL1ZGsdecint:
        res=(proctype==ProcessHandler::kZX ? 2 : 1);
        break;
      default:
        res=2;
        break;
      }
    }
    // Cross-category production modes should have large smaring at on-shell in L1
    if (
      hypo==ACHypothesisHelpers::kL1 && (
        ((proctype==ProcessHandler::kZH || proctype==ProcessHandler::kWH) && category==CategorizationHelpers::JJVBFTagged)
        ||
        (proctype==ProcessHandler::kVBF && category==CategorizationHelpers::HadVHTagged)
        )
      ) res *= 5;
  }
  else res=2;
  return res;
}

ProcessHandler const* TemplateHelpers::getOnshellProcessHandler(ProcessHandler::ProcessType type){
  switch (type){
  case ProcessHandler::kGG:
    return &TemplateHelpers::OnshellGGProcessHandle;
  case ProcessHandler::kVV:
    return &TemplateHelpers::OnshellVVProcessHandle;
  case ProcessHandler::kVBF:
    return &TemplateHelpers::OnshellVBFProcessHandle;
  case ProcessHandler::kZH:
    return &TemplateHelpers::OnshellZHProcessHandle;
  case ProcessHandler::kWH:
    return &TemplateHelpers::OnshellWHProcessHandle;
  case ProcessHandler::kTT:
    return &TemplateHelpers::OnshellTTProcessHandle;
  case ProcessHandler::kBB:
    return &TemplateHelpers::OnshellBBProcessHandle;
  case ProcessHandler::kQQBkg:
    return &TemplateHelpers::OnshellQQBkgProcessHandle;
  case ProcessHandler::kZX:
    return &TemplateHelpers::OnshellZXProcessHandle;
  default:
    return nullptr;
  }
}
ProcessHandler const* TemplateHelpers::getOffshellProcessHandler(ProcessHandler::ProcessType type){
  switch (type){
  case ProcessHandler::kGG:
    return &TemplateHelpers::OffshellGGProcessHandle;
  case ProcessHandler::kVV:
    return &TemplateHelpers::OffshellVVProcessHandle;
  case ProcessHandler::kVBF:
    return &TemplateHelpers::OffshellVBFProcessHandle;
  case ProcessHandler::kZH:
    return &TemplateHelpers::OffshellZHProcessHandle;
  case ProcessHandler::kWH:
    return &TemplateHelpers::OffshellWHProcessHandle;
  case ProcessHandler::kTT:
    return &TemplateHelpers::OffshellTTProcessHandle;
  case ProcessHandler::kBB:
    return &TemplateHelpers::OffshellBBProcessHandle;
  case ProcessHandler::kQQBkg:
    return &TemplateHelpers::OffshellQQBkgProcessHandle;
  case ProcessHandler::kZX:
    return &TemplateHelpers::OffshellZXProcessHandle;
  default:
    return nullptr;
  }
}
ProcessHandler const* TemplateHelpers::getProcessHandlerPerMassRegion(ProcessHandler::ProcessType type, CategorizationHelpers::MassRegion massregion){
  switch (massregion){
  case CategorizationHelpers::kOffshell: return getOffshellProcessHandler(type);
  case CategorizationHelpers::kOnshell: return getOnshellProcessHandler(type);
  default:
    MELAerr << "TemplateHelpers::getProcessHandlerPerMassRegion: Analysis region is not implemented." << endl;
    assert(0);
    return nullptr;
  }
}


/******************/
/* Z+X background */
/******************/
ZXFakeRateHandler* TemplateHelpers::getFakeRateHandler(ZXFakeRateHandler::FakeRateMethod FRMethod, SystematicsHelpers::SystematicVariationTypes syst){
  signed char useUpDn = 1*(syst==SystematicsHelpers::eZXStatsUp) - 1*(syst==SystematicsHelpers::eZXStatsDn);
  TString fname="../data/FakeRate_";
  fname += ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod) + "_";
  fname += theDataPeriod + ".root";
  return new ZXFakeRateHandler(fname, FRMethod, useUpDn);
}
