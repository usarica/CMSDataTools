#ifndef EXT_CATEGORY_H
#define EXT_CATEGORY_H

#include "external_Discriminants.h"

#include <iostream>
#include <cmath>
#include "TMath.h"
#include "TRandom3.h"

using namespace std;

#ifndef CAT_VERBOSE
#define CAT_VERBOSE true
#endif

//---------- RunI categorization 

enum CategoryLegacy {
  ZeroOneJet = 0,
  Dijet      = 1
};
int categoryLegacy(int nCleanedJetsPt30){
  if (CAT_VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryLegacy'" << endl;

  if (nCleanedJetsPt30>=2)
    return Dijet;
  else
    return ZeroOneJet;
}


//---------- Moriond 2016 categorization 

enum CategoryMor16 {
  UntaggedMor16  = 0,
  VBFTaggedMor16 = 1
};
int categoryMor16(
  int nCleanedJetsPt30,
  float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
  ){
  if (CAT_VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryMor16'" << endl;

  float vbfMela = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal / (p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal + p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);

  if (nCleanedJetsPt30>=2 && vbfMela>0.5)
    return VBFTaggedMor16;
  else
    return UntaggedMor16;
}



//---------- ICHEP 2016 categorization

enum CategoryIchep16 {
  UntaggedIchep16     = 0,
  VBF1jTaggedIchep16  = 1,
  VBF2jTaggedIchep16  = 2,
  VHLeptTaggedIchep16 = 3,
  VHHadrTaggedIchep16 = 4,
  ttHTaggedIchep16    = 5
};
int categoryIchep16(
  int nExtraLep,
  int nExtraZ,
  int nCleanedJetsPt30,
  int nCleanedJetsPt30BTagged_bTagSF,
  float* jetQGLikelihood,
  float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
  float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
  float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
  float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
  float* jetPhi,
  float ZZMass,
  bool useQGTagging
  ){
  if (CAT_VERBOSE) cout << "WARNING: using deprecated categorization function 'categoryIchep16'" << endl;

  float D_VBF2j = -2;
  float D_VBF1j = -2;
  float D_WHh   = -2;
  float D_ZHh   = -2;
  if (useQGTagging){
    if (nCleanedJetsPt30==1)
      D_VBF1j = DVBF1j_ME_QG(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    else if (nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }
  else{
    if (nCleanedJetsPt30==1)
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    else if (nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  if (nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j)
    return VBF1jTaggedIchep16;
  else if (nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j)
    return VBF2jTaggedIchep16;
  else if ((nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh))
    || (nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3) && nCleanedJetsPt30BTagged_bTagSF>=2))
    return VHHadrTaggedIchep16;
  else if ((nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1))
    || (nCleanedJetsPt30==0 && nExtraLep>=1))
    return VHLeptTaggedIchep16;
  else if ((nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1)
    || nExtraLep>=1)
    return ttHTaggedIchep16;
  else
    return UntaggedIchep16;
}


//---------- Moriond 2017 categorization

enum CategoryMor17{
  InclusiveMor17,
  UntaggedMor17,
  VBF2jTaggedMor17,
  VHHadrTaggedMor17,
  VHLeptTaggedMor17,
  ttHTaggedMor17,
  VHMETTaggedMor17,
  VBF1jTaggedMor17,
  nCategoriesMor17
};
TString nameCategoryMor17(int icat){
  if (icat==(int)CategoryMor17::InclusiveMor17) return TString("InclusiveMor17");
  else if (icat==(int)CategoryMor17::UntaggedMor17) return TString("UntaggedMor17");
  else if (icat==(int)CategoryMor17::VBF2jTaggedMor17) return TString("VBF2jTaggedMor17");
  else if (icat==(int)CategoryMor17::VHHadrTaggedMor17) return TString("VHHadrTaggedMor17");
  else if (icat==(int)CategoryMor17::VHLeptTaggedMor17) return TString("VHLeptTaggedMor17");
  else if (icat==(int)CategoryMor17::ttHTaggedMor17) return TString("ttHTaggedMor17");
  else if (icat==(int)CategoryMor17::VHMETTaggedMor17) return TString("VHMETTaggedMor17");
  else if (icat==(int)CategoryMor17::VBF1jTaggedMor17) return TString("VBF1jTaggedMor17");
  else return "";
}
CategoryMor17 categoryMor17(
  int nExtraLep,
  int nExtraZ,
  int nCleanedJetsPt30,
  int nCleanedJetsPt30BTagged_bTagSF,
  float* jetQGLikelihood,
  float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
  float p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
  float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
  float p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,
  float* jetPhi,
  float ZZMass,
  float PFMET,
  bool useVHMETTagged,
  bool useQGTagging
  ){

  float D_VBF2j = -2;
  float D_VBF1j = -2;
  float D_WHh   = -2;
  float D_ZHh   = -2;
  if (useQGTagging){
    if (nCleanedJetsPt30==1)
      D_VBF1j = DVBF1j_ME_QG(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    else if (nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_WHh   = DWHh_ME_QG(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
      D_ZHh   = DZHh_ME_QG(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
    }
  }
  else{
    if (nCleanedJetsPt30==1)
      D_VBF1j = DVBF1j_ME(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    else if (nCleanedJetsPt30>=2){
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_WHh   = DWHh_ME(p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
      D_ZHh   = DZHh_ME(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
    }
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);
  float WP_VBF1j = getDVBF1jetWP(ZZMass, useQGTagging);
  float WP_WHh = getDWHhWP(ZZMass, useQGTagging);
  float WP_ZHh = getDZHhWP(ZZMass, useQGTagging);

  if (nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j)
    return VBF2jTaggedMor17;
  else if (nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && (D_WHh>WP_WHh||D_ZHh>WP_ZHh))
    return VHHadrTaggedMor17;
  else if ((nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged_bTagSF==0 && (nExtraLep==1||nExtraZ>=1))
    || (nCleanedJetsPt30==0 && nExtraLep>=1))
    return VHLeptTaggedMor17;
  else if ((nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged_bTagSF>=1)
    || nExtraLep>=1)
    return ttHTaggedMor17;
  else if (useVHMETTagged && nExtraLep==0 && (nCleanedJetsPt30==0||nCleanedJetsPt30==1) && PFMET>100)
    return VHMETTaggedMor17;
  else if (nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j)
    return VBF1jTaggedMor17;
  else
    return UntaggedMor17;
}

int categoryMor17_LegacyStyle(
  int nExtraLep,
  int nCleanedJetsPt30,
  int nCleanedJetsPt30BTagged_bTagSF,
  float* jetQGLikelihood,
  float p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
  float p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
  float* jetPhi,
  float ZZMass,
  bool useQGTagging
  ){

  float D_VBF2j = -2;
  if (useQGTagging){
    if (nCleanedJetsPt30>=2)
      D_VBF2j = DVBF2j_ME_QG(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass, jetQGLikelihood, jetPhi);
  }
  else{
    if (nCleanedJetsPt30>=2)
      D_VBF2j = DVBF2j_ME(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, ZZMass);
  }

  float WP_VBF2j = getDVBF2jetsWP(ZZMass, useQGTagging);

  if (nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0)) && D_VBF2j>WP_VBF2j)
    return VBF2jTaggedMor17;
  else
    return UntaggedMor17;
}


#endif
