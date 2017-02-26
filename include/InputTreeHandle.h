#ifndef INPUTTREEHANDLE_H
#define INPUTTREEHANDLE_H

#include <cmath>
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"

class InputTreeHandle{
protected:
  int theType;

  TFile* theFile;
  TTree* theTree;
  TH1F* theCounters;

  void SetBranchAddresses();
  void InitVectors_FirstTime();

public:
  Int_t RunNumber;
  Long64_t EventNumber;
  Int_t LumiNumber;
  Short_t NRecoMu;
  Short_t NRecoEle;
  Short_t Nvtx;
  Short_t NObsInt;
  Float_t NTrueInt;
  Float_t PUWeight;
  Float_t PUWeight_Up;
  Float_t PUWeight_Dn;

  Float_t KFactor_QCD_ggZZ_Nominal;
  Float_t KFactor_QCD_ggZZ_PDFScaleDn;
  Float_t KFactor_QCD_ggZZ_PDFScaleUp;
  Float_t KFactor_QCD_ggZZ_QCDScaleDn;
  Float_t KFactor_QCD_ggZZ_QCDScaleUp;
  Float_t KFactor_QCD_ggZZ_AsDn;
  Float_t KFactor_QCD_ggZZ_AsUp;
  Float_t KFactor_QCD_ggZZ_PDFReplicaDn;
  Float_t KFactor_QCD_ggZZ_PDFReplicaUp;
  Float_t KFactor_EW_qqZZ;
  Float_t KFactor_EW_qqZZ_unc;
  Float_t KFactor_QCD_qqZZ_dPhi;
  Float_t KFactor_QCD_qqZZ_M;
  Float_t KFactor_QCD_qqZZ_Pt;
  Float_t PFMET;
  Float_t PFMET_jesUp;
  Float_t PFMET_jesDn;
  Float_t PFMETPhi;
  Float_t PFMETNoHF;
  Float_t PFMETNoHFPhi;
  Short_t nCleanedJets;
  Short_t nCleanedJetsPt30;
  Short_t nCleanedJetsPt30_jecUp;
  Short_t nCleanedJetsPt30_jecDn;
  Short_t nCleanedJetsPt30BTagged;
  Short_t nCleanedJetsPt30BTagged_bTagSF;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn;
  Short_t trigWord;
  Float_t ZZMass;
  Float_t ZZMassErr;
  Float_t ZZMassErrCorr;
  Float_t ZZMassPreFSR;
  Short_t ZZsel;
  Float_t ZZPt;
  Float_t ZZEta;
  Float_t ZZPhi;
  Int_t CRflag;
  Float_t Z1Mass;
  Float_t Z1Pt;
  Short_t Z1Flav;
  Float_t ZZMassRefit;
  Float_t ZZMassRefitErr;
  Float_t ZZMassUnrefitErr;
  Float_t ZZMassCFit;
  Float_t ZZChi2CFit;
  Float_t Z2Mass;
  Float_t Z2Pt;
  Short_t Z2Flav;
  Float_t costhetastar;
  Float_t helphi;
  Float_t helcosthetaZ1;
  Float_t helcosthetaZ2;
  Float_t phistarZ1;
  Float_t phistarZ2;
  Float_t xi;
  Float_t xistar;
  Float_t TLE_dR_Z; // Delta-R between a TLE and the Z it does not belong to.
  Float_t TLE_min_dR_3l; // Minimum DR between a TLE and any of the other leptons

  std::vector<float>* LepPt;
  std::vector<float>* LepEta;
  std::vector<float>* LepPhi;
  std::vector<short>* LepLepId;
  std::vector<float>* LepSIP;
  std::vector<float>* LepTime;
  std::vector<bool>* LepisID;
  std::vector<float>* LepBDT;
  std::vector<char>* LepMissingHit;
  std::vector<float>* LepCombRelIsoPF;
  std::vector<short>* LepisLoose;
  std::vector<float>* LepRecoSF;
  std::vector<float>* LepRecoSF_Unc;
  std::vector<float>* LepSelSF;
  std::vector<float>* LepSelSF_Unc;


  std::vector<float>* fsrPt;
  std::vector<float>* fsrEta;
  std::vector<float>* fsrPhi;
  std::vector<float>* fsrDR;
  std::vector<short>* fsrLept;
  std::vector<short>* fsrLeptID;
  std::vector<float>* fsrGenPt;
  Bool_t passIsoPreFSR;

  std::vector<float>* JetPt;
  std::vector<float>* JetEta;
  std::vector<float>* JetPhi;
  std::vector<float>* JetMass;
  std::vector<float>* JetBTagger;
  std::vector<float>* JetIsBtagged;
  std::vector<float>* JetIsBtaggedWithSF;
  std::vector<float>* JetIsBtaggedWithSFUp;
  std::vector<float>* JetIsBtaggedWithSFDn;
  std::vector<float>* JetQGLikelihood;
  std::vector<float>* JetAxis2;
  std::vector<float>* JetMult;
  std::vector<float>* JetPtD;
  std::vector<float>* JetSigma;
  std::vector<short>* JetHadronFlavour;
  std::vector<short>* JetPartonFlavour;

  std::vector<float>* JetPUValue;
  std::vector<short>* JetPUID;

  std::vector<float>* JetJERUp;
  std::vector<float>* JetJERDown;

  Float_t DiJetMass;
  Float_t DiJetDEta;
  Float_t DiJetFisher;
  Short_t nExtraLep;
  Short_t nExtraZ;
  std::vector<float>* ExtraLepPt;
  std::vector<float>* ExtraLepEta;
  std::vector<float>* ExtraLepPhi;
  std::vector<short>* ExtraLepLepId;
  Short_t genFinalState;
  Int_t genProcessId;
  Float_t genHEPMCweight;

  std::vector<float>* LHEMotherPz;
  std::vector<float>* LHEMotherE;
  std::vector<short>* LHEMotherId;
  std::vector<float>* LHEDaughterPt;
  std::vector<float>* LHEDaughterEta;
  std::vector<float>* LHEDaughterPhi;
  std::vector<float>* LHEDaughterMass;
  std::vector<short>* LHEDaughterId;
  std::vector<float>* LHEAssociatedParticlePt;
  std::vector<float>* LHEAssociatedParticleEta;
  std::vector<float>* LHEAssociatedParticlePhi;
  std::vector<float>* LHEAssociatedParticleMass;
  std::vector<short>* LHEAssociatedParticleId;

  Float_t LHEPDFScale;
  Float_t LHEweight_QCDscale_muR1_muF1;
  Float_t LHEweight_QCDscale_muR1_muF2;
  Float_t LHEweight_QCDscale_muR1_muF0p5;
  Float_t LHEweight_QCDscale_muR2_muF1;
  Float_t LHEweight_QCDscale_muR2_muF2;
  Float_t LHEweight_QCDscale_muR2_muF0p5;
  Float_t LHEweight_QCDscale_muR0p5_muF1;
  Float_t LHEweight_QCDscale_muR0p5_muF2;
  Float_t LHEweight_QCDscale_muR0p5_muF0p5;
  Float_t LHEweight_PDFVariation_Up;
  Float_t LHEweight_PDFVariation_Dn;
  Float_t LHEweight_AsMZ_Up;
  Float_t LHEweight_AsMZ_Dn;

  Short_t genExtInfo;
  Float_t xsec;
  Float_t dataMCWeight;
  Float_t trigEffWeight;
  Float_t HqTMCweight;
  Float_t ZXFakeweight;
  Float_t overallEventWeight;
  Float_t GenHMass;
  Float_t GenHPt;
  Float_t GenHRapidity;
  Float_t GenZ1Mass;
  Float_t GenZ1Eta;
  Float_t GenZ1Pt;
  Float_t GenZ1Phi;
  Float_t GenZ1Flav;
  Float_t GenZ2Mass;
  Float_t GenZ2Eta;
  Float_t GenZ2Pt;
  Float_t GenZ2Phi;
  Float_t GenZ2Flav;
  Float_t GenLep1Pt;
  Float_t GenLep1Eta;
  Float_t GenLep1Phi;
  Short_t GenLep1Id;
  Float_t GenLep2Pt;
  Float_t GenLep2Eta;
  Float_t GenLep2Phi;
  Short_t GenLep2Id;
  Float_t GenLep3Pt;
  Float_t GenLep3Eta;
  Float_t GenLep3Phi;
  Short_t GenLep3Id;
  Float_t GenLep4Pt;
  Float_t GenLep4Eta;
  Float_t GenLep4Phi;
  Short_t GenLep4Id;
  Float_t GenAssocLep1Pt;
  Float_t GenAssocLep1Eta;
  Float_t GenAssocLep1Phi;
  Short_t GenAssocLep1Id;
  Float_t GenAssocLep2Pt;
  Float_t GenAssocLep2Eta;
  Float_t GenAssocLep2Phi;
  Short_t GenAssocLep2Id;

  Float_t p_GG_SIG_ghg2_1_ghz1_1_JHUGen;
  Float_t p_QQB_BKG_MCFM;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECNominal;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECNominal;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECNominal;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECUp;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECUp;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECUp;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECUp;
  Float_t p_JJVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t p_JJQCD_SIG_ghg2_1_JHUGen_JECDn;
  Float_t p_JVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn;
  Float_t p_JQCD_SIG_ghg2_1_JHUGen_JECDn;
  Float_t p_HadWH_SIG_ghw1_1_JHUGen_JECDn;
  Float_t p_HadZH_SIG_ghz1_1_JHUGen_JECDn;

  Float_t p_Gen_CPStoBWPropRewgt;
  Float_t p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM;
  Float_t p_Gen_GG_BKG_MCFM;
  Float_t p_Gen_QQB_BKG_MCFM;

  InputTreeHandle(const TString inputFileName, const TString treeName, const TString countersName, int type=-1);
  virtual ~InputTreeHandle();

  void InitDefaults();

  float GetNGenEvents();
  int GetSampleType()const;
  int GetEntries()const;
  void GetEntry(int ev);

};

InputTreeHandle::InputTreeHandle(const TString inputFileName, const TString treeName, const TString countersName, int type) : theType(type), theFile(0), theTree(0), theCounters(0){
  InitVectors_FirstTime();
  InitDefaults();
  if (inputFileName!="" && treeName!=""){
    theFile = TFile::Open(inputFileName, "read");
    if (theFile!=0 && theFile->IsZombie()){
      if (theFile->IsOpen()) theFile->Close();
      delete theFile;
      theFile=0;
      std::cerr << "InputTreeHandle::InputTreeHandle: File " << inputFileName << " could not be opened." << std::endl;
    }
    else{
      theTree = (TTree*)theFile->Get(treeName);
      if (theTree==0) std::cerr << "InputTreeHandle::InputTreeHandle: Tree " << treeName << " DNE." << std::endl;
      if (countersName!=""){
        theCounters = (TH1F*)theFile->Get(countersName);
        if (theCounters==0) std::cerr << "InputTreeHandle::InputTreeHandle: Counters " << countersName << " DNE." << std::endl;
      }
      SetBranchAddresses();
    }
  }
  else std::cerr << "InputTreeHandle::InputTreeHandle: File name " << inputFileName << " or tree name " << treeName << " is empty." << std::endl;
}

InputTreeHandle::~InputTreeHandle(){
  if (theFile->IsOpen()) theFile->Close();
}

void InputTreeHandle::InitDefaults(){
  RunNumber=0;
  EventNumber=0;
  LumiNumber=0;
  NRecoMu=0;
  NRecoEle=0;
  Nvtx=0;
  NObsInt=0;
  NTrueInt=0;
  PUWeight=1;
  PUWeight_Up=1;
  PUWeight_Dn=1;

  KFactor_QCD_ggZZ_Nominal=1;
  KFactor_QCD_ggZZ_PDFScaleDn=1;
  KFactor_QCD_ggZZ_PDFScaleUp=1;
  KFactor_QCD_ggZZ_QCDScaleDn=1;
  KFactor_QCD_ggZZ_QCDScaleUp=1;
  KFactor_QCD_ggZZ_AsDn=1;
  KFactor_QCD_ggZZ_AsUp=1;
  KFactor_QCD_ggZZ_PDFReplicaDn=1;
  KFactor_QCD_ggZZ_PDFReplicaUp=1;
  KFactor_EW_qqZZ=1;
  KFactor_EW_qqZZ_unc=0;
  KFactor_QCD_qqZZ_dPhi=1;
  KFactor_QCD_qqZZ_M=1;
  KFactor_QCD_qqZZ_Pt=1;

  PFMET=-99;
  PFMET_jesUp=-99;
  PFMET_jesDn=-99;
  PFMETPhi=-99;
  PFMETNoHF=-99;
  PFMETNoHFPhi=-99;

  nCleanedJets=0;
  nCleanedJetsPt30=0;
  nCleanedJetsPt30_jecUp=0;
  nCleanedJetsPt30_jecDn=0;
  nCleanedJetsPt30BTagged=0;
  nCleanedJetsPt30BTagged_bTagSF=0;
  nCleanedJetsPt30BTagged_bTagSFUp=0;
  nCleanedJetsPt30BTagged_bTagSFDn=0;

  trigWord=0;

  ZZMass=0;
  ZZMassErr=0;
  ZZMassErrCorr=0;
  ZZMassPreFSR=0;
  ZZsel=0;
  ZZPt=0;
  ZZEta=0;
  ZZPhi=0;
  CRflag=0;
  Z1Mass=0;
  Z1Pt=0;
  Z1Flav=0;
  ZZMassRefit=0;
  ZZMassRefitErr=0;
  ZZMassUnrefitErr=0;
  ZZMassCFit=0;
  ZZChi2CFit=0;
  Z2Mass=0;
  Z2Pt=0;
  Z2Flav=0;
  costhetastar=0;
  helphi=0;
  helcosthetaZ1=0;
  helcosthetaZ2=0;
  phistarZ1=0;
  phistarZ2=0;
  xi=0;
  xistar=0;
  TLE_dR_Z=-1;
  TLE_min_dR_3l=999;

  passIsoPreFSR=0;

  DiJetMass=-99;
  //DiJetMassPlus=-99;
  //DiJetMassMinus=-99;
  DiJetDEta=-99;
  DiJetFisher=-99;
  nExtraLep=0;
  nExtraZ=0;
  genFinalState=0;
  genProcessId=0;
  genHEPMCweight=0;

  LHEPDFScale=1;
  LHEweight_QCDscale_muR1_muF1=1;
  LHEweight_QCDscale_muR1_muF2=1;
  LHEweight_QCDscale_muR1_muF0p5=1;
  LHEweight_QCDscale_muR2_muF1=1;
  LHEweight_QCDscale_muR2_muF2=1;
  LHEweight_QCDscale_muR2_muF0p5=1;
  LHEweight_QCDscale_muR0p5_muF1=1;
  LHEweight_QCDscale_muR0p5_muF2=1;
  LHEweight_QCDscale_muR0p5_muF0p5=1;
  LHEweight_PDFVariation_Up=1;
  LHEweight_PDFVariation_Dn=1;
  LHEweight_AsMZ_Up=1;
  LHEweight_AsMZ_Dn=1;

  genExtInfo=0;
  xsec=1;
  dataMCWeight=1;
  trigEffWeight=1;
  HqTMCweight=1;
  ZXFakeweight=1;
  overallEventWeight=1;
  GenHMass=-1;
  GenHPt=0;
  GenHRapidity=0;
  GenZ1Mass=0;
  GenZ1Eta=0;
  GenZ1Pt=0;
  GenZ1Phi=0;
  GenZ1Flav=0;
  GenZ2Mass=0;
  GenZ2Eta=0;
  GenZ2Pt=0;
  GenZ2Phi=0;
  GenZ2Flav=0;
  GenLep1Pt=0;
  GenLep1Eta=0;
  GenLep1Phi=0;
  GenLep1Id=0;
  GenLep2Pt=0;
  GenLep2Eta=0;
  GenLep2Phi=0;
  GenLep2Id=0;
  GenLep3Pt=0;
  GenLep3Eta=0;
  GenLep3Phi=0;
  GenLep3Id=0;
  GenLep4Pt=0;
  GenLep4Eta=0;
  GenLep4Phi=0;
  GenLep4Id=0;
  GenAssocLep1Pt=0;
  GenAssocLep1Eta=0;
  GenAssocLep1Phi=0;
  GenAssocLep1Id=0;
  GenAssocLep2Pt=0;
  GenAssocLep2Eta=0;
  GenAssocLep2Phi=0;
  GenAssocLep2Id=0;

  p_GG_SIG_ghg2_1_ghz1_1_JHUGen=-1;
  p_QQB_BKG_MCFM=-1;
  p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal=-1;
  p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal=-1;
  p_JVBF_SIG_ghv1_1_JHUGen_JECNominal=-1;
  pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal=-1;
  p_JQCD_SIG_ghg2_1_JHUGen_JECNominal=-1;
  p_HadWH_SIG_ghw1_1_JHUGen_JECNominal=-1;
  p_HadZH_SIG_ghz1_1_JHUGen_JECNominal=-1;
  p_JJVBF_SIG_ghv1_1_JHUGen_JECUp=-1;
  p_JJQCD_SIG_ghg2_1_JHUGen_JECUp=-1;
  p_JVBF_SIG_ghv1_1_JHUGen_JECUp=-1;
  pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp=-1;
  p_JQCD_SIG_ghg2_1_JHUGen_JECUp=-1;
  p_HadWH_SIG_ghw1_1_JHUGen_JECUp=-1;
  p_HadZH_SIG_ghz1_1_JHUGen_JECUp=-1;
  p_JJVBF_SIG_ghv1_1_JHUGen_JECDn=-1;
  p_JJQCD_SIG_ghg2_1_JHUGen_JECDn=-1;
  p_JVBF_SIG_ghv1_1_JHUGen_JECDn=-1;
  pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn=-1;
  p_JQCD_SIG_ghg2_1_JHUGen_JECDn=-1;
  p_HadWH_SIG_ghw1_1_JHUGen_JECDn=-1;
  p_HadZH_SIG_ghz1_1_JHUGen_JECDn=-1;

  p_Gen_CPStoBWPropRewgt=1;
  p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM=1;
  p_Gen_GG_BKG_MCFM=1;
  p_Gen_QQB_BKG_MCFM=1;

}

void InputTreeHandle::InitVectors_FirstTime(){
  LepPt=0;
  LepEta=0;
  LepPhi=0;
  LepLepId=0;
  LepSIP=0;
  LepTime=0;
  LepisID=0;
  LepBDT=0;
  LepMissingHit=0;
  LepCombRelIsoPF=0;
  LepisLoose=0;
  LepRecoSF=0;
  LepRecoSF_Unc=0;
  LepSelSF=0;
  LepSelSF_Unc=0;
  fsrPt=0;
  fsrEta=0;
  fsrPhi=0;
  fsrDR=0;
  fsrLept=0;
  fsrLeptID=0;
  fsrGenPt=0;
  JetPt=0;
  JetEta=0;
  JetPhi=0;
  JetMass=0;
  JetBTagger=0;
  JetIsBtagged=0;
  JetIsBtaggedWithSF=0;
  JetIsBtaggedWithSFUp=0;
  JetIsBtaggedWithSFDn=0;
  JetQGLikelihood=0;
  JetAxis2=0;
  JetMult=0;
  JetPtD=0;
  JetSigma=0;
  JetHadronFlavour=0;
  JetPartonFlavour=0;
  JetPUValue=0;
  JetPUID=0;
  JetJERUp=0;
  JetJERDown=0;
  ExtraLepPt=0;
  ExtraLepEta=0;
  ExtraLepPhi=0;
  ExtraLepLepId=0;
  LHEMotherPz=0;
  LHEMotherE=0;
  LHEMotherId=0;
  LHEDaughterPt=0;
  LHEDaughterEta=0;
  LHEDaughterPhi=0;
  LHEDaughterMass=0;
  LHEDaughterId=0;
  LHEAssociatedParticlePt=0;
  LHEAssociatedParticleEta=0;
  LHEAssociatedParticlePhi=0;
  LHEAssociatedParticleMass=0;
  LHEAssociatedParticleId=0;
}

void InputTreeHandle::SetBranchAddresses(){
  if (theTree==0) return;

  if (theTree->GetBranchStatus("RunNumber")) theTree->SetBranchAddress("RunNumber", &RunNumber);
  if (theTree->GetBranchStatus("EventNumber")) theTree->SetBranchAddress("EventNumber", &EventNumber);
  if (theTree->GetBranchStatus("LumiNumber")) theTree->SetBranchAddress("LumiNumber", &LumiNumber);
  if (theTree->GetBranchStatus("NRecoMu")) theTree->SetBranchAddress("NRecoMu", &NRecoMu);
  if (theTree->GetBranchStatus("NRecoEle")) theTree->SetBranchAddress("NRecoEle", &NRecoEle);
  if (theTree->GetBranchStatus("Nvtx")) theTree->SetBranchAddress("Nvtx", &Nvtx);
  if (theTree->GetBranchStatus("NObsInt")) theTree->SetBranchAddress("NObsInt", &NObsInt);
  if (theTree->GetBranchStatus("NTrueInt")) theTree->SetBranchAddress("NTrueInt", &NTrueInt);
  if (theTree->GetBranchStatus("PUWeight")) theTree->SetBranchAddress("PUWeight", &PUWeight);
  if (theTree->GetBranchStatus("PUWeight_Up")) theTree->SetBranchAddress("PUWeight_Up", &PUWeight_Up);
  if (theTree->GetBranchStatus("PUWeight_Dn")) theTree->SetBranchAddress("PUWeight_Dn", &PUWeight_Dn);

  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_Nominal")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_PDFScaleDn")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleDn", &KFactor_QCD_ggZZ_PDFScaleDn);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_PDFScaleUp")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleUp", &KFactor_QCD_ggZZ_PDFScaleUp);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_QCDScaleDn")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleDn", &KFactor_QCD_ggZZ_QCDScaleDn);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_QCDScaleUp")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleUp", &KFactor_QCD_ggZZ_QCDScaleUp);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_AsDn")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_AsDn", &KFactor_QCD_ggZZ_AsDn);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_AsUp")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_AsUp", &KFactor_QCD_ggZZ_AsUp);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_PDFReplicaDn")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaDn", &KFactor_QCD_ggZZ_PDFReplicaDn);
  if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_PDFReplicaUp")) theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaUp", &KFactor_QCD_ggZZ_PDFReplicaUp);
  if (theTree->GetBranchStatus("KFactor_EW_qqZZ")) theTree->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
  if (theTree->GetBranchStatus("KFactor_EW_qqZZ_unc")) theTree->SetBranchAddress("KFactor_EW_qqZZ_unc", &KFactor_EW_qqZZ_unc);
  if (theTree->GetBranchStatus("KFactor_QCD_qqZZ_dPhi")) theTree->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi);
  if (theTree->GetBranchStatus("KFactor_QCD_qqZZ_M")) theTree->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
  if (theTree->GetBranchStatus("KFactor_QCD_qqZZ_Pt")) theTree->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt);
  if (theTree->GetBranchStatus("PFMET")) theTree->SetBranchAddress("PFMET", &PFMET);
  if (theTree->GetBranchStatus("PFMET_jesUp")) theTree->SetBranchAddress("PFMET_jesUp", &PFMET_jesUp);
  if (theTree->GetBranchStatus("PFMET_jesDn")) theTree->SetBranchAddress("PFMET_jesDn", &PFMET_jesDn);
  if (theTree->GetBranchStatus("PFMETPhi")) theTree->SetBranchAddress("PFMETPhi", &PFMETPhi);
  if (theTree->GetBranchStatus("PFMETNoHF")) theTree->SetBranchAddress("PFMETNoHF", &PFMETNoHF);
  if (theTree->GetBranchStatus("PFMETNoHFPhi")) theTree->SetBranchAddress("PFMETNoHFPhi", &PFMETNoHFPhi);
  if (theTree->GetBranchStatus("nCleanedJets")) theTree->SetBranchAddress("nCleanedJets", &nCleanedJets);
  if (theTree->GetBranchStatus("nCleanedJetsPt30")) theTree->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30);
  if (theTree->GetBranchStatus("nCleanedJetsPt30_jecUp")) theTree->SetBranchAddress("nCleanedJetsPt30_jecUp", &nCleanedJetsPt30_jecUp);
  if (theTree->GetBranchStatus("nCleanedJetsPt30_jecDn")) theTree->SetBranchAddress("nCleanedJetsPt30_jecDn", &nCleanedJetsPt30_jecDn);
  if (theTree->GetBranchStatus("nCleanedJetsPt30BTagged")) theTree->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged);
  if (theTree->GetBranchStatus("nCleanedJetsPt30BTagged_bTagSF")) theTree->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nCleanedJetsPt30BTagged_bTagSF);
  if (theTree->GetBranchStatus("nCleanedJetsPt30BTagged_bTagSFUp")) theTree->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFUp", &nCleanedJetsPt30BTagged_bTagSFUp);
  if (theTree->GetBranchStatus("nCleanedJetsPt30BTagged_bTagSFDn")) theTree->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSFDn", &nCleanedJetsPt30BTagged_bTagSFDn);
  if (theTree->GetBranchStatus("trigWord")) theTree->SetBranchAddress("trigWord", &trigWord);
  if (theTree->GetBranchStatus("ZZMass")) theTree->SetBranchAddress("ZZMass", &ZZMass);
  if (theTree->GetBranchStatus("ZZMassErr")) theTree->SetBranchAddress("ZZMassErr", &ZZMassErr);
  if (theTree->GetBranchStatus("ZZMassErrCorr")) theTree->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr);
  if (theTree->GetBranchStatus("ZZMassPreFSR")) theTree->SetBranchAddress("ZZMassPreFSR", &ZZMassPreFSR);
  if (theTree->GetBranchStatus("ZZsel")) theTree->SetBranchAddress("ZZsel", &ZZsel);
  if (theTree->GetBranchStatus("ZZPt")) theTree->SetBranchAddress("ZZPt", &ZZPt);
  if (theTree->GetBranchStatus("ZZEta")) theTree->SetBranchAddress("ZZEta", &ZZEta);
  if (theTree->GetBranchStatus("ZZPhi")) theTree->SetBranchAddress("ZZPhi", &ZZPhi);
  if (theTree->GetBranchStatus("CRflag")) theTree->SetBranchAddress("CRflag", &CRflag);
  if (theTree->GetBranchStatus("Z1Mass")) theTree->SetBranchAddress("Z1Mass", &Z1Mass);
  if (theTree->GetBranchStatus("Z1Pt")) theTree->SetBranchAddress("Z1Pt", &Z1Pt);
  if (theTree->GetBranchStatus("Z1Flav")) theTree->SetBranchAddress("Z1Flav", &Z1Flav);
  if (theTree->GetBranchStatus("ZZMassRefit")) theTree->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
  if (theTree->GetBranchStatus("ZZMassRefitErr")) theTree->SetBranchAddress("ZZMassRefitErr", &ZZMassRefitErr);
  if (theTree->GetBranchStatus("ZZMassUnrefitErr")) theTree->SetBranchAddress("ZZMassUnrefitErr", &ZZMassUnrefitErr);
  if (theTree->GetBranchStatus("ZZMassCFit")) theTree->SetBranchAddress("ZZMassCFit", &ZZMassCFit);
  if (theTree->GetBranchStatus("ZZChi2CFit")) theTree->SetBranchAddress("ZZChi2CFit", &ZZChi2CFit);
  if (theTree->GetBranchStatus("Z2Mass")) theTree->SetBranchAddress("Z2Mass", &Z2Mass);
  if (theTree->GetBranchStatus("Z2Pt")) theTree->SetBranchAddress("Z2Pt", &Z2Pt);
  if (theTree->GetBranchStatus("Z2Flav")) theTree->SetBranchAddress("Z2Flav", &Z2Flav);
  if (theTree->GetBranchStatus("costhetastar")) theTree->SetBranchAddress("costhetastar", &costhetastar);
  if (theTree->GetBranchStatus("helphi")) theTree->SetBranchAddress("helphi", &helphi);
  if (theTree->GetBranchStatus("helcosthetaZ1")) theTree->SetBranchAddress("helcosthetaZ1", &helcosthetaZ1);
  if (theTree->GetBranchStatus("helcosthetaZ2")) theTree->SetBranchAddress("helcosthetaZ2", &helcosthetaZ2);
  if (theTree->GetBranchStatus("phistarZ1")) theTree->SetBranchAddress("phistarZ1", &phistarZ1);
  if (theTree->GetBranchStatus("phistarZ2")) theTree->SetBranchAddress("phistarZ2", &phistarZ2);
  if (theTree->GetBranchStatus("xi")) theTree->SetBranchAddress("xi", &xi);
  if (theTree->GetBranchStatus("xistar")) theTree->SetBranchAddress("xistar", &xistar);
  if (theTree->GetBranchStatus("TLE_dR_Z")) theTree->SetBranchAddress("TLE_dR_Z", &TLE_dR_Z);
  if (theTree->GetBranchStatus("TLE_min_dR_3l")) theTree->SetBranchAddress("TLE_min_dR_3l", &TLE_min_dR_3l);

  if (theTree->GetBranchStatus("LepPt")) theTree->SetBranchAddress("LepPt", &LepPt);
  if (theTree->GetBranchStatus("LepEta")) theTree->SetBranchAddress("LepEta", &LepEta);
  if (theTree->GetBranchStatus("LepPhi")) theTree->SetBranchAddress("LepPhi", &LepPhi);
  if (theTree->GetBranchStatus("LepLepId")) theTree->SetBranchAddress("LepLepId", &LepLepId);
  if (theTree->GetBranchStatus("LepSIP")) theTree->SetBranchAddress("LepSIP", &LepSIP);
  if (theTree->GetBranchStatus("LepTime")) theTree->SetBranchAddress("LepTime", &LepTime);
  if (theTree->GetBranchStatus("LepisID")) theTree->SetBranchAddress("LepisID", &LepisID);
  if (theTree->GetBranchStatus("LepBDT")) theTree->SetBranchAddress("LepBDT", &LepBDT);
  if (theTree->GetBranchStatus("LepMissingHit")) theTree->SetBranchAddress("LepMissingHit", &LepMissingHit);
  //if (theTree->GetBranchStatus("LepChargedHadIso")) theTree->SetBranchAddress("LepChargedHadIso", &LepChargedHadIso);
  //if (theTree->GetBranchStatus("LepNeutralHadIso")) theTree->SetBranchAddress("LepNeutralHadIso", &LepNeutralHadIso);
  //if (theTree->GetBranchStatus("LepPhotonIso")) theTree->SetBranchAddress("LepPhotonIso", &LepPhotonIso);
  if (theTree->GetBranchStatus("LepCombRelIsoPF")) theTree->SetBranchAddress("LepCombRelIsoPF", &LepCombRelIsoPF);
  if (theTree->GetBranchStatus("LepisLoose")) theTree->SetBranchAddress("LepisLoose", &LepisLoose);
  if (theTree->GetBranchStatus("LepRecoSF")) theTree->SetBranchAddress("LepRecoSF", &LepRecoSF);
  if (theTree->GetBranchStatus("LepRecoSF_Unc")) theTree->SetBranchAddress("LepRecoSF_Unc", &LepRecoSF_Unc);
  if (theTree->GetBranchStatus("LepSelSF")) theTree->SetBranchAddress("LepSelSF", &LepSelSF);
  if (theTree->GetBranchStatus("LepSelSF_Unc")) theTree->SetBranchAddress("LepSelSF_Unc", &LepSelSF_Unc);


  if (theTree->GetBranchStatus("fsrPt")) theTree->SetBranchAddress("fsrPt", &fsrPt);
  if (theTree->GetBranchStatus("fsrEta")) theTree->SetBranchAddress("fsrEta", &fsrEta);
  if (theTree->GetBranchStatus("fsrPhi")) theTree->SetBranchAddress("fsrPhi", &fsrPhi);
  if (theTree->GetBranchStatus("fsrDR")) theTree->SetBranchAddress("fsrDR", &fsrDR);
  if (theTree->GetBranchStatus("fsrLept")) theTree->SetBranchAddress("fsrLept", &fsrLept);
  if (theTree->GetBranchStatus("fsrLeptID")) theTree->SetBranchAddress("fsrLeptID", &fsrLeptID);
  if (theTree->GetBranchStatus("fsrGenPt")) theTree->SetBranchAddress("fsrGenPt", &fsrGenPt);
  if (theTree->GetBranchStatus("passIsoPreFSR")) theTree->SetBranchAddress("passIsoPreFSR", &passIsoPreFSR);

  if (theTree->GetBranchStatus("JetPt")) theTree->SetBranchAddress("JetPt", &JetPt);
  if (theTree->GetBranchStatus("JetEta")) theTree->SetBranchAddress("JetEta", &JetEta);
  if (theTree->GetBranchStatus("JetPhi")) theTree->SetBranchAddress("JetPhi", &JetPhi);
  if (theTree->GetBranchStatus("JetMass")) theTree->SetBranchAddress("JetMass", &JetMass);
  if (theTree->GetBranchStatus("JetBTagger")) theTree->SetBranchAddress("JetBTagger", &JetBTagger);
  if (theTree->GetBranchStatus("JetIsBtagged")) theTree->SetBranchAddress("JetIsBtagged", &JetIsBtagged);
  if (theTree->GetBranchStatus("JetIsBtaggedWithSF")) theTree->SetBranchAddress("JetIsBtaggedWithSF", &JetIsBtaggedWithSF);
  if (theTree->GetBranchStatus("JetIsBtaggedWithSFUp")) theTree->SetBranchAddress("JetIsBtaggedWithSFUp", &JetIsBtaggedWithSFUp);
  if (theTree->GetBranchStatus("JetIsBtaggedWithSFDn")) theTree->SetBranchAddress("JetIsBtaggedWithSFDn", &JetIsBtaggedWithSFDn);
  if (theTree->GetBranchStatus("JetQGLikelihood")) theTree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
  if (theTree->GetBranchStatus("JetAxis2")) theTree->SetBranchAddress("JetAxis2", &JetAxis2);
  if (theTree->GetBranchStatus("JetMult")) theTree->SetBranchAddress("JetMult", &JetMult);
  if (theTree->GetBranchStatus("JetPtD")) theTree->SetBranchAddress("JetPtD", &JetPtD);
  if (theTree->GetBranchStatus("JetSigma")) theTree->SetBranchAddress("JetSigma", &JetSigma);
  if (theTree->GetBranchStatus("JetHadronFlavour")) theTree->SetBranchAddress("JetHadronFlavour", &JetHadronFlavour);
  if (theTree->GetBranchStatus("JetPartonFlavour")) theTree->SetBranchAddress("JetPartonFlavour", &JetPartonFlavour);

  if (theTree->GetBranchStatus("JetPUValue")) theTree->SetBranchAddress("JetPUValue", &JetPUValue);
  if (theTree->GetBranchStatus("JetPUID")) theTree->SetBranchAddress("JetPUID", &JetPUID);

  if (theTree->GetBranchStatus("JetJERUp")) theTree->SetBranchAddress("JetJERUp", &JetJERUp);
  if (theTree->GetBranchStatus("JetJERDown")) theTree->SetBranchAddress("JetJERDown", &JetJERDown);

  if (theTree->GetBranchStatus("DiJetMass")) theTree->SetBranchAddress("DiJetMass", &DiJetMass);
  //if (theTree->GetBranchStatus("DiJetMassPlus")) theTree->SetBranchAddress("DiJetMassPlus", &DiJetMassPlus);
  //if (theTree->GetBranchStatus("DiJetMassMinus")) theTree->SetBranchAddress("DiJetMassMinus", &DiJetMassMinus);
  if (theTree->GetBranchStatus("DiJetDEta")) theTree->SetBranchAddress("DiJetDEta", &DiJetDEta);
  if (theTree->GetBranchStatus("DiJetFisher")) theTree->SetBranchAddress("DiJetFisher", &DiJetFisher);
  if (theTree->GetBranchStatus("nExtraLep")) theTree->SetBranchAddress("nExtraLep", &nExtraLep);
  if (theTree->GetBranchStatus("nExtraZ")) theTree->SetBranchAddress("nExtraZ", &nExtraZ);
  if (theTree->GetBranchStatus("ExtraLepPt")) theTree->SetBranchAddress("ExtraLepPt", &ExtraLepPt);
  if (theTree->GetBranchStatus("ExtraLepEta")) theTree->SetBranchAddress("ExtraLepEta", &ExtraLepEta);
  if (theTree->GetBranchStatus("ExtraLepPhi")) theTree->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi);
  if (theTree->GetBranchStatus("ExtraLepLepId")) theTree->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId);
  if (theTree->GetBranchStatus("genFinalState")) theTree->SetBranchAddress("genFinalState", &genFinalState);
  if (theTree->GetBranchStatus("genProcessId")) theTree->SetBranchAddress("genProcessId", &genProcessId);
  if (theTree->GetBranchStatus("genHEPMCweight")) theTree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);

  if (theTree->GetBranchStatus("LHEMotherPz")) theTree->SetBranchAddress("LHEMotherPz", &LHEMotherPz);
  if (theTree->GetBranchStatus("LHEMotherE")) theTree->SetBranchAddress("LHEMotherE", &LHEMotherE);
  if (theTree->GetBranchStatus("LHEMotherId")) theTree->SetBranchAddress("LHEMotherId", &LHEMotherId);
  if (theTree->GetBranchStatus("LHEDaughterPt")) theTree->SetBranchAddress("LHEDaughterPt", &LHEDaughterPt);
  if (theTree->GetBranchStatus("LHEDaughterEta")) theTree->SetBranchAddress("LHEDaughterEta", &LHEDaughterEta);
  if (theTree->GetBranchStatus("LHEDaughterPhi")) theTree->SetBranchAddress("LHEDaughterPhi", &LHEDaughterPhi);
  if (theTree->GetBranchStatus("LHEDaughterMass")) theTree->SetBranchAddress("LHEDaughterMass", &LHEDaughterMass);
  if (theTree->GetBranchStatus("LHEDaughterId")) theTree->SetBranchAddress("LHEDaughterId", &LHEDaughterId);
  if (theTree->GetBranchStatus("LHEAssociatedParticlePt")) theTree->SetBranchAddress("LHEAssociatedParticlePt", &LHEAssociatedParticlePt);
  if (theTree->GetBranchStatus("LHEAssociatedParticleEta")) theTree->SetBranchAddress("LHEAssociatedParticleEta", &LHEAssociatedParticleEta);
  if (theTree->GetBranchStatus("LHEAssociatedParticlePhi")) theTree->SetBranchAddress("LHEAssociatedParticlePhi", &LHEAssociatedParticlePhi);
  if (theTree->GetBranchStatus("LHEAssociatedParticleMass")) theTree->SetBranchAddress("LHEAssociatedParticleMass", &LHEAssociatedParticleMass);
  if (theTree->GetBranchStatus("LHEAssociatedParticleId")) theTree->SetBranchAddress("LHEAssociatedParticleId", &LHEAssociatedParticleId);

  if (theTree->GetBranchStatus("LHEPDFScale")) theTree->SetBranchAddress("LHEPDFScale", &LHEPDFScale);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR1_muF1")) theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF1", &LHEweight_QCDscale_muR1_muF1);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR1_muF2")) theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF2", &LHEweight_QCDscale_muR1_muF2);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR1_muF0p5")) theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5", &LHEweight_QCDscale_muR1_muF0p5);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR2_muF1")) theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF1", &LHEweight_QCDscale_muR2_muF1);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR2_muF2")) theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF2", &LHEweight_QCDscale_muR2_muF2);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR2_muF0p5")) theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5", &LHEweight_QCDscale_muR2_muF0p5);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR0p5_muF1")) theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1", &LHEweight_QCDscale_muR0p5_muF1);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR0p5_muF2")) theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2", &LHEweight_QCDscale_muR0p5_muF2);
  if (theTree->GetBranchStatus("LHEweight_QCDscale_muR0p5_muF0p5")) theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5);
  if (theTree->GetBranchStatus("LHEweight_PDFVariation_Up")) theTree->SetBranchAddress("LHEweight_PDFVariation_Up", &LHEweight_PDFVariation_Up);
  if (theTree->GetBranchStatus("LHEweight_PDFVariation_Dn")) theTree->SetBranchAddress("LHEweight_PDFVariation_Dn", &LHEweight_PDFVariation_Dn);
  if (theTree->GetBranchStatus("LHEweight_AsMZ_Up")) theTree->SetBranchAddress("LHEweight_AsMZ_Up", &LHEweight_AsMZ_Up);
  if (theTree->GetBranchStatus("LHEweight_AsMZ_Dn")) theTree->SetBranchAddress("LHEweight_AsMZ_Dn", &LHEweight_AsMZ_Dn);
  if (theTree->GetBranchStatus("genExtInfo")) theTree->SetBranchAddress("genExtInfo", &genExtInfo);
  if (theTree->GetBranchStatus("xsec")) theTree->SetBranchAddress("xsec", &xsec);
  if (theTree->GetBranchStatus("dataMCWeight")) theTree->SetBranchAddress("dataMCWeight", &dataMCWeight);
  if (theTree->GetBranchStatus("trigEffWeight")) theTree->SetBranchAddress("trigEffWeight", &trigEffWeight);
  if (theTree->GetBranchStatus("HqTMCweight")) theTree->SetBranchAddress("HqTMCweight", &HqTMCweight);
  if (theTree->GetBranchStatus("ZXFakeweight")) theTree->SetBranchAddress("ZXFakeweight", &ZXFakeweight);
  if (theTree->GetBranchStatus("overallEventWeight")) theTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
  if (theTree->GetBranchStatus("GenHMass")) theTree->SetBranchAddress("GenHMass", &GenHMass);
  if (theTree->GetBranchStatus("GenHPt")) theTree->SetBranchAddress("GenHPt", &GenHPt);
  if (theTree->GetBranchStatus("GenHRapidity")) theTree->SetBranchAddress("GenHRapidity", &GenHRapidity);
  if (theTree->GetBranchStatus("GenZ1Mass")) theTree->SetBranchAddress("GenZ1Mass", &GenZ1Mass);
  if (theTree->GetBranchStatus("GenZ1Eta")) theTree->SetBranchAddress("GenZ1Eta", &GenZ1Eta);
  if (theTree->GetBranchStatus("GenZ1Pt")) theTree->SetBranchAddress("GenZ1Pt", &GenZ1Pt);
  if (theTree->GetBranchStatus("GenZ1Phi")) theTree->SetBranchAddress("GenZ1Phi", &GenZ1Phi);
  if (theTree->GetBranchStatus("GenZ1Flav")) theTree->SetBranchAddress("GenZ1Flav", &GenZ1Flav);
  if (theTree->GetBranchStatus("GenZ2Mass")) theTree->SetBranchAddress("GenZ2Mass", &GenZ2Mass);
  if (theTree->GetBranchStatus("GenZ2Eta")) theTree->SetBranchAddress("GenZ2Eta", &GenZ2Eta);
  if (theTree->GetBranchStatus("GenZ2Pt")) theTree->SetBranchAddress("GenZ2Pt", &GenZ2Pt);
  if (theTree->GetBranchStatus("GenZ2Phi")) theTree->SetBranchAddress("GenZ2Phi", &GenZ2Phi);
  if (theTree->GetBranchStatus("GenZ2Flav")) theTree->SetBranchAddress("GenZ2Flav", &GenZ2Flav);
  if (theTree->GetBranchStatus("GenLep1Pt")) theTree->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
  if (theTree->GetBranchStatus("GenLep1Eta")) theTree->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
  if (theTree->GetBranchStatus("GenLep1Phi")) theTree->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
  if (theTree->GetBranchStatus("GenLep1Id")) theTree->SetBranchAddress("GenLep1Id", &GenLep1Id);
  if (theTree->GetBranchStatus("GenLep2Pt")) theTree->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
  if (theTree->GetBranchStatus("GenLep2Eta")) theTree->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
  if (theTree->GetBranchStatus("GenLep2Phi")) theTree->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
  if (theTree->GetBranchStatus("GenLep2Id")) theTree->SetBranchAddress("GenLep2Id", &GenLep2Id);
  if (theTree->GetBranchStatus("GenLep3Pt")) theTree->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
  if (theTree->GetBranchStatus("GenLep3Eta")) theTree->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
  if (theTree->GetBranchStatus("GenLep3Phi")) theTree->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
  if (theTree->GetBranchStatus("GenLep3Id")) theTree->SetBranchAddress("GenLep3Id", &GenLep3Id);
  if (theTree->GetBranchStatus("GenLep4Pt")) theTree->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
  if (theTree->GetBranchStatus("GenLep4Eta")) theTree->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
  if (theTree->GetBranchStatus("GenLep4Phi")) theTree->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
  if (theTree->GetBranchStatus("GenLep4Id")) theTree->SetBranchAddress("GenLep4Id", &GenLep4Id);
  if (theTree->GetBranchStatus("GenAssocLep1Pt")) theTree->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt);
  if (theTree->GetBranchStatus("GenAssocLep1Eta")) theTree->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta);
  if (theTree->GetBranchStatus("GenAssocLep1Phi")) theTree->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi);
  if (theTree->GetBranchStatus("GenAssocLep1Id")) theTree->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id);
  if (theTree->GetBranchStatus("GenAssocLep2Pt")) theTree->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt);
  if (theTree->GetBranchStatus("GenAssocLep2Eta")) theTree->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta);
  if (theTree->GetBranchStatus("GenAssocLep2Phi")) theTree->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi);
  if (theTree->GetBranchStatus("GenAssocLep2Id")) theTree->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id);

  if (theTree->GetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen")) theTree->SetBranchAddress("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", &p_GG_SIG_ghg2_1_ghz1_1_JHUGen);
  if (theTree->GetBranchStatus("p_QQB_BKG_MCFM")) theTree->SetBranchAddress("p_QQB_BKG_MCFM", &p_QQB_BKG_MCFM);
  if (theTree->GetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &p_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal")) theTree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &p_JQCD_SIG_ghg2_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal", &p_HadWH_SIG_ghw1_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal")) theTree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal", &p_HadZH_SIG_ghz1_1_JHUGen_JECNominal);
  if (theTree->GetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECUp", &p_JJVBF_SIG_ghv1_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECUp", &p_JJQCD_SIG_ghg2_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_JVBF_SIG_ghv1_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECUp", &p_JVBF_SIG_ghv1_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp")) theTree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECUp", &p_JQCD_SIG_ghg2_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_HadWH_SIG_ghw1_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECUp", &p_HadWH_SIG_ghw1_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_HadZH_SIG_ghz1_1_JHUGen_JECUp")) theTree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECUp", &p_HadZH_SIG_ghz1_1_JHUGen_JECUp);
  if (theTree->GetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECDn", &p_JJVBF_SIG_ghv1_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECDn", &p_JJQCD_SIG_ghg2_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("p_JVBF_SIG_ghv1_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_JVBF_SIG_ghv1_1_JHUGen_JECDn", &p_JVBF_SIG_ghv1_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn")) theTree->SetBranchAddress("pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn", &pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECDn", &p_JQCD_SIG_ghg2_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("p_HadWH_SIG_ghw1_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_HadWH_SIG_ghw1_1_JHUGen_JECDn", &p_HadWH_SIG_ghw1_1_JHUGen_JECDn);
  if (theTree->GetBranchStatus("p_HadZH_SIG_ghz1_1_JHUGen_JECDn")) theTree->SetBranchAddress("p_HadZH_SIG_ghz1_1_JHUGen_JECDn", &p_HadZH_SIG_ghz1_1_JHUGen_JECDn);

  if (theTree->GetBranchStatus("p_Gen_CPStoBWPropRewgt")) theTree->SetBranchAddress("p_Gen_CPStoBWPropRewgt", &p_Gen_CPStoBWPropRewgt);
  if (theTree->GetBranchStatus("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM")) theTree->SetBranchAddress("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM", &p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM);
  if (theTree->GetBranchStatus("p_Gen_GG_BKG_MCFM")) theTree->SetBranchAddress("p_Gen_GG_BKG_MCFM", &p_Gen_GG_BKG_MCFM);
  if (theTree->GetBranchStatus("p_Gen_QQB_BKG_MCFM")) theTree->SetBranchAddress("p_Gen_QQB_BKG_MCFM", &p_Gen_QQB_BKG_MCFM);

}

int InputTreeHandle::GetEntries()const{
  if (theTree!=0) return theTree->GetEntries();
  return 0;
}
void InputTreeHandle::GetEntry(int ev){ if (theTree!=0) theTree->GetEntry(ev); }
int InputTreeHandle::GetSampleType()const{ return theType; }

float InputTreeHandle::GetNGenEvents(){
  if (theCounters!=0) return theCounters->GetBinContent(40);
  else return 0;
}

#endif
