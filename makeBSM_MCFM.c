#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TColor.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"
#include "./data/HZZ4l_LeptonInterference.h"
#include "./data/Uncertainty_Tables.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>

using namespace std;

//Initializers
void makeBSM_MCFM_single(int folder, int erg_tev);

//Make lepton interference graph
TGraph* make_HZZ_LeptonInterferenceGraph(){
  float x[leptonInterf_YR3_Size];
  float y[leptonInterf_YR3_Size];
  for (int a=0; a<leptonInterf_YR3_Size; a++){
    x[a] = leptonInterf_YR3[a][0];
    y[a] = leptonInterf_YR3[a][1];
  }
  TGraph* tg = new TGraph(leptonInterf_YR3_Size, x, y);
  tg->SetName("tgHZZ_LeptonInterference");
  tg->SetTitle("H#rightarrowZZ#rightarrow4l Lepton Interference Weight on 4e, 4#mu wrt. 2e2#mu");

  return tg;
}

float findLeptonMass(int id){
  const float M_muon = 0.105658389;
  const float M_electron = 0.00051099907;
  const float M_tau = 1.777;

  if (abs(id)==11) return M_electron;
  else if (abs(id)==13) return M_muon;
  else if (abs(id)==15) return M_tau;
  else return 0;
}

bool passAcceptance(float pt, float eta, int id){
  if (abs(id)==13 && pt>5 && fabs(eta)<2.4) return true;
  else if (abs(id)==11 && pt>7 && fabs(eta)<2.5) return true;
  else return false;
}

//Main Function, runs over all desired iterations
void makeBSM_MCFM(){
  for (int CoM=7; CoM<9; ++CoM){
    for (int channel=0; channel<3; ++channel){
      makeBSM_MCFM_single(channel, CoM);
    }
  }
}

float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[SIZE_HVV][2]){
  float myprob=1.0;
  int myflavor=-1;
  if (abs(lepId[0])==abs(lepId[1]) &&
    abs(lepId[0])==abs(lepId[2]) &&
    abs(lepId[0])==abs(lepId[3])){
    if (abs(lepId[0])==11) myflavor=1;
    else myflavor=2;
  }
  else myflavor=3;

  if (myflavor>=0) myMela.computeP(angularOrdered[0], angularOrdered[1], angularOrdered[2], angularOrdered[3],
    angularOrdered[4], angularOrdered[5], angularOrdered[6], angularOrdered[7],
    myflavor,
    selfDHvvcoupl,
    myprob
    );
  return myprob;
};


void makeBSM_MCFM_single(int folder, int erg_tev){
  //  if (erg_tev==7) return;
  //  char TREE_NAME[] = "SelectedTree";
  char TREE_NAME[] = "GenTree";
  TString INPUT_NAME = "HZZ4l-125_6-";
  int tFitD = 0;
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV/", erg_tev);
  INPUT_NAME.Append(comstring);
  INPUT_NAME.Append("-");

  double mPOLE = 125.6;
  Mela mela(erg_tev, mPOLE);

  TString cinput_common = "/scratch0/hep/usarical/WidthSamples/" + erg_dir + "/";
  if (erg_tev!=7) cinput_common = cinput_common + user_folder[folder] + "/";

  int EnergyIndex = 1;
  if (erg_tev == 7) EnergyIndex = 0;

  float templateWeight = 1;
  float MC_weight;
  float ZZMass = 0;
  float Z1Mass = 0;
  float Z2Mass = 0;

  float GenHMass = 0;
  float GenZ1Mass = 0;
  float GenZ2Mass = 0;
  float GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;
  float GenLep1Mass, GenLep2Mass, GenLep3Mass, GenLep4Mass;
  float GenLep1Pt, GenLep2Pt, GenLep3Pt, GenLep4Pt;
  float GenLep1Phi, GenLep2Phi, GenLep3Phi, GenLep4Phi;
  float GenLep1Eta, GenLep2Eta, GenLep3Eta, GenLep4Eta;
  float GenLep1E, GenLep2E, GenLep3E, GenLep4E;

  float p0plus_VAJHU;
  float p0hplus_VAJHU;
  float p0minus_VAJHU;
  float p0_g1prime2_VAJHU;
  float p0_g1prime4_VAJHU;

  float pg1g2_VAJHU;
  float pg1g4_VAJHU;
  float pg1g1prime2_VAJHU;
  float pg1g1prime4_VAJHU;

  TString coutput_common = user_TemplateswithTrees_dir + "../BSMReweight_GenLevel/" + erg_dir;
  coutput_common += user_folder[folder] + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
  if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
  cinput_KDFactor = cinput_KDFactor + ".root";
  TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
  TString tgkfname = "KDScale_";
  tgkfname = tgkfname + "AllFlavors_UnNormalized";
  TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);
  TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

  cout<<endl;
  cout<<"==============================="<<endl;
  cout<<"CoM Energy: "<<erg_tev<<" TeV"<<endl;
  cout<<"Decay Channel: "<<user_folder[folder]<<endl;
  cout<<"==============================="<<endl;
  cout<<endl;


  const int kNumTemplates = 4;
  TString templatenames[kNumTemplates] ={ "ggF Sig", "gg Bkg", "ggF BSI", "ggF BSI25" };

  int nEntries;
  float fitYval = 0;
  MC_weight = 1;
  ZZMass = 0;

  double selfDHvvcoupl[SIZE_HVV][2] ={ { 0 } };
  double ggvvcoupl[2]={ 0, 0 };

  //Template and tree filler
  for (int t = 0; t < kNumTemplates; t++){
    TString OUTPUT_NAME;
    OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_BSMTrees_";
    OUTPUT_NAME += sample_suffix[t] + ".root";
    TString coutput = coutput_common + OUTPUT_NAME;
    TFile* foutput = new TFile(coutput, "recreate");

    foutput->cd();

    TChain* tree = new TChain(TREE_NAME);
    TTree* templateTree = new TTree(TREE_NAME, TREE_NAME);

    float wPOLE=4.15e-3;
    if (t==3) wPOLE*=25;

    TString cinput = cinput_common + sample_suffix[t] + "/" + INPUT_NAME + sample_suffix[t] + ".root";
    tree->Add(cinput);

    templateTree->Branch("MC_weight", &templateWeight);
    templateTree->Branch("ZZMass", &ZZMass);
    templateTree->Branch("Z1Mass", &Z1Mass);
    templateTree->Branch("Z2Mass", &Z2Mass);
    //    templateTree->Branch("GenHMass", &GenHMass);
    templateTree->Branch("p0plus_VAJHU", &p0plus_VAJHU);
    templateTree->Branch("p0hplus_VAJHU", &p0hplus_VAJHU);
    templateTree->Branch("p0minus_VAJHU", &p0minus_VAJHU);
    templateTree->Branch("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
    templateTree->Branch("p0_g1prime4_VAJHU", &p0_g1prime4_VAJHU);
    templateTree->Branch("pg1g2_VAJHU", &pg1g2_VAJHU);
    templateTree->Branch("pg1g4_VAJHU", &pg1g4_VAJHU);
    templateTree->Branch("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);
    templateTree->Branch("pg1g1prime4_VAJHU", &pg1g1prime4_VAJHU);

    //Making templates using appropriate weights
    //    tree->SetBranchAddress("ZZMass", &ZZMass);
    //    tree->SetBranchAddress("MC_weight", &MC_weight);
    tree->SetBranchAddress("GenZZMass", &GenHMass);
    tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
    tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
    tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
    tree->SetBranchAddress("GenLep4Id", &GenLep4Id);
    tree->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
    tree->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
    tree->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
    tree->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
    tree->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
    tree->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
    tree->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
    tree->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
    tree->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
    tree->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
    tree->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
    tree->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
    tree->SetBranchAddress("GenLep1E", &GenLep1E);
    tree->SetBranchAddress("GenLep2E", &GenLep2E);
    tree->SetBranchAddress("GenLep3E", &GenLep3E);
    tree->SetBranchAddress("GenLep4E", &GenLep4E);

    if (tFitD != 0){
      if (tree->GetBranchStatus(strFitDim[tFitD])){
        tree->SetBranchAddress("Dgg10_VAMCFM", &fitYval);
        templateTree->Branch(strFitDim[tFitD], &fitYval);
      }
      else{
        cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
        fitYval = 0;
      }
    }

    nEntries = tree->GetEntries();
    double xsec = 1./nEntries;
    if (t==0) xsec *= xsec_ggHZZ_MCFM[EnergyIndex];
    else if (t==1) xsec *= xsec_ggZZ_MCFM[EnergyIndex];
    else if (t==2) xsec *= xsec_ggHZZ_BSI_MCFM[EnergyIndex];
    else if (t==3) xsec *= xsec_ggHZZ_BSI25_MCFM[EnergyIndex];
    else xsec = 1;
    if (folder<2) xsec *= 0.25;
    else xsec *= 0.5;

    //    for (int ev = 0; ev < 100; ev++){

    for (int ev = 0; ev < nEntries; ev++){
      tree->GetEntry(ev);
      ZZMass = GenHMass;
      MC_weight = xsec;

      if (fitYval != fitYval) continue;
      // Protect against any KD exceeding boundaries
      if (tFitD!=0 && fitYval>=1){
        cout << "Found fitYval == " << fitYval;
        fitYval = 0.999;
        cout << ". Fixed to " << fitYval << endl;
      }
      if (tFitD!=0 && fitYval<0) fitYval = 0;

      if (t < 4 && erg_tev==7){
        if (folder==1){
          GenLep1Id=11;
          GenLep2Id=-11;
          GenLep3Id=11;
          GenLep4Id=-11;
        }
        else if (folder==0){
          GenLep1Id=13;
          GenLep2Id=-13;
          GenLep3Id=13;
          GenLep4Id=-13;
        }
      }

      double weight = MC_weight;
      if (t < 4) weight *= tgkf->Eval(GenHMass);
      if (t < 4
        && abs(GenLep1Id)==abs(GenLep2Id)
        && abs(GenLep1Id)==abs(GenLep3Id)
        && abs(GenLep1Id)==abs(GenLep4Id)
        ) weight *= tg_interf->Eval(GenHMass);
      templateWeight = weight;

      float mzz = mPOLE;
      float m1 = 0;
      float m2 = 0;
      float h1 = 0;
      float h2 = 0;
      float phi = 0;
      float hs = 0;
      float phi1 = 0;

      GenLep1Mass = findLeptonMass(GenLep1Id);
      GenLep2Mass = findLeptonMass(GenLep2Id);
      GenLep3Mass = findLeptonMass(GenLep3Id);
      GenLep4Mass = findLeptonMass(GenLep4Id);

      TLorentzVector p1, p2, p3, p4;
      p1.SetPtEtaPhiE(GenLep1Pt, GenLep1Eta, GenLep1Phi, GenLep1E);
      p2.SetPtEtaPhiE(GenLep2Pt, GenLep2Eta, GenLep2Phi, GenLep2E);
      p3.SetPtEtaPhiE(GenLep3Pt, GenLep3Eta, GenLep3Phi, GenLep3E);
      p4.SetPtEtaPhiE(GenLep4Pt, GenLep4Eta, GenLep4Phi, GenLep4E);
      TLorentzVector pZ1 = p1+p2;
      TLorentzVector pZ2 = p3+p4;
      TLorentzVector pZ1p = p1+p4;
      TLorentzVector pZ2p = p3+p2;
      TLorentzVector pZZ = pZ1+pZ2;
      mzz=pZZ.M();
      m1=pZ1.M();
      m2=pZ2.M();
      Z1Mass = m1;
      Z2Mass=m2;
      if (Z1Mass>=120 || Z1Mass<=40) continue;
      if (Z2Mass>=120 || Z2Mass<=12) continue;
      if (!passAcceptance(GenLep1Pt, GenLep1Eta, (int)GenLep1Id)) continue;
      if (!passAcceptance(GenLep2Pt, GenLep2Eta, (int)GenLep2Id)) continue;
      if (!passAcceptance(GenLep3Pt, GenLep3Eta, (int)GenLep3Id)) continue;
      if (!passAcceptance(GenLep4Pt, GenLep4Eta, (int)GenLep4Id)) continue;
      if (pZ1.M()<=4 ||pZ2.M()<=4 ||pZ1p.M()<=4 ||pZ2p.M()<=4) continue;
      TLorentzVector pp[4] ={ p1, p2, p3, p4 };
      int pass20=-1;
      int pass10=-1;
      for (int vv=0; vv<4; vv++){
        if (pp[vv].Pt()>20){
          pass20=vv; break;
        }
      }
      for (int vv=0; vv<4; vv++){
        if (vv==pass20) continue;
        if (pp[vv].Pt()>10){
          pass10=vv; break;
        }
      }
      if (pass20<0 || pass10<0) continue;

      mela::computeAngles(p1, GenLep1Id,
        p2, GenLep2Id,
        p3, GenLep3Id,
        p4, GenLep4Id,
        hs,
        h1,
        h2,
        phi,
        phi1);
      int lepIdOrdered[4] ={ (int)GenLep1Id, (int)GenLep2Id, (int)GenLep3Id, (int)GenLep4Id };
      float angularOrdered[8] ={ mzz, m1, m2, hs, h1, h2, phi, phi1 };

      if (t==0) mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
      if (t==2 || t==3) mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
      if (t!=1){
        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[0][0]=1;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[1][0]=1.638;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        p0hplus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[3][0]=2.521;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        p0minus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[11][0]=-12046.01;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[13][0]=-pow(10000.0/mPOLE, 2);
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        p0_g1prime4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);


        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[0][0]=1;
        selfDHvvcoupl[1][0]=1.638;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        pg1g2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[0][0]=1;
        selfDHvvcoupl[3][0]=2.521;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        pg1g4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[0][0]=1;
        selfDHvvcoupl[11][0]=-12046.01;
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
        selfDHvvcoupl[0][0]=1;
        selfDHvvcoupl[13][0]=-pow(10000.0/mPOLE, 2);
        mela.setMelaLeptonInterference(TVar::InterfOn);
        mela.setMelaHiggsWidth(wPOLE);
        pg1g1prime4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

        if (p0plus_VAJHU!=0){
          p0hplus_VAJHU/=p0plus_VAJHU;
          p0minus_VAJHU/=p0plus_VAJHU;
          p0_g1prime2_VAJHU/=p0plus_VAJHU;
          p0_g1prime4_VAJHU/=p0plus_VAJHU;

          pg1g2_VAJHU/=p0plus_VAJHU;
          pg1g4_VAJHU/=p0plus_VAJHU;
          pg1g1prime2_VAJHU/=p0plus_VAJHU;
          pg1g1prime4_VAJHU/=p0plus_VAJHU;
          p0plus_VAJHU=1;
        }
        else{
          p0plus_VAJHU=0;
          p0hplus_VAJHU=0;
          p0minus_VAJHU=0;
          p0_g1prime2_VAJHU=0;
          p0_g1prime4_VAJHU=0;

          pg1g2_VAJHU=0;
          pg1g4_VAJHU=0;
          pg1g1prime2_VAJHU=0;
          pg1g1prime4_VAJHU=0;
        }
      }
      else{
        p0plus_VAJHU=1;
        p0hplus_VAJHU=1;
        p0minus_VAJHU=1;
        p0_g1prime2_VAJHU=1;
        p0_g1prime4_VAJHU=1;

        pg1g2_VAJHU=1;
        pg1g4_VAJHU=1;
        pg1g1prime2_VAJHU=1;
        pg1g1prime4_VAJHU=1;
      }

      p0plus_VAJHU*=templateWeight;
      p0hplus_VAJHU*=templateWeight;
      p0minus_VAJHU*=templateWeight;
      p0_g1prime2_VAJHU*=templateWeight;
      p0_g1prime4_VAJHU*=templateWeight;

      pg1g2_VAJHU*=templateWeight;
      pg1g4_VAJHU*=templateWeight;
      pg1g1prime2_VAJHU*=templateWeight;
      pg1g1prime4_VAJHU*=templateWeight;

      templateTree->Fill();
    }

    foutput->WriteTObject(templateTree);
    delete tree;
    delete templateTree;
    foutput->Close();
  }

  delete tg_interf;
  delete tgkf;
  finput_KDFactor->Close();
}

void plot_BSM_MCFM(int SignalOnly=0){
  gROOT->ProcessLine(".x tdrstyle.cc");
  double mPOLE = 125.6;

  TString OUTPUT_NAME;
  OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_BSMPlots";
  if (SignalOnly==0) OUTPUT_NAME.Append(".root");
  else OUTPUT_NAME.Append("_SignalOnly.root");
  TString coutput_common = user_TemplateswithTrees_dir + "../BSMReweight_GenLevel/Plots/";
  gSystem->Exec("mkdir -p " + coutput_common);
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  foutput->cd();

  float ZZMass = 0;
  float p0plus_VAJHU;
  float p0hplus_VAJHU;
  float p0minus_VAJHU;
  float p0_g1prime2_VAJHU;
  float p0_g1prime4_VAJHU;
  float pg1g2_VAJHU;
  float pg1g4_VAJHU;
  float pg1g1prime2_VAJHU;
  float pg1g1prime4_VAJHU;

  TChain* tc[2][3][4];
  TH1F* hfill[4][9];
  int nbinsx = 73;
  double xlimits[2]={ 160, 1620 };
  if (SignalOnly==1){
    xlimits[0]=100;
    nbinsx = 76;
  }
  for (int t=0; t<4; t++){
    for (int ac=0; ac<9; ac++){
      hfill[t][ac]= new TH1F(Form("hSum_BSI%i_Hypo%i", t, ac), "", nbinsx, xlimits[0], xlimits[1]);
      hfill[t][ac]->SetXTitle("m_{4l} (GeV)");
      hfill[t][ac]->SetYTitle(Form("Events / %.0f GeV", (xlimits[1]-xlimits[0])/nbinsx));
    }
  }
  double nCounted[2][3][9]={ { { 0 } } };
  double nCountedScaled[2][3][9]={ { { 0 } } };
  for (int erg_tev=7; erg_tev<9; erg_tev++){
    for (int folder=0; folder<3; folder++){
      TString comstring;
      comstring.Form("%iTeV", erg_tev);
      TString erg_dir;
      erg_dir.Form("LHC_%iTeV/", erg_tev);

      int EnergyIndex = 1;
      if (erg_tev == 7) EnergyIndex = 0;
      TString cinput_common = user_TemplateswithTrees_dir + "../BSMReweight_GenLevel/";
      cinput_common.Append(+erg_dir);
      cinput_common += user_folder[folder] + "/";
      cout << cinput_common << endl;
      for (int t=0; t<4; t++){
        TString INPUT_NAME;
        INPUT_NAME = "HtoZZ4l_MCFM_125p6_BSMTrees_";
        INPUT_NAME += sample_suffix[t] + ".root";
        INPUT_NAME.Prepend(cinput_common);

        tc[EnergyIndex][folder][t] = new TChain("GenTree");
        if (t!=3) tc[EnergyIndex][folder][t]->Add(INPUT_NAME);
        tc[EnergyIndex][folder][t]->SetBranchAddress("ZZMass", &ZZMass);
        tc[EnergyIndex][folder][t]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("p0hplus_VAJHU", &p0hplus_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("p0minus_VAJHU", &p0minus_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("p0_g1prime4_VAJHU", &p0_g1prime4_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("pg1g2_VAJHU", &pg1g2_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("pg1g4_VAJHU", &pg1g4_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);
        tc[EnergyIndex][folder][t]->SetBranchAddress("pg1g1prime4_VAJHU", &pg1g1prime4_VAJHU);
      }

      double nsig_counted[9] ={ 0 };
      for (int ev=0; ev<tc[EnergyIndex][folder][0]->GetEntries(); ev++){
        tc[EnergyIndex][folder][0]->GetEntry(ev);
        if (fabs(ZZMass-mPOLE)<0.1){
          nsig_counted[0] += p0plus_VAJHU;
          nsig_counted[1] += p0hplus_VAJHU;
          nsig_counted[2] += p0minus_VAJHU;
          nsig_counted[3] += p0_g1prime2_VAJHU;
          nsig_counted[4] += p0_g1prime4_VAJHU;
          nsig_counted[5] += (2.*(p0plus_VAJHU+p0hplus_VAJHU)-pg1g2_VAJHU);
          nsig_counted[6] += pg1g4_VAJHU;
          nsig_counted[7] += pg1g1prime2_VAJHU;
          nsig_counted[8] += pg1g1prime4_VAJHU;
        }
      }

      for (int ac=0; ac<9; ac++) nCounted[EnergyIndex][folder][ac] = nsig_counted[ac];

    }
  }


  for (int erg_tev=7; erg_tev<9; erg_tev++){
    for (int folder=0; folder<3; folder++){
      int EnergyIndex = 1;
      if (erg_tev == 7) EnergyIndex = 0;

      double nsig_SM = nSM_ScaledPeak[EnergyIndex][folder];
      double SMscale = nsig_SM/nCounted[EnergyIndex][folder][0];
      double scale=1;
      for (int t=0; t<4; t++){
        for (int ev=0; ev<tc[EnergyIndex][folder][t]->GetEntries(); ev++){
          tc[EnergyIndex][folder][t]->GetEntry(ev);
          if (ZZMass<xlimits[0]) continue;
          if (ZZMass>=xlimits[1]) ZZMass=xlimits[1]*0.999;

          scale = SMscale;
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][0] = nCounted[EnergyIndex][folder][0]*scale;
          hfill[t][0]->Fill(ZZMass, p0plus_VAJHU*scale);
          scale = SMscale;
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][1] = nCounted[EnergyIndex][folder][1]*scale;
          hfill[t][1]->Fill(ZZMass, p0hplus_VAJHU*scale);
          scale = SMscale;
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][2] = nCounted[EnergyIndex][folder][2]*scale;
          hfill[t][2]->Fill(ZZMass, p0minus_VAJHU*scale);
          scale = SMscale;
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][3] = nCounted[EnergyIndex][folder][3]*scale;
          hfill[t][3]->Fill(ZZMass, p0_g1prime2_VAJHU*scale);
          scale = SMscale;
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][4] = nCounted[EnergyIndex][folder][4]*scale;
          hfill[t][4]->Fill(ZZMass, p0_g1prime4_VAJHU*scale);
          scale = SMscale;
          if (t>0){
            hfill[t][5]->Fill(ZZMass, pg1g2_VAJHU*scale);
          }
          else{
            scale *= (nCounted[0][2][0]+nCounted[1][2][0])/(nCounted[0][2][5]+nCounted[1][2][5]);
            hfill[t][5]->Fill(ZZMass, (2.*(p0plus_VAJHU+p0hplus_VAJHU)-pg1g2_VAJHU)*scale);
          }
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][5] = nCounted[EnergyIndex][folder][5]*scale;
          scale = SMscale;
          if (t==0) scale *= (nCounted[0][2][0]+nCounted[1][2][0])/(nCounted[0][2][6]+nCounted[1][2][6]);
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][6] = nCounted[EnergyIndex][folder][6]*scale;
          hfill[t][6]->Fill(ZZMass, pg1g4_VAJHU*scale);
          scale = SMscale;
          if (t==0) scale *= (nCounted[0][2][0]+nCounted[1][2][0])/(nCounted[0][2][7]+nCounted[1][2][7]);
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][7] = nCounted[EnergyIndex][folder][7]*scale;
          hfill[t][7]->Fill(ZZMass, pg1g1prime2_VAJHU*scale);
          scale = SMscale;
          //          if (t==0) scale *= (nCounted[0][2][0]+nCounted[1][2][0])/(nCounted[0][2][8]+nCounted[1][2][8]);
          if (t==0 && ev==0) nCountedScaled[EnergyIndex][folder][8] = nCounted[EnergyIndex][folder][8]*scale;
          hfill[t][8]->Fill(ZZMass, pg1g1prime4_VAJHU*scale);
        }
        delete tc[EnergyIndex][folder][t];
      }
    }
  }
  for (int ac=1; ac<9; ac++){
    double nTotal[2]={ 0 };
    for (int erg_tev=7; erg_tev<9; erg_tev++){
      for (int folder=0; folder<3; folder++){
        int EnergyIndex = 1;
        if (erg_tev == 7) EnergyIndex = 0;

        nTotal[0] += nCountedScaled[EnergyIndex][folder][0];
        nTotal[1] += nCountedScaled[EnergyIndex][folder][ac];
      }
    }
    if (ac!=8) hfill[0][ac]->Scale(nTotal[0]/nTotal[1]);
    else hfill[0][ac]->Scale(0.5);
  }


  double maxplot=0;
  for (int t=0; t<4; t++){
    for (int ac=0; ac<9; ac++){
      if (SignalOnly==0 && ac<5) maxplot = max(maxplot, hfill[t][ac]->GetMaximum());
      else if (SignalOnly==1 && !(ac<5 && ac>0) && t==0) maxplot = max(maxplot, hfill[t][ac]->GetMaximum());
      hfill[t][ac]->SetLineWidth(2);
      if (t==0 && ac>=5){
        hfill[t][ac]->SetLineStyle(7);
        //        hfill[t][ac]->Add(hfill[1][ac]);
      }
      if (t==1) hfill[t][ac]->SetLineStyle(3);
      if (t==3) hfill[t][ac]->SetLineStyle(9);
      hfill[t][ac]->GetXaxis()->SetLabelFont(42);
      hfill[t][ac]->GetXaxis()->SetLabelOffset(0.007);
      hfill[t][ac]->GetXaxis()->SetLabelSize(0.04);
      hfill[t][ac]->GetXaxis()->SetTitleSize(0.06);
      hfill[t][ac]->GetXaxis()->SetTitleOffset(0.9);
      hfill[t][ac]->GetXaxis()->SetTitleFont(42);
      hfill[t][ac]->GetYaxis()->SetNdivisions(505);
      hfill[t][ac]->GetYaxis()->SetLabelFont(42);
      hfill[t][ac]->GetYaxis()->SetLabelOffset(0.007);
      hfill[t][ac]->GetYaxis()->SetLabelSize(0.04);
      hfill[t][ac]->GetYaxis()->SetTitleSize(0.06);
      hfill[t][ac]->GetYaxis()->SetTitleOffset(1.1);
      hfill[t][ac]->GetYaxis()->SetTitleFont(42);
    }
  }

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  text->SetTextSize(0.0315);
  TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
  text = pt->AddText(0.537, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  float lxmin = 0.22;
  float lxwidth = 0.38;
  float lymax = 0.9;
  float lywidth = 0.3;
  float lxmax = lxmin + lxwidth;
  float lymin = lymax;
  if (SignalOnly==0) lymin -= lywidth*4./5.;
  else lymin -= lywidth;

  float lxmin2 = 0.22+0.39;
  float lymax2 = lymax;
  float lxmax2 = lxmin2 + lxwidth;
  float lymin2 = lymax2;
  if (SignalOnly==0) lymin2 -= lywidth*2./5.;
  else lymin2 -= lywidth*4./5.;

  if (SignalOnly==1){
    float lxmin3 = lxmin2;
    float lymax3 = lymax2;
    float lxmax3 = lxmax2;
    float lymin3 = lymin2;

    lxmin2 = lxmin;
    lxmax2 = lxmax;
    lymin2 = lymin;
    lymax2 = lymax;
    lxmin = lxmin3;
    lxmax = lxmax3;
    lymin = lymin3;
    lymax = lymax3;
  }

  float pxmin = 0.756;
  float pymin = 0.76;
  float pxmax = 0.85;
  if (SignalOnly==1){
    pymin -= 0.12;
  }
  float pymax = pymin+0.05;
  TPaveText* ptx = new TPaveText(pxmin, pymin, pxmax, pymax, "brNDC");
  ptx->SetBorderSize(0);
  ptx->SetTextFont(42);
  ptx->SetTextSize(0.04);
  ptx->SetLineColor(1);
  ptx->SetLineStyle(1);
  ptx->SetLineWidth(1);
  ptx->SetFillColor(0);
  ptx->SetFillStyle(0);
  text = ptx->AddText(0.01, 0.01, "gg#rightarrow4l");
  text->SetTextSize(0.04);

  TString canvasname = "cCanvas_MCFMBSM_GenLevel";
  if (SignalOnly==1) canvasname.Append("_SignalOnly");
  TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
  cc->cd();
  gStyle->SetOptStat(0);
  cc->SetFillColor(0);
  cc->SetBorderMode(0);
  cc->SetBorderSize(2);
  cc->SetTickx(1);
  cc->SetTicky(1);
  cc->SetLeftMargin(0.17);
  cc->SetRightMargin(0.05);
  cc->SetTopMargin(0.07);
  cc->SetBottomMargin(0.13);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);
  cc->SetLogy();

  TLegend* ll;
  TLegend* ll2;

  ll = new TLegend(lxmin2, lymin2, lxmax2, lymax2);
  ll2 = new TLegend(lxmin, lymin, lxmax, lymax);

  ll->SetBorderSize(0);
  ll->SetTextFont(42);
  ll->SetTextSize(0.04);
  ll->SetLineColor(1);
  ll->SetLineStyle(1);
  ll->SetLineWidth(1);
  ll->SetFillColor(0);
  ll->SetFillStyle(0);
  ll2->SetBorderSize(0);
  ll2->SetTextFont(42);
  ll2->SetTextSize(0.04);
  ll2->SetLineColor(1);
  ll2->SetLineStyle(1);
  ll2->SetLineWidth(1);
  ll2->SetFillColor(0);
  ll2->SetFillStyle(0);

  TString strACtitle[9]={ "",
    "f_{a2}=1", "f_{a3}=1", "f_{#Lambda1}=1", "f_{#LambdaQ}=1",
    "f_{a2}=0.5, #phi_{#lower[-0.2]{a2}}=#pi", "f_{a3}=0.5", "f_{#Lambda1}=0.5", "f_{#LambdaQ}=0.5"
  };

  int iDraw = 2 - 2*SignalOnly;

  if (SignalOnly==0) hfill[iDraw][0]->GetYaxis()->SetRangeUser(7e-3, maxplot*15.);
  else{
    double histmin = 7e-3;
    if (hfill[iDraw][0]->GetMinimum()>0) histmin = hfill[iDraw][0]->GetMinimum();
    hfill[iDraw][0]->GetYaxis()->SetRangeUser(histmin, maxplot*2000.);
  }
  hfill[iDraw][0]->GetXaxis()->SetRangeUser(xlimits[0], 800.);

  hfill[iDraw][0]->SetLineColor(kBlack);
  if (SignalOnly==0){
    hfill[iDraw][0]->SetFillColor(kAzure-2);
    hfill[iDraw][0]->SetFillStyle(1001);
  }
  hfill[iDraw][0]->Draw("hist");

  hfill[iDraw][1]->SetLineColor(kBlue);
  hfill[iDraw][1]->Draw("histsame");

  hfill[iDraw][2]->SetLineColor(kRed);
  hfill[iDraw][2]->Draw("histsame");

  hfill[iDraw][3]->SetLineColor(kViolet);
  hfill[iDraw][3]->Draw("histsame");

  hfill[iDraw][4]->SetLineColor(kGreen+2);
  hfill[iDraw][4]->Draw("histsame");

  if (SignalOnly==1){
    hfill[iDraw][5]->SetLineColor(kBlue);
    hfill[iDraw][5]->Draw("histsame");

    hfill[iDraw][6]->SetLineColor(kRed);
    hfill[iDraw][6]->Draw("histsame");

    hfill[iDraw][7]->SetLineColor(kViolet);
    hfill[iDraw][7]->Draw("histsame");

    hfill[iDraw][8]->SetLineColor(kGreen+2);
    hfill[iDraw][8]->Draw("histsame");
  }

  if (SignalOnly==0){
    hfill[1][0]->SetLineColor(kBlack);
    hfill[1][0]->SetLineStyle(3);
    hfill[1][0]->Draw("histsame");
  }

  hfill[iDraw][0]->Draw("histsame");

  TLegendEntry* legendtext;
  if (SignalOnly==0){
    legendtext = ll->AddEntry(hfill[iDraw][0], "SM total", "f");
    legendtext = ll->AddEntry(hfill[1][0], "SM bkg.", "f");
    legendtext->SetFillStyle(1001);
    legendtext->SetFillColor(hfill[1][0]->GetFillColor());
  }
  else{
    legendtext = ll->AddEntry(hfill[iDraw][0], "SM signal", "f");
    legendtext->SetFillStyle(3001);
  }

  if (SignalOnly==0){
    legendtext = ll2->AddEntry(hfill[iDraw][4], Form("#Gamma_{H}=#Gamma^{SM}_{H}, %s", strACtitle[4].Data()), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][4]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][2], Form("#Gamma_{H}=#Gamma^{SM}_{H}, %s", strACtitle[2].Data()), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][2]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][1], Form("#Gamma_{H}=#Gamma^{SM}_{H}, %s", strACtitle[1].Data()), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][1]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][3], Form("#Gamma_{H}=#Gamma^{SM}_{H}, %s", strACtitle[3].Data()), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][3]->GetFillColor());
  }
  else{
    legendtext = ll->AddEntry(hfill[iDraw][4], strACtitle[4].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][4]->GetFillColor());
    legendtext = ll->AddEntry(hfill[iDraw][2], strACtitle[2].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][2]->GetFillColor());
    legendtext = ll->AddEntry(hfill[iDraw][1], strACtitle[1].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][1]->GetFillColor());
    legendtext = ll->AddEntry(hfill[iDraw][3], strACtitle[3].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][3]->GetFillColor());
  }
  if (SignalOnly==1){
    legendtext = ll2->AddEntry(hfill[iDraw][8], strACtitle[8].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][8]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][5], strACtitle[5].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][5]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][7], strACtitle[7].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][7]->GetFillColor());
    legendtext = ll2->AddEntry(hfill[iDraw][6], strACtitle[6].Data(), "f");
    legendtext->SetFillStyle(3001);
    legendtext->SetFillColor(hfill[iDraw][6]->GetFillColor());
  }

  ll->Draw("same");
  ll2->Draw("same");
  ptx->Draw();
  pt->Draw();
  cc->RedrawAxis();
  cc->Update();

  canvasname.Prepend(coutput_common);
  TString canvasname_pdf = canvasname;
  TString canvasname_eps = canvasname;
  TString canvasname_png = canvasname;
  TString canvasname_root = canvasname;
  TString canvasname_c = canvasname;
  canvasname_pdf.Append(".pdf");
  canvasname_eps.Append(".eps");
  canvasname_png.Append(".png");
  canvasname_root.Append(".root");
  canvasname_c.Append(".C");
  cc->SaveAs(canvasname_pdf);
  cc->SaveAs(canvasname_eps);
  cc->SaveAs(canvasname_png);
  cc->SaveAs(canvasname_root);
  cc->SaveAs(canvasname_c);

  foutput->WriteTObject(cc);
  delete ll2;
  delete ll;
  cc->Close();
  delete ptx;
  delete pt;
  for (int t=0; t<4; t++){
    for (int ac=0; ac<5; ac++){
      foutput->WriteTObject(hfill[t][ac]);
      delete hfill[t][ac];
    }
  }
  foutput->Close();
}
