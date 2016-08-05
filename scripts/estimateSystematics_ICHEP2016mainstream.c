#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include "TFile.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "../data/ZZ4l_Samples.h"
#include "../include/external_cConstants.h"

using namespace std;

// Categories
enum Category{
  Inclusive = 0,
  Untagged = 1,
  VBF2jTagged = 2,
  VBF1jTagged = 3,
  LepVHTagged = 4,
  HadVHTagged = 5,
  ttHTagged = 6
};

// Constants to affect the template code
const TString fixedDate="";

// Initializers
void estimateSystematics_ICHEP2016mainstream_one(int channel, int erg_tev, int categorizationType, float ZZMass_low=-1, float ZZMass_high=-1);
void progressbar(int val, int tot);
vector<pair<string, double>> getFileList(string strinput);
TString todaysdate();
void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray);
void addByLowest(float weight, float val, std::vector<std::pair<float, float>>& valArray);
void addByLowest(float val, std::vector<float>& valArray);
int getCategory(
  float ZZMass,
  short nExtraLep,
  short nExtraZ,
  std::vector<float> JetPt,
  std::vector<float> JetPhi,
  std::vector<float> JetSigma,
  std::vector<float> JetIsBtagged,
  std::vector<float> JetQGLikelihood,
  float phjj_VAJHU_highestPTJets,
  float phj_VAJHU,
  float pvbf_VAJHU_highestPTJets,
  float pAux_vbf_VAJHU,
  float pwh_hadronic_VAJHU,
  float pzh_hadronic_VAJHU,
  int categorizationType,
  int JESvariation=0
  );

// Main Function, runs over all desired iterations
void estimateSystematics_ICHEP2016mainstream_one(int folder, int erg_tev, int categorizationType, float ZZMass_low, float ZZMass_high){
  //int EnergyIndex=(erg_tev==13 ? 0 : 0);

  float ZZMass_cut[2]={ 70., 3000. };
  float xwidth = 5;
  if (ZZMass_low>0. && ZZMass_high>0.){
    ZZMass_cut[0]=ZZMass_low;
    ZZMass_cut[1]=ZZMass_high;
    xwidth=ZZMass_cut[1]-ZZMass_cut[0];
  }
  const int nbinsx = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth;
  float kDXarray[nbinsx+1];
  for (int bin=0; bin<nbinsx+1; bin++) kDXarray[bin] = ZZMass_cut[0] + xwidth*bin;

  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strChannel = (folder>=0 ? user_folder[folder].Data() : "4l");
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString OUTPUT_LOG_NAME = Form("HtoZZ%s_Systematics", strChannel.Data());
  if (categorizationType==-2) OUTPUT_LOG_NAME.Append("_DjetwQG");
  else if (categorizationType==-1) OUTPUT_LOG_NAME.Append("_DjetnoQG");
  else if (categorizationType==0) OUTPUT_LOG_NAME.Append("_noCat");
  else if (categorizationType==1) OUTPUT_LOG_NAME.Append("_6CatnoQG");
  else if (categorizationType==2) OUTPUT_LOG_NAME.Append("_6CatwQG");
  if (ZZMass_low>0. && ZZMass_high>0.) OUTPUT_LOG_NAME.Append(Form("_mZZlow_%.1f_mZZhigh_%.1f", ZZMass_cut[0], ZZMass_cut[1]));
  TString OUTPUT_NAME = OUTPUT_LOG_NAME;
  OUTPUT_LOG_NAME += ".log";
  OUTPUT_NAME += ".root";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = user_dir + erg_dir + "Systematics/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  float overallEventWeight=1;
  float xsec=1;

  float KFactor_QCD_ggZZ_Nominal = 1;
  float KFactor_QCD_ggZZ_PDFScaleDn = 1;
  float KFactor_QCD_ggZZ_PDFScaleUp = 1;
  float KFactor_QCD_ggZZ_QCDScaleDn = 1;
  float KFactor_QCD_ggZZ_QCDScaleUp = 1;
  float KFactor_QCD_ggZZ_AsDn = 1;
  float KFactor_QCD_ggZZ_AsUp = 1;
  float KFactor_QCD_ggZZ_PDFReplicaDn = 1;
  float KFactor_QCD_ggZZ_PDFReplicaUp = 1;

  float KFactor_EW_qqZZ = 1;
  float KFactor_EW_qqZZ_unc = 1;
  float KFactor_QCD_qqZZ_M = 1;
  float LHEweight_QCDscale_muR1_muF1  = 1;
  float LHEweight_QCDscale_muR1_muF2  = 1;
  float LHEweight_QCDscale_muR1_muF0p5  = 1;
  float LHEweight_QCDscale_muR2_muF1  = 1;
  float LHEweight_QCDscale_muR2_muF2  = 1;
  float LHEweight_QCDscale_muR2_muF0p5  = 1;
  float LHEweight_QCDscale_muR0p5_muF1  = 1;
  float LHEweight_QCDscale_muR0p5_muF2  = 1;
  float LHEweight_QCDscale_muR0p5_muF0p5  = 1;

  float GenHMass;
  short Z1Flav, Z2Flav;
  short nExtraLep;
  short nExtraZ;
  int ZZFlav;
  float ZZMass;
  float phj_VAJHU;
  float pAux_vbf_VAJHU;
  float phjj_VAJHU_highestPTJets;
  float pvbf_VAJHU_highestPTJets;
  float pwh_hadronic_VAJHU;
  float pzh_hadronic_VAJHU;
  float phj_VAJHU_up;
  float pAux_vbf_VAJHU_up;
  float phjj_VAJHU_highestPTJets_up;
  float pvbf_VAJHU_highestPTJets_up;
  float pwh_hadronic_VAJHU_up;
  float pzh_hadronic_VAJHU_up;
  float phj_VAJHU_dn;
  float pAux_vbf_VAJHU_dn;
  float phjj_VAJHU_highestPTJets_dn;
  float pvbf_VAJHU_highestPTJets_dn;
  float pwh_hadronic_VAJHU_dn;
  float pzh_hadronic_VAJHU_dn;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float>* JetBTagger=0;
  std::vector<float>* JetIsBtagged=0;
  std::vector<float>* JetQGLikelihood=0;
  std::vector<float>* JetAxis2=0;
  std::vector<float>* JetMult=0;
  std::vector<float>* JetPtD=0;
  std::vector<float>* JetSigma=0;

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << strChannel.Data() << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << strChannel.Data() << endl;
  tout << "===============================" << endl;
  tout << endl;

  // Get list of samples
  const unsigned int nProcesses = 7;
  const TString strProcess[nProcesses][2]={
    { "Samples_gg_Sig_POWHEG.txt", "gg_Sig" },
    { "Samples_VBF_Sig_POWHEG.txt", "VBF_Sig" },
    { "Samples_WH_Sig_POWHEG.txt", "WH_Sig" },
    { "Samples_ZH_Sig_POWHEG.txt", "ZH_Sig" },
    { "Samples_tt_Sig_POWHEG.txt", "tt_Sig" },
    { "Samples_gg_Bkg_MCFM.txt", "gg_Bkg" },
    { "Samples_qq_Bkg_POWHEG.txt", "qq_Bkg" }
  };
  vector<pair<string, double>> fileNameList[nProcesses];
  for (unsigned int is=0; is<nProcesses; is++) fileNameList[is] = getFileList(string(strProcess[is][0].Data()));
  vector<TFile*> fileList[nProcesses];
  vector<pair<TTree*, TH1F*>> treeList[nProcesses];
  for (unsigned int is=0; is<nProcesses; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_gg2VV_location + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;
      tout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      TFile* finput = TFile::Open(cinput.c_str(), "read"); finput->cd();
      TTree* theTree = (TTree*)finput->Get(user_treename);
      TH1F* theCounters = (TH1F*)finput->Get(user_countersname);

      fileList[is].push_back(finput);
      foutput->cd();

      theTree->SetBranchAddress("ZZMass", &ZZMass);
      theTree->SetBranchAddress("Z1Flav", &Z1Flav);
      theTree->SetBranchAddress("Z2Flav", &Z2Flav);
      theTree->SetBranchAddress("Z2Flav", &Z2Flav);
      theTree->SetBranchAddress("phj_VAJHU", &phj_VAJHU);
      theTree->SetBranchAddress("pAux_vbf_VAJHU", &pAux_vbf_VAJHU);
      theTree->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets);
      theTree->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets);
      theTree->SetBranchAddress("pwh_hadronic_VAJHU", &pwh_hadronic_VAJHU);
      theTree->SetBranchAddress("pzh_hadronic_VAJHU", &pzh_hadronic_VAJHU);
      theTree->SetBranchAddress("phj_VAJHU_up", &phj_VAJHU_up);
      theTree->SetBranchAddress("pAux_vbf_VAJHU_up", &pAux_vbf_VAJHU_up);
      theTree->SetBranchAddress("phjj_VAJHU_highestPTJets_up", &phjj_VAJHU_highestPTJets_up);
      theTree->SetBranchAddress("pvbf_VAJHU_highestPTJets_up", &pvbf_VAJHU_highestPTJets_up);
      theTree->SetBranchAddress("pwh_hadronic_VAJHU_up", &pwh_hadronic_VAJHU_up);
      theTree->SetBranchAddress("pzh_hadronic_VAJHU_up", &pzh_hadronic_VAJHU_up);
      theTree->SetBranchAddress("phj_VAJHU_dn", &phj_VAJHU_dn);
      theTree->SetBranchAddress("pAux_vbf_VAJHU_dn", &pAux_vbf_VAJHU_dn);
      theTree->SetBranchAddress("phjj_VAJHU_highestPTJets_dn", &phjj_VAJHU_highestPTJets_dn);
      theTree->SetBranchAddress("pvbf_VAJHU_highestPTJets_dn", &pvbf_VAJHU_highestPTJets_dn);
      theTree->SetBranchAddress("pwh_hadronic_VAJHU_dn", &pwh_hadronic_VAJHU_dn);
      theTree->SetBranchAddress("pzh_hadronic_VAJHU_dn", &pzh_hadronic_VAJHU_dn);

      theTree->SetBranchAddress("JetPt", &JetPt);
      theTree->SetBranchAddress("JetEta", &JetEta);
      theTree->SetBranchAddress("JetPhi", &JetPhi);
      theTree->SetBranchAddress("JetMass", &JetMass);
      theTree->SetBranchAddress("JetSigma", &JetSigma);
      theTree->SetBranchAddress("JetBTagger", &JetBTagger);
      theTree->SetBranchAddress("JetIsBtagged", &JetIsBtagged);
      theTree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
      theTree->SetBranchAddress("JetAxis2", &JetAxis2);
      theTree->SetBranchAddress("JetMult", &JetMult);
      theTree->SetBranchAddress("JetPtD", &JetPtD);

      theTree->SetBranchAddress("nExtraLep", &nExtraLep);
      theTree->SetBranchAddress("nExtraZ", &nExtraZ);

      if (theTree->GetBranchStatus("GenHMass")) theTree->SetBranchAddress("GenHMass", &GenHMass);
      if (theTree->GetBranchStatus("overallEventWeight")) theTree->SetBranchAddress("overallEventWeight", &overallEventWeight);
      if (theTree->GetBranchStatus("xsec")) theTree->SetBranchAddress("xsec", &xsec);
      if (theTree->GetBranchStatus("KFactor_QCD_ggZZ_Nominal")){
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleDn", &KFactor_QCD_ggZZ_PDFScaleDn);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleUp", &KFactor_QCD_ggZZ_PDFScaleUp);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleDn", &KFactor_QCD_ggZZ_QCDScaleDn);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleUp", &KFactor_QCD_ggZZ_QCDScaleUp);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_AsDn", &KFactor_QCD_ggZZ_AsDn);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_AsUp", &KFactor_QCD_ggZZ_AsUp);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaDn", &KFactor_QCD_ggZZ_PDFReplicaDn);
        theTree->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaUp", &KFactor_QCD_ggZZ_PDFReplicaUp);
      }
      if (theTree->GetBranchStatus("KFactor_QCD_qqZZ_M")) theTree->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
      if (theTree->GetBranchStatus("KFactor_EW_qqZZ")){
        theTree->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
        theTree->SetBranchAddress("KFactor_EW_qqZZ_unc", &KFactor_EW_qqZZ_unc);
      }
      if (theTree->GetBranchStatus("LHEweight_QCDscale_muR1_muF1")){
        theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF1", &LHEweight_QCDscale_muR1_muF1);

        theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF2", &LHEweight_QCDscale_muR1_muF2);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5", &LHEweight_QCDscale_muR1_muF0p5);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF1", &LHEweight_QCDscale_muR2_muF1);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1", &LHEweight_QCDscale_muR0p5_muF1);

        theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF2", &LHEweight_QCDscale_muR2_muF2);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5", &LHEweight_QCDscale_muR2_muF0p5);
        theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2", &LHEweight_QCDscale_muR0p5_muF2);
      }
      treeList[is].push_back(pair<TTree*, TH1F*>(theTree, theCounters));
    }
  }

  const unsigned int nCategories=7; // Categories
  TString strCategory[nCategories] ={
    "Inclusive",
    "UnTagged",
    "VBF2jTagged",
    "VBF1jTagged",
    "LepVHTagged",
    "HadVHTagged",
    "ttHTagged"
  };
  const unsigned int nSyst = 13;
  TString strSyst[nSyst]={
    "Nominal",
    "QCDScaleUp", "QCDScaleDown",
    "PDFScaleUp", "PDFScaleDown",
    "AsUp", "AsDn",
    "PDFReplicaUp", "PDFReplicaDn",
    "EWUp", "EWDown",
    "JESUp", "JESDn"
  };

  tout << endl;
  tout << "Process" << " " << "Category" << " ";
  for (unsigned int isyst=0; isyst<nSyst; isyst++) tout << strSyst[isyst] << " ";
  tout << endl;

  for (unsigned int is=0; is<nProcesses; is++){
    cout << "Process " << is << " under scrutiny" << endl;

    TH1F* hSyst[nCategories][nSyst];
    for (unsigned int icat=0; icat<nCategories; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        hSyst[icat][isyst] = new TH1F(Form("Process_%s_Category_%s_Syst_%s", strProcess[is][1].Data(), strCategory[icat].Data(), strSyst[isyst].Data()), "", nbinsx, kDXarray);
        hSyst[icat][isyst]->Sumw2();
      }
    }

    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      TTree* theTree = treeList[is].at(it).first;
      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;

      float nGenEvents = treeList[is].at(it).second->GetBinContent(40);
      for (int ev=0; ev<nEntries; ev++){
        overallEventWeight=1;
        xsec=1;

        KFactor_QCD_ggZZ_Nominal = 1;
        KFactor_QCD_ggZZ_PDFScaleDn = 1;
        KFactor_QCD_ggZZ_PDFScaleUp = 1;
        KFactor_QCD_ggZZ_QCDScaleDn = 1;
        KFactor_QCD_ggZZ_QCDScaleUp = 1;
        KFactor_QCD_ggZZ_AsDn = 1;
        KFactor_QCD_ggZZ_AsUp = 1;
        KFactor_QCD_ggZZ_PDFReplicaDn = 1;
        KFactor_QCD_ggZZ_PDFReplicaUp = 1;

        KFactor_EW_qqZZ = 1;
        KFactor_EW_qqZZ_unc = 0;
        KFactor_QCD_qqZZ_M = 1;

        LHEweight_QCDscale_muR1_muF1  = 1;
        LHEweight_QCDscale_muR1_muF2  = 1;
        LHEweight_QCDscale_muR1_muF0p5  = 1;
        LHEweight_QCDscale_muR2_muF1  = 1;
        LHEweight_QCDscale_muR2_muF2  = 1;
        LHEweight_QCDscale_muR2_muF0p5  = 1;
        LHEweight_QCDscale_muR0p5_muF1  = 1;
        LHEweight_QCDscale_muR0p5_muF2  = 1;
        LHEweight_QCDscale_muR0p5_muF0p5  = 1;

        theTree->GetEntry(ev);
        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Fix for JHUGen v6 vegas bug in low mass samples
        if (strProcess[is][0].Contains("POWHEG") && mass>0 && mass<300 && GenHMass>mass+5.) continue;

        if (ZZMass<ZZMass_cut[0] || ZZMass>=ZZMass_cut[1]) continue;

        ZZFlav=Z1Flav*Z2Flav;
        if (folder==0 && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (folder==1 && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (folder==2 && abs(ZZFlav)!=pow(13*11, 2)) continue;

        int icat=getCategory(
          ZZMass,
          nExtraLep,
          nExtraZ,
          *JetPt,
          *JetPhi,
          *JetSigma,
          *JetIsBtagged,
          *JetQGLikelihood,
          phjj_VAJHU_highestPTJets,
          phj_VAJHU,
          pvbf_VAJHU_highestPTJets,
          pAux_vbf_VAJHU,
          pwh_hadronic_VAJHU,
          pzh_hadronic_VAJHU,
          categorizationType,
          0
          );
        int icat_jesup=getCategory(
          ZZMass,
          nExtraLep,
          nExtraZ,
          *JetPt,
          *JetPhi,
          *JetSigma,
          *JetIsBtagged,
          *JetQGLikelihood,
          phjj_VAJHU_highestPTJets_up,
          phj_VAJHU_up,
          pvbf_VAJHU_highestPTJets_up,
          pAux_vbf_VAJHU_up,
          pwh_hadronic_VAJHU_up,
          pzh_hadronic_VAJHU_up,
          categorizationType,
          1
          );
        int icat_jesdn=getCategory(
          ZZMass,
          nExtraLep,
          nExtraZ,
          *JetPt,
          *JetPhi,
          *JetSigma,
          *JetIsBtagged,
          *JetQGLikelihood,
          phjj_VAJHU_highestPTJets_dn,
          phj_VAJHU_dn,
          pvbf_VAJHU_highestPTJets_dn,
          pAux_vbf_VAJHU_dn,
          pwh_hadronic_VAJHU_dn,
          pzh_hadronic_VAJHU_dn,
          categorizationType,
          -1
          );

        float weight[nSyst]={ 0 };
        /*
        "Nominal",
        "QCDScaleUp", "QCDScaleDown",
        "PDFScaleUp", "PDFScaleDown",
        "AsUp", "AsDn",
        "PDFReplicaUp", "PDFReplicaDn",
        "EWUp", "EWDown"
        */
        if (is==0 || is==5){ // gg samples
          weight[0] = KFactor_QCD_ggZZ_Nominal;
          weight[1] = KFactor_QCD_ggZZ_QCDScaleUp;
          weight[2] = KFactor_QCD_ggZZ_QCDScaleDn;
          weight[3] = KFactor_QCD_ggZZ_PDFScaleUp;
          weight[4] = KFactor_QCD_ggZZ_PDFScaleDn;
          weight[5] = KFactor_QCD_ggZZ_AsUp;
          weight[6] = KFactor_QCD_ggZZ_AsDn;
          weight[7] = KFactor_QCD_ggZZ_PDFReplicaUp;
          weight[8] = KFactor_QCD_ggZZ_PDFReplicaDn;
          weight[9] = KFactor_QCD_ggZZ_Nominal;
          weight[10] = KFactor_QCD_ggZZ_Nominal;
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[0][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);

          if (is==0){ // LHE weights are meaningless for gg bkg samples at LO
            weight[1] = KFactor_QCD_ggZZ_Nominal * LHEweight_QCDscale_muR2_muF1/LHEweight_QCDscale_muR1_muF1;
            weight[2] = KFactor_QCD_ggZZ_Nominal * LHEweight_QCDscale_muR0p5_muF1/LHEweight_QCDscale_muR1_muF1;
            weight[3] = KFactor_QCD_ggZZ_Nominal * LHEweight_QCDscale_muR1_muF2/LHEweight_QCDscale_muR1_muF1;
            weight[4] = KFactor_QCD_ggZZ_Nominal * LHEweight_QCDscale_muR1_muF0p5/LHEweight_QCDscale_muR1_muF1;
          }

          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)Inclusive) hSyst[icat][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)Inclusive) hSyst[icat_jesup][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)Inclusive) hSyst[icat_jesdn][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
          }
        }
        else if(is!=6){ // VBF, VH and ttH samples
          weight[0] = 1;
          weight[1] = LHEweight_QCDscale_muR2_muF1/LHEweight_QCDscale_muR1_muF1;
          weight[2] = LHEweight_QCDscale_muR0p5_muF1/LHEweight_QCDscale_muR1_muF1;
          weight[3] = LHEweight_QCDscale_muR1_muF2/LHEweight_QCDscale_muR1_muF1;
          weight[4] = LHEweight_QCDscale_muR1_muF0p5/LHEweight_QCDscale_muR1_muF1;
          weight[5] = 1;
          weight[6] = 1;
          weight[7] = 1;
          weight[8] = 1;
          weight[9] = 1;
          weight[10] = 1;
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[0][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)Inclusive) hSyst[icat][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)Inclusive) hSyst[icat_jesup][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)Inclusive) hSyst[icat_jesdn][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
          }

        }
        else{
          weight[0] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ;
          weight[1] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ*LHEweight_QCDscale_muR2_muF1/LHEweight_QCDscale_muR1_muF1;
          weight[2] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ*LHEweight_QCDscale_muR0p5_muF1/LHEweight_QCDscale_muR1_muF1;
          weight[3] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ*LHEweight_QCDscale_muR1_muF2/LHEweight_QCDscale_muR1_muF1;
          weight[4] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ*LHEweight_QCDscale_muR1_muF0p5/LHEweight_QCDscale_muR1_muF1;
          weight[5] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ;
          weight[6] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ;
          weight[7] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ;
          weight[8] = KFactor_QCD_qqZZ_M*KFactor_EW_qqZZ;
          weight[9] = KFactor_QCD_qqZZ_M*(KFactor_EW_qqZZ+KFactor_EW_qqZZ_unc);
          weight[10] = KFactor_QCD_qqZZ_M*(KFactor_EW_qqZZ-KFactor_EW_qqZZ_unc);
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[0][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)Inclusive) hSyst[icat][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)Inclusive) hSyst[icat_jesup][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)Inclusive) hSyst[icat_jesdn][isyst]->Fill(ZZMass, weight[isyst]*overallEventWeight*xsec/nGenEvents);
          }

        }
      }
    }

    // Rescale sum of categories in each systematic scenario to inclusive
    for (unsigned int isyst=0; isyst<nSyst; isyst++){
      for (int ix=0; ix<=hSyst[0][isyst]->GetNbinsX()+1; ix++){
        double sum_inc=hSyst[0][isyst]->GetBinContent(ix);
        double sum_cat=0;
        for (unsigned int icat=1; icat<nCategories; icat++) sum_cat += hSyst[icat][isyst]->GetBinContent(ix);
        double rescale=1;
        cout << "isyst=" << isyst << " inc=" << sum_inc << " sum(cat)=" << sum_cat << endl;
        if (sum_cat!=0.) rescale = sum_inc / sum_cat;
        for (unsigned int icat=1; icat<nCategories; icat++){
          hSyst[icat][isyst]->SetBinContent(ix, hSyst[icat][isyst]->GetBinContent(ix)*rescale);
          hSyst[icat][isyst]->SetBinError(ix, hSyst[icat][isyst]->GetBinError(ix)*rescale);
        }
      }
    }

    // Print the ratio of each integral to inclusive nominal into the output file
    float integral_nom = hSyst[0][0]->Integral(1, hSyst[0][0]->GetNbinsX());
    for (unsigned int icat=0; icat<nCategories; icat++){
      tout << strProcess[is][1] << " " << strCategory[icat] << " ";
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        double integral = hSyst[icat][isyst]->Integral(1, hSyst[icat][isyst]->GetNbinsX());
        if (integral_nom!=0.) integral /= integral_nom;
        else integral=0;
        tout << integral << " ";
      }
      tout << endl;
    }
    tout << endl;

    for (unsigned int icat=0; icat<nCategories; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        if (icat==0 && isyst==0) continue; // Do not rescale inclusive nominal since you need this to divide everything by.
        for (int ix=0; ix<=hSyst[icat][isyst]->GetNbinsX()+1; ix++){
          double val = hSyst[icat][isyst]->GetBinContent(ix);
          double err = hSyst[icat][isyst]->GetBinError(ix);
          double val_nom = hSyst[0][0]->GetBinContent(ix);
          if (val_nom!=0.){ val /= val_nom; err /= val_nom; }
          else{ val=0; err=0; }
          hSyst[icat][isyst]->SetBinContent(ix, val);
          hSyst[icat][isyst]->SetBinError(ix, err);
        }
      }
    }

    for (unsigned int icat=0; icat<nCategories; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        foutput->WriteTObject(hSyst[icat][isyst]);
        delete hSyst[icat][isyst];
      }
    }

  }


  for (unsigned int is=0; is<nProcesses; is++){
    for (unsigned int ifile=0; ifile<fileList[is].size(); ifile++){
      if (fileList[is].at(ifile)!=0 && fileList[is].at(ifile)->IsOpen()) fileList[is].at(ifile)->Close();
    }
  }

  foutput->Close();
  tout.close();
}


void progressbar(int val, int tot){
  int percent=floor(0.01*tot);
  if (percent==0) percent=1;

  if (val%percent==0 && val!=tot){
    cout<<"[ "<<setw(3)<<val/percent<<"% |";
    for (int k=1; k<val/percent; k++) cout<<"=";
    if (val%percent!=100) cout<<">";
    for (int k=val/percent+1; k<100; k++) cout<<" ";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
  else if (val==tot){
    cout<<"[ 100% |";
    for (int k=0; k<100; k++) cout<<"=";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
}

TString todaysdate(){
  TString result;
  time_t now = time(0);
  tm* ltm = localtime(&now);
  int year = 1900 + ltm->tm_year-2000;
  int month = 1 + ltm->tm_mon;
  int day = ltm->tm_mday;
  if(month>=10) result = Form("%i%i%i", day, month, year);
  else result = Form("%i%i%i%i", day, 0, month, year);
  return result;
}


void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray){
  bool inserted = false;
  for (std::vector<std::pair<unsigned int, float>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).second>=val){
      inserted=true;
      valArray.insert(it, std::pair<unsigned int, float>(ev, val));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<unsigned int, float>(ev, val));
}

void addByLowest(float weight, float val, std::vector<std::pair<float, float>>& valArray){
  bool inserted = false;
  for (std::vector<std::pair<float, float>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).second>=val){
      inserted=true;
      valArray.insert(it, std::pair<float, float>(weight, val));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<float, float>(weight, val));
}

void addByLowest(float val, std::vector<float>& valArray){
  bool inserted = false;
  for (std::vector<float>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it)>=val){
      inserted=true;
      valArray.insert(it, val);
      break;
    }
  }
  if (!inserted) valArray.push_back(val);
}


vector<pair<string, double>> getFileList(string strinput){
  ifstream fin;
  string filename = user_gg2VV_location.Data();
  filename += "/"; filename += strinput;
  fin.open(filename);
  cout << "Checking list from " << filename << endl;
  vector<pair<string, double>> result;
  if (fin.good()){
    if (filename.find("Bkg")==string::npos){ // Signal or BSI samples
      while (!fin.eof()){
        string strtmp;
        double mass=-1;
        fin >> strtmp >> mass;
        if (strtmp!="") result.push_back(pair<string, double>(strtmp, mass));
        else break;
        cout << "Signal file " << strtmp << " with mass " << mass << endl;
      }
    }
    else{
      while (!fin.eof()){
        string strtmp;
        double mass=-1;
        fin >> strtmp;
        if (strtmp!="") result.push_back(pair<string, double>(strtmp, mass));
        else break;
        cout << "Bkg. file " << strtmp << " with mass " << mass << endl;
      }
    }
  }
  fin.close();
  return result;
}

int getCategory(
  float ZZMass,
  short nExtraLep,
  short nExtraZ,
  std::vector<float> JetPt,
  std::vector<float> JetPhi,
  std::vector<float> JetSigma,
  std::vector<float> JetIsBtagged,
  std::vector<float> JetQGLikelihood,
  float phjj_VAJHU_highestPTJets,
  float phj_VAJHU,
  float pvbf_VAJHU_highestPTJets,
  float pAux_vbf_VAJHU,
  float pwh_hadronic_VAJHU,
  float pzh_hadronic_VAJHU,
  int categorizationType,
  int JESvariation
  ){

  if (categorizationType==0) return Inclusive;

  float WP_VBF2j=0, WP_VBF1j=0, WP_hadWH=0, WP_hadZH=0;
  if (categorizationType==1){ // No QG tagging, 6 categories
    WP_VBF2j = 1.043-460./(ZZMass+634.);
    WP_VBF1j = 0.699;
    WP_hadWH = 0.959;
    WP_hadZH = 0.9946;
  }
  else if (categorizationType==2){ // Use QG tagging, 6 categories
    WP_VBF2j = 0.391;
    WP_VBF1j = 0.720;
    WP_hadWH = 0.973;
    WP_hadZH = 0.996;
  }
  else if (categorizationType==-1){ // No QG tagging, 2 categories
    WP_VBF2j = 1.043-460./(ZZMass+634.);
  }
  else if (categorizationType==-2){ // Use QG tagging, 2 categories
    WP_VBF2j = 0.391;
  }

  float c_Mela2j = getDVBF2jetsConstant(ZZMass);
  float c_Mela1j = getDVBF1jetConstant(ZZMass);
  float c_MelaWH = 100000.;
  float c_MelaZH = 10000.;

  vector<float> jetPgOverPq;
  unsigned int nCleanedJetsPt30=0;
  unsigned int nCleanedJetsPt30BTagged=0;
  for (unsigned int j=0; j<JetPt.size(); j++){
    float jetpt = JetPt.at(j) * (1.+JetSigma.at(j)*((float)JESvariation));
    if (jetpt>=30.){
      nCleanedJetsPt30++;
      if (JetIsBtagged.at(j)>0.) nCleanedJetsPt30BTagged++;
      if (JetQGLikelihood.at(j)<0.){
        TRandom3 rand;
        rand.SetSeed(abs(static_cast<int>(sin(JetPhi.at(j))*100000)));
        jetPgOverPq.push_back(1./rand.Uniform() - 1.);
      }
      else jetPgOverPq.push_back(1./JetQGLikelihood.at(j) - 1.);
    }
  }

  float D_VBF2j = (nCleanedJetsPt30>=2) ? 1/(1+ c_Mela2j*phjj_VAJHU_highestPTJets/pvbf_VAJHU_highestPTJets * (abs(categorizationType)==2 ? TMath::Power(jetPgOverPq.at(0)*jetPgOverPq.at(1), 1/3.) : 1.)) : -2;
  float D_VBF1j = (nCleanedJetsPt30>=1) ? 1/(1+ (c_Mela1j*phj_VAJHU)/(pvbf_VAJHU_highestPTJets*pAux_vbf_VAJHU) * (abs(categorizationType)==2 ? TMath::Power(jetPgOverPq.at(0), 1/3.) : 1.)) : -2;
  float D_hadWH = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaWH*phjj_VAJHU_highestPTJets/pwh_hadronic_VAJHU * (abs(categorizationType)==2 ? TMath::Power(jetPgOverPq.at(0)*jetPgOverPq.at(1), 1/3.) : 1.)) : -2;
  float D_hadZH = (nCleanedJetsPt30>=2) ? 1/(1+ c_MelaZH*phjj_VAJHU_highestPTJets/pzh_hadronic_VAJHU * (abs(categorizationType)==2 ? TMath::Power(jetPgOverPq.at(0)*jetPgOverPq.at(1), 1/3.) : 1.)) : -2;

  if (categorizationType>0){
    if (nExtraLep==0 && nCleanedJetsPt30==1 && D_VBF1j>WP_VBF1j) return VBF1jTagged;
    else if (nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged==0)) && D_VBF2j>WP_VBF2j) return VBF2jTagged;
    else if ((nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged==0)) && (D_hadWH>WP_hadWH||D_hadZH>WP_hadZH)) || (nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3) && nCleanedJetsPt30BTagged>=2)) return HadVHTagged;
    else if ((nCleanedJetsPt30<=3 && nCleanedJetsPt30BTagged==0 && (nExtraLep==1||nExtraZ>=1)) || (nCleanedJetsPt30==0 && nExtraLep>=1)) return LepVHTagged;
    else if ((nCleanedJetsPt30>=4 && nCleanedJetsPt30BTagged>=1) || nExtraLep>=1) return ttHTagged;
    else return Untagged;
  }
  else if (categorizationType<0){
    if (nCleanedJetsPt30>=2 && D_VBF2j>WP_VBF2j) return VBF2jTagged;
    else return Untagged;
  }
  else return Inclusive;
}



