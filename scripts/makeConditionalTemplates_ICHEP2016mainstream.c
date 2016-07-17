#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
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
#include <ZZMatrixElement/MELA/interface/Mela.h>

using namespace std;

// Constants to affect the template code
const TString fixedDate="";

// Initializers
void makeConditionalTemplates_ICHEP2016mainstream_one(int channel, int erg_tev, int Systematics);
void progressbar(int val, int tot);
vector<pair<string, double>> getFileList(string strinput);
void convertTGraphErrorsToTH1F(TGraphErrors* tg, TH1F* histo);
void convertTGraphAsymmErrorsToTH1F(TGraphAsymmErrors* tg, TH1F* histo);
vector<pair<TH1F*, TH1F*>> getZXFR_SS();
TString todaysdate();
SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPtEtaPhi[3][4]);
void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray);
bool test_bit(int mask, unsigned int iBit);

// Main Function, runs over all desired iterations
void makeConditionalTemplates_ICHEP2016mainstream(){
  for (int CoM=13; CoM<=13; CoM++){ for (int channel=0; channel<3; channel++){ for (int syst=0; syst<=0; syst++) makeConditionalTemplates_ICHEP2016mainstream_one(channel, CoM, syst); } }
}

// Function to build one template
// folder = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// erg_tev = 13 (CoM energy)
void makeConditionalTemplates_ICHEP2016mainstream_one(int folder, int erg_tev, int Systematics){
  float mPOLE=125.;
  float wPOLE=4.07e-3;
  int EnergyIndex=(erg_tev==13 ? 0 : 0);

  const float ZZMass_cut[2]={ 100., 3000. };
  const float xwidth = 2;
  const int nbinsx = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth;
  float kDXarray[nbinsx+1];
  for (int bin=0; bin<nbinsx+1; bin++) kDXarray[bin] = ZZMass_cut[0] + xwidth*bin;
  const int nbinsy=30;
  float kDYarray[nbinsy+1];
  const float kDY_bounds[2]={ 0, 1 };
  for (int bin=0; bin<nbinsy+1; bin++){
    float binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
    kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
  }


  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString INPUT_NAME = user_treefile;
  TString OUTPUT_NAME = Form("HtoZZ%s_ConditionalTemplatesForCombine_", user_folder[folder].Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else if (Systematics == 1) OUTPUT_NAME += "SysUp_QCD";
  else if (Systematics == -1) OUTPUT_NAME += "SysDown_QCD";
  else if (Systematics == 2) OUTPUT_NAME += "SysUp_PDF";
  else if (Systematics == -2) OUTPUT_NAME += "SysDown_PDF";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = user_dir + erg_dir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  float varKD=0;
  float weight=1;
  float overallEventWeight=1;
  float xsec=1;
  float reweight=1;

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

  short Z1Flav, Z2Flav;
  int ZZFlav;
  float ZZMass, p0plus_VAJHU, bkg_VAMCFM;
  float GenHMass;
  short GenLepId[4];
  float GenLepPtEtaPhi[3][4];

  TVar::VerbosityLevel verbosity = TVar::SILENT;
  Mela* mela = new Mela(erg_tev, mPOLE, verbosity);

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << user_folder[folder] << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << user_folder[folder] << endl;
  tout << "===============================" << endl;
  tout << endl;

  // Get list of samples
  const unsigned int nProcesses = 8;
  const TString strProcess[nProcesses][2]={
    { "Samples_gg_Sig_POWHEG.txt", "gg_Sig" },
    { "Samples_VBF_Sig_POWHEG.txt", "VBF_Sig" },
    { "Samples_VBF_Sig_Phantom.txt", "VBF_Sig" },
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
      theTree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      theTree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      if (theTree->GetBranchStatus("GenHMass")){
        theTree->SetBranchAddress("GenHMass", &GenHMass);
        for (unsigned int ilep=0; ilep<4; ilep++){
          theTree->SetBranchAddress(Form("GenLep%iId", ilep+1), &(GenLepId[ilep]));
          theTree->SetBranchAddress(Form("GenLep%iPt", ilep+1), &(GenLepPtEtaPhi[0][ilep]));
          theTree->SetBranchAddress(Form("GenLep%iEta", ilep+1), &(GenLepPtEtaPhi[1][ilep]));
          theTree->SetBranchAddress(Form("GenLep%iPhi", ilep+1), &(GenLepPtEtaPhi[2][ilep]));
        }
      }
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

      treeList[is].push_back(pair<TTree*, TH1F*>(theTree, theCounters));
    }
  }

  foutput->cd();
  const unsigned int nTemplates = 3;
  TString strTemplates[nTemplates]={
    "Sig",
    "ggBkg",
    "qqBkg"
  };
  TTree* rawtree[nTemplates];
  TTree* rewtree[nTemplates];
  TTree* newtree[nTemplates];

  TH1F* D_temp_1D_raw[nTemplates];
  TH1F* D_temp_1D_rew[nTemplates];
  TH1F* N_temp_1D_raw[nTemplates];
  TH1F* N_temp_1D_rew[nTemplates];
  TH2F* D_temp_2D[nTemplates];
  for (unsigned int t=0; t<nTemplates; t++){
    rawtree[t] = new TTree(Form("T_2D_%s_rwgt", strTemplates[t].Data()), "");
    rawtree[t]->Branch("ZZMass", &ZZMass);
    //rawtree[t]->Branch("GenHMass", &GenHMass);
    rawtree[t]->Branch("KD", &varKD);
    rawtree[t]->Branch("weight", &weight);
    D_temp_1D_raw[t] = new TH1F(Form("h%s_raw", strTemplates[t].Data()), "", nbinsx, kDXarray);
    N_temp_1D_raw[t] = new TH1F(Form("hn%s_raw", strTemplates[t].Data()), "", nbinsx, kDXarray);

    rewtree[t] = new TTree(Form("T_2D_%s_rwgt", strTemplates[t].Data()), "");
    rewtree[t]->Branch("ZZMass", &ZZMass);
    //rewtree[t]->Branch("GenHMass", &GenHMass);
    rewtree[t]->Branch("KD", &varKD);
    rewtree[t]->Branch("weight", &weight);
    rewtree[t]->Branch("reweight", &reweight);
    D_temp_1D_rew[t] = new TH1F(Form("h%s_rew", strTemplates[t].Data()), "", nbinsx, kDXarray);
    N_temp_1D_rew[t] = new TH1F(Form("hn%s_rew", strTemplates[t].Data()), "", nbinsx, kDXarray);

    newtree[t] = new TTree(Form("T_2D_%s", strTemplates[t].Data()), "");
    newtree[t]->Branch("ZZMass", &ZZMass);
    //newtree[t]->Branch("GenHMass", &GenHMass);
    newtree[t]->Branch("KD", &varKD);
    newtree[t]->Branch("weight", &weight);
    newtree[t]->Branch("reweight", &reweight);

    D_temp_2D[t] = new TH2F(Form("h%s", strTemplates[t].Data()), "", nbinsx, kDXarray, nbinsy, kDYarray);
    D_temp_2D[t]->Sumw2();
    D_temp_2D[t]->SetOption("colz");
    D_temp_2D[t]->GetYaxis()->SetTitle("D^{kin}_{bkg}");
    D_temp_2D[t]->GetXaxis()->SetTitle("m_{4l} (GeV)");
  }
  for (unsigned int is=0; is<nProcesses; is++){
    cout << "Process " << is << " under scrutiny" << endl;
    tout << "Process " << is << " under scrutiny" << endl;
    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      TTree* theTree = treeList[is].at(it).first;
      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;
      tout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;

      for (int ev=0; ev<nEntries; ev++){
        weight=1;
        overallEventWeight=1;
        xsec=1;
        reweight=1;

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
        KFactor_EW_qqZZ_unc = 1;
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
        if (strProcess[is][0].Contains("POWHEG") && mass>0 && mass<300 && GenHMass>mass+5.){
          //cout << "Process " << is << " sample " << it << " has mH = " << mass << " and mH* = " << GenHMass << ">mH+5 GeV. Ignoring event..." << endl;
          continue;
        }

        ZZFlav=Z1Flav*Z2Flav;
        if (folder==0 && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (folder==1 && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (folder==2 && abs(ZZFlav)!=pow(13*11, 2)) continue;

        varKD = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM);
        if (varKD!=varKD) continue;

        // weight is product of everything. Only some are non-1.
        weight = KFactor_QCD_ggZZ_Nominal*KFactor_EW_qqZZ*KFactor_QCD_qqZZ_M*overallEventWeight/**xsec*/;
        if (weight<0.) continue; // No negative weights!

        // Protect against any KD exceeding boundaries
        if (varKD>=kDYarray[nbinsy]){
          cout << "Found varKD == " << varKD;
          varKD = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
          cout << ". Corrected to have " << varKD << endl;
        }
        if (varKD<kDYarray[0]) varKD = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*ev/nEntries;
        if (varKD>=kDYarray[nbinsy] || varKD<kDYarray[0]) cout << "Fix has been numerically unsuccessful." << endl;

        unsigned int index=0;
        if (strProcess[is][1]=="gg_Bkg") index = 1;
        else if (strProcess[is][1]=="qq_Bkg") index = 2;

        rawtree[index]->Fill();
        D_temp_1D_raw[index]->Fill(ZZMass, weight);
        N_temp_1D_raw[index]->Fill(ZZMass, 1);

        if (strProcess[is][1]=="gg_Sig"){
          SimpleParticleCollection_t daughters = constructFourVectors(GenLepId, GenLepPtEtaPhi);
          mela->setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

          float sampleWeight=1;
          mela->setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
          mela->setMelaHiggsMassWidth(GenHMass, 1./GenHMass, 0); // Set prop=1
          mela->computeP(sampleWeight, false);

          float targetWeight=1;
          mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
          mela->computeP(targetWeight, false);

          mela->resetInputEvent();

          if (sampleWeight>0.){
            reweight=targetWeight/sampleWeight;
            index=1;
            rewtree[index]->Fill();
            D_temp_1D_rew[index]->Fill(ZZMass, weight*reweight);
            N_temp_1D_rew[index]->Fill(ZZMass, 1);
          }
        }

      }
    }
  }

  for (unsigned int t=0; t<nTemplates; t++){
    // Merge raw and reweighted trees by sum of weights and fill the histograms
    for (int ev=0; ev<rawtree[t]->GetEntries(); ev++){
      rawtree[t]->GetEntry(ev);
      reweight=1;
      int bin = D_temp_1D_raw[t]->GetXaxis()->FindBin(ZZMass);
      double rescale = 1.;
      if (D_temp_1D_raw[t]->GetBinContent(bin)>0.) rescale = N_temp_1D_raw[t]->GetBinContent(bin)/D_temp_1D_raw[t]->GetBinContent(bin);
      reweight *= rescale;
      newtree[t]->Fill();
      D_temp_2D[t]->Fill(ZZMass, varKD, weight*reweight);
    }
    for (int ev=0; ev<rewtree[t]->GetEntries(); ev++){
      rewtree[t]->GetEntry(ev);
      int bin = D_temp_1D_raw[t]->GetXaxis()->FindBin(ZZMass);
      double rescale = 1.;
      if (D_temp_1D_rew[t]->GetBinContent(bin)>0.) rescale = N_temp_1D_rew[t]->GetBinContent(bin)/D_temp_1D_rew[t]->GetBinContent(bin);
      reweight *= rescale;
      newtree[t]->Fill();
      D_temp_2D[t]->Fill(ZZMass, varKD, weight*reweight);
    }

    // Conditionalize the histograms
    for (int binx=0; binx<=D_temp_2D[t]->GetNbinsX()+1; binx++){
      double integral = D_temp_2D[t]->Integral(binx, binx, 0, D_temp_2D[t]->GetNbinsY()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int biny=1; biny<=D_temp_2D[t]->GetNbinsY(); biny++){
        D_temp_2D[t]->SetBinContent(binx, biny, D_temp_2D[t]->GetBinContent(binx, biny)/integral);
        D_temp_2D[t]->SetBinError(binx, biny, D_temp_2D[t]->GetBinError(binx, biny)/integral);
      }
    }
  }

  for (unsigned int t=0; t<nTemplates; t++){
    foutput->WriteTObject(newtree[t]);
    foutput->WriteTObject(D_temp_2D[t]);
    delete N_temp_1D_rew[t];
    delete N_temp_1D_raw[t];
    delete D_temp_1D_rew[t];
    delete D_temp_1D_raw[t];
    delete D_temp_2D[t];
    delete newtree[t];
    delete rewtree[t];
    delete rawtree[t];
  }

  for (unsigned int is=0; is<nProcesses; is++){
    for (unsigned int ifile=0; ifile<fileList[is].size(); ifile++){
      if (fileList[is].at(ifile)!=0 && fileList[is].at(ifile)->IsOpen()) fileList[is].at(ifile)->Close();
    }
  }

  foutput->Close();
  tout.close();
  if (mela!=0) delete mela;
}

// Function to build one template for Z+X
// folder = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// erg_tev = 13 (CoM energy)
void makeConditionalTemplates_ICHEP2016mainstream_ZX_one(int folder, int erg_tev, int Systematics){
  int EnergyIndex=(erg_tev==13 ? 0 : 0);

  const float ZZMass_cut[2]={ 100., 3000. };
  const float xwidth = 2;
  const int nbinsx = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth;
  float kDXarray[nbinsx+1];
  for (int bin=0; bin<nbinsx+1; bin++) kDXarray[bin] = ZZMass_cut[0] + xwidth*bin;
  const int nbinsy=30;
  float kDYarray[nbinsy+1];
  const float kDY_bounds[2]={ 0, 1 };
  for (int bin=0; bin<nbinsy+1; bin++){
    float binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
    kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
  }


  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString INPUT_NAME = user_treefile;
  TString OUTPUT_NAME = Form("HtoZZ%s_ConditionalTemplatesForCombine_ZX_", user_folder[folder].Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = user_dir + erg_dir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  float varKD=0;
  float weight=1;

  int CRflag=666;
  short Z1Flav, Z2Flav;
  int ZZFlav;
  float ZZMass, p0plus_VAJHU, bkg_VAMCFM;
  std::vector<short>* LepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << user_folder[folder] << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << user_folder[folder] << endl;
  tout << "===============================" << endl;
  tout << endl;

  vector<pair<TH1F*, TH1F*>> hFR_ZX_SS = getZXFR_SS();
  if (hFR_ZX_SS.size()!=2){
    cerr << "ERROR: ZX FR size " << hFR_ZX_SS.size() << "!=2. Aborting..." << endl;
    tout << "ERROR: ZX FR size " << hFR_ZX_SS.size() << "!=2. Aborting..." << endl;
    for (unsigned int ifr=0; ifr<hFR_ZX_SS.size(); ifr++){ delete hFR_ZX_SS.at(ifr).first; delete hFR_ZX_SS.at(ifr).second; }
    foutput->Close();
    tout.close();
    return;
  }
  cout << "Retrieved ZX SS FRs" << endl;
  tout << "Retrieved ZX SS FRs" << endl;

  // Get list of samples
  const unsigned int nProcesses = 1;
  const TString strProcess[nProcesses][2]={
    { "AllData", "ZX" }
  };
  vector<pair<string, double>> fileNameList[nProcesses];
  for (unsigned int is=0; is<nProcesses; is++) fileNameList[is].push_back(pair<string, double>(string(strProcess[is][0].Data()), -1));
  vector<TFile*> fileList[nProcesses];
  vector<pair<TTree*, TH1F*>> treeList[nProcesses];
  for (unsigned int is=0; is<nProcesses; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_gg2VV_location + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;
      tout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      TFile* finput = TFile::Open(cinput.c_str(), "read"); finput->cd();
      TTree* theTree = (TTree*)finput->Get(user_CRtreename);
      TH1F* theCounters = (TH1F*)finput->Get(user_CRcountersname);

      fileList[is].push_back(finput);
      foutput->cd();

      theTree->SetBranchAddress("CRflag", &CRflag);
      theTree->SetBranchAddress("ZZMass", &ZZMass);
      theTree->SetBranchAddress("Z1Flav", &Z1Flav);
      theTree->SetBranchAddress("Z2Flav", &Z2Flav);
      theTree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      theTree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      theTree->SetBranchAddress("LepLepId", &LepId);
      theTree->SetBranchAddress("LepPt", &LepPt);
      theTree->SetBranchAddress("LepEta", &LepEta);

      treeList[is].push_back(pair<TTree*, TH1F*>(theTree, theCounters));
    }
  }

  foutput->cd();
  const unsigned int nTemplates = 1;
  TString strTemplates[nTemplates]={
    strProcess[0][1]
  };
  TTree* newtree[nTemplates];
  TH2F* D_temp_2D[nTemplates];
  for (unsigned int t=0; t<nTemplates; t++){
    newtree[t] = new TTree(Form("T_2D_%s", strTemplates[t].Data()), "");
    newtree[t]->Branch("ZZMass", &ZZMass);
    newtree[t]->Branch("KD", &varKD);
    newtree[t]->Branch("weight", &weight);

    D_temp_2D[t] = new TH2F(Form("h%s", strTemplates[t].Data()), "", nbinsx, kDXarray, nbinsy, kDYarray);
    D_temp_2D[t]->Sumw2();
    D_temp_2D[t]->SetOption("colz");
    D_temp_2D[t]->GetYaxis()->SetTitle("D^{kin}_{bkg}");
    D_temp_2D[t]->GetXaxis()->SetTitle("m_{4l} (GeV)");
  }
  for (unsigned int is=0; is<nProcesses; is++){
    cout << "Process " << is << " under scrutiny" << endl;
    tout << "Process " << is << " under scrutiny" << endl;
    for (unsigned int it=0; it<treeList[is].size(); it++){
      TTree* theTree = treeList[is].at(it).first;
      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. Tree name: " << theTree->GetName() << ", #Events = " << nEntries << endl;
      tout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. Tree name: " << theTree->GetName() << ", #Events = " << nEntries << endl;

      for (int ev=0; ev<nEntries; ev++){
        weight=1;

        theTree->GetEntry(ev);
        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Best CR selection
        if (!test_bit(CRflag, (unsigned int)CJLST_FinalState::CRZLLss)) continue;

        ZZFlav=Z1Flav*Z2Flav;
        if (folder==0 && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (folder==1 && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (folder==2 && abs(ZZFlav)!=pow(13*11, 2)) continue;

        varKD = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM);
        if (varKD!=varKD) continue;

        // Apply FRs in SS
        if (CRflag==CJLST_FinalState::CRZLLss){
          for (unsigned int ilep=2; ilep<=3; ilep++){
            unsigned int whichLep = 1;
            if (abs(LepId->at(ilep))==11) whichLep=0;
            else if (abs(LepId->at(ilep))==13) whichLep=1;
            if (
              (fabs(LepEta->at(ilep))<1.2 && whichLep==0)
              ||
              (fabs(LepEta->at(ilep))<1.449 && whichLep==1)
              ){
              int bin = hFR_ZX_SS.at(whichLep).first->GetXaxis()->FindBin(LepPt->at(ilep));
              if (bin>hFR_ZX_SS.at(whichLep).first->GetNbinsX()) bin=hFR_ZX_SS.at(whichLep).first->GetNbinsX();
              weight *= hFR_ZX_SS.at(whichLep).first->GetBinContent(bin);
            }
            else{
              int bin = hFR_ZX_SS.at(whichLep).second->GetXaxis()->FindBin(LepPt->at(ilep));
              if (bin>hFR_ZX_SS.at(whichLep).second->GetNbinsX()) bin=hFR_ZX_SS.at(whichLep).second->GetNbinsX();
              weight *= hFR_ZX_SS.at(whichLep).second->GetBinContent(bin);
            }
          }
        }

        // Protect against any KD exceeding boundaries
        if (varKD>=kDYarray[nbinsy]){
          cout << "Found varKD == " << varKD;
          varKD = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
          cout << ". Corrected to have " << varKD << endl;
        }
        if (varKD<kDYarray[0]) varKD = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*ev/nEntries;
        if (varKD>=kDYarray[nbinsy] || varKD<kDYarray[0]) cout << "Fix has been numerically unsuccessful." << endl;

        const unsigned int index=0;
        newtree[index]->Fill();
        D_temp_2D[index]->Fill(ZZMass, varKD, weight);
      }
    }
  }

  cout << "Filled the trees and templates" << endl;
  tout << "Filled the trees and templates" << endl;

  // Not enough statistics for the moment, so use projection on Y
  for (unsigned int t=0; t<nTemplates; t++){
    for (int biny=0; biny<=D_temp_2D[t]->GetNbinsY()+1; biny++){
      double integral = D_temp_2D[t]->Integral(0, D_temp_2D[t]->GetNbinsX()+1, biny, biny);
      for (int binx=0; binx<=D_temp_2D[t]->GetNbinsX()+1; binx++){
        D_temp_2D[t]->SetBinContent(binx, biny, integral);
        D_temp_2D[t]->SetBinError(binx, biny, integral);
      }
    }
  }

  // Conditionalize the histograms
  for (unsigned int t=0; t<nTemplates; t++){
    for (int binx=0; binx<=D_temp_2D[t]->GetNbinsX()+1; binx++){
      double integral = D_temp_2D[t]->Integral(binx, binx, 0, D_temp_2D[t]->GetNbinsY()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int biny=1; biny<=D_temp_2D[t]->GetNbinsY(); biny++){
        D_temp_2D[t]->SetBinContent(binx, biny, D_temp_2D[t]->GetBinContent(binx, biny)/integral);
        D_temp_2D[t]->SetBinError(binx, biny, D_temp_2D[t]->GetBinError(binx, biny)/integral);
      }
    }
  }

  cout << "Templates conditionally normalized" << endl;
  tout << "Templates conditionally normalized" << endl;

  for (unsigned int ifr=0; ifr<hFR_ZX_SS.size(); ifr++){
    cout << "Deleting hFR_ZX_SS[" << ifr << "]" << endl;
    gROOT->cd();
    TH1F* hfirst = hFR_ZX_SS.at(ifr).first;
    //if (hfirst!=0) delete hfirst;
    cout << "Deleted hFR_ZX_SS[" << ifr << "] first" << endl;
    TH1F* hsecond = hFR_ZX_SS.at(ifr).second;
    //if (hsecond!=0) delete hsecond;
    cout << "Deleted hFR_ZX_SS[" << ifr << "] second" << endl;
  }

  cout << "Deleted hFR_ZX" << endl;
  tout << "Deleted hFR_ZX" << endl;

  foutput->cd();
  for (unsigned int t=0; t<nTemplates; t++){
    foutput->WriteTObject(newtree[t]);
    foutput->WriteTObject(D_temp_2D[t]);
    delete D_temp_2D[t];
    delete newtree[t];
  }

  cout << "Deleted objects" << endl;
  tout << "Deleted objects" << endl;

  for (unsigned int is=0; is<nProcesses; is++){
    for (unsigned int ifile=0; ifile<fileList[is].size(); ifile++){
      if (fileList[is].at(ifile)!=0 && fileList[is].at(ifile)->IsOpen()) fileList[is].at(ifile)->Close();
    }
  }

  cout << "Deleted file list" << endl;
  tout << "Deleted file list" << endl;

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


SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPtEtaPhi[3][4]){
  SimpleParticleCollection_t result;
  for (int idau=0; idau<4; idau++){
    float mass=0;
    if (abs(GenLepId[idau])==11) mass=0.000511;
    else if (abs(GenLepId[idau])==13) mass=0.10566;
    else if (abs(GenLepId[idau])==15) mass=1.7768;
    else if (abs(GenLepId[idau])==5) mass=4.75;
    else if (abs(GenLepId[idau])==4) mass=1.275;
    TLorentzVector mom; mom.SetPtEtaPhiM((double)GenLepPtEtaPhi[0][idau], (double)GenLepPtEtaPhi[1][idau], (double)GenLepPtEtaPhi[2][idau], (double)mass);
    result.push_back(SimpleParticle_t((int)GenLepId[idau], mom));
  }
  return result;
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


void convertTGraphErrorsToTH1F(TGraphErrors* tg, TH1F* histo){
  if (tg==0){
    cerr << "convertTGraphErrorsToTH1F: TGraph is 0!" << endl;
    histo=0;
    return;
  }

  double* xx = tg->GetX();
  double* yy = tg->GetY();
  double* ex = tg->GetEX();
  double* ey = tg->GetEY();
  // This is just stupid.
  /*
  double* ex_up = tg->GetEXhigh();
  double* ex_dn = tg->GetEXlow();
  double* ey_up = tg->GetEYhigh();
  double* ey_dn = tg->GetEYlow();
  */
  const int nbins = tg->GetN();
  double xarray[nbins+1];
  /*
  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex_dn[ix];
  xarray[nbins] = xx[nbins-1]+ex_up[nbins-1];
  */
  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex[ix];
  xarray[nbins] = xx[nbins-1]+ex[nbins-1];
  gROOT->cd();
  histo = new TH1F(Form("h_%s", tg->GetName()), "", nbins, xarray);
  for (int ix=1; ix<=nbins; ix++){
    cout << "x, y = " << xarray[ix-1] << " " << yy[ix-1] << endl;
    histo->SetBinContent(ix, yy[ix-1]);
    //histo->SetBinError(ix, sqrt((pow(ey_dn[ix-1], 2)+pow(ey_up[ix-1], 2))*0.5));
    histo->SetBinError(ix, ey[ix-1]);
  }
}
void convertTGraphAsymmErrorsToTH1F(TGraphAsymmErrors* tg, TH1F* histo){
  if (tg==0){
    cerr << "convertTGraphAsymmErrorsToTH1F: TGraph is 0!" << endl;
    histo=0;
    return;
  }

  double* xx = tg->GetX();
  double* yy = tg->GetY();
  //double* ex = tg->GetEX();
  //double* ey = tg->GetEY();
  double* ex_up = tg->GetEXhigh();
  double* ex_dn = tg->GetEXlow();
  double* ey_up = tg->GetEYhigh();
  double* ey_dn = tg->GetEYlow();
  const int nbins = tg->GetN();
  double xarray[nbins+1];
  
  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex_dn[ix];
  xarray[nbins] = xx[nbins-1]+ex_up[nbins-1];
  
  //for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex[ix];
  //xarray[nbins] = xx[nbins-1]+ex[nbins-1];
  gROOT->cd();
  histo = new TH1F(Form("h_%s", tg->GetName()), "", nbins, xarray);
  for (int ix=1; ix<=nbins; ix++){
    cout << "x, y = " << xarray[ix-1] << " " << yy[ix-1] << endl;
    histo->SetBinContent(ix, yy[ix-1]);
    histo->SetBinError(ix, sqrt((pow(ey_dn[ix-1], 2)+pow(ey_up[ix-1], 2))*0.5));
    //histo->SetBinError(ix, ey[ix-1]);
  }
}

// Pairs are (e-EB, e-EE), (mu-EB, mu-EE)
vector<pair<TH1F*, TH1F*>> getZXFR_SS(){
  TString tgname[2][2]={
    { "CorrFR_TG_Lep3_pT_all_EB_afterMET_wzremoved_ALL_80XB", "CorrFR_TG_Lep3_pT_all_EE_afterMET_wzremoved_ALL_80XB" },
    { "TG_Lep3_pT_all_EB_afterMET_wzremoved_ALL_80XB", "TG_Lep3_pT_all_EE_afterMET_wzremoved_ALL_80XB" }
  };
  TFile* finput[2]={ 0 };
  finput[0] = TFile::Open("../data/computedfakerate_Z2e_80XB.root", "read");
  finput[1] = TFile::Open("../data/computedfakerate_Z2mu_80XB.root", "read");
  vector<pair<TH1F*, TH1F*>> result;
  for (int f=0; f<2; f++){
    if (!(finput[f]!=0 && finput[f]->IsOpen())){ cerr << "getZXFR_SS: File " << f << " is not open!" << endl; return result; }
    else cout << "getZXFR_SS: File " << f << " opened" << endl;
    TH1F* htmp[2];
    for (int t=0; t<2; t++){
      cout << "\tAttempting to acquire " << tgname[f][t] << "...";
      TGraphErrors* tg = 0;
      TGraphAsymmErrors* tga = 0;
      TObject* to = finput[f]->Get(tgname[f][t]);
      if (dynamic_cast<TGraphAsymmErrors*>(to)!=0){
        tga = dynamic_cast<TGraphAsymmErrors*>(to);
        tga->SetName(tgname[f][t]);
        cout << tga->GetName() << endl;
        convertTGraphAsymmErrorsToTH1F(tga, htmp[t]);
      }
      else if (dynamic_cast<TGraphErrors*>(to)!=0){
        tg = dynamic_cast<TGraphErrors*>(to);
        tg->SetName(tgname[f][t]);
        cout << tg->GetName() << endl;
        convertTGraphErrorsToTH1F(tg, htmp[t]);
      }
    }
    result.push_back(pair<TH1F*, TH1F*>(htmp[0], htmp[1]));
  }
  for (int f=0; f<2; f++){ if (finput[f]!=0 && finput[f]->IsOpen()) finput[f]->Close(); }
  return result;
}


bool test_bit(int mask, unsigned int iBit){ return (mask >> iBit) & 1; }