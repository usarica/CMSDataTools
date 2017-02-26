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
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include "InputTreeHandle.h"
#include "ZZ4l_Samples.h"
#include "external_Category.h"
#include "external_cConstants.h"

using namespace std;

// Constants to affect the template code
const TString fixedDate="";

// Initializers
void makeConditionalTemplates_Moriond2017mainstream_one(int channel, int erg_tev, int Systematics);
void makeConditionalTemplates_Moriond2017mainstream_ZX_one(int ichan, int erg_tev, int Systematics);
void collectConditionalTemplates_Moriond2017mainstream_one(int ichan, int erg_tev, int Systematics, bool doSmooth=true);
void progressbar(int val, int tot);
vector<pair<string, double>> getFileList(string strinput);
void convertTGraphErrorsToTH1F(TGraphErrors* tg, TH1F* histo);
void convertTGraphAsymmErrorsToTH1F(TGraphAsymmErrors* tg, TH1F* histo);
template <typename T> vector<pair<T*, T*>> getZXFR_SS();
template <typename T> double evaluateTObject(T* obj, float val);
template<> double evaluateTObject<TH1F>(TH1F* obj, float val);
template<> double evaluateTObject<TGraphErrors>(TGraphErrors* obj, float val);
TString todaysdate();
SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPt[4], float GenLepEta[4], float GenLepPhi[4]);
void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray);
void addByLowest(float weight, float val, std::vector<std::pair<float, float>>& valArray);
void addByLowest(float val, std::vector<float>& valArray);
void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);
bool test_bit(int mask, unsigned int iBit);
void regularizeSlice(TGraph* tgSlice, std::vector<double>* fixedX=0, double omitbelow=0., int nIter_=-1, double threshold_=-1);
void regularizeHistogram(TH2F* histo, int nIter_, double threshold_);
void regularizeHistogram(TH3F* histo, int nIter_, double threshold_);
void conditionalizeHistogram(TH2F* histo, unsigned int axis=0);
void conditionalizeHistogram(TH3F* histo, unsigned int axis=0);
void wipeOverUnderFlows(TH1F* hwipe);
void wipeOverUnderFlows(TH2F* hwipe);
void wipeOverUnderFlows(TH3F* hwipe);


// Main Function, runs over all desired iterations
void makeConditionalTemplates_Moriond2017mainstream(){
  for (int CoM=13; CoM<=13; CoM++){
    for (int channel=0; channel<3; channel++){
      for (int syst=0; syst<=0; syst++){
        makeConditionalTemplates_Moriond2017mainstream_one(channel, CoM, syst);
        makeConditionalTemplates_Moriond2017mainstream_ZX_one(channel, CoM, syst);
        collectConditionalTemplates_Moriond2017mainstream_one(channel, CoM, syst, true);
      }
    }
  }
}

// Function to build one template
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// erg_tev = 13 (CoM energy)
void makeConditionalTemplates_Moriond2017mainstream_one(int ichan, int erg_tev, int Systematics){
  if (ichan>(int)nChannels) return;
  const TString strChannel = channame(ichan);
  constructSampleTypeList();

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
  TString OUTPUT_NAME = Form("HtoZZ%s_ConditionalTemplatesForCombine_", strChannel.Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else if (Systematics == 1) OUTPUT_NAME += "SysUp_QCD";
  else if (Systematics == -1) OUTPUT_NAME += "SysDown_QCD";
  else if (Systematics == 2) OUTPUT_NAME += "SysUp_PDF";
  else if (Systematics == -2) OUTPUT_NAME += "SysDown_PDF";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString cinput_common = user_input_dir + "/";
  TString coutput_common = user_output_dir + erg_dir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  float varKD=0;
  float weight=1;
  float reweight=1;
  float ZZMass, GenHMass;

  TVar::VerbosityLevel verbosity = TVar::SILENT;

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << strChannel << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << strChannel << endl;
  tout << "===============================" << endl;
  tout << endl;

  const string path_HiggsWidthFile = "../../../ZZMatrixElement/MELA/data/HiggsTotalWidth_YR3.txt";
  cout << "Initializing MELAHXSWidth at path " << path_HiggsWidthFile.substr(0, path_HiggsWidthFile.length()-23) << endl;
  MELAHXSWidth xswidth(path_HiggsWidthFile.substr(0, path_HiggsWidthFile.length()-23));
  cout << "MELAHXSWidth initialized" << endl;

  // Get list of samples
  vector<pair<string, double>> fileNameList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++){
    cout << "Constructing file name list for " << sampleTypeList[is].first << endl;
    fileNameList[is] = getFileList(string(sampleTypeList[is].first.Data()));
  }
  vector<InputTreeHandle*> treeList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      //if (pickSampleWithMass>0. && !(fileNameList[is].at(ifile).second==-1. || fileNameList[is].at(ifile).second==pickSampleWithMass)) continue;

      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_input_dir + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;
      tout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      InputTreeHandle* theTree = new InputTreeHandle(TString(cinput.c_str()), user_treename, user_countersname, sampleTypeList[is].second);
      foutput->cd();
      treeList[is].push_back(theTree);
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

    newtree[t] = new TTree(Form("T_2D_%s_Tree", strTemplates[t].Data()), "");
    newtree[t]->Branch("ZZMass", &ZZMass);
    //newtree[t]->Branch("GenHMass", &GenHMass);
    newtree[t]->Branch("KD", &varKD);
    newtree[t]->Branch("weight", &weight);
    newtree[t]->Branch("reweight", &reweight);

    D_temp_2D[t] = new TH2F(Form("T_2D_%s", strTemplates[t].Data()), "", nbinsx, kDXarray, nbinsy, kDYarray);
    D_temp_2D[t]->Sumw2();
    D_temp_2D[t]->SetOption("colz");
    D_temp_2D[t]->GetYaxis()->SetTitle("D^{kin}_{bkg}");
    D_temp_2D[t]->GetXaxis()->SetTitle("m_{4l} (GeV)");
  }
  for (int is=0; is<nSampleTypes; is++){
    TString SampleTypeName = nameSampleType(is);
    cout << "Process " << SampleTypeName << "[" << is << "] under scrutiny" << endl;
    tout << "Process " << SampleTypeName << "[" << is << "] under scrutiny" << endl;

    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      double width=0;
      if (
        mass>0.
        &&
        is!=(int)SampleType::gg_Bkg_MCFM && is!=(int)SampleType::VBF_Sig_Phantom && is!=(int)SampleType::VBF_Bkg_Phantom
        ) width = xswidth.HiggsWidth(mass);
      else if (mass==125.) width = 4.07e-3; // Phantom and MCFM samples use 4.07 MeV.

      InputTreeHandle* theTree = treeList[is].at(it);

      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;
      tout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;
      cout << "\tIdentified (mass,width) = ( " << mass << " , " << width << " )" << endl;

      for (int ev=0; ev<nEntries; ev++){
        theTree->InitDefaults();
        theTree->GetEntry(ev);
        weight=1;

        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Fix for JHUGen v6 vegas bug in low mass samples
        if (SampleTypeName.Contains("POWHEG") && mass>0 && mass<300 && theTree->GenHMass>mass+5.) continue;

        if (theTree->ZZsel<70.) continue; // Disregard non-selected events (need this if skipEmptyEvents=false)

        if (theTree->ZZMass<ZZMass_cut[0] || theTree->ZZMass>=ZZMass_cut[1]) continue;

        int ZZFlav=theTree->Z1Flav*theTree->Z2Flav;
        if (ichan==(int)k4mu && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (ichan==(int)k4e && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (ichan==(int)k2e2mu && abs(ZZFlav)!=pow(13*11, 2)) continue;

        ZZMass = theTree->ZZMass;
        varKD = theTree->p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(theTree->p_GG_SIG_ghg2_1_ghz1_1_JHUGen + theTree->p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));
        if (varKD!=varKD) continue;

        // weight is product of everything. Only some are non-1.
        weight = theTree->KFactor_QCD_ggZZ_Nominal*theTree->KFactor_EW_qqZZ*theTree->KFactor_QCD_qqZZ_M*theTree->overallEventWeight/**xsec*/;
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
        if (is==(int)SampleType::gg_Bkg_MCFM) index = 1;
        else if (is==(int)SampleType::qq_Bkg) index = 2;

        rawtree[index]->Fill();
        D_temp_1D_raw[index]->Fill(ZZMass, weight);
        N_temp_1D_raw[index]->Fill(ZZMass, 1);

        if (
          (is>=(int)SampleType::gg_Sig && is<=(int)SampleType::gg_Sig_PythiaTuneUp_MiNLO)
          ){
          short GenLepId[4]={ theTree->GenLep1Id, theTree->GenLep2Id, theTree->GenLep3Id, theTree->GenLep4Id };
          float GenLepPt[4]={ theTree->GenLep1Pt, theTree->GenLep2Pt, theTree->GenLep3Pt, theTree->GenLep4Pt };
          float GenLepEta[4]={ theTree->GenLep1Eta, theTree->GenLep2Eta, theTree->GenLep3Eta, theTree->GenLep4Eta };
          float GenLepPhi[4]={ theTree->GenLep1Phi, theTree->GenLep2Phi, theTree->GenLep3Phi, theTree->GenLep4Phi };

          SimpleParticleCollection_t daughters = constructFourVectors(GenLepId, GenLepPt, GenLepEta, GenLepPhi);

          TLorentzVector sumP(0, 0, 0, 0); for (unsigned int idau=0; idau<daughters.size(); idau++) sumP = sumP + daughters.at(idau).second;

          float propagatorBW = 1./(pow(pow(sumP.M(), 2)-pow(mass, 2), 2)+pow(mass*width, 2));
          float sampleWeight=theTree->p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM/propagatorBW;

          float targetWeight=theTree->p_Gen_GG_BKG_MCFM;

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
    conditionalizeHistogram(D_temp_2D[t], 0);
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

  for (unsigned int is=0; is<nSampleTypes; is++){ for (unsigned int it=0; it<treeList[is].size(); it++) delete treeList[is].at(it); }

  foutput->Close();
  tout.close();
}

// Function to build one template for Z+X
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// erg_tev = 13 (CoM energy)
void makeConditionalTemplates_Moriond2017mainstream_ZX_one(int ichan, int erg_tev, int Systematics){
  if (ichan>(int)nChannels) return;
  const TString strChannel = channame(ichan);
  constructSampleTypeList();
  typedef TGraphErrors FRtype;

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

  TString INPUT_NAME = Form("HtoZZ%s_ConditionalTemplatesForCombine_ZX_", strChannel.Data());
  TString OUTPUT_NAME = Form("HtoZZ%s_ConditionalTemplatesForCombine_ZX_", strChannel.Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else return;

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString cinput_common = user_input_dir + "/";
  TString coutput_common = user_output_dir + erg_dir + "Templates/" + strdate + "/";
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
  cout << "Decay Channel: " << strChannel << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << strChannel << endl;
  tout << "===============================" << endl;
  tout << endl;

  vector<pair<FRtype*, FRtype*>> hFR_ZX_SS = getZXFR_SS<FRtype>();
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
  const unsigned int nSampleTypes = 1;
  const TString strProcess[nSampleTypes][2]={
    { "AllData", "ZX" }
  };

  string cinput = "/AllData";
  cinput = user_input_dir + cinput + "/" + user_treefile;
  InputTreeHandle* theTree = new InputTreeHandle(TString(cinput.c_str()), user_CRtreename, user_CRcountersname);
  foutput->cd();

  const unsigned int nTemplates = 1;
  TString strTemplates[nTemplates]={ "ZX" };
  TTree* newtree[nTemplates];
  TH2F* D_temp_2D[nTemplates];
  for (unsigned int t=0; t<nTemplates; t++){
    newtree[t] = new TTree(Form("T_2D_%s_Tree", strTemplates[t].Data()), "");
    newtree[t]->Branch("ZZMass", &ZZMass);
    newtree[t]->Branch("KD", &varKD);
    newtree[t]->Branch("weight", &weight);

    D_temp_2D[t] = new TH2F(Form("T_2D_%s", strTemplates[t].Data()), "", nbinsx, kDXarray, nbinsy, kDYarray);
    D_temp_2D[t]->Sumw2();
    D_temp_2D[t]->SetOption("colz");
    D_temp_2D[t]->GetYaxis()->SetTitle("D^{kin}_{bkg}");
    D_temp_2D[t]->GetXaxis()->SetTitle("m_{4l} (GeV)");
  }

  cout << "Process ZX under scrutiny" << endl;
  tout << "Process ZX under scrutiny" << endl;

  int nEntries = theTree->GetEntries();

  cout << "\tSample[" << 0 << "]=ZX under scrutiny. #Events = " << nEntries << endl;
  tout << "\tSample[" << 0 << "]=ZX under scrutiny. #Events = " << nEntries << endl;

  for (int ev=0; ev<nEntries; ev++){
    theTree->InitDefaults();
    theTree->GetEntry(ev);
    weight=1;

    if (ev%10000==0) cout << "Event " << ev << "..." << endl;

    // Best CR selection
    if (!test_bit(theTree->CRflag, (unsigned int)CJLST_FinalState::CRZLLss)) continue;

    if (theTree->ZZMass<ZZMass_cut[0] || theTree->ZZMass>=ZZMass_cut[1]) continue;

    int ZZFlav=theTree->Z1Flav*theTree->Z2Flav;
    if (ichan==(int)k4mu && abs(ZZFlav)!=pow(13, 4)) continue;
    else if (ichan==(int)k4e && abs(ZZFlav)!=pow(11, 4)) continue;
    else if (ichan==(int)k2e2mu && abs(ZZFlav)!=pow(13*11, 2)) continue;

    ZZMass = theTree->ZZMass;
    varKD = theTree->p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(theTree->p_GG_SIG_ghg2_1_ghz1_1_JHUGen + theTree->p_QQB_BKG_MCFM*getDbkgkinConstant(ZZFlav, ZZMass));
    if (varKD!=varKD) continue;

    // Apply FRs in SS
    for (unsigned int ilep=2; ilep<=3; ilep++){
      unsigned int whichLep = 1;
      if (abs(theTree->LepLepId->at(ilep))==11) whichLep=0;
      else if (abs(theTree->LepLepId->at(ilep))==13) whichLep=1;
      if (
        (fabs(theTree->LepEta->at(ilep))<1.2 && whichLep==0)
        ||
        (fabs(theTree->LepEta->at(ilep))<1.449 && whichLep==1)
        ) weight=evaluateTObject<FRtype>(hFR_ZX_SS.at(whichLep).first, (float)theTree->LepPt->at(ilep));
      else weight=evaluateTObject<FRtype>(hFR_ZX_SS.at(whichLep).second, (float)theTree->LepPt->at(ilep));
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
  for (unsigned int t=0; t<nTemplates; t++) conditionalizeHistogram(D_temp_2D[t], 0);
  cout << "Templates conditionally normalized" << endl;
  tout << "Templates conditionally normalized" << endl;

  for (unsigned int ifr=0; ifr<hFR_ZX_SS.size(); ifr++){
    cout << "Deleting hFR_ZX_SS[" << ifr << "]" << endl;
    gROOT->cd();
    FRtype* hfirst = hFR_ZX_SS.at(ifr).first;
    if (hfirst!=0) delete hfirst;
    cout << "Deleted hFR_ZX_SS[" << ifr << "] first" << endl;
    FRtype* hsecond = hFR_ZX_SS.at(ifr).second;
    if (hsecond!=0) delete hsecond;
    cout << "Deleted hFR_ZX_SS[" << ifr << "] second" << endl;
  }

  foutput->cd();
  for (unsigned int t=0; t<nTemplates; t++){
    foutput->WriteTObject(newtree[t]);
    foutput->WriteTObject(D_temp_2D[t]);
    delete D_temp_2D[t];
    delete newtree[t];
  }

  delete theTree;

  foutput->Close();
  tout.close();
}

void collectConditionalTemplates_Moriond2017mainstream_one(int ichan, int erg_tev, int Systematics, bool doSmooth){
  if (ichan>(int)nChannels) return;
  const TString strChannel = channame(ichan);
  constructSampleTypeList();

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

  const float xwidth_rebin = 1;
  const int nbinsx_rebin = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth_rebin;
  float kDXarray_rebin[nbinsx_rebin+1];
  for (int bin=0; bin<nbinsx_rebin+1; bin++) kDXarray_rebin[bin] = ZZMass_cut[0] + xwidth_rebin*bin;

  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString coutput_common = user_output_dir + erg_dir + "Templates/" + strdate + "/";
  TString cinput_common = coutput_common;
  gSystem->Exec("mkdir -p " + coutput_common);

  TString strDjet="";
  TString strSyst;
  if (Systematics == 0) strSyst = "Nominal";
  else if (Systematics == 1) strSyst = "SysUp_QCD";
  else if (Systematics == -1) strSyst = "SysDown_QCD";
  else if (Systematics == 2) strSyst = "SysUp_PDF";
  else if (Systematics == -2) strSyst = "SysDown_PDF";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }
  if (doSmooth){
    TString strCmd_dir = strdate;
    TString strCmd_ergtev = Form("%i", erg_tev);
    TString strCmd_channel = strChannel;
    TString strCmd_xbinning = Form("%i,%.0f,%.0f", nbinsx, kDXarray[0], kDXarray[nbinsx]);
    TString strCmd_syst = strSyst;
    TString strCmd_djet = strDjet;
    TString strCmd = Form("%s %s %s %s %s %s",
      strCmd_dir.Data(),
      strCmd_ergtev.Data(),
      strCmd_channel.Data(),
      strCmd_xbinning.Data(),
      strCmd_syst.Data(),
      strCmd_djet.Data()
      );
    cout << "strCmd = " << strCmd << endl;
    gSystem->Exec("pushd ../TemplateBuilderTemplates/; source ../TemplateBuilderTemplates/buildJson.sh " + strCmd + "; popd; ");
  }

  TString INPUT_NAME = Form("HtoZZ%s_ConditionalSmoothTemplatesForCombine_", strChannel.Data());
  TString OUTPUT_NAME = Form("HtoZZ%s_ConditionalSmoothMergedTemplatesForCombine_", strChannel.Data());
  INPUT_NAME += strSyst;
  OUTPUT_NAME += strSyst;

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  INPUT_NAME += ".root";
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << strChannel << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << strChannel << endl;
  tout << "===============================" << endl;
  tout << endl;

  const unsigned int nTemplates = 4;
  TString strTemplates[nTemplates]={
    "Sig",
    "ggBkg",
    "qqBkg",
    "ZX"
  };
  vector<TH2F*> D_temp_2D;

  TString cinput = cinput_common + INPUT_NAME;
  cout << "Opening file " << cinput << endl;
  tout << "Opening file " << cinput << endl;
  TFile* finput = TFile::Open(cinput, "read");

  foutput->cd();
  for (unsigned int t=0; t<nTemplates; t++){
    D_temp_2D.push_back((TH2F*)finput->Get(Form("T_2D_%s", strTemplates[t].Data())));
    D_temp_2D.at(t)->SetOption("colz");
    D_temp_2D.at(t)->GetYaxis()->SetTitle("D^{kin}_{bkg}");
    D_temp_2D.at(t)->GetXaxis()->SetTitle("m_{4l} (GeV)");

    // Floor
    for (int binx=0; binx<=D_temp_2D.at(t)->GetNbinsX()+1; binx++){
      for (int biny=0; biny<=D_temp_2D.at(t)->GetNbinsY()+1; biny++){
        double bincontent = D_temp_2D.at(t)->GetBinContent(binx, biny);
        if (bincontent<=0.) D_temp_2D.at(t)->SetBinContent(binx, biny, 1e-15);
      }
    }

    // Assign same bin content for qqZZ for m>250 in ZX and >600 in qqZZ.
    /*
    if (t>=nTemplates-2){
      int bin_int_start = D_temp_2D.at(t)->GetXaxis()->FindBin(250.);
      if (t==nTemplates-2) bin_int_start = D_temp_2D.at(t)->GetXaxis()->FindBin(600.);
      for (int biny=0; biny<=D_temp_2D.at(t)->GetNbinsY()+1; biny++){
        double integral = D_temp_2D.at(t)->Integral(bin_int_start, D_temp_2D.at(t)->GetNbinsX()+1, biny, biny);
        for (int binx=bin_int_start; binx<=D_temp_2D.at(t)->GetNbinsX()+1; binx++){
          D_temp_2D.at(t)->SetBinContent(binx, biny, integral);
          D_temp_2D.at(t)->SetBinError(binx, biny, integral);
        }
      }
    }
    */
    // Conditionalize
    wipeOverUnderFlows(D_temp_2D.at(t));
    conditionalizeHistogram(D_temp_2D.at(t), 0);
    //if (t!=nTemplates-1){
      cout << "Regularizing template " << t << endl;
      tout << "Regularizing template " << t << endl;

      const int nbinsy = D_temp_2D.at(t)->GetNbinsY();
      foutput->cd();
      TDirectory* xcheckdir = foutput->mkdir(Form("xcheck_%s", D_temp_2D.at(t)->GetName()));
      TCanvas* cSlice[nbinsy];
      TH1F* hSlice[nbinsy];
      TH1F* hSlice_new[nbinsy];

      // Plot pre-smoothing slices
      xcheckdir->cd();
      for (int ish=0; ish<nbinsy; ish++){
        hSlice[ish] = new TH1F(Form("h%s_Slice%i", D_temp_2D.at(t)->GetName(), ish), "", nbinsx, kDXarray[0], kDXarray[nbinsx]);
        for (int binx=1; binx<=D_temp_2D.at(t)->GetNbinsX(); binx++) hSlice[ish]->SetBinContent(binx, D_temp_2D.at(t)->GetBinContent(binx, ish));
        hSlice[ish]->SetLineColor(kRed);
        hSlice[ish]->SetLineWidth(2);
        hSlice[ish]->SetMarkerColor(kRed);

        cSlice[ish] = new TCanvas(Form("c_%s_Slice%i", D_temp_2D.at(t)->GetName(), ish), "", 8, 30, 800, 800);
        cSlice[ish]->cd();
        gStyle->SetOptStat(0);
        cSlice[ish]->SetFillColor(0);
        cSlice[ish]->SetBorderMode(0);
        cSlice[ish]->SetBorderSize(2);
        cSlice[ish]->SetTickx(1);
        cSlice[ish]->SetTicky(1);
        cSlice[ish]->SetLeftMargin(0.17);
        cSlice[ish]->SetRightMargin(0.05);
        cSlice[ish]->SetTopMargin(0.07);
        cSlice[ish]->SetBottomMargin(0.13);
        cSlice[ish]->SetFrameFillStyle(0);
        cSlice[ish]->SetFrameBorderMode(0);
        cSlice[ish]->SetFrameFillStyle(0);
        cSlice[ish]->SetFrameBorderMode(0);

        cSlice[ish]->cd();
        hSlice[ish]->Draw("hist");
      }

      foutput->cd();
      regularizeHistogram(D_temp_2D.at(t), 50, 0.01);
      conditionalizeHistogram(D_temp_2D.at(t), 0);

      // Plot post-smoothing slices
      xcheckdir->cd();
      for (int ish=0; ish<nbinsy; ish++){
        hSlice_new[ish] = new TH1F(Form("h%s_Slice%i_smoothened", D_temp_2D.at(t)->GetName(), ish), "", nbinsx, kDXarray[0], kDXarray[nbinsx]);
        for (int binx=1; binx<=D_temp_2D.at(t)->GetNbinsX(); binx++) hSlice_new[ish]->SetBinContent(binx, D_temp_2D.at(t)->GetBinContent(binx, ish));
        hSlice_new[ish]->SetLineColor(kBlue);
        hSlice_new[ish]->SetLineWidth(2);
        hSlice_new[ish]->SetMarkerColor(kBlue);

        cSlice[ish]->cd();
        hSlice_new[ish]->Draw("histsame");
        cSlice[ish]->RedrawAxis();
        cSlice[ish]->Modified();
        cSlice[ish]->Update();
        xcheckdir->WriteTObject(cSlice[ish]);
      }

      for (int ish=0; ish<nbinsy; ish++){
        cSlice[ish]->Close();
        delete hSlice_new[ish];
        delete hSlice[ish];
      }
      xcheckdir->Close();
      foutput->cd();

    //}
  }
  for (unsigned int izx=1; izx<=2; izx++){
    D_temp_2D.push_back(
      (TH2F*)D_temp_2D.at(nTemplates-1)->Clone(Form("%s%s", D_temp_2D.at(nTemplates-1)->GetName(), (izx==1 ? "_StatUp" : "_StatDn")))
      );
    // Mirroring from qqBkg
    for (int binx=0; binx<=D_temp_2D.at(nTemplates-1+izx)->GetNbinsX()+1; binx++){
      for (int biny=0; biny<=D_temp_2D.at(nTemplates-1+izx)->GetNbinsY()+1; biny++){
        double bincontent = D_temp_2D.at(nTemplates-1)->GetBinContent(binx, biny);
        double bincontent_qqzz = D_temp_2D.at(nTemplates-2)->GetBinContent(binx, biny);
        if (izx==1) bincontent = bincontent_qqzz;
        else bincontent = 2.*bincontent - bincontent_qqzz;
        if (bincontent<=0.) D_temp_2D.at(nTemplates-1+izx)->SetBinContent(binx, biny, 1e-15);
        else D_temp_2D.at(nTemplates-1+izx)->SetBinContent(binx, biny, bincontent);
      }
    }

    // Re-conditionalize
    wipeOverUnderFlows(D_temp_2D.at(nTemplates-1+izx));
    conditionalizeHistogram(D_temp_2D.at(nTemplates-1+izx), 0);
  }

  vector<TH2F*> D_temp_2D_rebin;
  int nloops = (int)(xwidth/xwidth_rebin+0.5);
  for (unsigned int t=0; t<D_temp_2D.size(); t++){
    TString tmpnamecore=D_temp_2D.at(t)->GetName();
    TString tmpname = tmpnamecore + "_tmp";
    TH2F* htmp = new TH2F(tmpname, D_temp_2D.at(t)->GetTitle(), nbinsx_rebin, kDXarray_rebin, nbinsy, kDYarray);
    for (int binx=1; binx<=D_temp_2D.at(t)->GetNbinsX(); binx++){
      for (int biny=1; biny<=D_temp_2D.at(t)->GetNbinsY(); biny++){
        double bincontent = D_temp_2D.at(t)->GetBinContent(binx, biny);
        for (int ist=0; ist<nloops; ist++){
          int newbinx = nloops*(binx-1)+ist+1;
          D_temp_2D.at(t)->SetBinContent(newbinx, biny, bincontent);
        }
      }
    }
    delete D_temp_2D.at(t);
    htmp->SetName(tmpnamecore);
    D_temp_2D_rebin.push_back(htmp);
  }
  D_temp_2D.clear();
  for (unsigned int t=0; t<D_temp_2D_rebin.size(); t++) D_temp_2D.push_back(D_temp_2D_rebin.at(t));
  D_temp_2D_rebin.clear();

  for (unsigned int t=0; t<D_temp_2D.size(); t++){
    foutput->WriteTObject(D_temp_2D.at(t));
    delete D_temp_2D.at(t);
  }
  finput->Close();
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


SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPt[4], float GenLepEta[4], float GenLepPhi[4]){
  SimpleParticleCollection_t result;
  for (int idau=0; idau<4; idau++){
    float mass=0;
    if (abs(GenLepId[idau])==11) mass=0.000511;
    else if (abs(GenLepId[idau])==13) mass=0.10566;
    else if (abs(GenLepId[idau])==15) mass=1.7768;
    else if (abs(GenLepId[idau])==5) mass=4.75;
    else if (abs(GenLepId[idau])==4) mass=1.275;
    TLorentzVector mom; mom.SetPtEtaPhiM((double)GenLepPt[idau], (double)GenLepEta[idau], (double)GenLepPhi[idau], (double)mass);
    result.push_back(SimpleParticle_t((int)GenLepId[idau], mom));
  }
  return result;
}

vector<pair<string, double>> getFileList(string strinput){
  cout << strinput << endl;
  ifstream fin;
  string filename = user_input_dir.Data();
  filename += "/"; filename += strinput;
  fin.open(filename);
  cout << "Checking list from " << filename << endl;
  vector<pair<string, double>> result;
  if (fin.good()){
    if (filename.find("Sig")!=string::npos && filename.find("POWHEG")!=string::npos){ // Signal POWHEG
      while (!fin.eof()){
        string strtmp;
        double mass=-1;
        fin >> strtmp;

        if (strtmp=="") continue;

        vector<string> strsplit;
        splitOptionRecursive(strtmp, strsplit, 'H');
        if (strsplit.size()>1){
          string strmass = strsplit.at(1);
          strsplit.clear();
          splitOptionRecursive(strmass, strsplit, '_');
          strmass = strsplit.at(0);
          mass = std::stod(strmass);
        }
        else{
          cout << "getFileList: POWHEG filename " << strtmp << " does not conform requirements" << endl;
          cout << '\t';
          for (unsigned int isp=0; isp<strsplit.size(); isp++) cout << strsplit.at(isp) << " ";
          cout << endl;
        }
        result.push_back(pair<string, double>(strtmp, mass));
        cout << "POWHEG file " << strtmp << " with mass " << mass << endl;
      }
    }
    else{
      double mass=-1;
      if (filename.find("Sig")!=string::npos || filename.find("BSI")!=string::npos) mass = 125; // Sig or BSI
      while (!fin.eof()){
        string strtmp;
        fin >> strtmp;
        if (strtmp!="") result.push_back(pair<string, double>(strtmp, mass));
        else break;
        cout << "Specialized file " << strtmp << " with mass " << mass << endl;
      }
    }
  }
  fin.close();
  return result;
}

int getCategory(
  float ZZMass,

  vector<float> JetPt,
  vector<float> JetSigma,
  vector<float> JetPhi,
  vector<float> JetQGLikelihood,
  vector<float> JetIsBtaggedWithSF,

  int JESvariation,
  int nExtraLep,
  int nExtraZ,
  float PFMET,

  float p_JJQCD_SIG_ghg2_1_JHUGen,
  float p_JQCD_SIG_ghg2_1_JHUGen,
  float p_JJVBF_SIG_ghv1_1_JHUGen,
  float p_JVBF_SIG_ghv1_1_JHUGen,
  float pAux_JVBF_SIG_ghv1_1_JHUGen,
  float p_HadWH_SIG_ghw1_1_JHUGen,
  float p_HadZH_SIG_ghz1_1_JHUGen,

  int categorizationType
  ){
  bool useQGTagging = categorizationType<0;
  bool useVHMETTagged=true;

  int nCleanedJetsPt30=0;
  int nCleanedJetsPt30BTagged_bTagSF=0;
  float jetQGLikelihood[2]={ -1, -1 };
  float jetPhi[2]={ 0, 0 };
  float jetPt[2]={ 0, 0 };
  for (unsigned int j=0; j<JetPt.size(); j++){
    float jetpt = JetPt.at(j) * (1.+JetSigma.at(j)*((float)JESvariation));
    if (jetpt>=30.){
      nCleanedJetsPt30++;
      if (JetIsBtaggedWithSF.at(j)>0.) nCleanedJetsPt30BTagged_bTagSF++;
      if (jetpt>=jetPt[0]){
        jetPt[0]=jetpt;
        jetPhi[0]=JetPhi.at(j);
        jetQGLikelihood[0]=JetQGLikelihood.at(j);
      }
      else if (jetpt<jetPt[0] && jetpt>=jetPt[1]){
        jetPt[1]=jetpt;
        jetPhi[1]=JetPhi.at(j);
        jetQGLikelihood[1]=JetQGLikelihood.at(j);
      }
    }
  }

  if (categorizationType==0) return 0;
  else if (abs(categorizationType)==1) // Main analysis categories
    return categoryMor17(
    nExtraLep,
    nExtraZ,
    nCleanedJetsPt30,
    nCleanedJetsPt30BTagged_bTagSF,
    jetQGLikelihood,
    p_JJQCD_SIG_ghg2_1_JHUGen,
    p_JQCD_SIG_ghg2_1_JHUGen,
    p_JJVBF_SIG_ghv1_1_JHUGen,
    p_JVBF_SIG_ghv1_1_JHUGen,
    pAux_JVBF_SIG_ghv1_1_JHUGen,
    p_HadWH_SIG_ghw1_1_JHUGen,
    p_HadZH_SIG_ghz1_1_JHUGen,
    jetPhi,
    ZZMass,
    PFMET,
    useVHMETTagged,
    useQGTagging
    );
  else
    return categoryMor17_LegacyStyle(
    nExtraLep,
    nCleanedJetsPt30,
    nCleanedJetsPt30BTagged_bTagSF,
    jetQGLikelihood,
    p_JJQCD_SIG_ghg2_1_JHUGen,
    p_JJVBF_SIG_ghv1_1_JHUGen,
    jetPhi,
    ZZMass,
    useQGTagging
    );
}

void splitOption(const string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
Bool_t checkListVariable(const vector<string>& list, const string& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
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
template <typename T> vector<pair<T*, T*>> getZXFR_SS(){
  /*
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
  */
  TString hname[2][2]={
    { "FR_SS_electron_EB", "FR_SS_electron_EE" },
    { "FR_SS_muon_EB", "FR_SS_muon_EE" }
  };
  TFile* finput;
  finput = TFile::Open("../data/FakeRate_SS_Moriond368.root", "read");
  vector<pair<T*, T*>> result;
  for (int f=0; f<2; f++){
    if (!(finput!=0 && finput->IsOpen())){ cerr << "getZXFR_SS: File is not open!" << endl; return result; }
    else cout << "getZXFR_SS: File opened" << endl;
    T* htmp[2];
    gROOT->cd();
    for (int t=0; t<2; t++) htmp[t] = (T*)finput->Get(hname[f][t]);
    result.push_back(pair<T*, T*>((T*)htmp[0]->Clone(Form("%s_copy", hname[f][0].Data())), (T*)htmp[1]->Clone(Form("%s_copy", hname[f][1].Data()))));
  }
  if (finput!=0 && finput->IsOpen()) finput->Close();
  return result;
}

template<> double evaluateTObject<TH1F>(TH1F* obj, float val){
  int bin = obj->GetXaxis()->FindBin(val);
  if (bin>obj->GetXaxis()->GetNbins()) bin=obj->GetXaxis()->GetNbins();
  return obj->GetBinContent(bin);
}
template<> double evaluateTObject<TGraphErrors>(TGraphErrors* obj, float val){
  double* xx = obj->GetX();
  double* yy = obj->GetY();
  int n = obj->GetN();
  if (val>xx[n-1]) val=xx[n-1];
  return obj->Eval(val);
}


bool test_bit(int mask, unsigned int iBit){ return (mask >> iBit) & 1; }

/*
void checkDbkgkin(int flavor){
  if (flavor>2) return;
  TChain* tggH = new TChain("T_2D_Sig_Tree");
  TChain* tqqBkg = new TChain("T_2D_qqBkg_Tree");
  if (flavor<=0){
    tggH->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ4mu_ConditionalTemplatesForCombine_Nominal.root");
    tqqBkg->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ4mu_ConditionalTemplatesForCombine_Nominal.root");
  }
  if (flavor<0 || flavor==1){
    tggH->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ4e_ConditionalTemplatesForCombine_Nominal.root");
    tqqBkg->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ4e_ConditionalTemplatesForCombine_Nominal.root");
  }
  if (flavor<0 || flavor==2){
    tggH->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ2e2mu_ConditionalTemplatesForCombine_Nominal.root");
    tqqBkg->Add("/scratch0/hep/usarical/CJLST/Analysis/Moriond2017_mainstream/LHC_13TeV/Templates/180716/HtoZZ2e2mu_ConditionalTemplatesForCombine_Nominal.root");
  }
  int ZZflav=143;
  if (flavor==0) ZZflav=169;
  else if (flavor==1) ZZflav=121;
  else if (flavor==2) ZZflav=143;

  TString coutput = "test";
  if (flavor>=0) coutput.Append(Form("_%s", user_folder[flavor].Data()));
  coutput.Append(".root");

  float weight, reweight;
  float ZZMass, KD;
  tggH->SetBranchAddress("ZZMass", &ZZMass);
  tggH->SetBranchAddress("KD", &KD);
  tggH->SetBranchAddress("weight", &weight);
  tggH->SetBranchAddress("reweight", &reweight);
  tqqBkg->SetBranchAddress("ZZMass", &ZZMass);
  tqqBkg->SetBranchAddress("KD", &KD);
  tqqBkg->SetBranchAddress("weight", &weight);
  tqqBkg->SetBranchAddress("reweight", &reweight);

  // 70-3000 in 5 GeV steps
  const int nbins = (3000.-70)/5;
  vector<float> KDarray_ggH[nbins];
  vector<float> KDarray_qqBkg[nbins];
  for (int bin=0; bin<nbins; bin++){
    KDarray_ggH[bin].clear();
    KDarray_qqBkg[bin].clear();
  }
  cout << "#ggH = " << tggH->GetEntries() << endl;
  for (int ev=0; ev<tggH->GetEntries(); ev++){
    tggH->GetEntry(ev);
    if (ev%10000==0) cout << "Event " << ev << endl;
    int bin = (ZZMass-70.)/5.;
    if (bin>=nbins) continue;
    if (KD!=0.) KD = 1./((1./KD-1.)/getDbkgkinConstant(ZZflav, ZZMass)+1.);
    float KDinv = 1.-KD;
    addByLowest(KDinv, KDarray_ggH[bin]);
  }
  cout << "Done ggH" << endl;
  cout << "#qqBkg = " << tqqBkg->GetEntries() << endl;
  for (int ev=0; ev<tqqBkg->GetEntries(); ev++){
    tqqBkg->GetEntry(ev);
    if (ev%10000==0) cout << "Event " << ev << endl;
    int bin = (ZZMass-70.)/5.;
    if (bin>=nbins) continue;
    if (KD!=0.) KD = 1./((1./KD-1.)/getDbkgkinConstant(ZZflav, ZZMass)+1.);
    addByLowest(KD, KDarray_qqBkg[bin]);
  }
  cout << "Done qqBkg" << endl;

  delete tggH;
  delete tqqBkg;

  TFile* foutput = new TFile(coutput, "recreate");
  TH1F* histo = new TH1F("test", "", nbins, 70, 3000);

  for (int bin=0; bin<nbins; bin++){
    cout << "Checking bin " << bin << endl;

    const float KDinc=0.001;
    float closestApproach=0;
    float diffFrac=2;

    cout << "ggH size = " << KDarray_ggH[bin].size() << endl;
    cout << "qqBkg size = " << KDarray_qqBkg[bin].size() << endl;

    for (float iKD=KDinc; iKD<1.; iKD+=KDinc){
      float gKD = 1.-iKD;
      int index_ggH=0;
      while (index_ggH<(int)KDarray_ggH[bin].size() && KDarray_ggH[bin].at(index_ggH)<gKD) index_ggH++;

      int index_qqBkg=0;
      while (index_qqBkg<(int)KDarray_qqBkg[bin].size() && KDarray_qqBkg[bin].at(index_qqBkg)<iKD) index_qqBkg++;

      double f_ggH = ((double)index_ggH)/((double)KDarray_ggH[bin].size());
      double f_qqBkg = ((double)index_qqBkg)/((double)KDarray_qqBkg[bin].size());
      if (fabs(f_ggH-f_qqBkg)<diffFrac){
        diffFrac = fabs(f_ggH-f_qqBkg);
        closestApproach = iKD;
      }
    }

    histo->SetBinContent(bin+1, closestApproach);
  }

  foutput->WriteTObject(histo);
  foutput->Close();
}

void checkDVBF2jets(int erg_tev){
  const float ZZMass_cut[2]={ 70., 5000. };
  const float xwidth = 10;
  const int nbins = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth;

  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString INPUT_NAME = user_treefile;
  TString OUTPUT_NAME = "testVBF2jets.root";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = "./";

  float GenHMass, ZZMass, pvbf_VAJHU_highestPTJets, phjj_VAJHU_highestPTJets;

  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  TH1F* histo = new TH1F("test", "", nbins, ZZMass_cut[0], ZZMass_cut[1]);

  // Get list of samples
  const unsigned int nSampleTypes = 2;
  const TString strProcess[nSampleTypes][2]={
    { "Samples_gg_Sig_POWHEG.txt", "gg_Sig" },
    { "Samples_VBF_Sig_POWHEG.txt", "VBF_Sig" }
  };
  vector<pair<string, double>> fileNameList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++) fileNameList[is] = getFileList(string(strProcess[is][0].Data()));
  vector<TFile*> fileList[nSampleTypes];
  vector<pair<TTree*, TH1F*>> treeList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_gg2VV_location + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      TFile* finput = TFile::Open(cinput.c_str(), "read"); finput->cd();
      TTree* theTree = (TTree*)finput->Get(user_treename);
      TH1F* theCounters = (TH1F*)finput->Get(user_countersname);

      fileList[is].push_back(finput);
      foutput->cd();

      theTree->SetBranchAddress("GenHMass", &GenHMass);
      theTree->SetBranchAddress("ZZMass", &ZZMass);
      theTree->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets);
      theTree->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets);

      treeList[is].push_back(pair<TTree*, TH1F*>(theTree, theCounters));
    }
  }

  foutput->cd();
  vector<float> KDarray_HJJ[nbins];
  vector<float> KDarray_VBF[nbins];
  for (int bin=0; bin<nbins; bin++){
    KDarray_HJJ[bin].clear();
    KDarray_VBF[bin].clear();
  }
  for (unsigned int is=0; is<nSampleTypes; is++){
    cout << "Process " << is << " under scrutiny" << endl;
    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      TTree* theTree = treeList[is].at(it).first;
      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;

      for (int ev=0; ev<nEntries; ev++){
        theTree->GetEntry(ev);
        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Fix for JHUGen v6 vegas bug in low mass samples
        if (strProcess[is][0].Contains("POWHEG") && mass>0 && mass<300 && GenHMass>mass+5.){
          //cout << "Process " << is << " sample " << it << " has mH = " << mass << " and mH* = " << GenHMass << ">mH+5 GeV. Ignoring event..." << endl;
          continue;
        }

        if (pvbf_VAJHU_highestPTJets<0. || phjj_VAJHU_highestPTJets<0.) continue;

        float varKD = pvbf_VAJHU_highestPTJets/(pvbf_VAJHU_highestPTJets + phjj_VAJHU_highestPTJets);
        if (varKD!=varKD) continue;

        int bin = (ZZMass-ZZMass_cut[0])/xwidth;
        if (bin>=nbins) continue;

        if (strProcess[is][1]=="VBF_Sig"){
          float KDinv = 1.-varKD;
          addByLowest(KDinv, KDarray_VBF[bin]);
        }
        else addByLowest(varKD, KDarray_HJJ[bin]);
      }
    }
  }

  for (int bin=0; bin<nbins; bin++){
    cout << "Checking bin " << bin << endl;

    const float KDinc=0.001;
    float closestApproach=0;
    float diffFrac=2;

    cout << "HJJ size = " << KDarray_HJJ[bin].size() << endl;
    cout << "VBF size = " << KDarray_VBF[bin].size() << endl;

    for (float iKD=KDinc; iKD<1.; iKD+=KDinc){
      float gKD = 1.-iKD;
      int index_VBF=0;
      while (index_VBF<(int)KDarray_VBF[bin].size() && KDarray_VBF[bin].at(index_VBF)<gKD) index_VBF++;

      int index_HJJ=0;
      while (index_HJJ<(int)KDarray_HJJ[bin].size() && KDarray_HJJ[bin].at(index_HJJ)<iKD) index_HJJ++;

      double f_HJJ = ((double)index_HJJ)/((double)KDarray_HJJ[bin].size());
      double f_VBF = ((double)index_VBF)/((double)KDarray_VBF[bin].size());
      if (fabs(f_HJJ-f_VBF)<diffFrac){
        diffFrac = fabs(f_HJJ-f_VBF);
        closestApproach = iKD;
      }
    }

    histo->SetBinContent(bin+1, closestApproach);
  }

  foutput->WriteTObject(histo);

  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileList[is].size(); ifile++){
      if (fileList[is].at(ifile)!=0 && fileList[is].at(ifile)->IsOpen()) fileList[is].at(ifile)->Close();
    }
  }

  foutput->Close();
}

void checkDVBF1jet(int erg_tev){
  const float ZZMass_cut[2]={ 70., 5000. };
  const float xwidth = 10;
  const int nbins = (ZZMass_cut[1]-ZZMass_cut[0])/xwidth;

  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString INPUT_NAME = user_treefile;
  TString OUTPUT_NAME = "testVBF1jet.root";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = "./";

  float GenHMass, ZZMass, pvbf_VAJHU_highestPTJets, phj_VAJHU, pAux_vbf_VAJHU;

  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  TH1F* histo = new TH1F("test", "", nbins, ZZMass_cut[0], ZZMass_cut[1]);

  // Get list of samples
  const unsigned int nSampleTypes = 2;
  const TString strProcess[nSampleTypes][2]={
    { "Samples_gg_Sig_POWHEG.txt", "gg_Sig" },
    { "Samples_VBF_Sig_POWHEG.txt", "VBF_Sig" }
  };
  vector<pair<string, double>> fileNameList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++) fileNameList[is] = getFileList(string(strProcess[is][0].Data()));
  vector<TFile*> fileList[nSampleTypes];
  vector<pair<TTree*, TH1F*>> treeList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_gg2VV_location + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      TFile* finput = TFile::Open(cinput.c_str(), "read"); finput->cd();
      TTree* theTree = (TTree*)finput->Get(user_treename);
      TH1F* theCounters = (TH1F*)finput->Get(user_countersname);

      fileList[is].push_back(finput);
      foutput->cd();

      theTree->SetBranchAddress("GenHMass", &GenHMass);
      theTree->SetBranchAddress("ZZMass", &ZZMass);
      theTree->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets);
      theTree->SetBranchAddress("pAux_vbf_VAJHU", &pAux_vbf_VAJHU);
      theTree->SetBranchAddress("phj_VAJHU", &phj_VAJHU);

      treeList[is].push_back(pair<TTree*, TH1F*>(theTree, theCounters));
    }
  }

  foutput->cd();
  vector<float> KDarray_HJJ[nbins];
  vector<float> KDarray_VBF[nbins];
  for (int bin=0; bin<nbins; bin++){
    KDarray_HJJ[bin].clear();
    KDarray_VBF[bin].clear();
  }
  for (unsigned int is=0; is<nSampleTypes; is++){
    cout << "Process " << is << " under scrutiny" << endl;
    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      TTree* theTree = treeList[is].at(it).first;
      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;

      for (int ev=0; ev<nEntries; ev++){
        theTree->GetEntry(ev);
        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Fix for JHUGen v6 vegas bug in low mass samples
        if (strProcess[is][0].Contains("POWHEG") && mass>0 && mass<300 && GenHMass>mass+5.){
          //cout << "Process " << is << " sample " << it << " has mH = " << mass << " and mH* = " << GenHMass << ">mH+5 GeV. Ignoring event..." << endl;
          continue;
        }
        if (pvbf_VAJHU_highestPTJets<0. || phj_VAJHU<0.) continue;

        float varKD = pvbf_VAJHU_highestPTJets*pAux_vbf_VAJHU/(pvbf_VAJHU_highestPTJets*pAux_vbf_VAJHU + phj_VAJHU);
        if (varKD!=varKD) continue;

        int bin = (ZZMass-ZZMass_cut[0])/xwidth;
        if (bin>=nbins) continue;

        if (strProcess[is][1]=="VBF_Sig"){
          float KDinv = 1.-varKD;
          addByLowest(KDinv, KDarray_VBF[bin]);
        }
        else addByLowest(varKD, KDarray_HJJ[bin]);
      }
    }
  }

  for (int bin=0; bin<nbins; bin++){
    cout << "Checking bin " << bin << endl;

    const float KDinc=0.001;
    float closestApproach=0;
    float diffFrac=2;

    cout << "HJJ size = " << KDarray_HJJ[bin].size() << endl;
    cout << "VBF size = " << KDarray_VBF[bin].size() << endl;

    for (float iKD=KDinc; iKD<1.; iKD+=KDinc){
      float gKD = 1.-iKD;
      int index_VBF=0;
      while (index_VBF<(int)KDarray_VBF[bin].size() && KDarray_VBF[bin].at(index_VBF)<gKD) index_VBF++;

      int index_HJJ=0;
      while (index_HJJ<(int)KDarray_HJJ[bin].size() && KDarray_HJJ[bin].at(index_HJJ)<iKD) index_HJJ++;

      double f_HJJ = ((double)index_HJJ)/((double)KDarray_HJJ[bin].size());
      double f_VBF = ((double)index_VBF)/((double)KDarray_VBF[bin].size());
      if (fabs(f_HJJ-f_VBF)<diffFrac){
        diffFrac = fabs(f_HJJ-f_VBF);
        closestApproach = iKD;
      }
    }

    histo->SetBinContent(bin+1, closestApproach);
  }

  foutput->WriteTObject(histo);

  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileList[is].size(); ifile++){
      if (fileList[is].at(ifile)!=0 && fileList[is].at(ifile)->IsOpen()) fileList[is].at(ifile)->Close();
    }
  }

  foutput->Close();
}
*/
void regularizeSlice(TGraph* tgSlice, std::vector<double>* fixedX, double omitbelow, int nIter_, double threshold_){
  unsigned int nbins_slice = tgSlice->GetN();
  double* xy_slice[2]={
    tgSlice->GetX(),
    tgSlice->GetY()
  };

  double* xy_mod[2];
  for (unsigned int ix=0; ix<2; ix++){
    xy_mod[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++) xy_mod[ix][iy] = xy_slice[ix][iy];
  }
  unsigned int bin_first = 2, bin_last = nbins_slice-1;

  std::vector<int> fixedBins;
  if (fixedX!=0){
    for (unsigned int ifx=0; ifx<fixedX->size(); ifx++){
      double requestedVal = fixedX->at(ifx);
      double distance=1e15;
      int bin_to_fix=-1;
      for (unsigned int bin=0; bin<nbins_slice; bin++){ if (distance>fabs(xy_mod[0][bin]-requestedVal)){ bin_to_fix = bin; distance = fabs(xy_mod[0][bin]-requestedVal); } }
      if (bin_to_fix>=0) fixedBins.push_back(bin_to_fix);
      cout << "Requested to fix bin " << bin_to_fix << endl;
    }
  }
  if (omitbelow>0.){
    for (unsigned int bin=0; bin<nbins_slice; bin++){
      if (xy_mod[0][bin]<omitbelow) fixedBins.push_back(bin);
      cout << "Requested to fix bin " << bin << endl;
    }
  }

  double* xx_second;
  double* yy_second;

  int nIter = (nIter_<0 ? 1000 : nIter_);
  for (int it=0; it<nIter; it++){
    double threshold = (threshold_<0. ? 0.01 : threshold_);
    for (unsigned int binIt = bin_first; binIt<=bin_last; binIt++){
      bool doFix=false;
      for (unsigned int ifx=0; ifx<fixedBins.size(); ifx++){
        if ((int)(binIt-1)==fixedBins.at(ifx)){ doFix=true; /*cout << "Iteration " << it << " is fixing bin " << (binIt-1) << endl; */break; }
      }
      if (doFix) continue;

      int ctr = 0;
      int nbins_second = nbins_slice-1;
      xx_second = new double[nbins_second];
      yy_second = new double[nbins_second];
      for (unsigned int bin = 1; bin<=nbins_slice; bin++){
        if (bin==binIt) continue;
        xx_second[ctr] = xy_mod[0][bin-1];
        yy_second[ctr] = xy_mod[1][bin-1];
        ctr++;
      }

      TGraph* interpolator = new TGraph(nbins_second, xx_second, yy_second);
      double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
      double derivative_last = (yy_second[nbins_second-1]-yy_second[nbins_second-2])/(xx_second[nbins_second-1]-xx_second[nbins_second-2]);
      TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double center = xy_mod[0][binIt-1];
      double val = spline->Eval(center);
      if (fabs(xy_mod[1][binIt-1]-val)>threshold*val && val>0.) xy_mod[1][binIt-1]=val;

      delete spline;
      delete interpolator;

      delete[] yy_second;
      delete[] xx_second;
    }
  }

  for (unsigned int iy=0; iy<nbins_slice; iy++) xy_slice[1][iy] = xy_mod[1][iy];
  for (unsigned int ix=0; ix<2; ix++) delete[] xy_mod[ix];
}

void regularizeHistogram(TH2F* histo, int nIter_, double threshold_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();

  for (int iy=1; iy<=nbinsy; iy++){
    cout << "regularizeHistogram::Bin " << iy << " being regularized..." << endl;
    double xy[2][nbinsx];
    for (int ix=1; ix<=nbinsx; ix++){
      xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
      xy[1][ix-1] = histo->GetBinContent(ix, iy);
    }

    TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
    tg->SetName("tg_tmp");
    regularizeSlice(tg, 0, 0, nIter_, threshold_);
    for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, tg->Eval(xy[0][ix-1]));
  }
  conditionalizeHistogram(histo, 0);
}

void regularizeHistogram(TH3F* histo, int nIter_, double threshold_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();
  const int nbinsz = histo->GetNbinsZ();

  for (int iy=1; iy<=nbinsy; iy++){
    for (int iz=1; iz<=nbinsz; iz++){
      cout << "regularizeHistogram::Bin " << iy << ", " << iz << " being regularized..." << endl;
      double xy[2][nbinsx];
      for (int ix=1; ix<=nbinsx; ix++){
        xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
        xy[1][ix-1] = histo->GetBinContent(ix, iy, iz);
      }

      TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
      tg->SetName("tg_tmp");
      regularizeSlice(tg, 0, 0, nIter_, threshold_);
      for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, iz, tg->Eval(xy[0][ix-1]));
    }
  }
  conditionalizeHistogram(histo, 0);
}

void conditionalizeHistogram(TH2F* histo, unsigned int axis){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
  else{
    for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
      double integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy);
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
}

void conditionalizeHistogram(TH3F* histo, unsigned int axis){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1, 0, histo->GetNbinsZ()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
        for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
  else if (axis==1){
    for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
      double integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy, 0, histo->GetNbinsZ()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
  else{
    for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
      double integral = histo->Integral(0, histo->GetNbinsX()+1, 0, histo->GetNbinsY()+1, iz, iz);
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
}

void wipeOverUnderFlows(TH1F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    if (hwipe->GetBinContent(binx)!=0){
      hwipe->SetBinContent(binx, 0);
      cout << hwipe->GetName() << " binX = " << binx << " non-zero." << endl;
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}
void wipeOverUnderFlows(TH2F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      if (hwipe->GetBinContent(binx, biny)!=0){
        hwipe->SetBinContent(binx, biny, 0);
        cout << hwipe->GetName() << " binX = " << binx << " binY = " << biny << " non-zero." << endl;
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}
void wipeOverUnderFlows(TH3F* hwipe){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1, 0, hwipe->GetNbinsZ()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      for (int binz=0; binz<=hwipe->GetNbinsZ()+1; binz++){
        if (binz>=1 && binz<=hwipe->GetNbinsZ()) continue;
        if (hwipe->GetBinContent(binx, biny, binz)!=0){
          hwipe->SetBinContent(binx, biny, binz, 0);
          cout << hwipe->GetName() << " binX = " << binx << " binY = " << biny << " binZ = " << binz << " non-zero." << endl;
        }
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  hwipe->Scale(wipeScale);
}

