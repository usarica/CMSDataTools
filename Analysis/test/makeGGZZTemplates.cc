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
#include "HelperFunctions.h"
#include "SampleHelpers.h"
#include "CJLSTSet.h"
#include "CategorizationHelpers.h"
#include "Mela.h"

using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace CategorizationHelpers;

// Constants to affect the template code
const TString fixedDate="";
const TString user_output_dir = "output/";

// Initializers
void makeGGZZTemplates_one(int channel, int icat, int Systematics);
void collectGGZZTemplates_one(int ichan, int icat, int Systematics, bool doSmooth=true);

// FIXME:
// 2) WRITE A REWEIGHTING CLASS THAT TAKES THE CJLSTTREE AND THE WEIGHT STRING TO CONTROL IF PROPAGATOR NEEDS TO BE TAKEN OUT
// 3) ADD A REWEIGHTING CLASSES FIELD TO GETEVENTS
// 4) IMPLEMENT A CATEGORIZATION PLAN TO GETEVENTS AS WELL
// Event scanner
void getEvents(
  CJLSTSet* theSet,

  const TString& strTrackVar,
  const vector<TString>& strExtraWeightsList,
  const vector<pair<Discriminant*, vector<TString>>>& KDVars,

  vector<SimpleEntry>& index,

  Channel& chan
  );

// Main Function, runs over all desired iterations
void makeGGZZTemplates(){
  for (int ichan=(int) k4mu; ichan<(int) NChannels; ichan++){
    for (int syst=0; syst<=0; syst++){
      //for (int icat=(int) Untagged; icat<(int) nCategoriesMor17; icat++){
      for (int icat=(int) Inclusive; icat<=(int) Inclusive; icat++){
        makeGGZZTemplates_one(ichan, icat, syst);
        collectGGZZTemplates_one(ichan, icat, syst, true);
      }
    }
  }
}

// Function to build one templates
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// theSqrts = 13 (CoM energy) is fixed in Samples.h
void makeGGZZTemplates_one(int ichan, int icat, int Systematics){
  if (ichan>=(int)NChannels) return;

  const Channel channel = (Channel) ichan;
  const TString strChannel = getChannelName(channel);

  const WidthCategory category = (WidthCategory) icat;
  const TString strCategory = getWidthCategoryName(category);

  enum{
    kBkg=0,
    kSig=1,
    kBSI=2,
    nTemplates=3
  };
  TString strTemplateName[nTemplates]={
    "Bkg",
    "Sig",
    "BSI"
  };
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

  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form("%s_HtoZZ%s_Stage1_", strCategory.Data(), strChannel.Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  float varKD=0;
  float weight=1;
  float ZZMass;

  TVar::VerbosityLevel verbosity = TVar::SILENT;

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << "Opened file " << coutput << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << theSqrts << " TeV" << endl;
  cout << "Decay Channel: " << strChannel << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << theSqrts << " TeV" << endl;
  tout << "Decay Channel: " << strChannel << endl;
  tout << "===============================" << endl;
  tout << endl;

  // Get list of samples
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("gg_Sig_POWHEG");

  vector<TString> strSamples;
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples);
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);

  vector<SimpleEntry> entryCollection;



  foutput->cd();

  TTree* theFinalTree[nTemplates];
  TH2F* D_temp_2D[nTemplates];
  for (unsigned int t=0; t<nTemplates; t++){
    theFinalTree[t] = new TTree(Form("T_2D_%s_Tree", strTemplateName[t].Data()), "");
    theFinalTree[t]->Branch("ZZMass", &ZZMass);
    theFinalTree[t]->Branch("KD", &varKD);
    theFinalTree[t]->Branch("weight", &weight);

    D_temp_2D[t] = new TH2F(Form("T_2D_%s", strTemplateName[t].Data()), "", nbinsx, kDXarray, nbinsy, kDYarray);
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
      theFinalTree[t]->Fill();
      D_temp_2D[t]->Fill(ZZMass, varKD, weight*reweight);
    }
    for (int ev=0; ev<rewtree[t]->GetEntries(); ev++){
      rewtree[t]->GetEntry(ev);
      int bin = D_temp_1D_raw[t]->GetXaxis()->FindBin(ZZMass);
      double rescale = 1.;
      if (D_temp_1D_rew[t]->GetBinContent(bin)>0.) rescale = N_temp_1D_rew[t]->GetBinContent(bin)/D_temp_1D_rew[t]->GetBinContent(bin);
      reweight *= rescale;
      theFinalTree[t]->Fill();
      D_temp_2D[t]->Fill(ZZMass, varKD, weight*reweight);
    }

    // Conditionalize the histograms
    conditionalizeHistogram(D_temp_2D[t], 0);
  }

  for (unsigned int t=0; t<nTemplates; t++){
    foutput->WriteTObject(theFinalTree[t]);
    foutput->WriteTObject(D_temp_2D[t]);
    delete N_temp_1D_rew[t];
    delete N_temp_1D_raw[t];
    delete D_temp_1D_rew[t];
    delete D_temp_1D_raw[t];
    delete D_temp_2D[t];
    delete theFinalTree[t];
    delete rewtree[t];
    delete rawtree[t];
  }

  for (unsigned int is=0; is<nSampleTypes; is++){ for (unsigned int it=0; it<treeList[is].size(); it++) delete treeList[is].at(it); }

  foutput->Close();
  tout.close();
}


void collectGGZZTemplates_one(int ichan,  int Systematics, bool doSmooth){
  if (ichan>(int)NChannels) return;
  const TString strChannel = getChannelName(ichan);
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

  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
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
    TString strCmd_ergtev = Form("%i", theSqrts);
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
  cout << "CoM Energy: " << theSqrts << " TeV" << endl;
  cout << "Decay Channel: " << strChannel << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "Opened file " << coutput << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << theSqrts << " TeV" << endl;
  tout << "Decay Channel: " << strChannel << endl;
  tout << "===============================" << endl;
  tout << endl;

  const unsigned int nTemplates = 4;
  TString strTemplateName[nTemplates]={
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
    D_temp_2D.push_back((TH2F*)finput->Get(Form("T_2D_%s", strTemplateName[t].Data())));
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


void getEvents(
  CJLSTSet* theSet,

  const TString& strTrackVar,
  const vector<TString>& strExtraWeightsList,
  const vector<pair<Discriminant*, vector<TString>>>& KDVars,

  vector<SimpleEntry>& index,

  Channel& chan
  ){
  short matchdecid=-1;
  if (chan==k2e2mu) matchdecid=121*169;
  else if (chan==k4mu) matchdecid=169*169;
  else if (chan==k4e) matchdecid=121*121;

  int ev=0, ev_acc=0;
  CJLSTTree* tree=nullptr;
  while ((tree = theSet->getSelectedEvent(ev))){
    ev++;

    float scale = theSet->getPermanentWeight(tree)*theSet->getOverallEventWgt(tree);
    vector<float> valReco;

    bool doProcess = true;
    // Flavor check
    if (matchdecid>0){
      short Z1Id, Z2Id;
      tree->getVal("Z1Flav", Z1Id);
      tree->getVal("Z2Flav", Z2Id);
      doProcess &= (matchdecid==Z1Id*Z2Id);
    }

    for (auto& KDVar:KDVars){
      Discriminant*& KDbuilder = KDVar.first;
      vector<TString>& strKDVarsList = KDVar.second;
      vector<float> KDBuildVals;
      for (auto const& s : strKDVarsList){
        float tmp;
        tree->getVal(s, tmp);
        KDBuildVals.push_back(tmp);
      }
      float KD = KDbuilder->update(KDBuildVals);
      valReco.push_back(KD);
      doProcess &= !(std::isnan(KD) || std::isinf(KD));
    }

    if (!doProcess) continue;

    float varTrack;
    tree->getVal(strTrackVar, varTrack);

    float wgt = scale;
    for (auto const& w : strExtraWeightsList){
      float wval;
      tree->getVal(w, wval);
      wgt *= wval;
    }
    if (std::isnan(wgt) || std::isinf(wgt)){
      // If weight is NaN, it is a big problem.
      if (std::isnan(wgt) || std::isinf(wgt)) cerr << "Invalid weight " << wgt << " is being discarded at mass " << varTrack << endl;
      continue;
    }

    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;

    SimpleEntry theEntry(ev, varTrack, valReco, wgt);
    addByLowest(index, theEntry, false);

    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;
}