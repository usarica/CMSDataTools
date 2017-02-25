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
#include "InputTreeHandle.h"
#include "ZZ4l_Samples.h"
#include "external_Category.h"

using namespace std;

// Constants to affect the template code
const TString fixedDate="";

// Initializers
void estimateSystematics_Moriond2017mainstream_one(int channel, int erg_tev, int categorizationType, float ZZMass_low=-1, float ZZMass_high=-1, float pickSampleWithMass=-1);
void progressbar(int val, int tot);
vector<pair<string, double>> getFileList(string strinput);
TString todaysdate();
void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray);
void addByLowest(float weight, float val, std::vector<std::pair<float, float>>& valArray);
void addByLowest(float val, std::vector<float>& valArray);

void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
Bool_t checkListVariable(const vector<string>& list, const string& var);

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
  );

// Main Function, runs over all desired iterations
void estimateSystematics_Moriond2017mainstream_one(int ichan, int erg_tev, int categorizationType, float ZZMass_low, float ZZMass_high, float pickSampleWithMass){
  if (ichan>(int)nChannels) return;

  constructSampleTypeList();

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
  const TString strChannel = channame(ichan);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString OUTPUT_LOG_NAME = Form("HtoZZ%s_Systematics", strChannel.Data());
  if (categorizationType==0) OUTPUT_LOG_NAME.Append("_noCat");
  else if (categorizationType==1) OUTPUT_LOG_NAME.Append("_MainAnalysis_noQG");
  else if (categorizationType==-1) OUTPUT_LOG_NAME.Append("_MainAnalysis_wQG");
  else if (categorizationType==-2) OUTPUT_LOG_NAME.Append("_MainAnalysisLegacyStyle_wQG");
  else if (categorizationType==2) OUTPUT_LOG_NAME.Append("_MainAnalysisLegacyStyle_noQG");

  if (ZZMass_low>0. && ZZMass_high>0.) OUTPUT_LOG_NAME.Append(Form("_mZZlow_%.1f_mZZhigh_%.1f", ZZMass_cut[0], ZZMass_cut[1]));
  TString OUTPUT_NAME = OUTPUT_LOG_NAME;
  OUTPUT_LOG_NAME += ".log";
  OUTPUT_NAME += ".root";

  TString cinput_common = user_input_dir + "/";
  TString coutput_common = user_output_dir + erg_dir + "Systematics/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

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
  vector<pair<string, double>> fileNameList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++) fileNameList[is] = getFileList(string(sampleTypeList[is].first.Data()));
  vector<InputTreeHandle*> treeList[nSampleTypes];
  for (unsigned int is=0; is<nSampleTypes; is++){
    for (unsigned int ifile=0; ifile<fileNameList[is].size(); ifile++){
      if (pickSampleWithMass>0. && !(fileNameList[is].at(ifile).second==-1. || fileNameList[is].at(ifile).second==pickSampleWithMass)) continue;

      string cinput = fileNameList[is].at(ifile).first;
      cinput = user_input_dir + cinput + "/" + user_treefile;

      cout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;
      tout << "Opening file " << cinput << " with mass "  << fileNameList[is].at(ifile).second << endl;

      InputTreeHandle* theTree = new InputTreeHandle(TString(cinput.c_str()), user_treename, user_countersname, sampleTypeList[is].second);
      foutput->cd();
      treeList[is].push_back(theTree);
    }
  }

  const unsigned int nSyst = 13;
  TString strSyst[nSyst]={
    "Nominal",
    "QCDScaleUp", "QCDScaleDn",
    "PDFScaleUp", "PDFScaleDn",
    "AsUp", "AsDn",
    "PDFReplicaUp", "PDFReplicaDn",

    "EWUp", "EWDown",
    "JESUp", "JESDn"
  };

  double ggKfactorScale[nSyst]; for (unsigned int isyst=0; isyst<nSyst; isyst++) ggKfactorScale[isyst]=1;
  tout << endl;
  tout << "Extracting gg NNLO K factors" << endl;
  TString datadir = package_dir + "data/";
  TFile* finput_ggKfactors = TFile::Open(datadir+"Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root", "read");
  for (unsigned int isyst=0; isyst<9; isyst++){
    TString strkfactor = "sp_kfactor_"; strkfactor += strSyst[isyst];
    TSpline3* spkfactor_ggzz_nnlo = (TSpline3*)finput_ggKfactors->Get(strkfactor);
    ggKfactorScale[isyst] = spkfactor_ggzz_nnlo->Eval(125.);
    tout << "gg NNLO K factor at 125 GeV for " << strSyst[isyst] << " = " << ggKfactorScale[isyst] << endl;
    delete spkfactor_ggzz_nnlo;
    if (isyst>0) ggKfactorScale[isyst] = ggKfactorScale[isyst]/ggKfactorScale[0]-1.;
  }
  ggKfactorScale[0]=1;
  {
    double YR4_ggScaleUp=1.046;
    double YR4_ggScaleDn=0.946;
    double YR4_ggPDFAsUp=1.033;
    double YR4_ggPDFAsDn=0.967;

    double NNLO_ggScaleUp=1.+sqrt(pow(max(ggKfactorScale[1], ggKfactorScale[2]), 2) + pow(max(ggKfactorScale[3], ggKfactorScale[4]), 2)); tout << "NNLO_ggScaleUp = " << NNLO_ggScaleUp << endl;
    double NNLO_ggScaleDn=1.-sqrt(pow(min(ggKfactorScale[1], ggKfactorScale[2]), 2) + pow(min(ggKfactorScale[3], ggKfactorScale[4]), 2)); tout << "NNLO_ggScaleDn = " << NNLO_ggScaleDn << endl;
    double NNLO_ggPDFAsUp=1.+sqrt(pow(max(ggKfactorScale[5], ggKfactorScale[6]), 2) + pow(max(ggKfactorScale[7], ggKfactorScale[8]), 2)); tout << "NNLO_ggPDFAsUp = " << NNLO_ggPDFAsUp << endl;
    double NNLO_ggPDFAsDn=1.-sqrt(pow(min(ggKfactorScale[5], ggKfactorScale[6]), 2) + pow(min(ggKfactorScale[7], ggKfactorScale[8]), 2)); tout << "NNLO_ggPDFAsDn = " << NNLO_ggPDFAsDn << endl;

    for (unsigned int isyst=1; isyst<5; isyst++){
      if (ggKfactorScale[isyst]>0.) ggKfactorScale[isyst] = YR4_ggScaleUp/NNLO_ggScaleUp;
      else ggKfactorScale[isyst] = YR4_ggScaleDn/NNLO_ggScaleDn;
    }
    for (unsigned int isyst=5; isyst<9; isyst++){
      if (ggKfactorScale[isyst]>0.) ggKfactorScale[isyst] = YR4_ggPDFAsUp/NNLO_ggPDFAsUp;
      else ggKfactorScale[isyst] = YR4_ggPDFAsDn/NNLO_ggPDFAsDn;
    }
    for (unsigned int isyst=0; isyst<9; isyst++) tout << "Will scale K factor uncertainty for " << strSyst[isyst] << " by " << ggKfactorScale[isyst] << endl;
  }
  finput_ggKfactors->Close();
  tout << endl;
  tout << "Process" << " " << "Category" << " ";
  for (unsigned int isyst=0; isyst<nSyst; isyst++) tout << strSyst[isyst] << " ";
  tout << endl;

  for (int is=0; is<nSampleTypes; is++){
    TString SampleTypeName = nameSampleType(is);
    cout << "Process " << SampleTypeName << "[" << is << "] under scrutiny" << endl;

    TH1F* hSyst[CategoryMor17::nCategoriesMor17][nSyst];
    for (int icat=(int)CategoryMor17::InclusiveMor17; icat<(int)CategoryMor17::nCategoriesMor17; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        hSyst[icat][isyst] = new TH1F(Form("Process_%s_Category_%s_Syst_%s", SampleTypeName.Data(), nameCategoryMor17(icat).Data(), strSyst[isyst].Data()), "", nbinsx, kDXarray);
        hSyst[icat][isyst]->Sumw2();
      }
    }

    for (unsigned int it=0; it<treeList[is].size(); it++){
      double mass = fileNameList[is].at(it).second;
      InputTreeHandle* theTree = treeList[is].at(it);

      int nEntries = theTree->GetEntries();

      cout << "\tSample[" << it << "]=" << fileNameList[is].at(it).first << " under scrutiny. #Events = " << nEntries << endl;

      float nGenEvents = theTree->GetNGenEvents();
      for (int ev=0; ev<nEntries; ev++){
        theTree->InitDefaults();
        theTree->GetEntry(ev);

        if (ev%10000==0) cout << "Event " << ev << "..." << endl;

        // Fix for JHUGen v6 vegas bug in low mass samples
        if (SampleTypeName.Contains("POWHEG") && mass>0 && mass<300 && theTree->GenHMass>mass+5.) continue;

        if (theTree->ZZsel<70.) continue; // Disregard non-selected events (need this if skipEmptyEvents=false)

        if (theTree->ZZMass<ZZMass_cut[0] || theTree->ZZMass>=ZZMass_cut[1]) continue;

        int ZZFlav=theTree->Z1Flav*theTree->Z2Flav;
        if (ichan==0 && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (ichan==1 && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (ichan==2 && abs(ZZFlav)!=pow(13*11, 2)) continue;

        int icat = getCategory(
          theTree->ZZMass,

          *(theTree->JetPt),
          *(theTree->JetSigma),
          *(theTree->JetPhi),
          *(theTree->JetQGLikelihood),
          *(theTree->JetIsBtaggedWithSF),

          0,
          theTree->nExtraLep,
          theTree->nExtraZ,
          theTree->PFMET,

          theTree->p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,
          theTree->p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,
          theTree->p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
          theTree->p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
          theTree->pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
          theTree->p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
          theTree->p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,

          categorizationType
          );
        int icat_jesup=getCategory(
          theTree->ZZMass,

          *(theTree->JetPt),
          *(theTree->JetSigma),
          *(theTree->JetPhi),
          *(theTree->JetQGLikelihood),
          *(theTree->JetIsBtaggedWithSF),

          1,
          theTree->nExtraLep,
          theTree->nExtraZ,
          theTree->PFMET_jesUp,

          theTree->p_JJQCD_SIG_ghg2_1_JHUGen_JECUp,
          theTree->p_JQCD_SIG_ghg2_1_JHUGen_JECUp,
          theTree->p_JJVBF_SIG_ghv1_1_JHUGen_JECUp,
          theTree->p_JVBF_SIG_ghv1_1_JHUGen_JECUp,
          theTree->pAux_JVBF_SIG_ghv1_1_JHUGen_JECUp,
          theTree->p_HadWH_SIG_ghw1_1_JHUGen_JECUp,
          theTree->p_HadZH_SIG_ghz1_1_JHUGen_JECUp,

          categorizationType
          );
        int icat_jesdn=getCategory(
          theTree->ZZMass,

          *(theTree->JetPt),
          *(theTree->JetSigma),
          *(theTree->JetPhi),
          *(theTree->JetQGLikelihood),
          *(theTree->JetIsBtaggedWithSF),

          -1,
          theTree->nExtraLep,
          theTree->nExtraZ,
          theTree->PFMET_jesDn,

          theTree->p_JJQCD_SIG_ghg2_1_JHUGen_JECDn,
          theTree->p_JQCD_SIG_ghg2_1_JHUGen_JECDn,
          theTree->p_JJVBF_SIG_ghv1_1_JHUGen_JECDn,
          theTree->p_JVBF_SIG_ghv1_1_JHUGen_JECDn,
          theTree->pAux_JVBF_SIG_ghv1_1_JHUGen_JECDn,
          theTree->p_HadWH_SIG_ghw1_1_JHUGen_JECDn,
          theTree->p_HadZH_SIG_ghz1_1_JHUGen_JECDn,

          categorizationType
          );

        float weight[nSyst]={ 0 };
        /*
        "Nominal",
        "QCDScaleUp", "QCDScaleDown",
        "PDFScaleUp", "PDFScaleDown",
        "AsUp", "AsDn",
        "PDFReplicaUp", "PDFReplicaDn",
        "EWUp", "EWDown",
        "JESUp", "JESDn"
        */
        if ((is>=(int)SampleType::gg_Sig && is<=(int)SampleType::gg_Sig_PythiaTuneUp_MiNLO) || is==(int)SampleType::gg_Bkg_MCFM){ // gg samples
          weight[0] = theTree->KFactor_QCD_ggZZ_Nominal;
          weight[1] = theTree->KFactor_QCD_ggZZ_QCDScaleUp;
          weight[2] = theTree->KFactor_QCD_ggZZ_QCDScaleDn;
          weight[3] = theTree->KFactor_QCD_ggZZ_PDFScaleUp;
          weight[4] = theTree->KFactor_QCD_ggZZ_PDFScaleDn;
          weight[5] = theTree->KFactor_QCD_ggZZ_AsUp;
          weight[6] = theTree->KFactor_QCD_ggZZ_AsDn;
          weight[7] = theTree->KFactor_QCD_ggZZ_PDFReplicaUp;
          weight[8] = theTree->KFactor_QCD_ggZZ_PDFReplicaDn;
          weight[9] = theTree->KFactor_QCD_ggZZ_Nominal;
          weight[10] = theTree->KFactor_QCD_ggZZ_Nominal;
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->Fill(theTree->ZZMass, weight[isyst]*ggKfactorScale[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);

          if (is>=(int)SampleType::gg_Sig && is<=(int)SampleType::gg_Sig_PythiaTuneUp_MiNLO){ // theTree->LHE weights are meaningless for gg bkg samples at LO
            weight[1] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_QCDscale_muR2_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[2] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_QCDscale_muR0p5_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[3] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_QCDscale_muR1_muF2/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[4] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_QCDscale_muR1_muF0p5/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[5] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_AsMZ_Up/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[6] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_AsMZ_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[7] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_PDFVariation_Up/theTree->LHEweight_QCDscale_muR1_muF1;
            weight[8] = theTree->KFactor_QCD_ggZZ_Nominal * theTree->LHEweight_PDFVariation_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
          }

          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)CategoryMor17::InclusiveMor17) hSyst[icat][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesup][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesdn][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          }
        }
        else if (
          (is>=(int)SampleType::VBF_Sig && is<=(int)SampleType::VBF_Sig_PythiaTuneUp)
          ||
          (is>=(int)SampleType::WH_Sig && is<=(int)SampleType::WH_Sig_PythiaTuneUp)
          ||
          (is>=(int)SampleType::ZH_Sig && is<=(int)SampleType::ZH_Sig_PythiaTuneUp)
          ||
          (is>=(int)SampleType::tt_Sig && is<=(int)SampleType::tt_Sig_PythiaTuneUp)
          ){ // VBF, VH and ttH samples
          weight[0] = 1;
          weight[1] = theTree->LHEweight_QCDscale_muR2_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[2] = theTree->LHEweight_QCDscale_muR0p5_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[3] = theTree->LHEweight_QCDscale_muR1_muF2/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[4] = theTree->LHEweight_QCDscale_muR1_muF0p5/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[5] = theTree->LHEweight_AsMZ_Up/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[6] = theTree->LHEweight_AsMZ_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[7] = theTree->LHEweight_PDFVariation_Up/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[8] = theTree->LHEweight_PDFVariation_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[9] = 1;
          weight[10] = 1;
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)CategoryMor17::InclusiveMor17) hSyst[icat][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesup][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesdn][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          }

        }
        else if (is==(int)SampleType::qq_Bkg){
          weight[0] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ;
          weight[1] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_QCDscale_muR2_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[2] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_QCDscale_muR0p5_muF1/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[3] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_QCDscale_muR1_muF2/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[4] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_QCDscale_muR1_muF0p5/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[5] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_AsMZ_Up/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[6] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_AsMZ_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[7] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_PDFVariation_Up/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[8] = theTree->KFactor_QCD_qqZZ_M*theTree->KFactor_EW_qqZZ*theTree->LHEweight_PDFVariation_Dn/theTree->LHEweight_QCDscale_muR1_muF1;
          weight[9] = theTree->KFactor_QCD_qqZZ_M*(theTree->KFactor_EW_qqZZ+theTree->KFactor_EW_qqZZ_unc);
          weight[10] = theTree->KFactor_QCD_qqZZ_M*(theTree->KFactor_EW_qqZZ-theTree->KFactor_EW_qqZZ_unc);
          // JES Up, Dn have same weight as nominal
          weight[11] = weight[0];
          weight[12] = weight[0];

          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)CategoryMor17::InclusiveMor17) hSyst[icat][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesup][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesdn][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          }
        }
        else{
          for (unsigned int isyst=0; isyst<nSyst; isyst++) weight[isyst] = 1;
          for (unsigned int isyst=0; isyst<nSyst; isyst++) hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          for (unsigned int isyst=0; isyst<nSyst; isyst++){
            if (isyst<11 && icat!=(int)CategoryMor17::InclusiveMor17) hSyst[icat][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==11 && icat_jesup!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesup][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
            else if (isyst==12 && icat_jesdn!=(int)CategoryMor17::InclusiveMor17) hSyst[icat_jesdn][isyst]->Fill(theTree->ZZMass, weight[isyst]*theTree->overallEventWeight*theTree->xsec/nGenEvents);
          }
        }

      }
    }


    // Rescale sum of categories in each systematic scenario to inclusive
    for (unsigned int isyst=0; isyst<nSyst; isyst++){
      for (int ix=0; ix<=hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->GetNbinsX()+1; ix++){
        double sum_inc=hSyst[(int)CategoryMor17::InclusiveMor17][isyst]->GetBinContent(ix);
        double sum_cat=0;
        for (int icat=(int)CategoryMor17::InclusiveMor17+1; icat<(int)CategoryMor17::nCategoriesMor17; icat++) sum_cat += hSyst[icat][isyst]->GetBinContent(ix);
        double rescale=1;
        cout << "isyst=" << isyst << " inc=" << sum_inc << " sum(cat)=" << sum_cat << endl;
        if (sum_cat!=0.) rescale = sum_inc / sum_cat;
        for (int icat=(int)CategoryMor17::InclusiveMor17+1; icat<(int)CategoryMor17::nCategoriesMor17; icat++){
          hSyst[icat][isyst]->SetBinContent(ix, hSyst[icat][isyst]->GetBinContent(ix)*rescale);
          hSyst[icat][isyst]->SetBinError(ix, hSyst[icat][isyst]->GetBinError(ix)*rescale);
        }
      }
    }

    // Print the ratio of each integral to inclusive nominal into the output file
    float integral_nom = hSyst[(int)CategoryMor17::InclusiveMor17][0]->Integral(1, hSyst[(int)CategoryMor17::InclusiveMor17][0]->GetNbinsX());
    for (int icat=(int)CategoryMor17::InclusiveMor17; icat<(int)CategoryMor17::nCategoriesMor17; icat++){
      if (
        (categorizationType==0 && icat>(int)CategoryMor17::InclusiveMor17)
        ||
        (abs(categorizationType)==2 && icat>(int)CategoryMor17::VBF2jTaggedMor17)
        ) continue;
      tout << SampleTypeName << " " << nameCategoryMor17(icat) << " ";
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        double integral = hSyst[icat][isyst]->Integral(1, hSyst[icat][isyst]->GetNbinsX());
        if (integral_nom!=0.) integral /= integral_nom;
        else integral=0;
        tout << integral << " ";
      }
      tout << endl;
    }
    tout << endl;

    for (int icat=(int)CategoryMor17::InclusiveMor17; icat<(int)CategoryMor17::nCategoriesMor17; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        if (icat==(int)CategoryMor17::InclusiveMor17 && isyst==0) continue; // Do not rescale inclusive nominal since you need this to divide everything by.
        for (int ix=0; ix<=hSyst[icat][isyst]->GetNbinsX()+1; ix++){
          double val = hSyst[icat][isyst]->GetBinContent(ix);
          double err = hSyst[icat][isyst]->GetBinError(ix);
          double val_nom = hSyst[(int)CategoryMor17::InclusiveMor17][0]->GetBinContent(ix);
          if (val_nom!=0.){ val /= val_nom; err /= val_nom; }
          else{ val=0; err=0; }
          hSyst[icat][isyst]->SetBinContent(ix, val);
          hSyst[icat][isyst]->SetBinError(ix, err);
        }
      }
    }

    for (int icat=(int)CategoryMor17::InclusiveMor17; icat<(int)CategoryMor17::nCategoriesMor17; icat++){
      for (unsigned int isyst=0; isyst<nSyst; isyst++){
        foutput->WriteTObject(hSyst[icat][isyst]);
        delete hSyst[icat][isyst];
      }
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

