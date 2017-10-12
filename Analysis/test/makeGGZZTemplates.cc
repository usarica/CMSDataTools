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
#include "ReweightingBuilder.h"
#include "ReweightingFunctions.h"
#include "ExtendedBinning.h"
#include "StringStreamer.h"
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

// FIXME:
// 2) WRITE A REWEIGHTING CLASS THAT TAKES THE CJLSTTREE AND THE WEIGHT STRING TO CONTROL IF PROPAGATOR NEEDS TO BE TAKEN OUT
// 3) ADD A REWEIGHTING CLASSES FIELD TO GETEVENTS
// 4) IMPLEMENT A CATEGORIZATION PLAN TO GETEVENTS AS WELL
// Event scanner
void getEvents(
  CJLSTSet* theSet,

  const TString& strTrackVar,
  const vector<pair<Discriminant*, vector<TString>>>& KDVars,
  const vector<ReweightingBuilder*>& rewgtBuilders,

  TTree* indexTree,

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

  // Setup binning
  ExtendedBinning ZZMassBinning(2900/2., 100., 3000., "mZZ");
  ExtendedBinning KD1Binning(30, 0, 1, "KD1");

  // Setup the output directories
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
  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  StringStreamer tout(coutput_log.Data());
  tout << "Opened log file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");
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

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);

  foutput->cd();

  TTree* theFinalTree[nTemplates];
  for (int t=kBkg; t<nTemplates; t++){
    // Build the possible reweightings
    ReweightingBuilder* rewgtBuilder=nullptr;
    TString strWeight;
    switch (t){
    case kBkg:
      strWeight = "p_Gen_GG_BKG_MCFM";
      break;
    case kSig:
      strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";
      break;
    case kBSI:
      strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM";
      break;
    };
    rewgtBuilder = new ReweightingBuilder(strWeight, getSimpleWeight);
    for (auto& tree:theSampleSet->getCJLSTTreeList()) rewgtBuilder->setupWeightVariables(tree, ZZMassBinning.getBoundaryPairsList<float>(), "ZZMass");

    theFinalTree[t] = new TTree(Form("T_2D_%s_Tree", strTemplateName[t].Data()), "");


  }

  delete theSampleSet;

  for (int t=(int) kBkg; t<(int)nTemplates; t++){
    foutput->WriteTObject(theFinalTree[t]);
    delete theFinalTree[t];
  }

  foutput->Close();
}


void getEvents(
  CJLSTSet* theSet,

  const TString& strTrackVar,
  const vector<pair<Discriminant*, vector<TString>>>& KDVars,
  const vector<ReweightingBuilder*>& rewgtBuilders,

  TTree* indexTree,

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