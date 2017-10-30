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
#include "BaseTreeLooper.h"
#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ReweightingBuilder.h"
#include "ReweightingFunctions.h"
#include "ExtendedBinning.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"


using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace CategorizationHelpers;
using namespace DiscriminantClasses;
using namespace ReweightingFunctions;
using namespace MELAStreamHelpers;


// Constants to affect the template code
const TString fixedDate="";
const TString user_output_dir = "output/";


// ggH analyzer
class GGHAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool runEvent(CJLSTTree* tree, SimpleEntry& product);

public:
  GGHAnalyzer(Channel channel_, Category category_) : BaseTreeLooper(), channel(channel_), category(category_){}
  GGHAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) : BaseTreeLooper(inTree), channel(channel_), category(category_){}
  GGHAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) : BaseTreeLooper(inTreeList), channel(channel_), category(category_){}
  GGHAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) : BaseTreeLooper(inTreeSet), channel(channel_), category(category_){}

};

// Function to build one templates
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// theSqrts = 13 (CoM energy) is fixed in Samples.h
void makeGGZZTemplates_one(const Channel channel, const Category category, TString strSystematics){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

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
  ExtendedBinning ZZMassBinning(2900/2., 100., 3000., "ZZMass");
  ExtendedBinning KD1Binning(30, 0, 1, "KD1");
  //ExtendedBinning KD2Binning(30, -1, 1, "KD2");

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form("%s_HtoZZ%s_Stage1_%s", strCategory.Data(), strChannel.Data(), strSystematics.Data());
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";
  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  // Get list of samples
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("gg_Sig_POWHEG");
  vector<TString> strSamples;
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples);

  // Register the discriminants
  Discriminant* KD1 = DiscriminantClasses::constructKDFromType(DiscriminantClasses::kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
  vector<TString> KD1vars = DiscriminantClasses::getKDVars(DiscriminantClasses::kDbkgkin);
  Discriminant* KD2=nullptr;
  vector<TString> KD2vars;
  //Discriminant* KD2 = DiscriminantClasses::constructKDFromType(DiscriminantClasses::kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
  //vector<TString> KD2vars = DiscriminantClasses::getKDVars(DiscriminantClasses::kDggint);

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);
  // Book common variables
  theSampleSet->bookXS();
  theSampleSet->bookOverallEventWgt();
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("GenHMass", 0);
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for SM reweighting
    tree->bookBranch<float>("p_Gen_GG_BKG_MCFM", 0);
    tree->bookBranch<float>("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM", 0);
    tree->bookBranch<float>("p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM", 0);
    // Variables for KDs
    for (auto& v:KD1vars) tree->bookBranch<float>(v, 0);
    for (auto& v:KD2vars) tree->bookBranch<float>(v, 0);
  }

  // Setup GenHMass binning
  // Binning for PUGenHepRewgt
  ExtendedBinning GenHMassInclusiveBinning("GenHMass");
  GenHMassInclusiveBinning.addBinBoundary(0);
  GenHMassInclusiveBinning.addBinBoundary(theSqrts*1000.);

  // Binning for MELARewgt
  ExtendedBinning GenHMassBinning("GenHMass");
  GenHMassBinning.addBinBoundary(0);
  for (unsigned int is=0; is<theSampleSet->getCJLSTTreeList().size()-1; is++) GenHMassBinning.addBinBoundary(
    0.5*(theSampleSet->getCJLSTTreeList().at(is)->MHVal + theSampleSet->getCJLSTTreeList().at(is+1)->MHVal)
    );
  GenHMassBinning.addBinBoundary(theSqrts*1000.);

  // Construct PUGenHepRewgt
  vector<TString> strReweightingWeigths;
  strReweightingWeigths.push_back("PUWeight");
  strReweightingWeigths.push_back("genHEPMCweight");
  ReweightingBuilder* pugenheprewgtBuilder = new ReweightingBuilder(strReweightingWeigths, getSimpleWeight);
  pugenheprewgtBuilder->setWeightBinning(GenHMassInclusiveBinning);
  for (auto& tree:theSampleSet->getCJLSTTreeList()) pugenheprewgtBuilder->setupWeightVariables(tree, 1.);



  for (int t=kBkg; t<nTemplates; t++){
    foutput->cd();
    TTree* theFinalTree = new TTree(Form("T_ggH_%s_Tree", strTemplateName[t].Data()), ""); // The tree to record into the ROOT file

    /************* Reweighting setup *************/
    // There are two builders:
    // 1) Rewighting from MELA (x) K factors, which adjust the cross section
    // 2) PU and GenHepMCWeight reweighting, which are supposed to keep the cross section the same
    // Total weight is (1)x(2)

    // Build the possible MELA reweightings
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
    strReweightingWeigths.clear();
    strReweightingWeigths.push_back(strWeight);
    strReweightingWeigths.push_back("xsec");
    ReweightingBuilder* melarewgtBuilder = new ReweightingBuilder(strReweightingWeigths, getSimpleWeight);
    melarewgtBuilder->setWeightBinning(GenHMassBinning);
    for (auto& tree:theSampleSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(tree);

    // Build the analyzer and loop over the events
    GGHAnalyzer theAnalyzer(theSampleSet, channel, category);
    // Book common variables needed for analysis
    theAnalyzer.addConsumed<float>("dataMCWeight");
    theAnalyzer.addConsumed<float>("trigEffWeight");
    theAnalyzer.addConsumed<float>("GenHMass");
    theAnalyzer.addConsumed<float>("ZZMass");
    theAnalyzer.addConsumed<short>("Z1Flav");
    theAnalyzer.addConsumed<short>("Z2Flav");
    // Add discriminant builders
    theAnalyzer.addDiscriminantBuilder("KD1", KD1, KD1vars);
    theAnalyzer.addDiscriminantBuilder("KD2", KD2, KD2vars);
    // Add reweighting builders
    theAnalyzer.addReweightingBuilder("MELARewgt", melarewgtBuilder);
    theAnalyzer.addReweightingBuilder("PUGenHEPRewgt", pugenheprewgtBuilder);
    // Loop
    theAnalyzer.loop(true, false, true);

    const std::vector<SimpleEntry>& products = theAnalyzer.getProducts();
    unordered_map<TString, float> theTreeFloats;
    MELAout << "There are " << products.size() << " products" << endl;
    if (products.at(0).namedfloats.find("ZZMass")==products.at(0).namedfloats.end()) MELAerr << "Uh-oh! ZZMass could not be found in the products!" << endl;
    for (auto it=products.at(0).namedfloats.begin(); it!=products.at(0).namedfloats.end(); it++){
      MELAout << "Booking branch " << it->first << " in " << theFinalTree->GetName() << endl;
      theTreeFloats[it->first]=0;
      theFinalTree->Branch(it->first, &(theTreeFloats[it->first]));
    }
    for (auto& product:products){
      for (auto it=product.namedfloats.begin(); it!=product.namedfloats.end(); it++) theTreeFloats[it->first] = it->second;
      theFinalTree->Fill();
    }

    delete melarewgtBuilder;

    foutput->WriteTObject(theFinalTree);
    delete theFinalTree;
  }

  delete pugenheprewgtBuilder;
  delete KD2;
  delete KD1;
  delete theSampleSet;
  foutput->Close();
  MELAout.close();
}


bool GGHAnalyzer::runEvent(CJLSTTree* tree, SimpleEntry& product){
  short matchdecid=-1;
  if (channel==k2e2mu) matchdecid=121*169;
  else if (channel==k4mu) matchdecid=169*169;
  else if (channel==k4e) matchdecid=121*121;

  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    float& GenHMass = *(valfloats["GenHMass"]);
    product.setNamedVal("ZZMass", ZZMass);
    product.setNamedVal("GenHMass", GenHMass);

    // Construct the weights
    float pure_reco_wgt = (*(valfloats["dataMCWeight"]))*(*(valfloats["trigEffWeight"]));

    auto& pugenheprewgtBuilder = Rewgtbuilders["PUGenHEPRewgt"];
    float pugenhep_wgt = pugenheprewgtBuilder->getPostThresholdWeight(tree)/pugenheprewgtBuilder->getSumPostThresholdWeights(tree); // Normalized to unit

    auto& melarewgtBuilder = Rewgtbuilders["MELARewgt"];
    float mela_wgt = melarewgtBuilder->getPostThresholdWeight(tree)/melarewgtBuilder->getSumPostThresholdWeights(tree); // Normalized to unit
    mela_wgt *= melarewgtBuilder->getWeightedSumAllPostThresholdWeights(tree);

    float wgt = pure_reco_wgt*pugenhep_wgt*mela_wgt;
    product.setNamedVal("weight", wgt);
    if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.){
      MELAerr << "Invalid weight " << wgt << " is being discarded at mass " << ZZMass << endl;
      validProducts=false;
    }

    // Flavor check
    if (matchdecid>0){
      short& Z1Id = *(valshorts["Z1Flav"]);
      short& Z2Id = *(valshorts["Z2Flav"]);
      validProducts &= (matchdecid==Z1Id*Z2Id);
    }

    // Compute the KDs
    for (auto it=KDbuilders.cbegin(); it!=KDbuilders.cend(); it++){
      auto& KDbuilderpair = it->second;
      auto& KDbuilder = KDbuilderpair.first;
      auto& strKDVarsList = KDbuilderpair.second;
      vector<float> KDBuildVals; KDBuildVals.reserve(strKDVarsList.size());
      for (auto const& s : strKDVarsList) KDBuildVals.push_back(*(valfloats[s]));
      float KD = KDbuilder->update(KDBuildVals);
      product.setNamedVal(it->first, KD);
      validProducts &= !(std::isnan(KD) || std::isinf(KD));
    }
  }

  return validProducts;
}