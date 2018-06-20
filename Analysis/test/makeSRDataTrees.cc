#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "TLeaf.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

void makeSRDataTrees_one(const Channel channel, const Category category, ACHypothesis const& hypo, MassRegion massregion, const TString fixedDate){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strMassRegion = getMassRegionName(massregion);
  const TString strSqrts = Form("%i", theSqrts);
  const TString strYear = theDataPeriod;

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  MELAout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Data/" + strdate + "/" + strACHypo + "/" + strMassRegion + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString coutput = "[outdir]hzz[channel]_[category]_[sqrts]TeV_[year].root";
  replaceString<TString, const TString>(coutput, TString("[outdir]"), coutput_common);
  replaceString<TString, const TString>(coutput, TString("[channel]"), strChannel);
  replaceString<TString, const TString>(coutput, TString("[category]"), strCategory);
  replaceString<TString, const TString>(coutput, TString("[sqrts]"), strSqrts);
  replaceString<TString, const TString>(coutput, TString("[year]"), strYear);
  TString coutput_log = coutput;
  replaceString(coutput_log, ".root", ".log");
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  // Register the discriminants
  vector<KDspecs> KDlist;
  vector<TString> strExtraCatVars_short;
  getLikelihoodDiscriminants(channel, category, sNominal, KDlist);
  if (category!=Inclusive){
    getCategorizationDiscriminants(sNominal, KDlist);
    getExtraCategorizationVariables<short>(globalCategorizationScheme, sNominal, strExtraCatVars_short);
  }

  // Get list of samples
  vector<TString> strSamples;
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("AllData");
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples, sNominal);

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    tree->bookBranch<int>("RunNumber", 0);
    tree->bookBranch<long long>("EventNumber", 0);
    // Variables for KDs
    for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
    // Extra categorization variables
    for (auto& s:strExtraCatVars_short) tree->bookBranch<short>(s, -1);
    tree->silenceUnused(); // Will no longer book another branch
  }
  theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_None, false, false);

  ExtendedBinning binning_mass = getDiscriminantFineBinning(channel, category, hypo, "ZZMass", massregion);

  // Construct reweighting variables vector
  foutput->cd();

  TString treename = "data_obs_tree";
  BaseTree* theFinalTree = new BaseTree(treename); // The tree to record into the ROOT file

  // Build the analyzer and loop over the events
  TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
  theAnalyzer.setRecordCategorizationKDs(true);
  //theAnalyzer.setRecordKDVariables(true);
  theAnalyzer.setExternalProductTree(theFinalTree);
  theAnalyzer.addMassWindow(pair<float, float>(binning_mass.getMin(), binning_mass.getMax()));
  // Book common variables needed for analysis
  theAnalyzer.addConsumed<float>("ZZMass");
  theAnalyzer.addConsumed<short>("Z1Flav");
  theAnalyzer.addConsumed<short>("Z2Flav");
  theAnalyzer.addConsumed<int>("RunNumber");
  theAnalyzer.addConsumed<long long>("EventNumber");
  // Add discriminant builders
  for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
  // Add extra categorization variables
  for (auto& s:strExtraCatVars_short) theAnalyzer.addConsumed<short>(s);
  // Loop
  theAnalyzer.setSampleIdStorageOption(BaseTreeLooper::kStoreByRunAndEventNumber);
  theAnalyzer.loop(true, false, true);
  MELAout << "There are " << theFinalTree->getNEvents() << " products" << endl;

  // Setup actual tree variables
  TString catFlagName="";
  if (category!=Inclusive) catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;
  vector<TString> KDset;
  {
    vector<TString> KDset2=getACHypothesisKDNameSet(hypo, category, massregion);
    if (massregion!=kOnshell || hypo==kSM) KDset.push_back("ZZMass"); // Only off-shell, or on-shell SM use ZZMass
    appendVector(KDset, KDset2);
  }
  unsigned int nKDs = KDset.size();
  MELAout << "Number of template dimensions = " << nKDs << ": " << KDset << endl;
  TTree* tin = theFinalTree->getSelectedTree();
  TTree* tout=nullptr;
  {
    unsigned int jKD=0;
    for (unsigned int iKD=0; iKD<KDset.size(); iKD++){
      TString newbranchname;
      if (KDset.at(iKD)=="ZZMass") newbranchname="mass";
      else{
        jKD++;
        newbranchname = Form("KD%i", jKD);
      }
      TBranch* br = tin->GetBranch(KDset.at(iKD));
      br->SetName(newbranchname);
      TLeaf* lf = tin->GetLeaf(KDset.at(iKD));
      lf->SetName(newbranchname);
    }
  }
  tout = tin->CopyTree(catFlagName);
  tout->SetName("data_obs");
  foutput->WriteTObject(tout);

  delete tout;
  delete theFinalTree;
  for (auto& KD:KDlist) delete KD.KD;
  delete theSampleSet;
  foutput->Close();
  MELAout.close();
}

void makeSRDataTrees(CategorizationHelpers::MassRegion massregion, const TString fixedDate=""){
  if (!CheckSetTemplatesCategoryScheme(Inclusive)) return;
  vector<Category> allowedCats = getAllowedCategories(globalCategorizationScheme);
  for (int ch=0; ch<(int) NChannels; ch++){
    Channel channel = (Channel) ch;
    if (channel==k4l || channel==k2l2l) continue;
    for (int ih=0; ih<nACHypotheses; ih++){
      ACHypothesis hypo = (ACHypothesis) ih;
      for (auto& category:allowedCats) makeSRDataTrees_one(channel, category, hypo, massregion, fixedDate);
    }
  }
}
