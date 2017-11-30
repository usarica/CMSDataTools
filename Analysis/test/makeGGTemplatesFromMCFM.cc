#include "common_includes.h"
#include "GGAnalyzer.h"


// Constants to affect the template code
const TString fixedDate="";
const TString user_output_dir = "output/";
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
TString getMELAHypothesisWeight(int const ihypo){
  TString strWeight;
  switch (ihypo){
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
  return strWeight;
}

TTree* fixTreeWeights(TTree* tree);

// Function to build one templates
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// theSqrts = 13 (CoM energy) is fixed in Samples.h
void makeGGTemplatesFromMCFM_one(const Channel channel, const Category category, TString strSystematics){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

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

  TString OUTPUT_NAME = Form("%s_ggHtoZZ%s_MCFM_Stage1_%s", strCategory.Data(), strChannel.Data(), strSystematics.Data());
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
  const unsigned int nMCFMChannels=6;
  TString strMCFMChannels[nMCFMChannels]={
    "4mu","4e","2e2mu",
    "2e2tau","2mu2tau","4tau"
  };
  vector<TString> strSamples[nMCFMChannels];
  for (unsigned int ich=0; ich<nMCFMChannels; ich++){
    vector<TString> strSampleIdentifiers;
    strSampleIdentifiers.push_back(Form("gg_Sig_MCFM_%s", strMCFMChannels[ich].Data()));
    getSamplesList(theSqrts, strSampleIdentifiers, strSamples[ich]);
  }

  // Kfactor variable names
  vector<TString> strKfactorVars;
  if (strSystematics == "Nominal") strKfactorVars.push_back("KFactor_QCD_ggZZ_Nominal");
  else if (strSystematics == "PDFScaleDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFScaleDn");
  else if (strSystematics == "PDFScaleUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFScaleUp");
  else if (strSystematics == "QCDScaleDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_QCDScaleDn");
  else if (strSystematics == "QCDScaleUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_QCDScaleUp");
  else if (strSystematics == "AsDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_AsDn");
  else if (strSystematics == "AsUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_AsUp");
  else if (strSystematics == "PDFReplicaDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFReplicaDn");
  else if (strSystematics == "PDFReplicaUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFReplicaUp");

  // Register the discriminants
  Discriminant* KD1 = DiscriminantClasses::constructKDFromType(DiscriminantClasses::kDbkgkin, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
  vector<TString> KD1vars = DiscriminantClasses::getKDVars(DiscriminantClasses::kDbkgkin);
  Discriminant* KD2=nullptr;
  vector<TString> KD2vars;
  //Discriminant* KD2 = DiscriminantClasses::constructKDFromType(DiscriminantClasses::kDggint, Form("%s%s%s", "../data/SmoothKDConstant_m4l_Dbkgkin_", strChannel.Data(), "13TeV.root"), "sp_gr_varTrue_Constant_Smooth");
  //vector<TString> KD2vars = DiscriminantClasses::getKDVars(DiscriminantClasses::kDggint);

  // Get the CJLST set
  //vector<TString> newlist; newlist.push_back(strSamples.back()); newlist.push_back(strSamples.front());
  //strSamples = newlist;
  vector<CJLSTSet*> theSets;
  for (unsigned int ich=0; ich<nMCFMChannels; ich++){
    CJLSTSet* theSampleSet = new CJLSTSet(strSamples[ich]);
    // Book common variables
    theSampleSet->bookXS(); // "xsec"
    theSampleSet->bookOverallEventWgt(); // Gen weigts "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
    for (auto& tree:theSampleSet->getCJLSTTreeList()){
      // Book common variables needed for analysis
      tree->bookBranch<float>("GenHMass", 0);
      tree->bookBranch<float>("ZZMass", -1);
      tree->bookBranch<short>("Z1Flav", 0);
      tree->bookBranch<short>("Z2Flav", 0);
      // Common variables for reweighting
      for (auto& s:strKfactorVars) tree->bookBranch<float>(s, 1);
      // Variables for SM reweighting
      tree->bookBranch<float>("p_Gen_GG_BKG_MCFM", 0);
      tree->bookBranch<float>("p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM", 0);
      tree->bookBranch<float>("p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM", 0);
      // Variables for KDs
      for (auto& v:KD1vars) tree->bookBranch<float>(v, 0);
      for (auto& v:KD2vars) tree->bookBranch<float>(v, 0);
      tree->silenceUnused(); // Will no longer book another branch
    }
    theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_NgenOverNgenWPU, false, true);
    theSets.push_back(theSampleSet);
  }

  // Setup GenHMass binning
  // Binning for inclusive reweighting
  ExtendedBinning GenHMassInclusiveBinning("GenHMass");
  // Binning for MELARewgt
  ExtendedBinning GenHMassBinning("GenHMass");
  GenHMassBinning.addBinBoundary(0);
  GenHMassBinning.addBinBoundary(124);
  GenHMassBinning.addBinBoundary(126);
  GenHMassBinning.addBinBoundary(160);
  GenHMassBinning.addBinBoundary(220);
  GenHMassBinning.addBinBoundary(1000);
  GenHMassBinning.addBinBoundary(theSqrts*1000.);

  // Construct reweighting variables vector
  vector<TString> strReweightingWeigths;
  for (int t=kBkg; t<nTemplates; t++){
    foutput->cd();

    TDirectory* controlsDir = foutput->mkdir(Form("controls_%s", strTemplateName[t].Data()), "");
    controlsDir->cd();
    TH1F* h_MELARewgtSumAllNonZeroWgtEventsPerBin = new TH1F("MELARewgtSumAllNonZeroWgtEventsPerBin", "", GenHMassBinning.getNbins(), GenHMassBinning.getBinning());
    h_MELARewgtSumAllNonZeroWgtEventsPerBin->GetXaxis()->SetTitle(GenHMassBinning.getLabel());
    TH1F* h_MELARewgtNormComponentPerBin = new TH1F("MELARewgtNormComponentPerBin", "", GenHMassBinning.getNbins(), GenHMassBinning.getBinning());
    h_MELARewgtNormComponentPerBin->GetXaxis()->SetTitle(GenHMassBinning.getLabel());

    foutput->cd();

    TTree* theFinalTree = new TTree(Form("T_ggH_%s_Tree", strTemplateName[t].Data()), ""); // The tree to record into the ROOT file

    /************* Reweighting setup *************/
    // There are two builders:
    // 1) Rewighting from MELA (x) K factors, which adjust the cross section
    // 2) PU and GenHepMCWeight reweighting, which are supposed to keep the cross section the same
    // Total weight is (1)x(2)

    // Build the possible MELA reweightings
    TString strWeight = getMELAHypothesisWeight(t);
    strReweightingWeigths.clear();
    strReweightingWeigths.push_back(strWeight);
    for (auto& s:strKfactorVars) strReweightingWeigths.push_back(s);
    strReweightingWeigths.push_back("xsec");
    ReweightingBuilder* melarewgtBuilder = new ReweightingBuilder(strReweightingWeigths, getSimpleWeight);
    melarewgtBuilder->rejectNegativeWeights(true);
    melarewgtBuilder->setDivideByNSample(true);
    melarewgtBuilder->setWeightBinning(GenHMassBinning);
    for (auto& theSampleSet:theSets){ for (auto& tree:theSampleSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(tree, 0.999, 250); }

    controlsDir->cd();
    for (unsigned int bin=0; bin<GenHMassBinning.getNbins(); bin++){
      unsigned int nSANZWEPB = melarewgtBuilder->getSumAllNonZeroWgtEvents(bin);
      float NCPB = melarewgtBuilder->getNormComponent(bin);
      MELAout
        << "GenHMass bin " << bin << " has a normalization contribution in melarewgtBuilder of " << NCPB
        << " and number of events with non-zero weight: " << nSANZWEPB
        << endl;

      h_MELARewgtSumAllNonZeroWgtEventsPerBin->SetBinContent(bin, nSANZWEPB);
      h_MELARewgtNormComponentPerBin->SetBinContent(bin, NCPB);
    }
    MELAout << "Overall normalization in melarewgtBuilder is " << melarewgtBuilder->getNorm() << endl;

    controlsDir->WriteTObject(h_MELARewgtSumAllNonZeroWgtEventsPerBin); delete h_MELARewgtSumAllNonZeroWgtEventsPerBin;
    controlsDir->WriteTObject(h_MELARewgtNormComponentPerBin); delete h_MELARewgtNormComponentPerBin;

    foutput->cd();

    std::vector<SimpleEntry> products;
    // Build the analyzer and loop over the events
    for (auto& theSampleSet:theSets){
      GGAnalyzer theAnalyzer(theSampleSet, channel, category);
      theAnalyzer.setExternalProductList(&products);
      // Book common variables needed for analysis
      theAnalyzer.addConsumed<float>("PUWeight");
      theAnalyzer.addConsumed<float>("genHEPMCweight");
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
      // Loop
      theAnalyzer.loop(true, false, true);
    }
    MELAout << "There are " << products.size() << " products" << endl;
    SimpleEntry::writeToTree(products.cbegin(), products.cend(), theFinalTree);

    delete melarewgtBuilder;

    controlsDir->Close();
    foutput->WriteTObject(theFinalTree);
    delete theFinalTree;
  }

  delete KD2;
  delete KD1;
  for (auto& theSampleSet:theSets) delete theSampleSet; theSets.clear();
  foutput->Close();
  MELAout.close();
}

void makeGGTemplatesFromMCFM_two(const Channel channel, const Category category, TString strSystematics){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";

  TString INPUT_NAME = Form("%s_ggHtoZZ%s_MCFM_Stage1_%s", strCategory.Data(), strChannel.Data(), strSystematics.Data());
  INPUT_NAME += ".root";
  TString cinput = coutput_common + INPUT_NAME;
  if (gSystem->AccessPathName(cinput)) return;

  gSystem->Exec("mkdir -p " + coutput_common);
  TString OUTPUT_NAME = Form("%s_ggHtoZZ%s_MCFM_Stage2_%s", strCategory.Data(), strChannel.Data(), strSystematics.Data());
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

  HelperFunctions::CopyFile(cinput, fixTreeWeights, nullptr);
  foutput->ls();

  foutput->Close();
  MELAout.close();
}

TTree* fixTreeWeights(TTree* tree){
  float ZZMass, weight;
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("weight", &weight);

  TTree* newtree = tree->CloneTree(0);
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);

    newtree->Fill();
  }
  return newtree;
}

