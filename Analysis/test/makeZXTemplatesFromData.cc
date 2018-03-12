#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"
#include "CheckSetTemplatesCategoryScheme.h"


// Process handle
typedef ZXProcessHandler ProcessHandleType;
const ProcessHandleType& theProcess = TemplateHelpers::OffshellZXProcessHandle;

// Process-specific functions
void makeZXTemplatesFromData_one(const Channel channel, const Category category, const SystematicVariationTypes syst, const ZXFakeRateHandler::FakeRateMethod FRMethod, const TString fixedDate="");
void makeZXTemplatesFromData_checkstage(const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst, const unsigned int istage, const TString fixedDate="");

// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif
#ifndef checkstage_def
#define checkstage_def
typedef void(*CheckStageFcn)(const Channel, const Category, const ACHypothesis, const SystematicVariationTypes, const unsigned int, const TString);
CheckStageFcn checkstagefcn = &makeZXTemplatesFromData_checkstage;
#endif

void plotProcessCheckStage(const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst, const unsigned int istage, const TString fixedDate="", ProcessHandler::ProcessType proctype=theProcess.getProcessType(), const TString strGenerator="Data");
void plotProcessCheckStage_SystPairs(const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst, const unsigned int istage, const TString fixedDate="", ProcessHandler::ProcessType proctype=theProcess.getProcessType(), const TString strGenerator="Data");

// Function to build one templates
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// theSqrts = 13 (CoM energy) is fixed in Samples.h
void makeZXTemplatesFromData_one(const Channel channel, const Category category, const SystematicVariationTypes syst, const ZXFakeRateHandler::FakeRateMethod FRMethod, const TString fixedDate){
  typedef TGraphErrors SSFRtype;

  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  if (!systematicAllowed(category, channel, theProcess.getProcessType(), syst)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strSystematics = getSystematicsName(syst);
  const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  MELAout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  coutput_common+="Stage1/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Data",
    strChannel.Data(), strCategory.Data(),
    theProcess.getProcessName().Data(), FRMethodName.Data(),
    strSystematics.Data()
  );
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

  // Get ZX fake rate estimator
  ZXFakeRateHandler* ZXFRHandler = getFakeRateHandler(FRMethod, syst);

  // Get the required discriminants
  vector<KDspecs> KDlist;
  getLikelihoodDiscriminants(channel, category, syst, KDlist);
  if (category!=Inclusive) getCategorizationDiscriminants(syst, KDlist);

  // Get list of samples
  vector<TString> strSamples;
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("AllData");
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples, syst);

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples, TREE_CRZLL_NAME, "", COUNTERS_CRZLL_NAME);
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book ZX variables
    ZXFRHandler->registerTree(tree);
    // Book common variables needed for analysis
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for KDs
    for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
    tree->silenceUnused(); // Will no longer book another branch
  }
  theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_None, false, false);

  std::vector<ReweightingBuilder*> extraEvaluators;
  SystematicsClass* systhandle = constructSystematic(category, channel, theProcess.getProcessType(), syst, theSampleSet->getCJLSTTreeList(), extraEvaluators);

  ExtendedBinning ZZMassInclusiveBinning("ZZMass");

  // Construct reweighting variables vector
  for (int t=ProcessHandleType::ZX; t<(int) ProcessHandleType::nZXTypes; t++){
    foutput->cd();

    TString treename = theProcess.getOutputTreeName();
    BaseTree* theFinalTree = new BaseTree(treename); // The tree to record into the ROOT file

    // Build the analyzer and loop over the events
    TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
    theAnalyzer.setExternalProductTree(theFinalTree);
    // Loosen channel check if needed
    theAnalyzer.setAllowSSChannel((FRMethod==ZXFakeRateHandler::mSS));
    // Book common variables needed for analysis
    theAnalyzer.addConsumed<float>("ZZMass");
    theAnalyzer.addConsumed<short>("Z1Flav");
    theAnalyzer.addConsumed<short>("Z2Flav");
    // Add discriminant builders
    for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
    // Add ZX FR handle
    theAnalyzer.addZXFakeRateHandler(Form("ZXFR%s", FRMethodName.Data()), ZXFRHandler);
    // Add systematics handle
    theAnalyzer.addSystematic(strSystematics, systhandle);
    // Loop
    theAnalyzer.loop(true, false, true);

    MELAout << "There are " << theFinalTree->getNEvents() << " products" << endl;
    theFinalTree->writeToFile(foutput);
    delete theFinalTree;
  }

  delete systhandle;
  for (auto& rb:extraEvaluators) delete rb;
  for (auto& KD:KDlist) delete KD.KD;
  delete theSampleSet;
  delete ZXFRHandler;
  foutput->Close();
  MELAout.close();
}

void makeZXTemplatesFromData_checkstage(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage,
  const TString fixedDate
){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  if (!systematicAllowed(category, channel, theProcess.getProcessType(), syst)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strSystematics = getSystematicsName(syst);
  const TString strStage = Form("Stage%i", istage);
  const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(ZXFakeRateHandler::mSS);
  std::vector<ProcessHandleType::HypothesisType> tplset; tplset.push_back(ProcessHandleType::ZX);
  const unsigned int ntpls = tplset.size();

  // Get the KDs needed for the AC hypothesis
  vector<TString> KDset = ACHypothesisHelpers::getACHypothesisKDNameSet(hypo, category);
  const unsigned int nKDs = KDset.size();
  unordered_map<TString, float> KDvars;
  for (auto& KDname:KDset) KDvars[KDname]=0;

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Check_" + strStage + "/";

  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Data.root",
    strChannel.Data(), strCategory.Data(),
    theProcess.getProcessName().Data(), FRMethodName.Data(),
    strSystematics.Data()
  );
  TString cinput = cinput_common + INPUT_NAME;
  if (gSystem->AccessPathName(cinput)){
    MELAerr << "File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
    return;
  }
  TFile* finput = TFile::Open(cinput, "read");

  gSystem->Exec("mkdir -p " + coutput_common);
  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_Data_Check%sDiscriminants",
    strChannel.Data(), strCategory.Data(),
    theProcess.getProcessName().Data(),
    strSystematics.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";
  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;
  gStyle->SetOptStat(0);
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  const int infMass=70;
  const int offshellMassBegin=220;
  const int offshellMassWidth=20;
  const int supMass=1000;
  vector<TH2F*> finalTemplates_2D; finalTemplates_2D.assign(ntpls, nullptr);
  vector<TH3F*> finalTemplates_3D; finalTemplates_3D.assign(ntpls, nullptr);
  vector<TH1F*> htpl_1D; htpl_1D.assign(ntpls, nullptr);
  vector<vector<TH2F*>> htpl_2D; htpl_2D.assign(ntpls, vector<TH2F*>());
  ExtendedBinning binning_KDpure(30, 0, 1, "KDpure");
  ExtendedBinning binning_KDint(30, -1, 1, "KDint");
  ExtendedBinning binning_mass_offshell((supMass-offshellMassBegin)/offshellMassWidth, offshellMassBegin, supMass, "m_{4l}");
  ExtendedBinning binning_mass("m_{4l}");
  binning_mass.addBinBoundary(70);
  for (unsigned int bin=0; bin<=(140-105)/1; bin++) binning_mass.addBinBoundary(105 + bin*1);
  binning_mass.addBinBoundary(180);
  for (unsigned int bin=0; bin<=(supMass-offshellMassBegin)/offshellMassWidth; bin++) binning_mass.addBinBoundary(offshellMassBegin + bin*offshellMassWidth);
  for (unsigned int t=0; t<ntpls; t++){
    TString templatename = theProcess.getTemplateName();
    TString treename = theProcess.getOutputTreeName();
    MELAout << "Setting up tree " << treename << " and template " << templatename << endl;

    finput->cd();
    TTree* tree = (TTree*)finput->Get(treename);
    foutput->cd();

    bool isCategory=(category==Inclusive);
    float ZZMass, weight;
    tree->SetBranchAddress("ZZMass", &ZZMass);
    tree->SetBranchAddress("weight", &weight);
    for(auto const& KDname:KDset){
      MELAout << "Setting up KD " << KDname << " tree variable" << endl;
      tree->SetBranchAddress(KDname, &(KDvars[KDname]));
    }
    if (!isCategory){
      TString catFlagName = TString("is_") + strCategory + TString("_") + getACHypothesisName(hypo);
      tree->SetBranchAddress(catFlagName, &isCategory);
    }

    htpl_1D[t] = new TH1F(Form("h1D_%s", treename.Data()), templatename, binning_mass.getNbins(), binning_mass.getBinning());
    htpl_1D[t]->GetXaxis()->SetTitle(Form("%s (GeV)", binning_mass.getLabel().Data()));
    htpl_1D[t]->SetOption("hist");
    for (unsigned int iKD=0; iKD<nKDs; iKD++){
      const TString& KDname = KDset.at(iKD);
      MELAout << "Setting up KD " << KDname << " histograms" << endl;
      TString hname = Form("h2D_%s_ZZMass_vs_%s", treename.Data(), KDname.Data());
      TH2F* htmp = new TH2F(
        hname, templatename,
        binning_mass.getNbins(), binning_mass.getBinning(),
        (KDname.Contains("int") ? binning_KDint.getNbins() : binning_KDpure.getNbins()),
        (KDname.Contains("int") ? binning_KDint.getBinning() : binning_KDpure.getBinning())
      );
      htmp->GetXaxis()->SetTitle(Form("%s (GeV)", binning_mass.getLabel().Data()));
      htmp->GetYaxis()->SetTitle(KDname);
      htmp->SetOption("colz");
      htpl_2D[t].push_back(htmp);
    }
    if (nKDs==1) finalTemplates_2D[t] = new TH2F(
      templatename, templatename,
      binning_mass_offshell.getNbins(), binning_mass_offshell.getBinning(),
      (KDset.at(0).Contains("int") ? binning_KDint.getNbins() : binning_KDpure.getNbins()),
      (KDset.at(0).Contains("int") ? binning_KDint.getBinning() : binning_KDpure.getBinning())
    );
    else finalTemplates_3D[t] = new TH3F(
      templatename, templatename,
      binning_mass_offshell.getNbins(), binning_mass_offshell.getBinning(),
      (KDset.at(0).Contains("int") ? binning_KDint.getNbins() : binning_KDpure.getNbins()),
      (KDset.at(0).Contains("int") ? binning_KDint.getBinning() : binning_KDpure.getBinning()),
      (KDset.at(1).Contains("int") ? binning_KDint.getNbins() : binning_KDpure.getNbins()),
      (KDset.at(1).Contains("int") ? binning_KDint.getBinning() : binning_KDpure.getBinning())
    );

    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);

      if (!isCategory) continue;

      htpl_1D[t]->Fill(ZZMass, weight);
      
      for (auto& KD:KDvars){ if (KD.second==1.) KD.second -= 0.001*float(ev)/float(tree->GetEntries()); }
      unsigned int iKD=0;
      for (auto& KDname:KDset){
        htpl_2D[t].at(iKD)->Fill(ZZMass, KDvars[KDname], weight);
        iKD++;
      }
      if (nKDs==1) finalTemplates_2D[t]->Fill(ZZMass, KDvars[KDset.at(0)], weight);
      else finalTemplates_3D[t]->Fill(ZZMass, KDvars[KDset.at(0)], KDvars[KDset.at(1)], weight);
    }
  }

  MELAout << "Extracting the 1D distributions of various components" << endl;
  theProcess.recombineHistogramsToTemplates(htpl_1D);
  MELAout << "Extracting the 2/3D templates" << endl;
  if (nKDs==1) theProcess.recombineHistogramsToTemplates(finalTemplates_2D);
  else theProcess.recombineHistogramsToTemplates(finalTemplates_3D);
  MELAout << "Extracting the 2D distributions of various components" << endl;
  for (unsigned int iKD=0; iKD<nKDs; iKD++){
    vector<TH2F*> htmp;
    for (unsigned int t=0; t<ntpls; t++) htmp.push_back(htpl_2D[t].at(iKD));
    theProcess.recombineHistogramsToTemplates(htmp);
  }
  MELAout << "Extracted all components" << endl;
  for (unsigned int t=0; t<ntpls; t++){
    if (nKDs==1) MELAout << "Template " << finalTemplates_2D[t]->GetName() << " integral: " << HelperFunctions::computeIntegral(finalTemplates_2D[t], true) << endl;
    else MELAout << "Template " << finalTemplates_3D[t]->GetName() << " integral: " << HelperFunctions::computeIntegral(finalTemplates_3D[t], true) << endl;
  }

  for (unsigned int iKD=0; iKD<nKDs; iKD++){
    vector<TH2F*> htmp; htmp.reserve(ntpls);
    for (unsigned int t=0; t<ntpls; t++) htmp.push_back(htpl_2D[t].at(iKD));
    theProcess.conditionalizeTemplates(htmp, hypo, 0);
  }
  for (int t=ProcessHandleType::ZXTpl; t<(int) ProcessHandleType::nZXTplTypes; t++){
    foutput->WriteTObject(htpl_1D[t]);
    delete htpl_1D[t];
    for (auto& htmp:htpl_2D[t]){
      foutput->WriteTObject(htmp);
      delete htmp;
    }
    if (finalTemplates_2D[t]){
      foutput->WriteTObject(finalTemplates_2D[t]);
      delete finalTemplates_2D[t];
    }
    else{
      foutput->WriteTObject(finalTemplates_3D[t]);
      delete finalTemplates_3D[t];
    }
  }
  foutput->Close();
  finput->Close();
  MELAout.close();
}

#include "plotProcessCheckStage.cc"

