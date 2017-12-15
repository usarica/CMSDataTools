#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"


// Constants to affect the template code
const TString user_output_dir = "output/";

TTree* fixTreeWeights(TTree* tree);

// Function to build one templates
// ichan = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// theSqrts = 13 (CoM energy) is fixed in Samples.h
void makeVVTemplatesFromPOWHEG_one(const Channel channel, const Category category, ACHypothesis hypo, TString strSystematics, const TString fixedDate=""){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const std::vector<TemplateHelpers::VVHypothesisType> tplset = getVVHypothesesForACHypothesis(hypo);
  std::vector<TString> melawgtvars; for (auto& hypotype:tplset) melawgtvars.push_back(getMELAVVHypothesisWeight(hypotype, hypo));

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_POWHEG_Stage1",
    strChannel.Data(), strCategory.Data(),
    getACHypothesisName(hypo).Data(), getVVProcessName(true).Data(),
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

  // Get list of samples
  vector<TString> strSampleIdentifiers;
  strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
  vector<TString> strSamples;
  getSamplesList(theSqrts, strSampleIdentifiers, strSamples);

  // Kfactor variable names
  vector<TString> strKfactorVars;
  strKfactorVars.push_back("p_Gen_CPStoBWPropRewgt");
  /*
  if (strSystematics == "Nominal") strKfactorVars.push_back("KFactor_QCD_ggZZ_Nominal");
  else if (strSystematics == "PDFScaleDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFScaleDn");
  else if (strSystematics == "PDFScaleUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFScaleUp");
  else if (strSystematics == "QCDScaleDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_QCDScaleDn");
  else if (strSystematics == "QCDScaleUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_QCDScaleUp");
  else if (strSystematics == "AsDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_AsDn");
  else if (strSystematics == "AsUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_AsUp");
  else if (strSystematics == "PDFReplicaDn") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFReplicaDn");
  else if (strSystematics == "PDFReplicaUp") strKfactorVars.push_back("KFactor_QCD_ggZZ_PDFReplicaUp");
  */

  // Register the discriminants
  vector<KDspecs> KDlist;
  getLikelihoodDiscriminants(channel, category, strSystematics, KDlist);
  if (category!=Inclusive) getCategorizationDiscriminants(strSystematics, KDlist);

  // Get the CJLST set
  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);
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
    // Variables for MELA reweighting
    for (TString const& wgtvar:melawgtvars) tree->bookBranch<float>(wgtvar, 0);
    // Variables for KDs
    for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
    tree->silenceUnused(); // Will no longer book another branch
  }
  theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_NgenOverNgenWPU, false, true);

  // Setup GenHMass binning
  // Binning for inclusive reweighting
  ExtendedBinning GenHMassInclusiveBinning("GenHMass");
  // Binning for MELARewgt
  ExtendedBinning GenHMassBinning("GenHMass");
  for (unsigned int is=0; is<theSampleSet->getCJLSTTreeList().size()-1; is++){
    if (theSampleSet->getCJLSTTreeList().at(is)->MHVal>0. && theSampleSet->getCJLSTTreeList().at(is+1)->MHVal>0.){
      float boundary = (theSampleSet->getCJLSTTreeList().at(is)->MHVal + theSampleSet->getCJLSTTreeList().at(is+1)->MHVal)/2.;
      GenHMassBinning.addBinBoundary(boundary);
    }
  }
  GenHMassBinning.addBinBoundary(0);
  GenHMassBinning.addBinBoundary(theSqrts*1000.);

  // Construct reweighting variables vector
  for (unsigned int t=0; t<tplset.size(); t++){
    auto& hypotype = tplset.at(t);
    foutput->cd();

    /************* Reweighting setup *************/
    // There are two builders:
    // 1) Rewighting from MELA (x) K factors, which adjust the cross section
    // 2) PU and GenHepMCWeight reweighting, which are supposed to keep the cross section the same
    // Total weight is (1)x(2)
    vector<TString> strReweightingWeights;
    strReweightingWeights.push_back(melawgtvars.at(t));
    for (auto& s:strKfactorVars) strReweightingWeights.push_back(s);
    strReweightingWeights.push_back("xsec");

    TString treename = getVVOutputTreeName(hypotype, true);
    BaseTree* theFinalTree = new BaseTree(treename); // The tree to record into the ROOT file

    // Build the analyzer and loop over the events
      // Build the MELA reweightings separately for each POWHEG 2f2f' sample set since NormComponent combines xsec
    ReweightingBuilder* melarewgtBuilder = new ReweightingBuilder(strReweightingWeights, getSimpleWeight);
    melarewgtBuilder->rejectNegativeWeights(true);
    melarewgtBuilder->setDivideByNSample(true);
    melarewgtBuilder->setWeightBinning(GenHMassBinning);
    for (auto& tree:theSampleSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(tree, 0.999, 250);

    // Make reweighting conrol plots
    TDirectory* controlsDir = foutput->mkdir(Form("controls_%s", treename.Data()), "");
    controlsDir->cd();
    TH1F* h_MELARewgtSumAllNonZeroWgtEventsPerBin = new TH1F("MELARewgtSumAllNonZeroWgtEventsPerBin", "", GenHMassBinning.getNbins(), GenHMassBinning.getBinning());
    h_MELARewgtSumAllNonZeroWgtEventsPerBin->GetXaxis()->SetTitle(GenHMassBinning.getLabel());
    TH1F* h_MELARewgtNormComponentPerBin = new TH1F("MELARewgtNormComponentPerBin", "", GenHMassBinning.getNbins(), GenHMassBinning.getBinning());
    h_MELARewgtNormComponentPerBin->GetXaxis()->SetTitle(GenHMassBinning.getLabel());
    for (unsigned int bin=0; bin<GenHMassBinning.getNbins(); bin++){
      unsigned int nSANZWEPB = melarewgtBuilder->getSumAllNonZeroWgtEvents(bin); // sum{t} N_tj
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

    controlsDir->Close();
    foutput->cd();

    TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
    theAnalyzer.setExternalProductTree(theFinalTree);
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
    for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
    // Add reweighting builders
    theAnalyzer.addReweightingBuilder("MELARewgt", melarewgtBuilder);
    // Loop
    theAnalyzer.loop(true, false, true);

    delete melarewgtBuilder;

    MELAout << "There are " << theFinalTree->getNEvents() << " products" << endl;
    theFinalTree->writeToFile(foutput);
    delete theFinalTree;
  }

  for (auto& KD:KDlist) delete KD.KD;
  delete theSampleSet;
  foutput->Close();
  MELAout.close();
}

void makeVVTemplatesFromPOWHEG_two(const Channel channel, const Category category, ACHypothesis hypo, TString strSystematics, const TString fixedDate=""){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";

  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_POWHEG_Stage1",
    strChannel.Data(), strCategory.Data(),
    getACHypothesisName(hypo).Data(), getVVProcessName(true).Data(),
    strSystematics.Data()
  );
  INPUT_NAME += ".root";
  TString cinput = coutput_common + INPUT_NAME;
  if (gSystem->AccessPathName(cinput)) makeVVTemplatesFromPOWHEG_one(channel, category, hypo, strSystematics, fixedDate);

  gSystem->Exec("mkdir -p " + coutput_common);
  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_POWHEG_Stage2",
    strChannel.Data(), strCategory.Data(),
    getACHypothesisName(hypo).Data(), getVVProcessName(true).Data(),
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

  HelperFunctions::CopyFile(cinput, fixTreeWeights, nullptr);
  foutput->ls();
  foutput->Close();

  MELAout << "===============================" << endl;
  MELAout << "Stage 1 file " << cinput << " is no longer needed. Removing the file..." << endl;
  gSystem->Exec("rm -f " + cinput);
  MELAout << "===============================" << endl;

  MELAout.close();
}

void makeVVTemplatesFromPOWHEG_checkstage(
  const Channel channel, const Category category, ACHypothesisHelpers::ACHypothesis hypo, TString strSystematics,
  const unsigned int istage,
  const TString fixedDate=""
){
  if (channel==NChannels) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  std::vector<TemplateHelpers::VVHypothesisType> tplset = getVVHypothesesForACHypothesis(kSM);
  if (hypo!=kSM){
    std::vector<TemplateHelpers::VVHypothesisType> tplset_tmp = getVVHypothesesForACHypothesis(hypo);
    for (TemplateHelpers::VVHypothesisType& v:tplset_tmp) tplset.push_back(v);
  }
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
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";

  vector<TFile*> finputList;
  for (unsigned int f=0; f<2; f++){
    if (f==1 && hypo==kSM) break;
    ACHypothesis fhypo = (f==1 ? hypo : kSM);
    TString INPUT_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_POWHEG_Stage%i",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(fhypo).Data(), getVVProcessName(true).Data(),
      strSystematics.Data(), istage
    );
    INPUT_NAME += ".root";
    TString cinput = coutput_common + INPUT_NAME;
    if (gSystem->AccessPathName(cinput)){
      if (istage==1) makeVVTemplatesFromPOWHEG_one(channel, category, fhypo, strSystematics, fixedDate);
      else if (istage==2) makeVVTemplatesFromPOWHEG_two(channel, category, fhypo, strSystematics, fixedDate);
      else return;
    }
    TFile* ftmp = TFile::Open(cinput, "read");
    finputList.push_back(ftmp);
  }

  gSystem->Exec("mkdir -p " + coutput_common);
  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_POWHEG_Check%sDiscriminants_Stage%i",
    strChannel.Data(), strCategory.Data(),
    getVVProcessName(true).Data(),
    strSystematics.Data(),
    getACHypothesisName(hypo).Data(), istage
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
    VVHypothesisType const& treetype = tplset.at(t);
    VVTemplateType tpltype = castIntToVVTemplateType(castVVHypothesisTypeToInt(treetype));
    TString templatename = getVVTemplateName(tpltype, true);
    TString treename = getVVOutputTreeName(treetype, true);
    MELAout << "Setting up tree " << treename << " and template " << templatename << endl;

    TFile*& finput = finputList.at(castVVHypothesisTypeToInt(treetype)>=castVVHypothesisTypeToInt(nVVSMTypes));
    finput->cd();
    TTree* tree = (TTree*) finput->Get(treename);
    foutput->cd();

    bool isCategory=(category==Inclusive);
    float ZZMass, weight;
    tree->SetBranchAddress("ZZMass", &ZZMass);
    tree->SetBranchAddress("weight", &weight);
    for (auto const& KDname:KDset){
      MELAout << "Setting up KD " << KDname << " tree variable" << endl;
      tree->SetBranchAddress(KDname, &(KDvars[KDname]));
    }
    if (!isCategory){
      TString catFlagName = TString("is_") + strCategory + TString("_") + getACHypothesisName(hypo);
      tree->SetBranchAddress(catFlagName, &isCategory);
    }

    htpl_1D[t] = new TH1F(Form("h1D_%s", templatename.Data()), templatename, binning_mass.getNbins(), binning_mass.getBinning());
    htpl_1D[t]->GetXaxis()->SetTitle(Form("%s (GeV)", binning_mass.getLabel().Data()));
    htpl_1D[t]->SetOption("hist");
    for (unsigned int iKD=0; iKD<nKDs; iKD++){
      const TString& KDname = KDset.at(iKD);
      MELAout << "Setting up KD " << KDname << " histograms" << endl;
      TString hname = Form("h2D_%s_ZZMass_vs_%s", templatename.Data(), KDname.Data());
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
  recombineVVHistogramsToTemplates(htpl_1D, hypo);
  MELAout << "Extracting the 2/3D templates" << endl;
  if (nKDs==1) recombineVVHistogramsToTemplates(finalTemplates_2D, hypo);
  else recombineVVHistogramsToTemplates(finalTemplates_3D, hypo);
  MELAout << "Extracting the 2D distributions of various components" << endl;
  for (unsigned int iKD=0;iKD<nKDs;iKD++){
    vector<TH2F*> htmp;
    for (unsigned int t=0; t<ntpls; t++) htmp.push_back(htpl_2D[t].at(iKD));
    recombineVVHistogramsToTemplates(htmp, hypo);
  }
  MELAout << "Extracted all components" << endl;
  for (unsigned int t=0; t<ntpls; t++){
    if (nKDs==1) MELAout << "Template " << finalTemplates_2D[t]->GetName() << " integral: " << HelperFunctions::computeIntegral(finalTemplates_2D[t], true) << endl;
    else MELAout << "Template " << finalTemplates_3D[t]->GetName() << " integral: " << HelperFunctions::computeIntegral(finalTemplates_3D[t], true) << endl;
  }

  for (unsigned int iKD=0; iKD<nKDs; iKD++){ for (unsigned int t=0; t<ntpls; t++) HelperFunctions::conditionalizeHistogram(htpl_2D[t].at(iKD), 0); }
  for (unsigned int t=0; t<ntpls; t++){
    foutput->WriteTObject(htpl_1D[t]);
    delete htpl_1D[t];
    for (auto& htmp:htpl_2D[t]){
      foutput->WriteTObject(htmp);
      delete htmp;
    }
    if (nKDs==1){
      foutput->WriteTObject(finalTemplates_2D[t]);
      delete finalTemplates_2D[t];
    }
    else{
      foutput->WriteTObject(finalTemplates_3D[t]);
      delete finalTemplates_3D[t];
    }
  }
  foutput->Close();
  for (TFile*& finput:finputList)finput->Close();
  MELAout.close();
}

TTree* fixTreeWeights(TTree* tree){
  const unsigned int nMarginalMax = 100;
  const unsigned int nMarginalMaxMult = 1000;
  const float nMarginalMaxFrac = 1./static_cast<float const>(nMarginalMaxMult);
  const unsigned int countThreshold=nMarginalMaxMult*nMarginalMax;

  const TString treename=tree->GetName();
  float ZZMass, weight;
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("weight", &weight);

  int const nbinsraw = (1000*theSqrts-70)/5;
  TH1F* hmass = new TH1F("hmass", "", nbinsraw, 70, 13000);
  TTree* newtree = tree->CloneTree(0);
  // Initial loop over the tree to count the events in each bin
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    hmass->Fill(ZZMass);
  }

  // Determine the final binning to set the weight thresholds
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Determining the final binning to set the weight thresholds"
    << endl;
  ExtendedBinning binning;
  binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(hmass->GetNbinsX()+1));
  vector<unsigned int> counts;
  unsigned int count=0;
  for (int bin=hmass->GetNbinsX(); bin>=0; bin--){
    count += hmass->GetBinContent(bin);
    if (count>countThreshold || bin==0){
      counts.push_back(count);
      binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(bin));
      count=0;
    }
  }
  delete hmass;
  std::reverse(counts.begin(), counts.end());
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;
  // These lines guarantee count>countThreshold in every bin
  if (counts.at(0)<countThreshold){
    counts.at(1) += counts.at(0);
    counts.erase(counts.begin());
    binning.removeBinLowEdge(1);
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;

  // Collect the count*nMarginalMaxFrac events with highest weights
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Collecting the count*" << nMarginalMaxFrac << " events with highest weights in " << binning.getNbins() << " bins"
    << endl;
  vector<vector<float>> wgtcollList;
  wgtcollList.assign(binning.getNbins(), vector<float>());
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int)binning.getNbins()){
      vector<float>& wgtcoll=wgtcollList.at(bin);
      const unsigned int maxPrunedSize = std::ceil(float(counts.at(bin))*nMarginalMaxFrac);
      if (wgtcoll.size()<maxPrunedSize) addByHighest(wgtcoll, fabs(weight), false);
      else if (wgtcoll.back()<fabs(weight)){
        addByHighest(wgtcoll, fabs(weight), false);
        wgtcoll.pop_back();
      }
    }
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Determining the weight thresholds"
    << endl;
  vector<float> wgtThresholds; wgtThresholds.reserve(binning.getNbins());
  for (auto const& wgtcoll:wgtcollList){
    unsigned int ns=wgtcoll.size();
    float threshold=0.5*(wgtcoll.at(ns-1)+wgtcoll.at(ns-2));
    if (wgtcoll.front()*5.<threshold) threshold=wgtcoll.front();
    else MELAout
      << "fixTreeWeights(" << treename << "): "
      << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
      << endl;
    wgtThresholds.push_back(threshold);
  }
  
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int)binning.getNbins() && fabs(weight)>wgtThresholds.at(bin)) weight = pow(wgtThresholds.at(bin), 2)/weight;
    newtree->Fill();
  }
  return newtree;
}

