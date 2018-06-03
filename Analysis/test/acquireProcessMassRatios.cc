#ifndef ACQUIREPROCESSMASSRATIOS_H
#define ACQUIREPROCESSMASSRATIOS_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "TemplatesEventAnalyzer.h"
#include "HistogramSmootherWithGaussianKernel.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

using namespace HistogramSmootherWithGaussianKernel;


ExtendedBinning getMassBinning(TTree* tree, float& ZZMass, bool& flag, bool separateZ4l=false);
void acquireMassRatio_ProcessSystToNominal_PythiaMINLO_one(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler::ProcessType proctype, const TString strGenerator
);

void acquireMassRatio_ProcessNominalToNominalInclusive_one(
  const Channel channel, const Category category, const ACHypothesis hypo,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler::ProcessType proctype, const TString strGenerator
){
  const SystematicVariationTypes syst=sNominal;

  if (channel==NChannels) return;
  if (category==Inclusive) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst, strGenerator)) return;
  if (proctype==ProcessHandler::kZX && strGenerator!="Data") return;

  TDirectory* rootdir = gDirectory;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strCategory_Inclusive = getCategoryName(Inclusive);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/";

  gSystem->Exec("mkdir -p " + coutput_common);
  vector<TString> INPUT_BASELINE_NAME;
  vector<TString> INPUT_NAME;
  if (proctype==ProcessHandler::kQQBkg){
    INPUT_BASELINE_NAME.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory_Inclusive.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    INPUT_NAME.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
  }
  else if (proctype==ProcessHandler::kZX){
    vector<ZXFakeRateHandler::FakeRateMethod> FRMethods;
    FRMethods.push_back(ZXFakeRateHandler::mSS);
    for (ZXFakeRateHandler::FakeRateMethod& FRMethod:FRMethods){
      const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod);
      INPUT_BASELINE_NAME.push_back(Form(
        "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
        strChannel.Data(), strCategory_Inclusive.Data(),
        thePerProcessHandle->getProcessName().Data(), FRMethodName.Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
      INPUT_NAME.push_back(Form(
        "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        thePerProcessHandle->getProcessName().Data(), FRMethodName.Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
    }
  }
  else{
    INPUT_BASELINE_NAME.push_back(Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory_Inclusive.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    INPUT_NAME.push_back(Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    if (hypo!=kSM){
      INPUT_BASELINE_NAME.push_back(Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory_Inclusive.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
      INPUT_NAME.push_back(Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
    }
  }

  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_MassRatiosToNominalInclusive_%s_%s_%s",
    strChannel.Data(), strCategory.Data(),
    strACHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  // Open the input files
  const unsigned int InputTypeSize=2; // Inclusive or dedicated category
  vector<vector<TString>> cinputList; cinputList.assign(InputTypeSize, vector<TString>());
  vector<vector<TFile*>> finputList; finputList.assign(InputTypeSize, vector<TFile*>());
  vector<vector<TTree*>> treeLists; treeLists.assign(InputTypeSize, vector<TTree*>());

  for (auto& s:INPUT_BASELINE_NAME) cinputList.at(0).push_back(cinput_common + s);
  for (auto& s:INPUT_NAME) cinputList.at(1).push_back(cinput_common + s);
  assert(cinputList.at(0).size()==cinputList.at(1).size());

  for (unsigned int it=0; it<InputTypeSize; it++){
    for (auto const& cinput:cinputList.at(it)){
      if (gSystem->AccessPathName(cinput)){
        MELAerr << "File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
        for (unsigned int it2=0; it2<InputTypeSize; it2++){ for (auto& finput:finputList.at(it2)) finput->Close(); }
        return;
      }
      if (cinput!="") finputList.at(it).push_back(TFile::Open(cinput, "read"));
    }
  }

  unsigned int ntreespertype=0;
  for (unsigned int it=0; it<InputTypeSize; it++){
    // Extract SM or AC trees
    for (auto& finput:finputList.at(it)){ if (finput) extractTreesFromDirectory(finput, treeLists.at(it)); }
    // Check number of trees
    unsigned int ntrees = treeLists.at(it).size();
    if (ntreespertype==0) ntreespertype=ntrees;
    else if (ntreespertype!=ntrees){
      MELAerr << "Number of trees in type " << it << " (" << ntrees << ") does not match Ntrees=" << ntreespertype << "! Aborting..." << endl;
      for (unsigned int it2=0; it2<InputTypeSize; it2++){ for (auto& finput2:finputList.at(it2)) finput2->Close(); }
      return;
    }
  }
  MELAout << "Acquired " << ntreespertype << " trees per type..." << endl;

  vector<vector<TTree*>> treeGroups;
  treeGroups.assign(ntreespertype, vector<TTree*>());
  for (auto& treeGroup:treeGroups) treeGroup.assign(InputTypeSize, nullptr);
  for (unsigned int it=0; it<InputTypeSize; it++){
    for (unsigned int ihypo=0; ihypo<ntreespertype; ihypo++){
      treeGroups.at(ihypo).at(it) = treeLists.at(it).at(ihypo);
    }
  }

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;

  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened output file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  rootdir->cd();
  vector<ExtendedHistogram_1D> hRatio; hRatio.reserve(treeGroups.size());
  for (auto& treeGroup:treeGroups){
    int nEntriesMin=-1;
    TTree* treeChosen=nullptr;
    for (TTree*& tree:treeGroup){
      MELAout << "Tree has " << tree->GetEntries() << " entries" << endl;
      if (tree->GetEntries()<=nEntriesMin || nEntriesMin<0){
        treeChosen=tree;
        nEntriesMin=tree->GetEntries();
      }
    }

    bool isCategory=false;
    float ZZMass, weight;
    for (unsigned int it=0; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      isCategory=(it==0);
      tree->SetBranchStatus("*", 0);
      bookBranch(tree, "weight", &weight);
      bookBranch(tree, "ZZMass", &ZZMass);
      if (!isCategory){
        TString catFlagName = TString("is_") + strCategory + TString("_") + getACHypothesisName(hypo);
        bookBranch(tree, catFlagName, &isCategory);
      }
    }
    isCategory=true;
    ExtendedBinning binning("ZZMass");
    if (proctype==ProcessHandler::kTT || proctype==ProcessHandler::kBB || strGenerator=="JHUGen"){ // On-shell - only samples
      binning.addBinBoundary(70.);
      binning.addBinBoundary(theSqrts*1000.);
    }
    else binning = getMassBinning(treeChosen, ZZMass, isCategory, proctype==ProcessHandler::kQQBkg);

    vector<ExtendedHistogram_1D*> hMass; hMass.assign(InputTypeSize, nullptr);
    std::vector<TreeHistogramAssociation_1D> treeAssocList; treeAssocList.reserve(InputTypeSize);
    for (unsigned int it=0; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      TString hname = Form("MassDistribution_%s_%i", tree->GetName(), it);
      treeAssocList.emplace_back(
        hname+"_Smooth", "",
        tree, ZZMass, weight, isCategory
      );
      ExtendedHistogram_1D*& hh = hMass.at(it); hh = new ExtendedHistogram_1D(hname, "", binning);
      isCategory=(it==0);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        if (!isCategory) continue;
        hh->fill(ZZMass, weight);
      }
    }
    vector<TH1F*> hSmoothList = getSimultaneousSmoothHistograms(binning, treeAssocList, 5 + 5*(proctype==ProcessHandler::kZX));
    for (unsigned int it=0; it<InputTypeSize; it++){
      ExtendedHistogram_1D*& hh = hMass.at(it);
      TString hname = hh->getName();
      TH1F*& hSmooth=hSmoothList.at(it);
      *(hh->getHistogram()) = *hSmooth;
      delete hSmooth;
      hh->getHistogram()->SetName(hname);
      for (int bin=1; bin<=hh->getHistogram()->GetNbinsX(); bin++) MELAout << "Bin " << bin << ": " << hh->getHistogram()->GetBinContent(bin) << " +- " << hh->getHistogram()->GetBinError(bin) << endl;
    }
    for (unsigned int it=1; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      ExtendedHistogram_1D*& hh = hMass.at(it);
      ExtendedHistogram_1D*& hi = hMass.at(0);
      hRatio.push_back(ExtendedHistogram_1D::divideHistograms(*hh, *hi, true, Form("MassRatio_%s", tree->GetName())));
    }
    for (auto& hh:hMass) delete hh;
  }

  vector<TGraph*> grRatioList; grRatioList.reserve(hRatio.size());
  bool useEfficiencyAtan=false;
  for (auto& hr:hRatio){
    foutput->WriteTObject(hr.getHistogram());
    TGraphErrors* gr = hr.getGraph(Form("gr_%s", hr.getName().Data()));
    foutput->WriteTObject(gr);

    if (proctype==ProcessHandler::kQQBkg){
      double omitbelow=220;
      regularizeSlice(gr, nullptr, omitbelow, gr->GetX()[gr->GetN()-1]-1e-4, 10000, 0.2, 0.01);
    }
    else{
      double omitbelow=1500;
      regularizeSlice(gr, nullptr, omitbelow, gr->GetX()[gr->GetN()-1]-1e-4, 100, 0.2, 0.01);
    }
    gr->SetName(Form("%s_Smooth", gr->GetName()));
    foutput->WriteTObject(gr);
    grRatioList.push_back(gr);

    TGraph* grPatched = genericPatcher(
      gr, Form("%s_Patched", gr->GetName()),
      70., theSqrts*1000.,
      getFcn_a0plusa1timesXN<1>, getFcn_a0plusa1overXN<1>,
      false, true,
      nullptr
    );
    if (grPatched->Eval(ZZMass_Supremum)<0.){
      delete grPatched;
      grPatched = genericPatcher(
        gr, Form("%s_Patched", gr->GetName()),
        70., theSqrts*1000.,
        getFcn_a0plusa1timesXN<1>, getFcn_EfficiencyAtan,
        false, false,
        nullptr
      );
      if (grPatched->Eval(ZZMass_Supremum)<0.) MELAerr << "acquireMassRatio_ProcessNominalToNominalInclusive_one: ERROR! grPatched " << grPatched->GetName() << " is still negative at sqrts!" << endl;
      else useEfficiencyAtan=true;
    }
    delete grPatched;
  }
  if (useEfficiencyAtan) MELAout << "acquireMassRatio_ProcessNominalToNominalInclusive_one: Will use the getFcn_EfficiencyAtan function for high mass..." << endl;
  for (auto*& gr:grRatioList){
    TGraph* grPatched = nullptr;
    if (!useEfficiencyAtan) grPatched = genericPatcher(
      gr, Form("%s_Patched", gr->GetName()),
      70., theSqrts*1000.,
      getFcn_a0plusa1timesXN<1>, getFcn_a0plusa1overXN<1>,
      false, true,
      nullptr
    );
    else grPatched = genericPatcher(
      gr, Form("%s_Patched", gr->GetName()),
      70., theSqrts*1000.,
      getFcn_a0plusa1timesXN<1>, getFcn_EfficiencyAtan,
      false, false,
      nullptr
    );

    foutput->WriteTObject(grPatched);
    TSpline3* spPatched = convertGraphToSpline3(grPatched, false, false);
    TString spname=grPatched->GetName(); replaceString(spname, "gr_", "");
    spPatched->SetName(spname); spPatched->SetTitle(gr->GetName());
    foutput->WriteTObject(spPatched);

    delete spPatched;
    delete grPatched;
    delete gr;
  }

  for (unsigned int it=0; it<InputTypeSize; it++){ for (auto& finput:finputList.at(it)) finput->Close(); }
  foutput->Close();
  MELAout.close();
}


void acquireMassRatio_ProcessSystToNominal_one(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler::ProcessType proctype, const TString strGenerator
){
  if (channel==NChannels) return;
  if (syst==sNominal) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst, strGenerator)) return;
  if (proctype==ProcessHandler::kZX && strGenerator!="Data") return;
  if (
    syst==tPythiaScaleDn || syst==tPythiaScaleUp
    ||
    syst==tPythiaTuneDn || syst==tPythiaTuneUp
    ||
    syst==tMINLODn || syst==tMINLOUp
    ){
    acquireMassRatio_ProcessSystToNominal_PythiaMINLO_one(channel, category, hypo, syst, istage, fixedDate, proctype, strGenerator);
    return;
  }

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/";

  gSystem->Exec("mkdir -p " + coutput_common);
  vector<TString> INPUT_BASELINE_NAME;
  vector<TString> INPUT_NAME;
  if (proctype==ProcessHandler::kQQBkg){
    INPUT_BASELINE_NAME.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    INPUT_NAME.push_back(Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
  }
  else if (proctype==ProcessHandler::kZX){
    vector<ZXFakeRateHandler::FakeRateMethod> FRMethods;
    FRMethods.push_back(ZXFakeRateHandler::mSS);
    for (ZXFakeRateHandler::FakeRateMethod& FRMethod:FRMethods){
      const TString FRMethodName = ZXFakeRateHandler::TranslateFakeRateMethodToString(FRMethod);
      INPUT_BASELINE_NAME.push_back(Form(
        "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        thePerProcessHandle->getProcessName().Data(), FRMethodName.Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
      INPUT_NAME.push_back(Form(
        "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        thePerProcessHandle->getProcessName().Data(), FRMethodName.Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
    }
  }
  else{
    INPUT_BASELINE_NAME.push_back(Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    INPUT_NAME.push_back(Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    )
    );
    if (hypo!=kSM){
      INPUT_BASELINE_NAME.push_back(Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
      INPUT_NAME.push_back(Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
      )
      );
    }
  }

  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_MassRatios_SystToNominal_%s_%s_%s",
    strChannel.Data(), strCategory.Data(),
    strACHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  // Open the input files
  const unsigned int InputTypeSize=2; // Inclusive or dedicated category
  vector<vector<TString>> cinputList; cinputList.assign(InputTypeSize, vector<TString>());
  vector<vector<TFile*>> finputList; finputList.assign(InputTypeSize, vector<TFile*>());
  vector<vector<TTree*>> treeLists; treeLists.assign(InputTypeSize, vector<TTree*>());

  for (auto& s:INPUT_BASELINE_NAME) cinputList.at(0).push_back(cinput_common + s);
  for (auto& s:INPUT_NAME) cinputList.at(1).push_back(cinput_common + s);
  assert(cinputList.at(0).size()==cinputList.at(1).size());

  for (unsigned int it=0; it<InputTypeSize; it++){
    for (auto const& cinput:cinputList.at(it)){
      if (gSystem->AccessPathName(cinput)){
        MELAerr << "File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
        for (unsigned int it2=0; it2<InputTypeSize; it2++){ for (auto& finput:finputList.at(it2)) finput->Close(); }
        return;
      }
      if (cinput!="") finputList.at(it).push_back(TFile::Open(cinput, "read"));
    }
  }

  unsigned int ntreespertype=0;
  for (unsigned int it=0; it<InputTypeSize; it++){
    // Extract SM or AC trees
    for (auto& finput:finputList.at(it)){ if (finput) extractTreesFromDirectory(finput, treeLists.at(it)); }
    // Check number of trees
    unsigned int ntrees = treeLists.at(it).size();
    if (ntreespertype==0) ntreespertype=ntrees;
    else if (ntreespertype!=ntrees){
      MELAerr << "Number of trees in type " << it << " (" << ntrees << ") does not match Ntrees=" << ntreespertype << "! Aborting..." << endl;
      for (unsigned int it2=0; it2<InputTypeSize; it2++){ for (auto& finput2:finputList.at(it2)) finput2->Close(); }
      return;
    }
  }
  MELAout << "Acquired " << ntreespertype << " trees per type..." << endl;

  vector<vector<TTree*>> treeGroups;
  treeGroups.assign(ntreespertype, vector<TTree*>());
  for (auto& treeGroup:treeGroups) treeGroup.assign(InputTypeSize, nullptr);
  for (unsigned int it=0; it<InputTypeSize; it++){
    for (unsigned int ihypo=0; ihypo<ntreespertype; ihypo++){
      treeGroups.at(ihypo).at(it) = treeLists.at(it).at(ihypo);
    }
  }

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;

  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened output file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  for (auto& treeGroup:treeGroups){
    int nEntriesMin=-1;
    TTree* treeChosen=nullptr;
    for (TTree*& tree:treeGroup){
      MELAout << "Tree has " << tree->GetEntries() << " entries" << endl;
      if (tree->GetEntries()<=nEntriesMin || nEntriesMin<0){
        treeChosen=tree;
        nEntriesMin=tree->GetEntries();
      }
    }

    bool isCategory=(category==Inclusive);
    float ZZMass, weight;
    for (unsigned int it=0; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      tree->SetBranchStatus("*", 0);
      bookBranch(tree, "weight", &weight);
      bookBranch(tree, "ZZMass", &ZZMass);
      if (!isCategory){
        TString catFlagName = TString("is_") + strCategory + TString("_") + getACHypothesisName(hypo);
        bookBranch(tree, catFlagName, &isCategory);
      }
    }
    ExtendedBinning binning("ZZMass");
    if (proctype==ProcessHandler::kTT || proctype==ProcessHandler::kBB || strGenerator=="JHUGen"){ // On-shell - only samples
      binning.addBinBoundary(70.);
      binning.addBinBoundary(theSqrts*1000.);
    }
    else binning = getMassBinning(treeChosen, ZZMass, isCategory, proctype==ProcessHandler::kQQBkg);

    vector<ExtendedHistogram_1D*> hMass; hMass.assign(InputTypeSize, nullptr);
    std::vector<TreeHistogramAssociation_1D> treeAssocList; treeAssocList.reserve(InputTypeSize);
    for (unsigned int it=0; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      TString hname = Form("MassDistribution_%s_%i", tree->GetName(), it);
      treeAssocList.emplace_back(
        hname+"_Smooth", "",
        tree, ZZMass, weight, isCategory
      );
      ExtendedHistogram_1D*& hh = hMass.at(it); hh = new ExtendedHistogram_1D(hname, "", binning);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        if (!isCategory) continue;
        hh->fill(ZZMass, weight);
      }
    }
    vector<TH1F*> hSmoothList = getSimultaneousSmoothHistograms(binning, treeAssocList, 5);
    for (unsigned int it=0; it<InputTypeSize; it++){
      ExtendedHistogram_1D*& hh = hMass.at(it);
      TString hname = hh->getName();
      TH1F*& hSmooth=hSmoothList.at(it);
      *(hh->getHistogram()) = *hSmooth;
      delete hSmooth;
      hh->getHistogram()->SetName(hname);
      for (int bin=1; bin<=hh->getHistogram()->GetNbinsX(); bin++) MELAout << "Bin " << bin << ": " << hh->getHistogram()->GetBinContent(bin) << " +- " << hh->getHistogram()->GetBinError(bin) << endl;
    }
    vector<ExtendedHistogram_1D> hRatio;
    for (unsigned int it=1; it<InputTypeSize; it++){
      TTree*& tree = treeGroup.at(it);
      ExtendedHistogram_1D*& hh = hMass.at(it);
      ExtendedHistogram_1D*& hi = hMass.at(0);
      hRatio.push_back(ExtendedHistogram_1D::divideHistograms(*hh, *hi, true, Form("MassRatio_%s", tree->GetName())));
      for (int bin=0; bin<=hRatio.back().getHistogram()->GetNbinsX()+1; bin++){
        if (hi->getHistogram()->GetBinContent(bin)>0.) hRatio.back().getHistogram()->SetBinError(bin, hh->getHistogram()->GetBinError(bin) / hi->getHistogram()->GetBinContent(bin));
        else hRatio.back().getHistogram()->SetBinError(bin, 0);
      }
    }
    for (auto& hh:hMass) delete hh;
    for (auto& hr:hRatio){
      foutput->WriteTObject(hr.getHistogram());
      TGraphErrors* gr = hr.getGraph(Form("gr_%s", hr.getName().Data()));
      foutput->WriteTObject(gr);
      {
        if (proctype==ProcessHandler::kQQBkg){
          double omitbelow=220;
          regularizeSlice(gr, nullptr, omitbelow, gr->GetX()[gr->GetN()-1]-1e-4, 100, 0.1, 0.01);
        }
        else{
          double omitbelow=1500;
          regularizeSlice(gr, nullptr, omitbelow, gr->GetX()[gr->GetN()-1]-1e-4, 100, 0.1, 0.01);
        }
        gr->SetName(Form("%s_Smooth", gr->GetName()));
        foutput->WriteTObject(gr);

        TGraph* grPatched = genericPatcher(
          gr, Form("%s_Patched", gr->GetName()),
          70., theSqrts*1000.,
          getFcn_a0plusa1timesXN<1>, getFcn_a0plusa1overXN<1>,
          false, true,
          nullptr
        );
        foutput->WriteTObject(grPatched);
        TSpline3* spPatched = convertGraphToSpline3(grPatched, false, false);
        TString spname=grPatched->GetName(); replaceString(spname, "gr_", "");
        spPatched->SetName(spname); spPatched->SetTitle(gr->GetName());
        foutput->WriteTObject(spPatched);
        delete spPatched;
        delete grPatched;
      }
      delete gr;
    }
  }

  for (unsigned int it=0; it<InputTypeSize; it++){ for (auto& finput:finputList.at(it)) finput->Close(); }
  foutput->Close();
  MELAout.close();
}


void acquireMassRatio_ProcessSystToNominal_PythiaMINLO_one(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage, const TString fixedDate,
  ProcessHandler::ProcessType proctype, const TString strGenerator
){
  if (channel==NChannels) return;
  if (syst==sNominal) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst, strGenerator)) return;
  if (proctype==ProcessHandler::kQQBkg || proctype==ProcessHandler::kZX) return;
  if (strGenerator!="POWHEG") return;
  bool doInvertRatio = (syst==tMINLODn);

  TDirectory* rootdir = gDirectory;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  const TString strSystematics = getSystematicsName(syst);
  const TString strSystematics_Nominal = getSystematicsName(sNominal);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/MassRatios/" + strStage + "/";

  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_MassRatios_SystToNominal_%s_%s_%s",
    strChannel.Data(), strCategory.Data(),
    strACHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;

  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened output file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  // Need to loop over each sample from scratch, so begin setting up samples

  vector<int> mHListGlobal; // For ZZMass binning
  BaseTree* theOutputTree[2]={ nullptr };
  for (unsigned int i=0; i<2; i++){
    TString treename = Form("OutputTree_%i", i);
    theOutputTree[i] = new BaseTree(treename);
  }

  // Loop over the samples from scratch
  vector<TString> strSampleIdentifiers;
  if (proctype==ProcessHandler::kGG) strSampleIdentifiers.push_back("gg_Sig_POWHEG");
  else if (proctype==ProcessHandler::kTT) strSampleIdentifiers.push_back("tt_Sig_POWHEG");
  else if (proctype==ProcessHandler::kBB) strSampleIdentifiers.push_back("bb_Sig_POWHEG");
  else if (proctype==ProcessHandler::kVBF) strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
  else if (proctype==ProcessHandler::kZH) strSampleIdentifiers.push_back("ZH_Sig_POWHEG");
  else if (proctype==ProcessHandler::kWH){
    strSampleIdentifiers.push_back("WminusH_Sig_POWHEG");
    strSampleIdentifiers.push_back("WplusH_Sig_POWHEG");
  }
  else assert(0);

  // Ignore any Kfactors
  // ...

  // Register the discriminants
  vector<KDspecs> KDlist;
  getLikelihoodDiscriminants(channel, category, syst, KDlist);
  getCategorizationDiscriminants(syst, KDlist);

  for (TString const& identifier:strSampleIdentifiers){
    // For the non-nominal tree
    std::vector<ReweightingBuilder*> extraEvaluators;
    SystematicsClass* systhandle = nullptr;

    vector<TString> strSamples[2];
    vector<TString> idvector; idvector.push_back(identifier);
    getSamplesList(theSqrts, idvector, strSamples[0], sNominal);
    getSamplesList(theSqrts, idvector, strSamples[1], syst);
    unordered_map<int, std::vector<TString>> mh_samplelist_map[2];
    vector<int> mHList;
    for (unsigned int i=0; i<2; i++){
      for (TString& strSample:strSamples[i]){
        int MHVal = SampleHelpers::findPoleMass(strSample);
        TString cinput = CJLSTTree::constructCJLSTSamplePath(strSample);
        if (MHVal>0 && !gSystem->AccessPathName(cinput)){
          auto it=mh_samplelist_map[i].find(MHVal);
          if (it!=mh_samplelist_map[i].end()) it->second.push_back(strSample);
          else{
            vector<TString> vtmp; vtmp.push_back(strSample);
            mh_samplelist_map[i][MHVal]=vtmp;
          }
        }
      }
      strSamples[i].clear();
    }
    for (auto it=mh_samplelist_map[0].begin(); it!=mh_samplelist_map[0].end(); it++){
      bool mhfound=false;
      for (auto jt=mh_samplelist_map[1].begin(); jt!=mh_samplelist_map[1].end(); jt++){
        if (it->first==jt->first){ mhfound=true; break; }
      }
      if (mhfound){
        addByLowest(mHList, it->first, true);
        addByLowest(mHListGlobal, it->first, true);
      }
    }

    for (int const& mh:mHList){
      for (unsigned int i=0; i<2; i++){
        for (auto& v:mh_samplelist_map[i][mh]) strSamples[i].push_back(v);
      }
    }

    for (unsigned int i=0; i<2; i++){
      MELAout << "Looping over tree set " << i << endl;

      CJLSTSet* theSampleSet=new CJLSTSet(strSamples[i]);
      // Book common variables
      theSampleSet->bookXS(); // "xsec"
      theSampleSet->bookOverallEventWgt(); // Gen weigts "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
      for (auto& tree:theSampleSet->getCJLSTTreeList()){
        // Book common variables needed for analysis
        tree->bookBranch<float>("GenHMass", 0);
        tree->bookBranch<float>("ZZMass", -1);
        tree->bookBranch<short>("Z1Flav", 0);
        tree->bookBranch<short>("Z2Flav", 0);
        // Variables for KDs
        for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
        tree->silenceUnused(); // Will no longer book another branch
      }
      theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_OneOverNgen, false, true); // One/Ngen is a better choice when we have a sparse set of samples

      if (i==1) systhandle = constructSystematic(category, channel, proctype, syst, theSampleSet->getCJLSTTreeList(), extraEvaluators, strGenerator);

      // Build the analyzer and loop over the events
      TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
      theAnalyzer.setVerbosity(true);
      theAnalyzer.setExternalProductTree(theOutputTree[i]);
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
      // Add systematics handle
      theAnalyzer.addSystematic(strSystematics, systhandle);
      // Loop
      theAnalyzer.loop(true, false, true);

      MELAout << "Looped over tree set " << i << endl;

      delete theSampleSet;
    }

    delete systhandle;
    for (auto& rb:extraEvaluators) delete rb;

    MELAout << "End Loop over tree sets" << endl;
  }
  for (auto& KD:KDlist) delete KD.KD;

  // The output tree now has all the information needed
  for (unsigned int i=0; i<2; i++){
    MELAout << "There are " << theOutputTree[i]->getNEvents() << " products" << endl;
    //theOutputTree[i]->writeToFile(foutput);
  }

  ExtendedBinning binning_mass("ZZMass");
  binning_mass.addBinBoundary(70.);
  binning_mass.addBinBoundary(theSqrts*1000.);
  if (mHListGlobal.size()>2){ // Approximate ZZMass binning with GenHMass binning
    for (unsigned int imh=0; imh<mHListGlobal.size()-1; imh++) binning_mass.addBinBoundary((mHListGlobal.at(imh) + mHListGlobal.at(imh+1)) / 2.);
  }
  else if (mHListGlobal.size()==2 && mHListGlobal.at(1)>210.) binning_mass.addBinBoundary(210);

  vector<TString> treenamelist;
  if (proctype==ProcessHandler::kGG){
    vector<GGProcessHandler::HypothesisType> tplset = ((GGProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(kSM);
    if (hypo!=kSM){
      vector<GGProcessHandler::HypothesisType> tplset_tmp = ((GGProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(hypo);
      for (GGProcessHandler::HypothesisType& v:tplset_tmp) tplset.push_back(v);
    }
    for (auto& t:tplset) treenamelist.push_back(((GGProcessHandler*) thePerProcessHandle)->getOutputTreeName(t));
  }
  else if (proctype==ProcessHandler::kTT){
    vector<TTProcessHandler::HypothesisType> tplset = ((TTProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(kSM);
    if (hypo!=kSM){
      vector<TTProcessHandler::HypothesisType> tplset_tmp = ((TTProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(hypo);
      for (TTProcessHandler::HypothesisType& v:tplset_tmp) tplset.push_back(v);
    }
    for (auto& t:tplset) treenamelist.push_back(((TTProcessHandler*) thePerProcessHandle)->getOutputTreeName(t));
  }
  else if (proctype==ProcessHandler::kBB){
    vector<BBProcessHandler::HypothesisType> tplset = ((BBProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(kSM);
    if (hypo!=kSM){
      vector<BBProcessHandler::HypothesisType> tplset_tmp = ((BBProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(hypo);
      for (BBProcessHandler::HypothesisType& v:tplset_tmp) tplset.push_back(v);
    }
    for (auto& t:tplset) treenamelist.push_back(((BBProcessHandler*) thePerProcessHandle)->getOutputTreeName(t));
  }
  else if (proctype==ProcessHandler::kVBF || proctype==ProcessHandler::kZH || proctype==ProcessHandler::kWH){
    vector<VVProcessHandler::HypothesisType> tplset = ((VVProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(kSM);
    if (hypo!=kSM){
      vector<VVProcessHandler::HypothesisType> tplset_tmp = ((VVProcessHandler*) thePerProcessHandle)->getHypothesesForACHypothesis(hypo);
      for (VVProcessHandler::HypothesisType& v:tplset_tmp) tplset.push_back(v);
    }
    for (auto& t:tplset) treenamelist.push_back(((VVProcessHandler*) thePerProcessHandle)->getOutputTreeName(t));
  }
  else if (proctype==ProcessHandler::kQQBkg){
    treenamelist.push_back(((QQBkgProcessHandler*) thePerProcessHandle)->getOutputTreeName());
  }
  else assert(0);

  {
    foutput->cd();
    vector<ExtendedHistogram_1D*> hMass; hMass.assign(2, nullptr);
    for (unsigned int i=0; i<2; i++){
      ExtendedHistogram_1D*& hh = hMass.at(i);
      hh = new ExtendedHistogram_1D(Form("MassDistribution_%i", i), "", binning_mass);
      TTree* tree = theOutputTree[i]->getSelectedTree();
      bool isCategory=(category==Inclusive);
      float ZZMass, weight;
      tree->ResetBranchAddresses();
      tree->SetBranchAddress("ZZMass", &ZZMass);
      tree->SetBranchAddress("weight", &weight);
      if (!isCategory){
        TString catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;
        tree->SetBranchAddress(catFlagName, &isCategory);
      }
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        if (!isCategory) continue;
        if (
          proctype==ProcessHandler::kZH && i==1 && CJLSTversion<=180224 && theDataPeriod=="2016"
          && (
            syst==tPythiaScaleDn || syst==tPythiaScaleUp
            ||
            syst==tPythiaTuneDn || syst==tPythiaTuneUp
            )
          ) weight /= 0.148; // FIXME: Rescale Pythia variation samples for missing filter efficiency; should be part of xsec in the future
        hh->fill(ZZMass, weight);
      }
      MELAout << "Mass integral of " << hh->getName() << ": " << hh->getHistogram()->Integral() << endl;
      //foutput->WriteTObject(hh->getHistogram());
      //foutput->WriteTObject(hh->getProfileX());
    }
    ExtendedHistogram_1D*& hh = hMass.at(1);
    ExtendedHistogram_1D*& hi = hMass.at(0);
    ExtendedHistogram_1D hRatio = (
      !doInvertRatio ?
      ExtendedHistogram_1D::divideHistograms(*hh, *hi, false, Form("MassRatio"))
      :
      ExtendedHistogram_1D::divideHistograms(*hi, *hh, false, Form("MassRatio"))
      );
    if ((syst==tMINLODn || syst==tMINLOUp) && theDataPeriod=="2016"){ // Mass ratio is set to 1 in 2016 MC
      for (int ix=1; ix<=hRatio.getHistogram()->GetNbinsX(); ix++){
        hRatio.getHistogram()->SetBinContent(ix, 1);
        hRatio.getHistogram()->SetBinError(ix, 0);
      }
    }
    foutput->WriteTObject(hRatio.getHistogram());
    foutput->WriteTObject(hRatio.getProfileX());
    for (auto& htmp:hMass) delete htmp;

    for (auto& name:treenamelist){
      TString hname=hRatio.getName() + "_" + name;
      TGraphErrors* gr = hRatio.getGraph(Form("gr_%s", hname.Data()));
      foutput->WriteTObject(gr);
      {
        gr->SetName(Form("%s_Smooth", gr->GetName()));
        foutput->WriteTObject(gr);

        double* yy=gr->GetY();
        double* xx=gr->GetX();
        int np=gr->GetN();
        std::vector<pair<double, double>> xy_new;
        for (int ip=0; ip<int(binning_mass.getMax()-binning_mass.getMin()); ip++){ // Add a point every GeV
          double xval=ip+binning_mass.getMin();
          if (xval<xx[0]){
            xy_new.emplace_back(xval, yy[0]);
          }
          else if (xval>=xx[np-1]){
            xy_new.emplace_back(xval, yy[np-1]);
          }
          else{
            for (int jp=0; jp<np-1; jp++){
              if (xval>=xx[jp] && xval<xx[jp+1]){
                double f=(xval-xx[jp])/(xx[jp+1]-xx[jp]);
                xy_new.emplace_back(xval, yy[jp]*(1.-f)+f*yy[jp+1]);
              }
            }
          }
        }

        TGraph* grPatched = makeGraphFromPair(xy_new, TString(gr->GetName())+"_Patched");
        foutput->WriteTObject(grPatched);
        TSpline3* spPatched = convertGraphToSpline3(grPatched, false, false);
        TString spname=grPatched->GetName(); replaceString(spname, "gr_", "");
        spPatched->SetName(spname); spPatched->SetTitle(gr->GetName());
        foutput->WriteTObject(spPatched);
        delete spPatched;
        delete grPatched;
      }
      delete gr;
    }
  }

  // Get KD binning
  for (int massregion=0; massregion<(int) NMassRegions; massregion++){
    TString massregionname = getMassRegionName((CategorizationHelpers::MassRegion) massregion);
    foutput->cd();
    TDirectory* savedir = foutput->mkdir(massregionname); savedir->cd();
    MELAout << "Making " << massregionname << " distributions" << endl;

    bool isCategory=(category==Inclusive);
    float weight=0;
    vector<TString> KDset;
    {
      vector<TString> KDset2=getACHypothesisKDNameSet(hypo, category, (CategorizationHelpers::MassRegion) massregion);
      if (massregion!=kOnshell || hypo==kSM) KDset.push_back("ZZMass"); // Only off-shell, or on-shell SM use ZZMass
      appendVector(KDset, KDset2);
    }
    MELAout << "KDs: " << KDset << endl;
    unordered_map<TString, float> KDvars;
    for (TString const& KDname:KDset) KDvars[KDname]=0;
    vector<ExtendedBinning> KDbinning;
    for (TString const& KDname:KDset){
      if (KDname!="ZZMass") KDbinning.push_back(getDiscriminantFineBinning(channel, category, KDname, (CategorizationHelpers::MassRegion) massregion));
      else KDbinning.push_back(binning_mass);
    }
    unsigned int nKDs = KDset.size();
    TString catFlagName="";
    if (!isCategory) catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;

    for (unsigned int i=0; i<2; i++){
      TTree* tree = theOutputTree[i]->getSelectedTree();
      tree->ResetBranchAddresses();
      bookBranch(tree, "weight", &weight);
      for (TString const& KDname:KDset) bookBranch(tree, KDname, &(KDvars.find(KDname)->second));
      if (catFlagName!="") bookBranch(tree, catFlagName, &isCategory);
    }

    if (nKDs==2){
      assert(
        KDvars.find(KDset.at(0))!=KDvars.cend()
        &&
        KDvars.find(KDset.at(1))!=KDvars.cend()
      );

      std::vector<TreeHistogramAssociation_2D> treeAssocList; treeAssocList.reserve(2);
      for (unsigned int i=0; i<2; i++){
        TTree* tree = theOutputTree[i]->getSelectedTree();
        treeAssocList.emplace_back(
          (i==0 ? strSystematics_Nominal.Data() : strSystematics.Data()), "",
          tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, weight, isCategory
        );
      }
      vector<TH2F*> hDistro = getSimultaneousSmoothHistograms(
        KDbinning.at(0), KDbinning.at(1), treeAssocList,
        1, 10
      );
      for (auto& hh:hDistro){
        MELAout << "Integral of distribution : " << hh->Integral() << endl;
        if (KDset.at(0)=="ZZMass") conditionalizeHistogram<TH2F>(hh, 0, nullptr, false, true);
        else{
          double hist_integral = getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), 1, hh->GetNbinsY(), false, nullptr);
          hh->Scale(1./hist_integral);
        }
      }
      TH2F* hRatio = new TH2F("RatioWithKD", "", KDbinning.at(0).getNbins(), KDbinning.at(0).getBinning(), KDbinning.at(1).getNbins(), KDbinning.at(1).getBinning());
      if (!doInvertRatio) divideHistograms(hDistro.at(1), hDistro.at(0), hRatio, false);
      else divideHistograms(hDistro.at(0), hDistro.at(1), hRatio, false);
      savedir->WriteTObject(hRatio); delete hRatio;
      for (auto& hh:hDistro) delete hh;
    }
    else if (nKDs==3){
      assert(
        KDvars.find(KDset.at(0))!=KDvars.cend()
        &&
        KDvars.find(KDset.at(1))!=KDvars.cend()
        &&
        KDvars.find(KDset.at(2))!=KDvars.cend()
      );

      std::vector<TreeHistogramAssociation_3D> treeAssocList; treeAssocList.reserve(2);
      for (unsigned int i=0; i<2; i++){
        TTree* tree = theOutputTree[i]->getSelectedTree();
        treeAssocList.emplace_back(
          (i==0 ? strSystematics_Nominal.Data() : strSystematics.Data()), "",
          tree, KDvars.find(KDset.at(0))->second, KDvars.find(KDset.at(1))->second, KDvars.find(KDset.at(2))->second, weight, isCategory
        );
      }
      vector<TH3F*> hDistro = getSimultaneousSmoothHistograms(
        KDbinning.at(0), KDbinning.at(1), KDbinning.at(2), treeAssocList,
        1, 10, 10
      );
      for (auto& hh:hDistro){
        MELAout << "Integral of distribution : " << hh->Integral() << endl;
        if (KDset.at(0)=="ZZMass") conditionalizeHistogram<TH3F>(hh, 0, nullptr, false, true);
        else{
          double hist_integral = getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), 1, hh->GetNbinsY(), 1, hh->GetNbinsZ(), false, nullptr);
          hh->Scale(1./hist_integral);
        }
      }
      TH3F* hRatio = new TH3F("RatioWithKD", "", KDbinning.at(0).getNbins(), KDbinning.at(0).getBinning(), KDbinning.at(1).getNbins(), KDbinning.at(1).getBinning(), KDbinning.at(2).getNbins(), KDbinning.at(2).getBinning());
      if (!doInvertRatio) divideHistograms(hDistro.at(1), hDistro.at(0), hRatio, false);
      else divideHistograms(hDistro.at(0), hDistro.at(1), hRatio, false);
      savedir->WriteTObject(hRatio); delete hRatio;
      for (auto& hh:hDistro) delete hh;
    }
  }
  foutput->cd();

  for (unsigned int i=0; i<2; i++) delete theOutputTree[i];
  foutput->Close();
  MELAout.close();
}


ExtendedBinning getMassBinning(TTree* tree, float& ZZMass, bool& flag, bool separateZ4l){
  ExtendedBinning binning("ZZMass");
  if (!tree) return binning;

  const int nEntries = tree->GetEntries();
  const int nEntriesWithCut = tree->GetEntries("ZZMass>=220");

  unsigned int countThreshold=20000;
  unsigned int nbins_est = nEntriesWithCut/countThreshold;
  if (nbins_est<15) nbins_est=15;
  else if (nbins_est>30) nbins_est=30;
  //if (nbins_est<10) nbins_est=10;
  //else if (nbins_est>25) nbins_est=25;
  countThreshold=nEntriesWithCut/nbins_est;

  const float massLow = 220.;
  const float massHigh = theSqrts*1000.;
  const float massWidth = 0.01;
  int const nbinsraw = (massHigh-massLow)/massWidth+0.5;
  TH1D* hmass = new TH1D("hmass", "", nbinsraw, massLow, massHigh);
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    if (!flag) continue;
    if (ZZMass>=220.) hmass->Fill(ZZMass);
  }

  // Get binning from the histogram
  binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(hmass->GetNbinsX()+1));
  vector<unsigned int> counts;
  unsigned int count=0;
  for (int bin=hmass->GetNbinsX(); bin>=1; bin--){
    count += hmass->GetBinContent(bin);
    if (count>countThreshold || bin==1){
      counts.push_back(count);
      binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(bin));
      count=0;
    }
  }
  delete hmass;
  std::reverse(counts.begin(), counts.end());
  MELAout << counts << endl;
  MELAout << binning.getBinningVector() << endl;
  // These lines guarantee count>countThreshold in every bin
  if (counts.at(0)<countThreshold && counts.size()>1){
    counts.at(1) += counts.at(0);
    counts.erase(counts.begin());
    binning.removeBinLowEdge(1);
  }
  MELAout << counts << endl;
  MELAout << binning.getBinningVector() << endl;
  if (counts.back()<countThreshold && counts.size()>1){
    counts.at(counts.size()-2) += counts.at(counts.size()-1);
    counts.erase(counts.begin()+counts.size()-1);
    binning.removeBinLowEdge(binning.getNbins()-1);
  }
  MELAout << counts << endl;
  MELAout << binning.getBinningVector() << endl;
  // Merge every two bins except the last one
  {
    unsigned int nbins_old=binning.getNbins();
    for (unsigned int bin=0; bin<nbins_old-2; bin++){
      if (bin%2==0) continue;
      binning.removeBinLowEdge(nbins_old-bin-1);
    }
  }
  MELAout << counts << endl;
  MELAout << binning.getBinningVector() << endl;
  // Finally, add the low mass boundaries
  binning.addBinBoundary(70.);
  if (separateZ4l) binning.addBinBoundary(105.);
  binning.addBinBoundary(140.);
  //binning.addBinBoundary(185.);
  MELAout << "Final binning for tree " << tree->GetName() << ": nbins=" << binning.getNbins() << " [ " << binning.getBinningVector() << " ] ( " << countThreshold << " / " << nEntries << " )" << endl;
  return binning;
}


#endif
