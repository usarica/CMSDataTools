#ifndef ACQUIREPROCESSMASSRATIOS_H
#define ACQUIREPROCESSMASSRATIOS_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


ExtendedBinning getMassBinning(TTree* tree);

void acquireMassRatio_ProcessNominalToNominalInclusive(
  const Channel channel, const Category category, const ACHypothesis hypo,
  const unsigned int istage,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  const SystematicVariationTypes syst=sNominal;

  if (channel==NChannels) return;
  if (category==Inclusive) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=nullptr;
  switch (proctype){
  case ProcessHandler::kGG:
    thePerProcessHandle = &TemplateHelpers::OffshellGGProcessHandle;
    break;
  case ProcessHandler::kVV:
    thePerProcessHandle = &TemplateHelpers::OffshellVVProcessHandle;
    break;
  case ProcessHandler::kQQBkg:
    thePerProcessHandle = &TemplateHelpers::OffshellQQBkgProcessHandle;
    break;
  case ProcessHandler::kZX:
    thePerProcessHandle = &TemplateHelpers::OffshellZXProcessHandle;
    break;
  default:
    break;
  };
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst)) return;

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
  TString INPUT_INCLUSIVE_NAME="", INPUT_INCLUSIVE_AC_NAME="", INPUT_NAME="", INPUT_AC_NAME="";
  if (proctype==ProcessHandler::kQQBkg){
    INPUT_INCLUSIVE_NAME = Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory_Inclusive.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    );
    INPUT_NAME = Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    );
  }
  else{
    INPUT_INCLUSIVE_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory_Inclusive.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    );
    INPUT_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    );
    if (hypo!=kSM){
      INPUT_INCLUSIVE_AC_NAME = Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory_Inclusive.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      );
      INPUT_AC_NAME = Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
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
  vector<pair<TString, TString>> cinputList;
  cinputList.push_back(pair<TString, TString>(cinput_common + INPUT_INCLUSIVE_NAME, cinput_common + INPUT_INCLUSIVE_AC_NAME));
  cinputList.push_back(pair<TString, TString>(cinput_common + INPUT_NAME, cinput_common + INPUT_AC_NAME));
  vector<pair<TFile*, TFile*>> finputList;
  for (auto const& cinput:cinputList){
    if (gSystem->AccessPathName(cinput.first)){
      MELAerr << "File " << cinput.first << " is not found! Run " << strStage << " functions first." << endl;
      for (auto& finput:finputList){
        if (finput.first) finput.first->Close();
        if (finput.second) finput.second->Close();
      }
      return;
    }
    pair<TFile*, TFile*> fnew(nullptr, nullptr);
    if (cinput.first!="") fnew.first = TFile::Open(cinput.first, "read");
    if (cinput.second!="") fnew.second = TFile::Open(cinput.second, "read");
    finputList.push_back(fnew);
  }
  const unsigned int nfiles=cinputList.size();
  vector<vector<TTree*>> treeLists; treeLists.assign(nfiles, vector<TTree*>());
  unsigned int ntreesperfile=0;
  for (unsigned int ifile=0; ifile<nfiles; ifile++){
    TFile* finput;
    // Extract SM trees
    finput=finputList.at(ifile).first;
    if (finput) extractTreesFromDirectory(finput, treeLists.at(ifile));
    // Extract AC trees
    finput=finputList.at(ifile).second;
    if (finput) extractTreesFromDirectory(finput, treeLists.at(ifile));
    if (ntreesperfile==0) ntreesperfile=treeLists.at(ifile).size();
    else if (ntreesperfile!=treeLists.at(ifile).size()){
      MELAerr << "Number of trees in files " << cinputList.at(ifile).first << " and " << cinputList.at(ifile).second << " does not match Ntrees=" << ntreesperfile << "! Aborting..." << endl;
      for (auto& ftmp:finputList){
        if (ftmp.first) ftmp.first->Close();
        if (ftmp.second) ftmp.second->Close();
      }
      return;
    }
    MELAout << "Acquired " << ntreesperfile << " trees..." << endl;
  }
  vector<vector<TTree*>> treeGroups;
  treeGroups.assign(ntreesperfile, vector<TTree*>());
  for (auto& treeGroup:treeGroups) treeGroup.assign(nfiles, nullptr);
  for (unsigned int ifile=0; ifile<nfiles; ifile++){
    for (unsigned int ihypo=0; ihypo<ntreesperfile; ihypo++){
      treeGroups.at(ihypo).at(ifile) = treeLists.at(ifile).at(ihypo);
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
    ExtendedBinning binning=getMassBinning(treeChosen);

    vector<ExtendedHistogram_1D*> hMass; hMass.assign(nfiles, nullptr);
    for (unsigned int ifile=0; ifile<nfiles; ifile++){
      ExtendedHistogram_1D*& hh = hMass.at(ifile);
      TTree*& tree = treeGroup.at(ifile);
      hh = new ExtendedHistogram_1D(Form("MassDistribution_%s_%i", tree->GetName(), ifile), "", binning);
      float ZZMass, weight;
      tree->SetBranchAddress("ZZMass", &ZZMass);
      tree->SetBranchAddress("weight", &weight);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        hh->fill(ZZMass, weight);
      }
    }
    vector<ExtendedHistogram_1D> hRatio;
    for (unsigned int ifile=1; ifile<nfiles; ifile++){
      TTree*& tree = treeGroup.at(ifile);
      ExtendedHistogram_1D*& hh = hMass.at(ifile);
      ExtendedHistogram_1D*& hi = hMass.at(0);
      hRatio.push_back(ExtendedHistogram_1D::divideHistograms(*hh, *hi, true, Form("MassRatio_%s", tree->GetName())));
    }
    for (auto& hh:hMass) delete hh;
    for (auto& hr:hRatio){
      foutput->WriteTObject(hr.getHistogram());
      TGraphErrors* gr = hr.getGraph(Form("gr_%s", hr.getName().Data()));
      foutput->WriteTObject(gr);
      {
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
      }
      delete gr;
    }
  }

  for (auto& finput:finputList){
    if (finput.first) finput.first->Close();
    if (finput.second) finput.second->Close();
  }
  foutput->Close();
  MELAout.close();
}


void acquireMassRatio_ProcessSystToNominal(
  const Channel channel, const Category category, const ACHypothesis hypo,
  const SystematicVariationTypes syst,
  const unsigned int istage,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=nullptr;
  switch (proctype){
  case ProcessHandler::kGG:
    thePerProcessHandle = &TemplateHelpers::OffshellGGProcessHandle;
    break;
  case ProcessHandler::kVV:
    thePerProcessHandle = &TemplateHelpers::OffshellVVProcessHandle;
    break;
  case ProcessHandler::kQQBkg:
    thePerProcessHandle = &TemplateHelpers::OffshellQQBkgProcessHandle;
    break;
  case ProcessHandler::kZX:
    thePerProcessHandle = &TemplateHelpers::OffshellZXProcessHandle;
    break;
  default:
    break;
  };
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst)) return;

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
  TString INPUT_BASELINE_NAME="", INPUT_BASELINE_AC_NAME="", INPUT_NAME="", INPUT_AC_NAME="";
  if (proctype==ProcessHandler::kQQBkg){
    INPUT_BASELINE_NAME = Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    );
    INPUT_NAME = Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    );
  }
  else{
    INPUT_BASELINE_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    );
    INPUT_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      getACHypothesisName(kSM).Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    );
    if (hypo!=kSM){
      INPUT_BASELINE_AC_NAME = Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics_Nominal.Data(),
        strGenerator.Data(),
        ".root"
      );
      INPUT_AC_NAME = Form(
        "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
        strChannel.Data(), strCategory.Data(),
        strACHypo.Data(),
        thePerProcessHandle->getProcessName().Data(),
        strSystematics.Data(),
        strGenerator.Data(),
        ".root"
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
  vector<pair<TString, TString>> cinputList;
  cinputList.push_back(pair<TString, TString>(cinput_common + INPUT_BASELINE_NAME, cinput_common + INPUT_BASELINE_AC_NAME));
  cinputList.push_back(pair<TString, TString>(cinput_common + INPUT_NAME, cinput_common + INPUT_AC_NAME));
  vector<pair<TFile*, TFile*>> finputList;
  for (auto const& cinput:cinputList){
    if (gSystem->AccessPathName(cinput.first)){
      MELAerr << "File " << cinput.first << " is not found! Run " << strStage << " functions first." << endl;
      for (auto& finput:finputList){
        if (finput.first) finput.first->Close();
        if (finput.second) finput.second->Close();
      }
      return;
    }
    pair<TFile*, TFile*> fnew(nullptr, nullptr);
    if (cinput.first!="") fnew.first = TFile::Open(cinput.first, "read");
    if (cinput.second!="") fnew.second = TFile::Open(cinput.second, "read");
    finputList.push_back(fnew);
  }
  const unsigned int nfiles=cinputList.size();
  vector<vector<TTree*>> treeLists; treeLists.assign(nfiles, vector<TTree*>());
  unsigned int ntreesperfile=0;
  for (unsigned int ifile=0; ifile<nfiles; ifile++){
    TFile* finput;
    // Extract SM trees
    finput=finputList.at(ifile).first;
    if (finput) extractTreesFromDirectory(finput, treeLists.at(ifile));
    // Extract AC trees
    finput=finputList.at(ifile).second;
    if (finput) extractTreesFromDirectory(finput, treeLists.at(ifile));
    if (ntreesperfile==0) ntreesperfile=treeLists.at(ifile).size();
    else if (ntreesperfile!=treeLists.at(ifile).size()){
      MELAerr << "Number of trees in files " << cinputList.at(ifile).first << " and " << cinputList.at(ifile).second << " does not match Ntrees=" << ntreesperfile << "! Aborting..." << endl;
      for (auto& ftmp:finputList){
        if (ftmp.first) ftmp.first->Close();
        if (ftmp.second) ftmp.second->Close();
      }
      return;
    }
    MELAout << "Acquired " << ntreesperfile << " trees..." << endl;
  }
  vector<vector<TTree*>> treeGroups;
  treeGroups.assign(ntreesperfile, vector<TTree*>());
  for (auto& treeGroup:treeGroups) treeGroup.assign(nfiles, nullptr);
  for (unsigned int ifile=0; ifile<nfiles; ifile++){
    for (unsigned int ihypo=0; ihypo<ntreesperfile; ihypo++){
      treeGroups.at(ihypo).at(ifile) = treeLists.at(ifile).at(ihypo);
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
    ExtendedBinning binning=getMassBinning(treeChosen);

    vector<ExtendedHistogram_1D*> hMass; hMass.assign(nfiles, nullptr);
    for (unsigned int ifile=0; ifile<nfiles; ifile++){
      ExtendedHistogram_1D*& hh = hMass.at(ifile);
      TTree*& tree = treeGroup.at(ifile);
      hh = new ExtendedHistogram_1D(Form("MassDistribution_%s_%i", tree->GetName(), ifile), "", binning);
      float ZZMass, weight;
      tree->SetBranchAddress("ZZMass", &ZZMass);
      tree->SetBranchAddress("weight", &weight);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        hh->fill(ZZMass, weight);
      }
    }
    vector<ExtendedHistogram_1D> hRatio;
    for (unsigned int ifile=1; ifile<nfiles; ifile++){
      TTree*& tree = treeGroup.at(ifile);
      ExtendedHistogram_1D*& hh = hMass.at(ifile);
      ExtendedHistogram_1D*& hi = hMass.at(0);
      hRatio.push_back(ExtendedHistogram_1D::divideHistograms(*hh, *hi, true, Form("MassRatio_%s", tree->GetName())));
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

  for (auto& finput:finputList){
    if (finput.first) finput.first->Close();
    if (finput.second) finput.second->Close();
  }
  foutput->Close();
  MELAout.close();
}


ExtendedBinning getMassBinning(TTree* tree){
  ExtendedBinning binning("ZZMass");
  if (!tree) return binning;

  float ZZMass;
  tree->SetBranchAddress("ZZMass", &ZZMass);
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
  // These lines guarantee count>countThreshold in every bin
  if (counts.at(0)<countThreshold && counts.size()>1){
    counts.at(1) += counts.at(0);
    counts.erase(counts.begin());
    binning.removeBinLowEdge(1);
  }
  // Merge every wo bins except the last one
  {
    unsigned int nbins_old=binning.getNbins();
    for (unsigned int bin=0; bin<nbins_old-2; bin++){
      if (bin%2==0) continue;
      binning.removeBinLowEdge(nbins_old-bin-1);
    }
  }
  // Finally, add the low mass boundaries
  binning.addBinBoundary(70.);
  binning.addBinBoundary(105.);
  binning.addBinBoundary(140.);
  //binning.addBinBoundary(185.);
  MELAout << "Final binning for tree " << tree->GetName() << ": nbins=" << binning.getNbins() << " [ " << binning.getBinningVector() << " ] ( " << countThreshold << " / " << nEntries << " )" << endl;
  return binning;
}


#endif
