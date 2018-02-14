#ifndef PLOTPROCESSCHECKSTAGE_H
#define PLOTPROCESSCHECKSTAGE_H

#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "smoothenHistograms.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


ExtendedBinning getMassBinning(TTree* tree);

void acquireProcessMassRatioToNominalInclusive(
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
  TString INPUT_INCLUSIVE_NAME, INPUT_NAME;
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
      strACHypo.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics_Nominal.Data(),
      strGenerator.Data(),
      ".root"
    );
    INPUT_NAME = Form(
      "HtoZZ%s_%s_%s_FinalTemplates_%s_%s_%s%s",
      strChannel.Data(), strCategory.Data(),
      strACHypo.Data(),
      thePerProcessHandle->getProcessName().Data(),
      strSystematics.Data(),
      strGenerator.Data(),
      ".root"
    );
  }
  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_%s_MassRatiosToNominalInclusive_%s_%s_%s%s",
    strChannel.Data(), strCategory.Data(),
    strACHypo.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    ".root"
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += "_plot.root";
  OUTPUT_LOG_NAME += "_plot.log";

  // Open the input files
  vector<TString> cinputList;
  cinputList.push_back(cinput_common + INPUT_INCLUSIVE_NAME);
  cinputList.push_back(cinput_common + INPUT_NAME);
  vector<TFile*> finputList;
  for (TString const& cinput:cinputList){
    if (gSystem->AccessPathName(cinput)){
      MELAerr << "File " << cinput << " is not found! Run " << strStage << " functions first." << endl;
      for (TFile*& finput:finputList) finput->Close();
      return;
    }
    TFile* ftmp = TFile::Open(cinput, "read");
    finputList.push_back(ftmp);
  }
  const unsigned int nfiles=cinputList.size();
  vector<vector<TTree*>> treeLists; treeLists.assign(nfiles, vector<TTree*>());
  unsigned int ntreesperfile=0;
  for (unsigned int ifile=0; ifile<nfiles; ifile++){
    TFile*& finput=finputList.at(ifile);
    extractTreesFromDirectory(finput, treeLists.at(ifile));
    if (ntreesperfile==0) ntreesperfile=treeLists.at(ifile).size();
    else if (ntreesperfile!=treeLists.at(ifile).size()){
      MELAerr << "Number of trees in file " << cinputList.at(ifile) << " does not match Ntrees=" << ntreesperfile << "! Aborting..." << endl;
      for (TFile*& ftmp:finputList) ftmp->Close();
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

    vector<TH1F*> hMass; hMass.assign(nfiles, nullptr);
    for (unsigned int ifile=0; ifile<nfiles; ifile++){
      TH1F*& hh = hMass.at(ifile);
      TTree*& tree = treeGroup.at(ifile);
      hh = new TH1F(Form("MassDistribution_%s_%i", tree->GetName(), ifile), "", binning.getNbins(), binning.getBinning());
      hh->Sumw2();
      float ZZMass, weight;
      tree->SetBranchAddress("ZZMass", &ZZMass);
      tree->SetBranchAddress("weight", &weight);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        hh->Fill(ZZMass, weight);
      }
    }
    vector<TH1F*> hRatio; hRatio.assign(nfiles-1, nullptr);
    for (unsigned int ifile=1; ifile<nfiles; ifile++){
      TTree*& tree = treeGroup.at(ifile);
      TH1F*& hr = hRatio.at(ifile-1); hr = new TH1F(Form("MassRatio_%s_%i", tree->GetName(), ifile), "", binning.getNbins(), binning.getBinning());
      TH1F*& hh = hMass.at(ifile);
      TH1F*& hi = hMass.at(0);
      for (int bin=1; bin<=hh->GetNbinsX(); bin++){
        hr->SetBinContent(bin, hh->GetBinContent(bin)/hi->GetBinContent(bin));
      }
    }
    for (TH1F*& hh:hMass) delete hh;
    for (TH1F*& hr:hRatio){
      foutput->WriteTObject(hr);
      delete hr;
    }
  }

  for (TFile*& finput:finputList) finput->Close();
  foutput->Close();
  MELAout.close();
}

ExtendedBinning getMassBinning(TTree* tree){
  ExtendedBinning binning("ZZMass");
  if (!tree) return binning;

  float ZZMass;
  tree->SetBranchAddress("ZZMass", &ZZMass);
  const int nEntries = tree->GetEntries();

  unsigned int countThreshold=100;
  if (nEntries>1000000) countThreshold=20000;
  else if (nEntries>500000) countThreshold=10000;
  else if (nEntries>100000) countThreshold=5000;
  else if (nEntries>50000) countThreshold=1000;
  else if (nEntries>10000) countThreshold=500;

  const float massLow = 70.;
  const float massHigh = theSqrts*1000.;
  const float massWidth = 0.01;
  int const nbinsraw = (massHigh-massLow)/massWidth+0.5;
  TH1D* hmass = new TH1D("hmass", "", nbinsraw, massLow, massHigh);
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    hmass->Fill(ZZMass);
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
  MELAout << "Final binning for tree " << tree->GetName() << ": [ " << binning.getBinningVector() << " ] ( " << countThreshold << " / " << nEntries << " )" << endl;
  return binning;
}


#endif
