#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooBinning.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "TIterator.h"
#include "Math/Minimizer.h"
#include "HiggsAnalysis/CombinedLimit/interface/AsymPow.h"
#include "HiggsAnalysis/CombinedLimit/interface/AsymQuad.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooPiecewisePolynomial.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooRealFlooredSumPdf.h"
#include "common_includes.h"
#include "CheckSetTemplatesCategoryScheme.h"
#include "TemplatesEventAnalyzer.h"


using namespace RooFit;

// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif


ExtendedBinning getIntermediateBinning(TH1F const* hwgt, TH1F* const hunwgt){
  ExtendedBinning res;
  double integralerror=0;
  double integral=hwgt->IntegralAndError(1, hwgt->GetNbinsX(), integralerror);
  double sumn=hunwgt->Integral(1, hunwgt->GetNbinsX());
  unsigned int nbinsmin=20;
  unsigned int nevtsmin=10000;
  unsigned int nevtsmax=100000;
  unsigned int nbins = sumn/nevtsmax;
  if (nbins<nbinsmin) nbins=nbinsmin;
  if (sumn/nbins<nevtsmin) nbins=sumn/nevtsmin;
  double sumNthr=sumn/double(nbins);

  res.addBinBoundary(hwgt->GetXaxis()->GetBinLowEdge(1));
  vector<float> sumNList;
  float sumN=0;
  float sumW=0, sumWsq=0;
  for (int bin=1; bin<=hwgt->GetNbinsX(); bin++){
    sumN += hunwgt->GetBinContent(bin);
    sumW += hwgt->GetBinContent(bin);
    sumWsq += pow(hwgt->GetBinError(bin), 2);
    if (sumN>sumNthr || bin==hwgt->GetNbinsX()){
      res.addBinBoundary(hwgt->GetXaxis()->GetBinLowEdge(bin+1));
      sumNList.push_back(sumN);
      sumN=0;
      sumW=0;
      sumWsq=0;
    }
  }
  // These lines guarantee sum>sumThreshold in every bin
  if (sumNList.size()>1 && sumNList.back()<sumNthr){
    res.removeBinLowEdge(res.getNbins()-1);
  }
  return res;
}


void acquireResolution_one(const Channel channel, const Category category, const TString fixedDate, ProcessHandler::ProcessType proctype, const TString strGenerator){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (proctype==ProcessHandler::kZX) return;
  if (strGenerator!="POWHEG") return;

  TDirectory* curdir = gDirectory;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(kSM);
  const TString strSqrts = Form("%i", theSqrts);
  const TString strYear = theDataPeriod;
  const TString strSqrtsYear = strSqrts + "TeV_" + strYear;

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Resolution/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME_CORE = Form(
    "HtoZZ%s_%s_Resolution_%s",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data()
  );
  TString OUTPUT_NAME=OUTPUT_NAME_CORE;
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  TString canvasnamecore = coutput_common + "c_" + OUTPUT_NAME + "_" + strGenerator;
  TString coutput = coutput_common + OUTPUT_NAME + ".root";
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME + ".log";
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  // Need to loop over each sample from scratch, so begin setting up samples
  vector<int> mHListGlobal; // For ZZMass binning
  BaseTree* theOutputTree = new BaseTree("OutputTree");

  // Loop over the samples from scratch
  vector<TString> strSampleIdentifiers;
  if (proctype==ProcessHandler::kGG) strSampleIdentifiers.push_back("gg_Sig_POWHEG");
  else if (proctype==ProcessHandler::kVBF) strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
  else if (proctype==ProcessHandler::kZH) strSampleIdentifiers.push_back("ZH_Sig_POWHEG");
  else if (proctype==ProcessHandler::kWH) strSampleIdentifiers.push_back("WH_Sig_POWHEG");
  else if (proctype==ProcessHandler::kQQBkg) strSampleIdentifiers.push_back("qq_Bkg_Combined");
  else assert(0);

  // Ignore any Kfactors
  // ...

  // Kfactor variable names
  vector<TString> strKfactorVars;
  if (proctype==ProcessHandler::kQQBkg){
    //strKfactorVars.push_back("KFactor_QCD_qqZZ_M"); // No need for this K factor
    strKfactorVars.push_back("KFactor_EW_qqZZ");
  }

  // Register only the categorization the discriminants
  vector<KDspecs> KDlist;
  getCategorizationDiscriminants(sNominal, KDlist);

  for (TString const& identifier:strSampleIdentifiers){
    // For the non-nominal tree
    std::vector<ReweightingBuilder*> extraEvaluators;
    SystematicsClass* systhandle = nullptr;

    vector<TString> strSamples;
    vector<TString> idvector; idvector.push_back(identifier);
    getSamplesList(theSqrts, idvector, strSamples, sNominal);
    unordered_map<int, std::vector<TString>> mh_samplelist_map;
    for (TString& strSample:strSamples){
      int MHVal = SampleHelpers::findPoleMass(strSample);
      TString cinput = CJLSTTree::constructCJLSTSamplePath(strSample);
      if (MHVal>0 && !gSystem->AccessPathName(cinput)){
        auto it=mh_samplelist_map.find(MHVal);
        if (it!=mh_samplelist_map.end()) it->second.push_back(strSample);
        else{
          vector<TString> vtmp; vtmp.push_back(strSample);
          mh_samplelist_map[MHVal]=vtmp;
        }
      }
    }
    strSamples.clear();
    for (auto it=mh_samplelist_map.begin(); it!=mh_samplelist_map.end(); it++) addByLowest(mHListGlobal, it->first, true);
    for (int const& mh:mHListGlobal){
      for (auto& v:mh_samplelist_map[mh]) strSamples.push_back(v);
    }

    CJLSTSet* theSampleSet=new CJLSTSet(strSamples);
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
      // Variables for KDs
      for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
      tree->silenceUnused(); // Will no longer book another branch
    }
    theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_XsecOverNgen_RelRenormToSumNgen, true, true); // We don't care about total xsec, but we care abou realative normalization in W- vs W+ H

    //systhandle = constructSystematic(category, channel, proctype, syst, theSampleSet->getCJLSTTreeList(), extraEvaluators, strGenerator);

    // Additional EW K factor reweighting in QQBkg
    ReweightingBuilder* regularewgtBuilder=nullptr;
    if (!strKfactorVars.empty()){
      ExtendedBinning GenHMassInclusiveBinning("GenHMass");
      regularewgtBuilder = new ReweightingBuilder(strKfactorVars, getSimpleWeight);
      regularewgtBuilder->rejectNegativeWeights(true);
      regularewgtBuilder->setDivideByNSample(false);
      regularewgtBuilder->setWeightBinning(GenHMassInclusiveBinning);
      for (auto& tree:theSampleSet->getCJLSTTreeList()) regularewgtBuilder->setupWeightVariables(tree, -1, 0);
    }

    // Build the analyzer and loop over the events
    TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
    theAnalyzer.setExternalProductTree(theOutputTree);
    // Book common variables needed for analysis
    theAnalyzer.addConsumed<float>("PUWeight");
    theAnalyzer.addConsumed<float>("genHEPMCweight");
    theAnalyzer.addConsumed<float>("dataMCWeight");
    theAnalyzer.addConsumed<float>("trigEffWeight");
    theAnalyzer.addConsumed<float>("GenHMass");
    theAnalyzer.addConsumed<float>("ZZMass");
    theAnalyzer.addConsumed<short>("Z1Flav");
    theAnalyzer.addConsumed<short>("Z2Flav");
    // Add reweighting builders
    if (regularewgtBuilder) theAnalyzer.addReweightingBuilder("RegularRewgt", regularewgtBuilder);
    // Add discriminant builders
    for (auto& KD:KDlist){ theAnalyzer.addDiscriminantBuilder(KD.KDname, KD.KD, KD.KDvars); }
    // Add systematics handle
    //theAnalyzer.addSystematic(strSystematics, systhandle);
    // Loop
    theAnalyzer.loop(true, false, true);

    delete regularewgtBuilder;
    delete systhandle;
    for (auto& rb:extraEvaluators) delete rb;
    delete theSampleSet;

    MELAout << "End Loop over tree sets" << endl;
  }
  for (auto& KD:KDlist) delete KD.KD;

  // The output tree now has all the information needed
  MELAout << "There are " << theOutputTree->getNEvents() << " products" << endl;

  // Setup the binning
  ExtendedBinning binning_mass("ZZMass");
  binning_mass.addBinBoundary(70.);
  binning_mass.addBinBoundary(theSqrts*1000.);
  if (mHListGlobal.size()>2){ // Approximate ZZMass binning with GenHMass binning
    for (unsigned int imh=0; imh<mHListGlobal.size()-1; imh++){
      binning_mass.addBinBoundary((mHListGlobal.at(imh) + mHListGlobal.at(imh+1)) / 2.);
    }
  }
  else{
    binning_mass.addBinBoundary(105);
    binning_mass.addBinBoundary(124);
    binning_mass.addBinBoundary(140);
    binning_mass.addBinBoundary(160);
    binning_mass.addBinBoundary(220);
    binning_mass.addBinBoundary(1000);
  }

  RooRealVar var_mreco("var_mreco", "m^{reco}_{4l} (GeV)", 125, 70, 13000);
  RooRealVar var_mtrue("var_mtrue", "m^{true}_{4l} (GeV)", 125, 0, 13000);
  RooRealVar var_mdiff("var_mdiff", "m^{reco}_{4l}-m^{true}_{4l} (GeV)", 0, -2000, 2000);
  RooRealVar var_weight("var_weight", "Event weight", 1, -10, 10); var_weight.removeMin(); var_weight.removeMax();
  RooArgSet treevars; treevars.add(var_mdiff); treevars.add(var_weight);
  RooArgSet incltreevars; incltreevars.add(var_mtrue); incltreevars.add(var_mreco); incltreevars.add(var_weight);
  RooArgSet conditionals; conditionals.add(var_mtrue);

  TString prefix = "CMS_zz4l_";
  vector<RooRealVar*> CB_parameter_list; CB_parameter_list.reserve(6);
  RooRealVar CB_mean(prefix + "CB_mean", "", 0, -3, 3); CB_parameter_list.push_back(&CB_mean);
  RooRealVar CB_width(prefix + "CB_width", "", 1, 0.3, 15); CB_parameter_list.push_back(&CB_width);
  RooRealVar CB_alpha1(prefix + "CB_alpha1", "", 2, 0, 4); CB_parameter_list.push_back(&CB_alpha1);
  RooRealVar CB_alpha2(prefix + "CB_alpha2", "", 3, 0, 10); CB_parameter_list.push_back(&CB_alpha2);
  RooRealVar CB_n1(prefix + "CB_n1", "", 1, 0, 10); CB_parameter_list.push_back(&CB_n1);
  RooRealVar CB_n2(prefix + "CB_n2", "", 1.5, 0, 40); CB_parameter_list.push_back(&CB_n2);
  vector<double> CB_parameter_init; CB_parameter_init.reserve(6);
  for (auto*& par:CB_parameter_list) CB_parameter_init.push_back(par->getVal());
  vector<vector<RooRealVar>> CB_piecewisepolypars_list; CB_piecewisepolypars_list.assign(6, vector<RooRealVar>());
  vector<RooArgList> CB_piecewisepolypars_args; CB_piecewisepolypars_args.assign(6, RooArgList());
  vector<RooPiecewisePolynomial> CB_piecewisepoly_list; CB_piecewisepoly_list.reserve(6);

  RooDoubleCB pdf(
    "pdf", "",
    var_mdiff,
    CB_mean, CB_width,
    CB_alpha1, CB_n1,
    CB_alpha2, CB_n2
  );

  TTree* tree = theOutputTree->getSelectedTree();
  float mtrue, mreco, wgt;
  bool isCategory=(category==Inclusive);
  tree->SetBranchAddress("ZZMass", &mreco);
  tree->SetBranchAddress("GenHMass", &mtrue);
  tree->SetBranchAddress("weight", &wgt);
  if (!isCategory){
    TString catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;
    tree->SetBranchAddress(catFlagName, &isCategory);
  }

  // Setup inclusive data
  // Add maually because the inclusive tree might be too big to simply copy
  RooDataSet incldata("incldata", "incldata", incltreevars, var_weight.GetName());
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (!isCategory) continue;
    var_mreco.setVal(mreco);
    var_mtrue.setVal(mtrue);
    var_weight.setVal(wgt);
    incldata.add(incltreevars);
  }
  var_mreco.setVal(125);
  var_mtrue.setVal(125);
  var_weight.setVal(1);

  vector<vector<pair<float, float>>> CB_parameter_val_list; CB_parameter_val_list.assign(CB_parameter_list.size(), vector<pair<float, float>>());
  vector<vector<pair<float, float>>> CB_parameter_errdn_list; CB_parameter_errdn_list.assign(CB_parameter_list.size(), vector<pair<float, float>>());
  vector<vector<pair<float, float>>> CB_parameter_errup_list; CB_parameter_errup_list.assign(CB_parameter_list.size(), vector<pair<float, float>>());

  for (unsigned int bin=0; bin<binning_mass.getNbins(); bin++){
    MELAout << "Fitting bin " << bin << " / " << binning_mass.getNbins() << endl;
    TTree newtree("DataTree", "");
    float mdiff;
    newtree.Branch("var_weight", &wgt);
    newtree.Branch("var_mdiff", &mdiff);
    float sum_wgt=0;
    float sum_mtrueXwgt=0;
    float sum_mtruesqXwgt=0;
    float mdiffmin=13000, mdiffmax=-1;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      if (mtrue<binning_mass.getBinLowEdge(bin) || mtrue>=binning_mass.getBinHighEdge(bin)) continue;
      if (!isCategory) continue;
      mdiff=mreco-mtrue;
      mdiffmin=std::min(mdiffmin, mdiff);
      mdiffmax=std::max(mdiffmax, mdiff);
    }
    var_mdiff.setRange(mdiffmin, mdiffmax);
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      if (mtrue<binning_mass.getBinLowEdge(bin) || mtrue>=binning_mass.getBinHighEdge(bin)) continue;
      if (!isCategory) continue;
      mdiff=mreco-mtrue;
      sum_mtrueXwgt += mtrue*wgt;
      sum_mtruesqXwgt += mtrue*mtrue*wgt;
      sum_wgt += wgt;
      newtree.Fill();
    }
    float avg_mtrue = 0, avg_mtrue_err=0; if (sum_wgt!=0.){ avg_mtrue = sum_mtrueXwgt/sum_wgt; avg_mtrue_err = sqrt(sum_mtruesqXwgt/sum_wgt-avg_mtrue*avg_mtrue); }
    var_mtrue.setVal(avg_mtrue);
    cout << "Number of entries in tree " << newtree.GetName() << ": " << newtree.GetEntries() << endl;
    cout << "Average mtrue = " << avg_mtrue << endl;

    // Prepare the dataset
    RooDataSet data("data", "data", &newtree, treevars, nullptr, var_weight.GetName());
    CB_mean.setRange(-var_mtrue.getVal()/50., var_mtrue.getVal()/50.);
    CB_width.setMax(var_mtrue.getVal()/10.);
    for (unsigned int ipar=0; ipar<CB_parameter_list.size(); ipar++) CB_parameter_list.at(ipar)->setVal(CB_parameter_init.at(ipar));
    const unsigned int nit=3;
    bool fitSuccessful=false;
    for (unsigned int it=0; it<nit; it++){
      RooLinkedList cmdList;
      RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
      //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
      RooCmdArg sumw2Arg = RooFit::SumW2Error(true); cmdList.Add((TObject*) &sumw2Arg);
      RooCmdArg hesseArg = RooFit::Hesse(it>0); cmdList.Add((TObject*) &hesseArg);
      //RooCmdArg initialhesseArg = RooFit::InitialHesse(it==nit-1); cmdList.Add((TObject*) &initialhesseArg);
      //RooCmdArg minosArg = RooFit::Minos(it==nit-1); cmdList.Add((TObject*) &minosArg);
      //RooCmdArg minimizerArg = RooFit::Minimizer("Minuit", "migrad"); cmdList.Add((TObject*)&minimizerArg);
      RooCmdArg minimizerStrategyArg = RooFit::Strategy(2); cmdList.Add((TObject*) &minimizerStrategyArg);
      // Misc. options
      RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
      //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

      RooFitResult* fitResult=pdf.fitTo(data, cmdList);
      int fitStatus=-1;
      if (fitResult){
        fitStatus = fitResult->status();
        cout << "Fit status is " << fitStatus << endl;
        //cout << "Fit properties:" << endl;
        //fitResult->Print("v");
      }
      delete fitResult;
      if (it==nit-1 && fitStatus%10==4){
        fitResult=pdf.fitTo(data, cmdList);
        if (fitResult){
          fitStatus = fitResult->status();
          cout << "Retried with status " << fitStatus << endl;
        }
        delete fitResult;
      }
      fitSuccessful=(fitStatus==0);
    }
    //for (unsigned int ipar=0; ipar<CB_parameter_list.size(); ipar++){
    //  RooRealVar* par = CB_parameter_list.at(ipar);
    //  double parerrorlo=par->getAsymErrorLo(); if (parerrorlo==0.) fitSuccessful=false;
    //  double parerrorhi=par->getAsymErrorHi(); if (parerrorhi==0.) fitSuccessful=false;
    //}
    for (unsigned int ipar=0; ipar<CB_parameter_list.size(); ipar++){
      RooRealVar* par = CB_parameter_list.at(ipar);
      CB_parameter_val_list.at(ipar).emplace_back(var_mtrue.getVal(), par->getVal());
      bool hitThr=false;
      double parerrorlo=par->getAsymErrorLo(); if (parerrorlo==0.) parerrorlo=par->getError(); if (parerrorlo==0. || !fitSuccessful){ parerrorlo=par->getVal()-par->getMin(); hitThr=true; }
      double parerrorhi=par->getAsymErrorHi(); if (parerrorhi==0.) parerrorhi=par->getError(); if (parerrorhi==0. || !fitSuccessful){ parerrorhi=-par->getVal()+par->getMax(); hitThr=true; }
      //double parerrorlo=par->getAsymErrorLo(); if (!fitSuccessful) parerrorlo=par->getVal()-par->getMin();
      //double parerrorhi=par->getAsymErrorHi(); if (!fitSuccessful) parerrorhi=-par->getVal()+par->getMax();
      if (parerrorlo<0.) parerrorlo *= -1.;
      if (fabs(parerrorlo-parerrorhi)>5.*std::min(parerrorhi, parerrorlo) || hitThr){
        parerrorlo=par->getVal()-par->getMin();
        parerrorhi=-par->getVal()+par->getMax();
      }
      if (par->getVal()-parerrorlo<par->getMin()) parerrorlo=par->getVal()-par->getMin();
      if (par->getVal()+parerrorhi>par->getMax()) parerrorhi=-par->getVal()+par->getMax();
      CB_parameter_errdn_list.at(ipar).emplace_back(avg_mtrue_err, parerrorlo);
      CB_parameter_errup_list.at(ipar).emplace_back(avg_mtrue_err, parerrorhi);
      par->setAsymError(-parerrorlo, parerrorhi);
      cout << "Fitted " << par->GetName() << " = " << par->getVal() << " +" << parerrorhi << "/-" << parerrorlo << endl;
    }
  }

  delete theOutputTree;

  const unsigned int polyndof=3;
  const unsigned int nnodes=2;
  const unsigned int nfcn=nnodes+1;
  const unsigned int npolypars = nnodes + (polyndof-2)*(nfcn-2) + (polyndof-1)*2;
  vector<RooRealVar> nodeVars;
  for (unsigned int inode=0; inode<nnodes; inode++){
    double mean=(inode==0 ? 200 : 1000);
    double low=mean*0.8;
    double high=mean*1.2;
    nodeVars.emplace_back(prefix+Form("CB_node_%i", inode), "", mean, low, high);
  }

  vector<TGraphAsymmErrors*> grlist; grlist.reserve(CB_parameter_list.size());
  for (unsigned int ipar=0; ipar<CB_parameter_list.size(); ipar++){
    RooRealVar* par = CB_parameter_list.at(ipar);

    // Construct the TGraph first
    TGraphAsymmErrors* gr = makeGraphAsymErrFromPair(
      CB_parameter_val_list.at(ipar),
      CB_parameter_errdn_list.at(ipar),
      CB_parameter_errup_list.at(ipar),
      "gr_" + TString(par->GetName())
    );
    gr->SetLineWidth(2);
    gr->SetLineColor(kBlack);
    gr->SetMarkerStyle(30);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerColor(kBlack);
    foutput->WriteTObject(gr);
    grlist.push_back(gr);

    RooPiecewisePolynomial qfitpoly("qfit", "", nfcn, polyndof);
    RooPiecewisePolynomial* qfitpolyptr=&qfitpoly;
    typedef double (RooPiecewisePolynomial::*RPPEvalFcn)(double*, double*) const;
    RPPEvalFcn evalFcn=&RooPiecewisePolynomial::evaluate;
    TF1 quickFit("quickFit", qfitpolyptr, evalFcn, gr->GetX()[0], gr->GetX()[gr->GetN()-1], npolypars);
    quickFit.SetLineColor(kRed); quickFit.SetLineWidth(2);
    for (unsigned int inode=0; inode<nnodes; inode++){ quickFit.SetParLimits(inode, 0, 13000); quickFit.FixParameter(inode, nodeVars.at(inode).getVal()); }
    for (unsigned int ip=nnodes; ip<npolypars; ip++){ quickFit.SetParameter(ip, 0); quickFit.SetParLimits(ip, -100, 100); }
    const double* params;
    const double* parerrs;
    MELAout << "Fit for parameters of " << par->GetName() << endl;
    // Try the fit until there is no need to fix the parameters
    bool fitretry=true;
    while (fitretry){
      // Prefit 5 times
      int fitCtr=5;
      while (fitCtr>=0){
        gr->Fit(&quickFit, "EX0M0");
        fitCtr--;
      }
      fitretry=false;
      params=quickFit.GetParameters();
      parerrs=quickFit.GetParErrors();
      unsigned int parctr=nnodes;
      for (unsigned int ifcn=0; ifcn<nfcn; ifcn++){
        const unsigned int nfcnpars=(ifcn==0 || ifcn==nfcn-1 ? polyndof-1 : polyndof-2);
        MELAout << "- Checking the "<< nfcnpars << " parameters of function " << ifcn << " of the CB parameter " << par->GetName() << endl;
        for (unsigned int jpar=0; jpar<nfcnpars; jpar++){
          if (jpar>0){
            double const& val=params[parctr];
            double const& err=parerrs[parctr];
            MELAout << "\t- Checking parameter " << jpar << "/" << nfcnpars << ": (val, err) = (" << val << ", " << err << ")" << endl;
            if (err!=0. && fabs(val)<=fabs(err)){
              MELAout << "\t\t- Will retry the fit by fixing parameter " << parctr << endl;
              quickFit.SetParError(parctr, 0);
              quickFit.FixParameter(parctr, 0);
              fitretry=true;
            }
          }
          parctr++;
        }
      }
    }
    MELAout << "Fit for parameters is done." << endl;

    // Construct parameters for the piecewise polynomial
    unsigned int parctr=nnodes;
    for (unsigned int ifcn=0; ifcn<nfcn; ifcn++){
      const unsigned int nfcnpars=(ifcn==0 || ifcn==nfcn-1 ? polyndof-1 : polyndof-2);
      MELAout << "Constructing the "<< nfcnpars << " parameters of function " << ifcn << " of the CB parameter " << par->GetName() << endl;
      for (unsigned int jpar=0; jpar<nfcnpars; jpar++){
        MELAout << "\t- Constructing parameter " << jpar << "/" << nfcnpars << endl;
        double cval=params[parctr], lval, hval;
        double err=parerrs[parctr];
        if (err==0.){ lval=cval; hval=cval; }
        else{
          lval=cval-parerrs[parctr];
          hval=cval+parerrs[parctr];
        }
        CB_piecewisepolypars_list.at(ipar).emplace_back(
          Form("%s_fcn%i_par%i", par->GetName(), ifcn, jpar), "",
          cval, lval, hval
        );
        if (err==0.) CB_piecewisepolypars_list.at(ipar).back().setConstant(true);
        parctr++;
      }
    }

    MELAout << "Adding parameters to the relevant RooArgList" << endl;
    for (RooRealVar& nodeVar:nodeVars){
      MELAout << "\t- Adding " << nodeVar.GetName() << ": " << nodeVar.getVal() << " [ " << nodeVar.getMin() << ", " << nodeVar.getMax() << " ]" << endl;
      CB_piecewisepolypars_args.at(ipar).add(nodeVar);
    }
    for (RooRealVar& var:CB_piecewisepolypars_list.at(ipar)){
      MELAout << "\t- Adding " << var.GetName() << ": " << var.getVal() << " [ " << var.getMin() << ", " << var.getMax() << " ]" << endl;
      CB_piecewisepolypars_args.at(ipar).add(var);
    }
    MELAout << "Constructing " << Form("%s_piecewisepolynomial", par->GetName()) << endl;
    CB_piecewisepoly_list.emplace_back(
      Form("%s_piecewisepolynomial", par->GetName()), "", var_mtrue, CB_piecewisepolypars_args.at(ipar), nfcn, polyndof
    );

    TCanvas cgr(Form("c_%s", par->GetName()), "");
    gr->Draw("ap");
    quickFit.Draw("csame");
    cgr.SaveAs(Form("%s_%s.pdf", canvasnamecore.Data(), par->GetName()));
    cgr.SaveAs(Form("%s_%s.png", canvasnamecore.Data(), par->GetName()));
    cgr.Close();
  }


  // Inclusive pdf
  TString newprefix = "CMS_"+OUTPUT_NAME_CORE+"_"+strSqrtsYear+"_";

  RooConstVar scale_uncval_center(newprefix + "final_CB_CMS_scale_emcenter", "", 0);
  RooRealVar scale_uncvar_e("CMS_scale_e", "CMS_scale_e", 0, -7, 7);
  RooConstVar scale_uncval_e_up(newprefix + "final_CB_CMS_scale_eUp", "", 0.002);
  RooConstVar scale_uncval_e_dn(newprefix + "final_CB_CMS_scale_eDown", "", -0.002);
  RooRealVar scale_uncvar_mu("CMS_scale_m", "CMS_scale_m", 0, -7, 7);
  RooConstVar scale_uncval_mu_up(newprefix + "final_CB_CMS_scale_mUp", "", 0.001);
  RooConstVar scale_uncval_mu_dn(newprefix + "final_CB_CMS_scale_mDown", "", -0.001);
  RooArgList scalemeanthetaarglist; RooArgList scalemeanfcnarglist;
  scalemeanfcnarglist.add(scale_uncval_center);
  if (channel==k2e2mu){
    scalemeanthetaarglist.add(scale_uncvar_e);
    scalemeanfcnarglist.add(scale_uncval_e_up);
    scalemeanfcnarglist.add(scale_uncval_e_dn);
    scalemeanthetaarglist.add(scale_uncvar_mu);
    scalemeanfcnarglist.add(scale_uncval_mu_up);
    scalemeanfcnarglist.add(scale_uncval_mu_dn);
  }
  else if (channel==k4mu){
    scalemeanthetaarglist.add(scale_uncvar_mu);
    scalemeanfcnarglist.add(scale_uncval_mu_up);
    scalemeanfcnarglist.add(scale_uncval_mu_dn);
  }
  else if (channel==k4e){
    scalemeanthetaarglist.add(scale_uncvar_e);
    scalemeanfcnarglist.add(scale_uncval_e_up);
    scalemeanfcnarglist.add(scale_uncval_e_dn);
  }
  AsymQuad scale_uncval(newprefix + "final_CB_CMS_scale_em_AsymQuad", "", scalemeanfcnarglist, scalemeanthetaarglist, 1., 2);
  TString strscalemeanFormula; RooArgList scalemeanarglist;
  scalemeanarglist.add(CB_piecewisepoly_list.at(0));
  scalemeanarglist.add(var_mtrue);
  scalemeanarglist.add(scale_uncval);
  strscalemeanFormula="(@0+@1)*(1.+@2)-@1"; // Until a new procedure is found, keep var_mtrue as part of the scale unc. definition

  RooRealVar res_uncvar_e("CMS_res_e", "CMS_res_e", 0, -7, 7);
  RooConstVar res_uncval_e_up(newprefix + "final_CB_CMS_res_eUp", "", 1.2);
  RooConstVar res_uncval_e_dn(newprefix + "final_CB_CMS_res_eDown", "", 1./1.2);
  AsymPow res_uncval_e(newprefix + "final_CB_CMS_res_e_AsymPow", "", res_uncval_e_dn, res_uncval_e_up, res_uncvar_e);
  RooRealVar res_uncvar_mu("CMS_res_m", "CMS_res_m", 0, -7, 7);
  RooConstVar res_uncval_mu_up(newprefix + "final_CB_CMS_res_mUp", "", 1.2);
  RooConstVar res_uncval_mu_dn(newprefix + "final_CB_CMS_res_mDown", "", 1./1.2);
  AsymPow res_uncval_mu(newprefix + "final_CB_CMS_res_m_AsymPow", "", res_uncval_mu_dn, res_uncval_mu_up, res_uncvar_mu);
  TString strreswidthFormula; RooArgList reswidtharglist;
  reswidtharglist.add(CB_piecewisepoly_list.at(1));
  if (channel==k2e2mu){
    strreswidthFormula="@0*@1*@2";
    reswidtharglist.add(res_uncval_e);
    reswidtharglist.add(res_uncval_mu);
  }
  else if (channel==k4mu){
    strreswidthFormula="@0*@1";
    reswidtharglist.add(res_uncval_mu);
  }
  else if (channel==k4e){
    strreswidthFormula="@0*@1";
    reswidtharglist.add(res_uncval_e);
  }

  RooFormulaVar incl_CB_mean(newprefix + "final_CB_mean", strscalemeanFormula, scalemeanarglist);
  RooFormulaVar incl_CB_width(newprefix + "final_CB_width", strreswidthFormula, reswidtharglist);
  RooFormulaVar incl_CB_alpha1(newprefix + "final_CB_alpha1", "max(@0,0.1)", RooArgList(CB_piecewisepoly_list.at(2)));
  RooFormulaVar incl_CB_alpha2(newprefix + "final_CB_alpha2", "max(@0,0.1)", RooArgList(CB_piecewisepoly_list.at(3)));
  RooFormulaVar incl_CB_n1(newprefix + "final_CB_n1", "max(@0,0.1)", RooArgList(CB_piecewisepoly_list.at(4)));
  RooFormulaVar incl_CB_n2(newprefix + "final_CB_n2", "max(@0,0.1)", RooArgList(CB_piecewisepoly_list.at(5)));
  RooDoubleCB incl_pdf(
    newprefix + "final_CB", newprefix + "final_CB",
    var_mreco, var_mtrue,
    incl_CB_mean, incl_CB_width,
    incl_CB_alpha1, incl_CB_n1,
    incl_CB_alpha2, incl_CB_n2
  );

  /*
  // Fit all parameters simultaneously
  res_uncvar_e.setConstant(true);
  res_uncvar_mu.setConstant(true);
  scale_uncvar_e.setConstant(true);
  scale_uncvar_mu.setConstant(true);
  {
    RooLinkedList cmdList;
    RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
    RooCmdArg condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*) &condObsArg);
    RooCmdArg sumw2Arg = RooFit::SumW2Error(true); cmdList.Add((TObject*) &sumw2Arg);
    RooCmdArg hesseArg = RooFit::Hesse(false); cmdList.Add((TObject*) &hesseArg);
    RooCmdArg minimizerStrategyArg = RooFit::Strategy(0); cmdList.Add((TObject*) &minimizerStrategyArg);
    // Misc. options
    RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
    RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
    //RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

    RooFitResult* fitResult=incl_pdf.fitTo(incldata, cmdList);
    if (fitResult){
      int fitStatus = fitResult->status();
      cout << "Fit status is " << fitStatus << endl;
      cout << "Fit properties:" << endl;
      fitResult->Print("v");
    }
    delete fitResult;
  }
  res_uncvar_e.setConstant(false);
  res_uncvar_mu.setConstant(false);
  scale_uncvar_e.setConstant(false);
  scale_uncvar_mu.setConstant(false);
  */

  {
    RooAbsData* reducedData=incldata.reduce("var_mreco>=105 && var_mreco<=140 && var_mtrue>=124.5 && var_mtrue<=125.5");

    RooPlot incl_plot(var_mreco, 105, 140, 80);
    reducedData->plotOn(&incl_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2));
    incl_pdf.plotOn(&incl_plot, LineColor(kRed), LineWidth(2), ProjWData(RooArgSet(var_mtrue), *reducedData));

    TCanvas can("ZZMass_105_140", "");
    incl_plot.Draw();
    can.SaveAs(Form("%s_%s.pdf", canvasnamecore.Data(), can.GetName()));
    can.SaveAs(Form("%s_%s.png", canvasnamecore.Data(), can.GetName()));
    can.Close();

    delete reducedData;
  }

  // Rename all variables
  for (auto& var:nodeVars){
    TString varname=var.GetName();
    replaceString<TString, const TString>(varname, prefix, newprefix);
    var.SetName(varname);
    var.SetTitle(varname);
    var.setConstant(true);
  }
  for (auto& v:CB_piecewisepolypars_list){
    for (auto& var:v){
      TString varname=var.GetName();
      replaceString<TString, const TString>(varname, prefix, newprefix);
      var.SetName(varname);
      var.SetTitle(varname);
      var.setConstant(true);
    }
  }
  for (auto& var:CB_piecewisepoly_list){
    TString varname=var.GetName();
    replaceString<TString, const TString>(varname, prefix, newprefix);
    var.SetName(varname);
    var.SetTitle(varname);
  }
  var_mtrue.SetName("MH");
  var_mreco.SetName("mass");
  incl_pdf.SetName("ResolutionModel");
  var_mtrue.SetTitle("MH");
  var_mreco.SetTitle("mass");
  incl_pdf.SetTitle("ResolutionModel");

  RooWorkspace w("w", "");
  w.import(incl_pdf, RecycleConflictNodes());
  //w.import(var_mreco, RenameVariable(var_mreco.GetName(), "newmass"));
  //RooAbsArg* pdfnew=w.factory("EDIT::ResolutionModelCopy(ResolutionModel, mass=newmass)");
  //w.import(*pdfnew, RecycleConflictNodes());
  foutput->WriteTObject(&w);

  for (auto*& gr:grlist){
    delete gr;
  }
  foutput->Close();
  MELAout.close();
}

void acquireH125OnshellMassShape_one(const Channel channel, const Category category, const TString fixedDate, ProcessHandler::ProcessType proctype, const TString strGenerator){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (proctype==ProcessHandler::kZX || proctype==ProcessHandler::kQQBkg) return;
  if (strGenerator!="POWHEG") return;

  TDirectory* curdir = gDirectory;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(kSM);
  const TString strSqrts = Form("%i", theSqrts);
  const TString strYear = theDataPeriod;
  const TString strSqrtsYear = strSqrts + "TeV_" + strYear;

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Resolution/";
  gSystem->Exec("mkdir -p " + coutput_common);

  TString OUTPUT_NAME_CORE = Form(
    "HtoZZ%s_%s_FinalMassShape_%s",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data()
  );
  TString OUTPUT_NAME=OUTPUT_NAME_CORE;
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  TString canvasnamecore = coutput_common + "c_" + OUTPUT_NAME + "_" + strGenerator;
  TString coutput = coutput_common + OUTPUT_NAME + ".root";
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME + ".log";
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout << "Opened file " << coutput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  // Need to loop over each sample from scratch, so begin setting up samples
  vector<int> mHListGlobal; // For ZZMass binning
  BaseTree* theOutputTree = new BaseTree("OutputTree");

  // Loop over the samples from scratch
  vector<TString> strSampleIdentifiers;
  if (proctype==ProcessHandler::kGG) strSampleIdentifiers.push_back("gg_Sig_POWHEG");
  else if (proctype==ProcessHandler::kVBF) strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
  else if (proctype==ProcessHandler::kZH) strSampleIdentifiers.push_back("ZH_Sig_POWHEG");
  else if (proctype==ProcessHandler::kWH) strSampleIdentifiers.push_back("WH_Sig_POWHEG");
  else assert(0);

  // Ignore any Kfactors
  // ...

  // Kfactor variable names
  //vector<TString> strKfactorVars;

  // Register only the categorization the discriminants
  vector<KDspecs> KDlist;
  getCategorizationDiscriminants(sNominal, KDlist);

  for (TString const& identifier:strSampleIdentifiers){
    // For the non-nominal tree
    std::vector<ReweightingBuilder*> extraEvaluators;
    SystematicsClass* systhandle = nullptr;

    vector<TString> strSamples;
    vector<TString> idvector; idvector.push_back(identifier);
    getSamplesList(theSqrts, idvector, strSamples, sNominal);
    unordered_map<int, std::vector<TString>> mh_samplelist_map;
    for (TString& strSample:strSamples){
      int MHVal = SampleHelpers::findPoleMass(strSample);
      if (MHVal!=125) continue;
      TString cinput = CJLSTTree::constructCJLSTSamplePath(strSample);
      if (!gSystem->AccessPathName(cinput)){
        auto it=mh_samplelist_map.find(MHVal);
        if (it!=mh_samplelist_map.end()) it->second.push_back(strSample);
        else{
          vector<TString> vtmp; vtmp.push_back(strSample);
          mh_samplelist_map[MHVal]=vtmp;
        }
      }
    }
    strSamples.clear();
    for (auto it=mh_samplelist_map.begin(); it!=mh_samplelist_map.end(); it++) addByLowest(mHListGlobal, it->first, true);
    for (int const& mh:mHListGlobal){
      for (auto& v:mh_samplelist_map[mh]) strSamples.push_back(v);
    }

    CJLSTSet* theSampleSet=new CJLSTSet(strSamples);
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
      //for (auto& s:strKfactorVars) tree->bookBranch<float>(s, 1);
      // Variables for KDs
      for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
      tree->silenceUnused(); // Will no longer book another branch
    }
    theSampleSet->setPermanentWeights(CJLSTSet::NormScheme_XsecOverNgen_RelRenormToSumNgen, true, true); // We don't care about total xsec, but we care abou realative normalization in W- vs W+ H

    //systhandle = constructSystematic(category, channel, proctype, syst, theSampleSet->getCJLSTTreeList(), extraEvaluators, strGenerator);

    // Build the analyzer and loop over the events
    TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
    theAnalyzer.setExternalProductTree(theOutputTree);
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
    //theAnalyzer.addSystematic(strSystematics, systhandle);
    // Loop
    theAnalyzer.loop(true, false, true);

    delete systhandle;
    for (auto& rb:extraEvaluators) delete rb;
    delete theSampleSet;

    MELAout << "End Loop over tree sets" << endl;
  }
  for (auto& KD:KDlist) delete KD.KD;

  // The output tree now has all the information needed
  MELAout << "There are " << theOutputTree->getNEvents() << " products" << endl;

  // Setup the binning
  ExtendedBinning binning_mass=getDiscriminantFineBinning(channel, category, "ZZMass", kOnshell);

  RooRealVar var_mreco("ZZMass", "m^{reco}_{4l} (GeV)", 125, binning_mass.getMin(), binning_mass.getMax());
  { // Set binning of ZZMass
    RooBinning var_mreco_binning(binning_mass.getNbins(), binning_mass.getBinning());
    var_mreco.setBinning(var_mreco_binning);
  }
  RooConstVar var_mtrue("MH", "MH", 125);
  RooRealVar var_weight("weight", "Event weight", 1, -10, 10); var_weight.removeMin(); var_weight.removeMax();

  // Setup inclusive data
  RooArgSet treevars; treevars.add(var_mreco); treevars.add(var_weight);
  RooDataSet data("data", "data", treevars, var_weight.GetName());
  {
    TTree* tree = theOutputTree->getSelectedTree();
    float mreco, wgt;
    bool isCategory=(category==Inclusive);
    tree->SetBranchAddress("ZZMass", &mreco);
    tree->SetBranchAddress("weight", &wgt);
    if (!isCategory){
      TString catFlagName = TString("is_") + strCategory + TString("_") + strACHypo;
      tree->SetBranchAddress(catFlagName, &isCategory);
    }
    // Add maually because of the category flag
    float sumwgts=0;
    unsigned int ndata=0;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      if (!isCategory) continue;
      if (mreco<var_mreco.getMin() || mreco>var_mreco.getMax()) continue;
      var_mreco.setVal(mreco);
      var_weight.setVal(wgt);
      sumwgts += wgt;
      ndata++;
      data.add(treevars);
    }
    var_mreco.setVal(125);
    var_weight.setVal(1);
    data.Print("v");
    MELAout << "Average weight: " << sumwgts / float(ndata) << endl;
  }
  delete theOutputTree;

  TString prefix = "CMS_"+OUTPUT_NAME_CORE+"_"+strSqrtsYear+"_";
  vector<RooRealVar> CB_parameter_list; CB_parameter_list.reserve(6);
  CB_parameter_list.emplace_back(prefix + "CB_mean", "", -0.1, -3, 3);
  CB_parameter_list.emplace_back(prefix + "CB_width", "", 1, 0.3, 15);
  CB_parameter_list.emplace_back(prefix + "CB_alpha1", "", 1, 0, 10);
  CB_parameter_list.emplace_back(prefix + "CB_alpha2", "", 1, 0, 10);
  CB_parameter_list.emplace_back(prefix + "CB_n1", "", 1, 0, 10);
  CB_parameter_list.emplace_back(prefix + "CB_n2", "", 1, 0, 40);
  vector<double> CB_parameter_init; CB_parameter_init.reserve(6);
  for (auto& par:CB_parameter_list) CB_parameter_init.push_back(par.getVal());

  RooConstVar scale_uncval_center(prefix + "final_CB_CMS_scale_emcenter", "", 0);
  RooRealVar scale_uncvar_e("CMS_scale_e", "CMS_scale_e", 0, -7, 7);
  RooConstVar scale_uncval_e_up(prefix + "final_CB_CMS_scale_eUp", "", 0.002);
  RooConstVar scale_uncval_e_dn(prefix + "final_CB_CMS_scale_eDown", "", -0.002);
  RooRealVar scale_uncvar_mu("CMS_scale_m", "CMS_scale_m", 0, -7, 7);
  RooConstVar scale_uncval_mu_up(prefix + "final_CB_CMS_scale_mUp", "", 0.001);
  RooConstVar scale_uncval_mu_dn(prefix + "final_CB_CMS_scale_mDown", "", -0.001);
  RooArgList scalemeanthetaarglist; RooArgList scalemeanfcnarglist;
  scalemeanfcnarglist.add(scale_uncval_center);
  if (channel==k2e2mu){
    scalemeanthetaarglist.add(scale_uncvar_e);
    scalemeanfcnarglist.add(scale_uncval_e_up);
    scalemeanfcnarglist.add(scale_uncval_e_dn);
    scalemeanthetaarglist.add(scale_uncvar_mu);
    scalemeanfcnarglist.add(scale_uncval_mu_up);
    scalemeanfcnarglist.add(scale_uncval_mu_dn);
  }
  else if (channel==k4mu){
    scalemeanthetaarglist.add(scale_uncvar_mu);
    scalemeanfcnarglist.add(scale_uncval_mu_up);
    scalemeanfcnarglist.add(scale_uncval_mu_dn);
  }
  else if (channel==k4e){
    scalemeanthetaarglist.add(scale_uncvar_e);
    scalemeanfcnarglist.add(scale_uncval_e_up);
    scalemeanfcnarglist.add(scale_uncval_e_dn);
  }
  AsymQuad scale_uncval(prefix + "final_CB_CMS_scale_em_AsymQuad", "", scalemeanfcnarglist, scalemeanthetaarglist, 1., 2);
  TString strscalemeanFormula; RooArgList scalemeanarglist;
  scalemeanarglist.add(CB_parameter_list.at(0));
  scalemeanarglist.add(var_mtrue);
  scalemeanarglist.add(scale_uncval);
  strscalemeanFormula="(@0+@1)*(1.+@2)-@1"; // Until a new procedure is found, keep var_mtrue as part of the scale unc. definition

  RooRealVar res_uncvar_e("CMS_res_e", "CMS_res_e", 0, -7, 7);
  RooConstVar res_uncval_e_up(prefix + "final_CB_CMS_res_eUp", "", 1.2);
  RooConstVar res_uncval_e_dn(prefix + "final_CB_CMS_res_eDown", "", 1./1.2);
  AsymPow res_uncval_e(prefix + "final_CB_CMS_res_e_AsymPow", "", res_uncval_e_dn, res_uncval_e_up, res_uncvar_e);
  RooRealVar res_uncvar_mu("CMS_res_m", "CMS_res_m", 0, -7, 7);
  RooConstVar res_uncval_mu_up(prefix + "final_CB_CMS_res_mUp", "", 1.2);
  RooConstVar res_uncval_mu_dn(prefix + "final_CB_CMS_res_mDown", "", 1./1.2);
  AsymPow res_uncval_mu(prefix + "final_CB_CMS_res_m_AsymPow", "", res_uncval_mu_dn, res_uncval_mu_up, res_uncvar_mu);
  TString strreswidthFormula; RooArgList reswidtharglist;
  reswidtharglist.add(CB_parameter_list.at(1));
  if (channel==k2e2mu){
    strreswidthFormula="@0*@1*@2";
    reswidtharglist.add(res_uncval_e);
    reswidtharglist.add(res_uncval_mu);
  }
  else if (channel==k4mu){
    strreswidthFormula="@0*@1";
    reswidtharglist.add(res_uncval_mu);
  }
  else if (channel==k4e){
    strreswidthFormula="@0*@1";
    reswidtharglist.add(res_uncval_e);
  }

  RooFormulaVar incl_CB_mean(prefix + "final_CB_mean", strscalemeanFormula, scalemeanarglist);
  RooFormulaVar incl_CB_width(prefix + "final_CB_width", strreswidthFormula, reswidtharglist);
  RooFormulaVar incl_CB_alpha1(prefix + "final_CB_alpha1", "max(@0,0.001)", RooArgList(CB_parameter_list.at(2)));
  RooFormulaVar incl_CB_alpha2(prefix + "final_CB_alpha2", "max(@0,0.001)", RooArgList(CB_parameter_list.at(3)));
  RooFormulaVar incl_CB_n1(prefix + "final_CB_n1", "max(@0,0.001)", RooArgList(CB_parameter_list.at(4)));
  RooFormulaVar incl_CB_n2(prefix + "final_CB_n2", "max(@0,0.001)", RooArgList(CB_parameter_list.at(5)));
  RooDoubleCB incl_pdf(
    prefix + "final_CB", prefix + "final_CB",
    var_mreco, var_mtrue,
    incl_CB_mean, incl_CB_width,
    incl_CB_alpha1, incl_CB_n1,
    incl_CB_alpha2, incl_CB_n2
  );

  // Fit all parameters simultaneously
  res_uncvar_e.setConstant(true);
  res_uncvar_mu.setConstant(true);
  scale_uncvar_e.setConstant(true);
  scale_uncvar_mu.setConstant(true);
  unsigned int nfits=5;
  for (unsigned int ifit=0; ifit<nfits; ifit++){
    RooLinkedList cmdList;
    RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
    //RooCmdArg condObsArg = RooFit::ConditionalObservables(conditionals); cmdList.Add((TObject*) &condObsArg);
    RooCmdArg sumw2Arg = RooFit::SumW2Error(true); cmdList.Add((TObject*) &sumw2Arg);
    RooCmdArg hesseArg = RooFit::Hesse(true); cmdList.Add((TObject*) &hesseArg);
    RooCmdArg minimizerStrategyArg = RooFit::Strategy(2); cmdList.Add((TObject*) &minimizerStrategyArg);
    // Misc. options
    RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
    RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
    //RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

    RooFitResult* fitResult=incl_pdf.fitTo(data, cmdList);
    if (fitResult){
      int fitStatus = fitResult->status();
      cout << "Fit status is " << fitStatus << endl;
      cout << "Fit properties:" << endl;
      fitResult->Print("v");
    }
    delete fitResult;
  }
  res_uncvar_e.setConstant(false);
  res_uncvar_mu.setConstant(false);
  scale_uncvar_e.setConstant(false);
  scale_uncvar_mu.setConstant(false);

  {
    RooPlot incl_plot(var_mreco, var_mreco.getMin(), var_mreco.getMax(), 80);
    data.plotOn(&incl_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2));
    incl_pdf.plotOn(&incl_plot, LineColor(kRed), LineWidth(2));

    TCanvas can(Form("ZZMass_%.0f_%.0f", var_mreco.getMin(), var_mreco.getMax()), "");
    incl_plot.Draw();
    can.SaveAs(Form("%s_%s.pdf", canvasnamecore.Data(), can.GetName()));
    can.SaveAs(Form("%s_%s.png", canvasnamecore.Data(), can.GetName()));
    can.Close();
  }

  // Set all variables constant
  for (auto& var:CB_parameter_list) var.setConstant(true);
  var_mtrue.SetName("MH");
  var_mreco.SetName("mass");
  incl_pdf.SetName("ResolutionModel");
  var_mtrue.SetTitle("MH");
  var_mreco.SetTitle("mass");
  incl_pdf.SetTitle("ResolutionModel");

  RooWorkspace w("w", "");
  w.import(incl_pdf, RecycleConflictNodes());
  //w.import(var_mreco, RenameVariable(var_mreco.GetName(), "newmass"));
  //RooAbsArg* pdfnew=w.factory("EDIT::ResolutionModelCopy(ResolutionModel, mass=newmass)");
  //w.import(*pdfnew, RecycleConflictNodes());
  foutput->WriteTObject(&w);
  w.pdf(incl_pdf.GetName())->Print("v");
  foutput->Close();
  MELAout.close();
}
