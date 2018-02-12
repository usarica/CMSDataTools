#include "common_includes.h"
#include "TemplatesEventAnalyzer.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

void fixWeights(std::vector<SimpleEntry>& index);

void testGenHMassDistributions(
  const Channel channel, const Category category,
  const ProcessHandler& processA, unsigned int hypoA, unsigned int AChypoA,
  const ProcessHandler& processB, unsigned int hypoB, unsigned int AChypoB,
  const bool applyRewgt=true
){
  const SystematicVariationTypes syst = sNominal;
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  // Get list of samples
  vector<TString> KDnames; KDnames.push_back("GenHMass");
  const ProcessHandler* process[2]={ &processA, &processB };
  unsigned int hypo[2]={ hypoA, hypoB };
  vector<TH2F*> hKD[2]; for (unsigned int is=0; is<2; is++) hKD[is].assign(KDnames.size(), nullptr);
  vector<vector<TString>> strSamples[2];
  vector<TString> melawgtvars[2];

  TString canvasname = "TrueMasscomparison_";
  ExtendedBinning binning_ZZMass("ZZMass");
  binning_ZZMass.addBinBoundary(70);
  binning_ZZMass.addBinBoundary(100);
  binning_ZZMass.addBinBoundary(120);
  binning_ZZMass.addBinBoundary(135);
  binning_ZZMass.addBinBoundary(155);
  binning_ZZMass.addBinBoundary(180);
  binning_ZZMass.addBinBoundary(220);
  binning_ZZMass.addBinBoundary(350);
  binning_ZZMass.addBinBoundary(650);
  binning_ZZMass.addBinBoundary(1100);
  binning_ZZMass.addBinBoundary(theSqrts*1000.);
  TString hnamecore[2];
  TString proclabel[2];
  for (unsigned int is=0; is<2; is++){
    if (process[is]->getProcessType()==ProcessHandler::kGG){
      hnamecore[is]=((const GGProcessHandler*) process[is])->getOutputTreeName(((const GGProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]));
      proclabel[is]=((const GGProcessHandler*) process[is])->getProcessLabel(((const GGProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]), kSM);
    }
    else if (process[is]->getProcessType()==ProcessHandler::kVV){
      hnamecore[is]=((const VVProcessHandler*) process[is])->getOutputTreeName(((const VVProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]));
      proclabel[is]=((const VVProcessHandler*) process[is])->getProcessLabel(((const VVProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]), kSM);
    }
    else if (process[is]->getProcessType()==ProcessHandler::kQQBkg){
      hnamecore[is]=((const QQBkgProcessHandler*) process[is])->getOutputTreeName();
      proclabel[is]=((const QQBkgProcessHandler*) process[is])->getProcessLabel();
    }
    for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
      TString hname = hnamecore[is] + "_ " + KDnames.at(ikd);
      ExtendedBinning binning_KD("GenHMass");
      binning_KD.addBinBoundary(binning_ZZMass.getBinHighEdge(binning_ZZMass.getNbins()-1));
      for (unsigned int bin=0; bin<binning_ZZMass.getNbins(); bin++){
        double binwidth=binning_ZZMass.getBinHighEdge(bin)-binning_ZZMass.getBinLowEdge(bin);
        unsigned int ndiv=20;
        if (binwidth>1000.) ndiv=40;
        for (unsigned int i=0; i<ndiv; i++) binning_KD.addBinBoundary(binning_ZZMass.getBinLowEdge(bin) + double(i)*binwidth/double(ndiv));
      }
      hKD[is][ikd] = new TH2F(hname, proclabel[is], binning_ZZMass.getNbins(), binning_ZZMass.getBinning(), binning_KD.getNbins(), binning_KD.getBinning());
    }
  }

  if (process[0]->getProcessType()==ProcessHandler::kGG){
    vector<TString> strSampleIdentifiers;
    strSampleIdentifiers.push_back("gg_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[0].push_back(ss);
    }

    melawgtvars[0].push_back(((const GGProcessHandler*) process[0])->getMELAHypothesisWeight(((const GGProcessHandler*) process[0])->castIntToHypothesisType(hypo[0]), ACHypothesisHelpers::ACHypothesis(AChypoA)));
    melawgtvars[0].push_back("p_Gen_CPStoBWPropRewgt");
    melawgtvars[0].push_back("KFactor_QCD_ggZZ_Nominal");
  }
  else if (process[0]->getProcessType()==ProcessHandler::kVV){
    vector<TString> strSampleIdentifiers;
    if (category==HadVHTagged){
      strSampleIdentifiers.push_back("WminusH_Sig_POWHEG");
      strSampleIdentifiers.push_back("WplusH_Sig_POWHEG");
      strSampleIdentifiers.push_back("ZH_Sig_POWHEG");
    }
    else strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[0].push_back(ss);
    }

    melawgtvars[0].push_back(((const VVProcessHandler*) process[0])->getMELAHypothesisWeight(((const VVProcessHandler*) process[0])->castIntToHypothesisType(hypo[0]), ACHypothesisHelpers::ACHypothesis(AChypoA)));
    melawgtvars[0].push_back("p_Gen_CPStoBWPropRewgt");
  }
  else if (process[0]->getProcessType()==ProcessHandler::kQQBkg){
    vector<TString> strSampleIdentifiers;
    strSampleIdentifiers.push_back("qq_Bkg_Combined");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[0].push_back(ss);
    }

    melawgtvars[0].push_back("KFactor_QCD_qqZZ_M");
    melawgtvars[0].push_back("KFactor_EW_qqZZ");
  }

  if (process[1]->getProcessType()==ProcessHandler::kGG){
    vector<TString> strSampleIdentifiers;
    strSampleIdentifiers.push_back("gg_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[1].push_back(ss);
    }

    melawgtvars[1].push_back(((const GGProcessHandler*) process[1])->getMELAHypothesisWeight(((const GGProcessHandler*) process[1])->castIntToHypothesisType(hypo[1]), ACHypothesisHelpers::ACHypothesis(AChypoB)));
    melawgtvars[1].push_back("p_Gen_CPStoBWPropRewgt");
    melawgtvars[1].push_back("KFactor_QCD_ggZZ_Nominal");
  }
  else if (process[1]->getProcessType()==ProcessHandler::kVV){
    vector<TString> strSampleIdentifiers;
    if (category==HadVHTagged){
      strSampleIdentifiers.push_back("WminusH_Sig_POWHEG");
      strSampleIdentifiers.push_back("WplusH_Sig_POWHEG");
      strSampleIdentifiers.push_back("ZH_Sig_POWHEG");
    }
    else strSampleIdentifiers.push_back("VBF_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[1].push_back(ss);
    }

    melawgtvars[1].push_back(((const VVProcessHandler*) process[1])->getMELAHypothesisWeight(((const VVProcessHandler*) process[1])->castIntToHypothesisType(hypo[1]), ACHypothesisHelpers::ACHypothesis(AChypoB)));
    melawgtvars[1].push_back("p_Gen_CPStoBWPropRewgt");
  }
  else if (process[1]->getProcessType()==ProcessHandler::kQQBkg){
    vector<TString> strSampleIdentifiers;
    strSampleIdentifiers.push_back("qq_Bkg_Combined");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[1].push_back(ss);
    }

    melawgtvars[1].push_back("KFactor_QCD_qqZZ_M");
    melawgtvars[1].push_back("KFactor_EW_qqZZ");
  }

  // Register the discriminants
  vector<KDspecs> KDlist;
  if (category!=Inclusive) getCategorizationDiscriminants(syst, KDlist);

  // Open the output file
  TFile* foutput = TFile::Open("tmp.root", "recreate");

  // Get the CJLST set
  for (unsigned int is=0; is<2; is++){
    for (vector<TString> slist:strSamples[is]){
      CJLSTSet* theSampleSet = new CJLSTSet(slist);
      // Book common variables
      theSampleSet->bookXS(); // "xsec"
      theSampleSet->bookOverallEventWgt(); // Gen weights "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
      for (auto& tree:theSampleSet->getCJLSTTreeList()){
        // Book common variables needed for analysis
        tree->bookBranch<float>("GenHMass", 0);
        tree->bookBranch<float>("ZZMass", -1);
        tree->bookBranch<short>("Z1Flav", 0);
        tree->bookBranch<short>("Z2Flav", 0);
        // Variables for MELA reweighting
        for (TString const& wgtvar:melawgtvars[is]) tree->bookBranch<float>(wgtvar, 0);
        // Variables for KDs
        for (auto& KD:KDlist){ for (auto& v:KD.KDvars) tree->bookBranch<float>(v, 0); }
        tree->silenceUnused(); // Will no longer book another branch
      }
      ExtendedBinning GenHMassBinning = ReweightingBuilder::getTrueMassBinning(theSampleSet->getCJLSTTreeList());

      vector<TString> strReweightingWeights;
      for (TString const& s:melawgtvars[is]) strReweightingWeights.push_back(s);
      SampleHelpers::addXsecBranchNames(strReweightingWeights);
      ReweightingBuilder* melarewgtBuilder = nullptr;
      if (applyRewgt || process[is]->getProcessType()==ProcessHandler::kQQBkg){
        melarewgtBuilder = new ReweightingBuilder(strReweightingWeights, getSimpleWeight);
        melarewgtBuilder->rejectNegativeWeights(true);
        melarewgtBuilder->setDivideByNSample(true);
        melarewgtBuilder->setWeightBinning(GenHMassBinning);
        for (auto& tree:theSampleSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(
          tree,
          (process[is]->getProcessType()==ProcessHandler::kQQBkg ? -1 : 0.999),
          (process[is]->getProcessType()==ProcessHandler::kQQBkg ? 0 : 250)
        );
      }

      vector<SimpleEntry> events;
      TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
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
      theAnalyzer.addReweightingBuilder((process[is]->getProcessType()==ProcessHandler::kQQBkg ? "RegularRewgt" : "MELARewgt"), melarewgtBuilder);
      // Loop
      theAnalyzer.setExternalProductList(&events);
      theAnalyzer.loop(true, false, true);
      fixWeights(events);

      delete melarewgtBuilder;
      delete theSampleSet;

      MELAout << "There are " << events.size() << " products" << endl;

      for (SimpleEntry& ev:events){
        float ZZMass, weight;
        ev.getNamedVal("ZZMass", ZZMass);
        ev.getNamedVal("weight", weight);
        for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
          float KDval;
          ev.getNamedVal(KDnames.at(ikd), KDval);
          hKD[is][ikd]->Fill(ZZMass, KDval, weight);
        }
      }
    }

    // Conditionalize the histogram in ZZMass
    for (unsigned int ikd=0; ikd<KDnames.size(); ikd++) conditionalizeHistogram(hKD[is][ikd], 0);
  }

  // Plot
  gStyle->SetOptStat(0);
  {
    int colors[100];
    Double_t Red[]    ={ 0.3, 0.4, 1.0 };
    Double_t Green[]  ={ 0.0, 1.0, 0.8 };
    Double_t Blue[]   ={ 1.0, 0.0, 0.3 };
    Double_t Length[] ={ 0.00, 0.50, 1.00 };
    int FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    if (FI<0) MELAout << "Failed to set color palette" << endl;
    else{
      for (unsigned int ic=0; ic<100; ic++) colors[ic] = FI+ic;
      gStyle->SetPalette(100, colors);
    }
    MELAout << "Ncolors: " << ncolors << endl;
  }

  for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
    TH2F* histo2D[2]={ hKD[0][ikd], hKD[1][ikd] };
    for (unsigned int bin=0; bin<binning_ZZMass.getNbins(); bin++){
      TString strZZMassRange = Form("ZZMass_%.0f_%.0f", binning_ZZMass.getBinLowEdge(bin), binning_ZZMass.getBinHighEdge(bin));
      TString strZZMassLabel = Form("(%.0f-%.0f GeV)", binning_ZZMass.getBinLowEdge(bin), binning_ZZMass.getBinHighEdge(bin));
      TString canvasname =
        hnamecore[0] + "_vs_" + hnamecore[1]
        + Form("_%s_%s_%s", strChannel.Data(), strCategory.Data(), strZZMassRange.Data());
      for (unsigned int is=0; is<2; is++){
        replaceString(canvasname, "T_", "");
        replaceString(canvasname, "_Tree", "");
      }
      if (!applyRewgt) canvasname += "_NoRewgt";

      TH1D* histo[2];
      double maxY=-1;
      for (unsigned int is=0; is<2; is++){
        TString htitle = TString(histo2D[is]->GetTitle()) + " " + strZZMassLabel;
        histo[is]=(TH1D*) histo2D[is]->ProjectionY(Form("%s_%s", histo2D[is]->GetName(), strZZMassRange.Data()), bin+1, bin+1);
        histo[is]->SetTitle(htitle);
        double integral = histo[is]->Integral(1, histo[is]->GetNbinsX());
        if (integral>0.) histo[is]->Scale(1./integral);
        for (int bin=1; bin<=histo[is]->GetNbinsX(); bin++) maxY = std::max(maxY, histo[is]->GetBinContent(bin));
      }

      TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 900, 800);
      canvas->cd();
      gStyle->SetOptStat(0);
      canvas->SetFillColor(0);
      canvas->SetBorderMode(0);
      canvas->SetBorderSize(2);
      canvas->SetTickx(1);
      canvas->SetTicky(1);
      canvas->SetLeftMargin(0.14);
      canvas->SetRightMargin(0.12);
      canvas->SetTopMargin(0.07);
      canvas->SetBottomMargin(0.13);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);
      canvas->SetFrameFillStyle(0);
      canvas->SetFrameBorderMode(0);

      TLegend* legend = new TLegend(0.90-0.37, 0.70-0.15/3.*2., 0.90, 0.70);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.03);
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(1);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);

      TPaveText* pt = new TPaveText(0.12, 0.93, 0.85, 1, "brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.045);
      TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      text = pt->AddText(0.155, 0.42, "#font[52]{Simulation}");
      text->SetTextSize(0.0315);
      TString strTitle;
      strTitle = Form("#font[42]{%i TeV}", theSqrts);
      text = pt->AddText(0.92, 0.45, strTitle);
      text->SetTextSize(0.0315);
      {
        const unsigned int ncolors = gStyle->GetNumberOfColors();
        const unsigned int stepsize = ncolors/2;
        for (unsigned int is=0; is<2; is++){
          int colorToUse = gStyle->GetColorPalette(is*stepsize);
          histo[is]->SetMarkerColor(colorToUse);
          histo[is]->SetLineColor(colorToUse);
          histo[is]->SetLineWidth(2);
          histo[is]->GetXaxis()->SetNdivisions(505);
          histo[is]->GetXaxis()->SetLabelFont(42);
          histo[is]->GetXaxis()->SetLabelOffset(0.007);
          histo[is]->GetXaxis()->SetLabelSize(0.04);
          histo[is]->GetXaxis()->SetTitleSize(0.06);
          histo[is]->GetXaxis()->SetTitleOffset(0.9);
          histo[is]->GetXaxis()->SetTitleFont(42);
          histo[is]->GetYaxis()->SetNdivisions(505);
          histo[is]->GetYaxis()->SetLabelFont(42);
          histo[is]->GetYaxis()->SetLabelOffset(0.007);
          histo[is]->GetYaxis()->SetLabelSize(0.04);
          histo[is]->GetYaxis()->SetTitleSize(0.06);
          histo[is]->GetYaxis()->SetTitleOffset(1);
          histo[is]->GetYaxis()->SetTitleFont(42);
          histo[is]->GetYaxis()->SetRangeUser(0, maxY*1.2);
          if (binning_ZZMass.getBinLowEdge(bin)<350) histo[is]->GetXaxis()->SetRangeUser(binning_ZZMass.getBinLowEdge(0), 349.9);
          else histo[is]->GetXaxis()->SetRangeUser(350, 4000);
          legend->AddEntry(histo[is], histo[is]->GetTitle(), "l");
          histo[is]->SetTitle("");
          if (is==0) histo[is]->Draw("hist");
          else histo[is]->Draw("histsame");
        }
      }
      legend->Draw("same");
      pt->Draw();

      canvas->RedrawAxis();
      canvas->Modified();
      canvas->Update();
      canvas->SaveAs(Form("%s%s", canvasname.Data(), ".pdf"));
      foutput->WriteTObject(canvas);
      delete pt;
      delete legend;
      canvas->Close();
      for (unsigned int is=0; is<2; is++) delete histo[is];
    }
  }

  for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
    for (unsigned int is=0; is<2; is++){
      foutput->WriteTObject(hKD[is][ikd]);
      delete hKD[is][ikd];
    }
  }
  foutput->Close();

  for (auto& KD:KDlist) delete KD.KD;
}

void fixWeights(std::vector<SimpleEntry>& index){
  const unsigned int nMarginalMax = 100;
  const unsigned int nMarginalMaxMult = 1000;
  const float nMarginalMaxFrac = 1./static_cast<float const>(nMarginalMaxMult);
  const unsigned int countThreshold=nMarginalMaxMult*nMarginalMax;

  int const nbinsraw = (1000*theSqrts-70)/5;
  TH1F* hmass = new TH1F("hmass", "", nbinsraw, 70, 13000);
  // Initial loop over the tree to count the events in each bin
  for (SimpleEntry const& product:index){
    float ZZMass;
    product.getNamedVal("ZZMass", ZZMass);
    hmass->Fill(ZZMass); // Do not use weight; just count
  }
                                                                           // Determine the final binning to set the weight thresholds
  MELAout
    << "fixWeights: "
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
    << "fixWeights: "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;
  // These lines guarantee count>countThreshold in every bin
  if (counts.at(0)<countThreshold){
    counts.at(1) += counts.at(0);
    counts.erase(counts.begin());
    binning.removeBinLowEdge(1);
  }
  MELAout
    << "fixWeights: "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;

  // Collect the count*nMarginalMaxFrac events with highest weights
  MELAout
    << "fixWeights: "
    << "Collecting the count*" << nMarginalMaxFrac << " events with highest weights in " << binning.getNbins() << " bins"
    << endl;
  vector<vector<float>> wgtcollList;
  wgtcollList.assign(binning.getNbins(), vector<float>());
  for (SimpleEntry const& product:index){
    float ZZMass, weight;
    product.getNamedVal("ZZMass", ZZMass);
    product.getNamedVal("weight", weight);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int) binning.getNbins()){
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
    << "fixWeights: "
    << "Determining the weight thresholds"
    << endl;
  vector<float> wgtThresholds; wgtThresholds.reserve(binning.getNbins());
  for (auto const& wgtcoll:wgtcollList){
    unsigned int ns=wgtcoll.size();
    float threshold=0.5*(wgtcoll.at(ns-1)+wgtcoll.at(ns-2));
    if (wgtcoll.front()*5.<threshold) threshold=wgtcoll.front();
    else MELAout
      << "fixWeights: "
      << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
      << endl;
    wgtThresholds.push_back(threshold);
  }

  // Fix the weights
  for (SimpleEntry& product:index){
    float ZZMass, weight;
    product.getNamedVal("ZZMass", ZZMass);
    product.getNamedVal("weight", weight);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int) binning.getNbins() && fabs(weight)>wgtThresholds.at(bin)) weight = pow(wgtThresholds.at(bin), 2)/weight;
    product.setNamedVal("weight", weight);
  }
}
