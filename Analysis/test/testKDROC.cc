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

void testKDROC(
  const TString strKD,
  const Channel channel, const Category category,
  float ZZMassbegin, float ZZMassend,
  const ProcessHandler& processA, unsigned int hypoA,
  const ProcessHandler& processB, unsigned int hypoB
){
  const SystematicVariationTypes syst = sNominal;
  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  vector<TString> KDnames;
  splitOptionRecursive(strKD, KDnames, ',');
  MELAout << "Processing KDs " << KDnames << endl;

  // Get list of samples
  const ProcessHandler* process[2]={ &processA, &processB };
  unsigned int hypo[2]={ hypoA, hypoB };
  vector<TH1F*> hKD[2]; for (unsigned int is=0; is<2; is++) hKD[is].assign(KDnames.size(), nullptr);
  vector<vector<TString>> strSamples[2];
  vector<TString> melawgtvars[2];

  if (process[0]->getProcessType()==ProcessHandler::kGG){
    vector<TString> strSampleIdentifiers;
    if (strKD.Contains("DbkgjjEWQCD")) strSampleIdentifiers.push_back("gg_Sig_POWHEG_MINLO");
    else strSampleIdentifiers.push_back("gg_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[0].push_back(ss);
    }

    melawgtvars[0].push_back(((const GGProcessHandler*) process[0])->getMELAHypothesisWeight(((const GGProcessHandler*) process[0])->castIntToHypothesisType(hypo[0]), kSM));
    melawgtvars[0].push_back("p_Gen_CPStoBWPropRewgt");
    //melawgtvars[0].push_back("KFactor_QCD_ggZZ_Nominal");
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

    melawgtvars[0].push_back(((const VVProcessHandler*) process[0])->getMELAHypothesisWeight(((const VVProcessHandler*) process[0])->castIntToHypothesisType(hypo[0]), kSM));
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
    if (strKD.Contains("DbkgjjEWQCD")) strSampleIdentifiers.push_back("gg_Sig_POWHEG_MINLO");
    else strSampleIdentifiers.push_back("gg_Sig_POWHEG");
    for (TString& s:strSampleIdentifiers){
      vector<TString> ss;
      vector<TString> slist; slist.push_back(s);
      getSamplesList(theSqrts, slist, ss);
      strSamples[1].push_back(ss);
    }

    melawgtvars[1].push_back(((const GGProcessHandler*) process[1])->getMELAHypothesisWeight(((const GGProcessHandler*) process[1])->castIntToHypothesisType(hypo[1]), kSM));
    melawgtvars[1].push_back("p_Gen_CPStoBWPropRewgt");
    //melawgtvars[1].push_back("KFactor_QCD_ggZZ_Nominal");
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

    melawgtvars[1].push_back(((const VVProcessHandler*) process[1])->getMELAHypothesisWeight(((const VVProcessHandler*) process[1])->castIntToHypothesisType(hypo[1]), kSM));
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
  vector<TString> strExtraCatVars_short;
  getLikelihoodDiscriminants(channel, category, syst, KDlist);
  if (category!=Inclusive){
    getCategorizationDiscriminants(syst, KDlist);
    getExtraCategorizationVariables<short>(globalCategorizationScheme, syst, strExtraCatVars_short);
  }

  // Open the output file
  TFile* foutput = TFile::Open("tmp.root", "recreate");
  TString canvasname = "ROC_";

  // Construct histograms
  for (unsigned int is=0; is<2; is++){
    TString hnamecore;
    TString proclabel;
    if (process[is]->getProcessType()==ProcessHandler::kGG){
      hnamecore=((const GGProcessHandler*) process[is])->getOutputTreeName(((const GGProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]));
      proclabel=((const GGProcessHandler*) process[is])->getProcessLabel(((const GGProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]), kSM);
    }
    else if (process[is]->getProcessType()==ProcessHandler::kVV){
      hnamecore=((const VVProcessHandler*) process[is])->getOutputTreeName(((const VVProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]));
      proclabel=((const VVProcessHandler*) process[is])->getProcessLabel(((const VVProcessHandler*) process[is])->castIntToHypothesisType(hypo[is]), kSM);
    }
    else if (process[is]->getProcessType()==ProcessHandler::kQQBkg){
      hnamecore=((const QQBkgProcessHandler*) process[is])->getOutputTreeName();
      proclabel=((const QQBkgProcessHandler*) process[is])->getProcessLabel();
    }
    for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
      TString hname = hnamecore + "_ " + KDnames.at(ikd);
      hKD[is][ikd] = new TH1F(hname, proclabel, 100, 0, 1);
    }
    canvasname += hnamecore;
    if (is==0) canvasname += "_vs_";
  }

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
        // Extra categorization variables
        for (auto& s:strExtraCatVars_short) tree->bookBranch<short>(s, -1);
        tree->silenceUnused(); // Will no longer book another branch
      }
      ExtendedBinning GenHMassBinning("GenHMass");
      for (unsigned int is=0; is<theSampleSet->getCJLSTTreeList().size()-1; is++){
        if (theSampleSet->getCJLSTTreeList().at(is)->MHVal>0. && theSampleSet->getCJLSTTreeList().at(is+1)->MHVal>0.){
          float boundary = (theSampleSet->getCJLSTTreeList().at(is)->MHVal + theSampleSet->getCJLSTTreeList().at(is+1)->MHVal)/2.;
          GenHMassBinning.addBinBoundary(boundary);
        }
      }
      GenHMassBinning.addBinBoundary(0);
      GenHMassBinning.addBinBoundary(theSqrts*1000.);

      vector<TString> strReweightingWeights;
      for (TString const& s:melawgtvars[is]) strReweightingWeights.push_back(s);
      SampleHelpers::addXsecBranchNames(strReweightingWeights);
      ReweightingBuilder* melarewgtBuilder = new ReweightingBuilder(strReweightingWeights, getSimpleWeight);
      melarewgtBuilder->rejectNegativeWeights(true);
      melarewgtBuilder->setDivideByNSample(true);
      melarewgtBuilder->setWeightBinning(GenHMassBinning);
      for (auto& tree:theSampleSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(
        tree,
        (process[is]->getProcessType()==ProcessHandler::kQQBkg ? -1 : 0.999),
        (process[is]->getProcessType()==ProcessHandler::kQQBkg ? 0 : 250)
      );

      vector<SimpleEntry> events;
      TemplatesEventAnalyzer theAnalyzer(theSampleSet, channel, category);
      theAnalyzer.addMassWindow(std::pair<float, float>(ZZMassbegin, ZZMassend));
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
      // Add extra categorization variables
      for (auto& s:strExtraCatVars_short) theAnalyzer.addConsumed<short>(s);
      // Add reweighting builders
      theAnalyzer.addReweightingBuilder((process[is]->getProcessType()==ProcessHandler::kQQBkg ? "RegularRewgt" : "MELARewgt"), melarewgtBuilder);
      // Loop
      theAnalyzer.setExternalProductList(&events);
      theAnalyzer.loop(true, false, true);

      delete melarewgtBuilder;
      delete theSampleSet;

      MELAout << "There are " << events.size() << " products" << endl;

      for (SimpleEntry& ev:events){
        float weight;
        ev.getNamedVal("weight", weight);
        for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
          float KDval;
          ev.getNamedVal(KDnames.at(ikd), KDval);
          hKD[is][ikd]->Fill(KDval, weight);
        }
      }
    }
  }

  vector<TGraph*> grlist;
  for (unsigned int ikd=0; ikd<KDnames.size(); ikd++){
    TGraph* gr = createROCFromDistributions(hKD[0][ikd], hKD[1][ikd], TString("ROC_")+KDnames.at(ikd));
    gr->SetTitle(getKDLabel(KDnames.at(ikd)));
    foutput->WriteTObject(gr);
    grlist.push_back(gr);
    for (unsigned int is=0; is<2; is++){
      foutput->WriteTObject(hKD[is][ikd]);
      delete hKD[is][ikd];
    }
  }

  // Plot
  canvasname += Form("_%s_%s_ZZMass_%.0f_%.0f", strChannel.Data(), strCategory.Data(), ZZMassbegin, ZZMassend);
  for (unsigned int is=0; is<2; is++){
    replaceString(canvasname, "T_", "");
    replaceString(canvasname, "_Tree", "");
  }
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

  {
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

    TLegend* legend = new TLegend(0.90-0.37, 0.70-0.15/3.*float(grlist.size()), 0.90, 0.70);
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
      unsigned int igr=0;
      {
        const unsigned int ncolors = gStyle->GetNumberOfColors();
        const unsigned int stepsize = ncolors/grlist.size();
        for (; igr<grlist.size(); igr++){
          int colorToUse = gStyle->GetColorPalette(igr*stepsize);
          grlist.at(igr)->SetMarkerColor(colorToUse);
          grlist.at(igr)->SetLineColor(colorToUse);
          grlist.at(igr)->SetLineWidth(2);
        }
      }
      igr=0;
      for (TGraph*& gr:grlist){
        gr->GetXaxis()->SetNdivisions(505);
        gr->GetXaxis()->SetLabelFont(42);
        gr->GetXaxis()->SetLabelOffset(0.007);
        gr->GetXaxis()->SetLabelSize(0.04);
        gr->GetXaxis()->SetTitleSize(0.06);
        gr->GetXaxis()->SetTitleOffset(0.9);
        gr->GetXaxis()->SetTitleFont(42);
        gr->GetYaxis()->SetNdivisions(505);
        gr->GetYaxis()->SetLabelFont(42);
        gr->GetYaxis()->SetLabelOffset(0.007);
        gr->GetYaxis()->SetLabelSize(0.04);
        gr->GetYaxis()->SetTitleSize(0.06);
        gr->GetYaxis()->SetTitleOffset(1);
        gr->GetYaxis()->SetTitleFont(42);
        legend->AddEntry(gr, gr->GetTitle(), "l");
        gr->SetTitle("");

        if (igr==0) gr->Draw("al");
        else gr->Draw("lsame");
        igr++;
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
  }
  for (TGraph*& gr:grlist) delete gr;
  foutput->Close();

  for (auto& KD:KDlist) delete KD.KD;
}

