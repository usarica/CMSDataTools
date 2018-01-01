#ifndef PLOTPROCESSCHECKSTAGE_H
#define PLOTPROCESSCHECKSTAGE_H

#include "common_includes.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


// Constants to affect the template code
#ifndef outputdir_def
#define outputdir_def
const TString user_output_dir = "output/";
#endif

void plotTH1Fs(TString coutdir, TString canvasname, std::vector<TH1F*>& histolist);
void plotTH2Fs(TString coutdir, std::vector<TH2F*>& histolist);
 
void plotProcessCheckStage(
  const Channel channel, const Category category, ACHypothesisHelpers::ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  if (channel==NChannels) return;
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
  if (!systematicAllowed(category, thePerProcessHandle->getProcessType(), syst)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strSystematics = getSystematicsName(syst);
  const TString strACHypo = getACHypothesisName(hypo);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Plots/" + strACHypo + "/" + strChannel + "/" + strCategory + "/" + strGenerator + "/" + strSystematics + "/";

  gSystem->Exec("mkdir -p " + coutput_common);
  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants_Stage%i",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data(),
    istage
  );
  TString OUTPUT_LOG_NAME = INPUT_NAME;
  INPUT_NAME += ".root";
  OUTPUT_LOG_NAME += "_plot.log";
  TString cinput = cinput_common + INPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  MELAout.open(coutput_log.Data());
  MELAout << "Opened log file " << coutput_log << endl;

#ifdef checkstage_def
  if (gSystem->AccessPathName(cinput)) checkstagefcn(channel, category, hypo, syst, istage, fixedDate);
#endif
  TFile* finput = TFile::Open(cinput, "read");
  if (!finput) return;
  MELAout << "Opened file " << cinput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  TFile* finputNominal=nullptr;
  if (syst!=sNominal){
    TString INPUT_NOMINAL_NAME = Form(
      "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants_Stage%i%s",
      strChannel.Data(), strCategory.Data(),
      thePerProcessHandle->getProcessName().Data(),
      getSystematicsName(sNominal).Data(),
      strGenerator.Data(),
      getACHypothesisName(hypo).Data(),
      istage,
      ".root"
    );
    TString cinputNominal = cinput_common + INPUT_NOMINAL_NAME;
#ifdef checkstage_def
    if (gSystem->AccessPathName(cinputNominal)) checkstagefcn(channel, category, hypo, sNominal, istage, fixedDate);
#endif
    finputNominal = TFile::Open(cinputNominal, "read");
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

  finput->ls();

  vector<TH1F*> h1dlist;
  vector<TH2F*> h2dlist;
  MELAout << "Extracting 1D histograms..." << endl;
  HelperFunctions::extractHistogramsFromDirectory(finput, h1dlist);
  MELAout << "N(1D histograms) = " << h1dlist.size() << endl;
  MELAout << "Extracting 2D histograms..." << endl;
  HelperFunctions::extractHistogramsFromDirectory(finput, h2dlist);
  MELAout << "N(2D histograms) = " << h2dlist.size() << endl;
  plotTH2Fs(coutput_common, h2dlist);
  plotTH1Fs(coutput_common, Form("c_M4lDistribution_%s", thePerProcessHandle->getProcessName().Data()), h1dlist);

  if (finputNominal){
    vector<TH1F*> h1dnomlist;
    vector<TH2F*> h2dnomlist;
    MELAout << "Extracting nominal 1D histograms..." << endl;
    HelperFunctions::extractHistogramsFromDirectory(finput, h1dnomlist);
    MELAout << "N(Nom. 1D histograms) = " << h1dnomlist.size() << endl;
    MELAout << "Extracting nominal 2D histograms..." << endl;
    HelperFunctions::extractHistogramsFromDirectory(finput, h2dnomlist);
    MELAout << "N(Nom. 2D histograms) = " << h2dnomlist.size() << endl;

    if (h1dnomlist.size()==h1dlist.size()){
      for (unsigned int ih=0; ih<h1dlist.size(); ih++){
        h1dlist.at(ih)->Divide(h1dnomlist.at(ih));
        h1dlist.at(ih)->SetName(Form("%s%s", h1dlist.at(ih)->GetName(), "_RatioToNominal"));
      }
      for (unsigned int ih=0; ih<h2dlist.size(); ih++){
        h2dlist.at(ih)->Divide(h2dnomlist.at(ih));
        h2dlist.at(ih)->SetName(Form("%s%s", h2dlist.at(ih)->GetName(), "_RatioToNominal"));
      }

      plotTH2Fs(coutput_common, h2dlist);
      plotTH1Fs(coutput_common, Form("c_M4lDistribution_%s_RatioToNominal", thePerProcessHandle->getProcessName().Data()), h1dlist);
    }

    finputNominal->Close();
  }
  finput->Close();
  MELAout.close();
}

void plotTH2Fs(TString coutdir, std::vector<TH2F*>& histolist){
  for (TH2F*& hh:histolist){
    if (hh->GetName()[0]=='T') continue;

    TString canvasname = Form("c_%s", hh->GetName());
    HelperFunctions::replaceString(canvasname, "h2D_T_", "");
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

    TPaveText* pt = new TPaveText(0.12, 0.93, 0.85, 1, "brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);
    TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
    text->SetTextSize(0.044);
    text = pt->AddText(0.155, 0.42, Form("#font[52]{Simulation %s}", hh->GetTitle()));
    text->SetTextSize(0.0315);
    TString strTitle;
    strTitle = Form("#font[42]{%i TeV}", theSqrts);
    text = pt->AddText(0.92, 0.45, strTitle);
    text->SetTextSize(0.0315);

    hh->SetTitle("");
    hh->GetXaxis()->SetRangeUser(220, hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX()));
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetLabelFont(42);
    hh->GetXaxis()->SetLabelOffset(0.007);
    hh->GetXaxis()->SetLabelSize(0.04);
    hh->GetXaxis()->SetTitleSize(0.06);
    hh->GetXaxis()->SetTitleOffset(0.9);
    hh->GetXaxis()->SetTitleFont(42);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetLabelFont(42);
    hh->GetYaxis()->SetLabelOffset(0.007);
    hh->GetYaxis()->SetLabelSize(0.04);
    hh->GetYaxis()->SetTitleSize(0.06);
    hh->GetYaxis()->SetTitleOffset(1);
    hh->GetYaxis()->SetTitleFont(42);

    hh->Draw("colz");
    pt->Draw();

    canvas->RedrawAxis();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".png"));
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".pdf"));
    delete pt;
    canvas->Close();
  }

}

void plotTH1Fs(TString coutdir, TString canvasname,  std::vector<TH1F*>& histolist){
  HelperFunctions::replaceString(canvasname, "h1D_T_", "");
  {
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    const unsigned int stepsize = ncolors/histolist.size();
    MELAout << "plotTH1Fs: Ncolors: " << ncolors << ", step size " << stepsize << endl;
    for (unsigned int ih=0; ih<histolist.size(); ih++){
      int colorToUse = gStyle->GetColorPalette(ih*stepsize);
      histolist.at(ih)->SetMarkerColor(colorToUse);
      histolist.at(ih)->SetLineColor(colorToUse);
      histolist.at(ih)->SetLineWidth(2);
    }
  }

  TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
  canvas->cd();
  gStyle->SetOptStat(0);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.13);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);

  TLegend* legend = new TLegend(0.38, 0.90-0.10/3.*float(histolist.size()), 0.75, 0.90);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{1 fb^{-1} %i TeV}", theSqrts);
  text = pt->AddText(0.9, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  bool first=true;
  float minY=9e9, maxY=-9e9;
  for (TH1F*& hh:histolist){
    if (hh->GetName()[0]=='T') continue;

    hh->GetXaxis()->SetRangeUser(220, hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX()));
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetLabelFont(42);
    hh->GetXaxis()->SetLabelOffset(0.007);
    hh->GetXaxis()->SetLabelSize(0.04);
    hh->GetXaxis()->SetTitleSize(0.06);
    hh->GetXaxis()->SetTitleOffset(0.9);
    hh->GetXaxis()->SetTitleFont(42);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetLabelFont(42);
    hh->GetYaxis()->SetLabelOffset(0.007);
    hh->GetYaxis()->SetLabelSize(0.04);
    hh->GetYaxis()->SetTitleSize(0.06);
    hh->GetYaxis()->SetTitleOffset(1.1);
    hh->GetYaxis()->SetTitleFont(42);
    hh->GetYaxis()->SetTitle("# events");
    minY = std::min(float(hh->GetMinimum()), minY);
    maxY = std::max(float(hh->GetMaximum()), maxY);
    legend->AddEntry(hh, hh->GetTitle(), "l");
    hh->SetTitle("");
  }
  for (TH1F*& hh:histolist){
    if (hh->GetName()[0]=='T') continue;
    hh->GetYaxis()->SetRangeUser((minY<0. ? minY*1.2 : minY*0.8), maxY*1.2);
    if (first){
      hh->Draw("hist");
      first=false;
    }
    else{
      hh->Draw("histsame");
    }
  }
  legend->Draw("same");
  pt->Draw();
  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".png"));
  canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".pdf"));
  delete legend;
  canvas->Close();
}


#endif
