#include "common_includes.h"


// Constants to affect the template code
const TString user_output_dir = "output/";

void plotTH2Fs(TString coutdir, std::vector<TH2F*>& histolist);

void plotProcessCheckStage(
  ProcessHandler::ProcessType proctype,
  const Channel channel, const Category category, ACHypothesisHelpers::ACHypothesis hypo, TString strSystematics,
  const unsigned int istage,
  const TString strGenerator="POWHEG",
  const TString fixedDate=""
){
  if (channel==NChannels) return;
  ProcessHandler const* theProcess=nullptr;
  switch (proctype){
  case ProcessHandler::kGG:
    theProcess = &TemplateHelpers::OffshellGGProcessHandle;
    break;
  case ProcessHandler::kVV:
    theProcess = &TemplateHelpers::OffshellVVProcessHandle;
    break;
  case ProcessHandler::kQQBkg:
    theProcess = &TemplateHelpers::OffshellQQBkgProcessHandle;
    break;
  default:
    break;
  };
  if (!theProcess) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Plots/";

  gSystem->Exec("mkdir -p " + coutput_common);
  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants_Stage%i",
    strChannel.Data(), strCategory.Data(),
    theProcess->getProcessName().Data(),
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

  TFile* finput = TFile::Open(cinput, "read");
  if (!finput) return;
  MELAout << "Opened file " << cinput << endl;
  MELAout << "===============================" << endl;
  MELAout << "CoM Energy: " << theSqrts << " TeV" << endl;
  MELAout << "Decay Channel: " << strChannel << endl;
  MELAout << "===============================" << endl;
  MELAout << endl;

  gStyle->SetOptStat(0);

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

  finput->Close();
  MELAout.close();
}

void plotTH2Fs(TString coutdir, std::vector<TH2F*>& histolist){
  gStyle->SetOptStat(0);

  for (TH2F*& hh:histolist){
    TString canvasname = Form("c_%s", hh->GetName());
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

    /*
    TLegend *legend = new TLegend(0.20, 0.80, 0.58, 0.90);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColor(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    */

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

    hh->Draw("colz");

    canvas->RedrawAxis();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".png"));
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".pdf"));
    //delete legend;
    canvas->Close();
  }

}
