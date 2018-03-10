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

void plotTH1Fs(TFile* foutput, TString coutdir, TString canvasname, TString ytitle, std::vector<TH1F*>& histolist);
void plotSystRatioTH1Fs(TFile* foutput, TString coutdir, TString canvasname, TString strSystematics, std::vector<TH1F*>* histolist[2]);
void plotTH2Fs(TFile* foutput, TString coutdir, std::vector<TH2F*>& histolist);
 
void plotProcessCheckStage(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst)) return;

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strSystematics = getSystematicsName(syst);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Check_" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Plots/" + strStage + "/"
    + strACHypo + "/" + strChannel + "/" + strCategory + "/" + strGenerator + "/" + strSystematics + "/";

  gSystem->Exec("mkdir -p " + coutput_common);
  TString INPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString OUTPUT_NAME = INPUT_NAME;
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  INPUT_NAME += ".root";
  OUTPUT_NAME += "_plot.root";
  OUTPUT_LOG_NAME += "_plot.log";
  TString cinput = cinput_common + INPUT_NAME;
  TString coutput = coutput_common + OUTPUT_NAME;
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

  TFile* foutput = TFile::Open(coutput, "recreate");

  vector<TH1F*> h1dlist;
  vector<TH2F*> h2dlist;
  vector<TH1F*> h1dnomlist;
  vector<TH2F*> h2dnomlist;
  finput->cd();
  MELAout << "Extracting 1D histograms..." << endl;
  HelperFunctions::extractHistogramsFromDirectory(finput, h1dlist);
  MELAout << "N(1D histograms) = " << h1dlist.size() << endl;
  MELAout << "Extracting 2D histograms..." << endl;
  HelperFunctions::extractHistogramsFromDirectory(finput, h2dlist);
  MELAout << "N(2D histograms) = " << h2dlist.size() << endl;
  if (finputNominal){
    finputNominal->cd();
    MELAout << "Extracting nominal 1D histograms..." << endl;
    HelperFunctions::extractHistogramsFromDirectory(finputNominal, h1dnomlist);
    MELAout << "N(Nom. 1D histograms) = " << h1dnomlist.size() << endl;
    MELAout << "Extracting nominal 2D histograms..." << endl;
    HelperFunctions::extractHistogramsFromDirectory(finputNominal, h2dnomlist);
    MELAout << "N(Nom. 2D histograms) = " << h2dnomlist.size() << endl;
  }

  plotTH1Fs(foutput, coutput_common, Form("c_M4lDistribution_%s", thePerProcessHandle->getProcessName().Data()), "# events", h1dlist);
  plotTH2Fs(foutput, coutput_common, h2dlist);

  if (finputNominal){
    if (h1dnomlist.size()==h1dlist.size()){
      for (unsigned int ih=0; ih<h1dlist.size(); ih++){
        TString htitle = h1dlist.at(ih)->GetTitle();
        h1dlist.at(ih)->Divide(h1dnomlist.at(ih));
        h1dlist.at(ih)->SetName(Form("%s%s", h1dlist.at(ih)->GetName(), "_RatioToNominal"));
        h1dlist.at(ih)->SetTitle(htitle);
      }
      for (unsigned int ih=0; ih<h2dlist.size(); ih++){
        TString htitle = h2dlist.at(ih)->GetTitle();
        h2dlist.at(ih)->Divide(h2dnomlist.at(ih));
        h2dlist.at(ih)->SetName(Form("%s%s", h2dlist.at(ih)->GetName(), "_RatioToNominal"));
        h2dlist.at(ih)->SetTitle(htitle);
      }

      plotTH1Fs(foutput, coutput_common, Form("c_M4lDistribution_%s_RatioToNominal", thePerProcessHandle->getProcessName().Data()), "Ratio", h1dlist);
      plotTH2Fs(foutput, coutput_common, h2dlist);
    }
  }

  foutput->Close();
  if (finputNominal) finputNominal->Close();
  finput->Close();
  MELAout.close();
}


void plotTH1Fs(TFile* foutput, TString coutdir, TString canvasname, TString ytitle, std::vector<TH1F*>& histolist){
  TDirectory* curdir = gDirectory;
  if (foutput) foutput->cd();

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
  unordered_map<TH1F*, TString> htitles;
  for (TH1F*& hh:histolist){
    if (hh->GetName()[0]=='T') continue;

    float minX=220;
    float maxX=hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX());
    int binXlow = hh->GetXaxis()->FindBin(minX);
    int binXhigh = hh->GetXaxis()->FindBin(maxX);

    hh->GetXaxis()->SetRangeUser(minX, maxX);
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
    hh->GetYaxis()->SetTitle(ytitle);
    for (int ix=binXlow; ix<binXhigh; ix++){
      float bc = hh->GetBinContent(ix);
      if (bc!=0.){
        minY = std::min(bc, minY);
        maxY = std::max(bc, maxY);
      }
    }
    cout << "Min = " << minY << ", max = " << maxY << " after " << hh->GetName() << " (" << hh->GetTitle() << ")" << endl;
    legend->AddEntry(hh, hh->GetTitle(), "l");
    htitles[hh]=hh->GetTitle();
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
  foutput->WriteTObject(canvas);
  delete pt;
  delete legend;
  canvas->Close();

  // Restore titles
  for (TH1F*& hh:histolist) hh->SetTitle(htitles[hh]);

  curdir->cd();
}
void plotSystRatioTH1Fs(TFile* foutput, TString coutdir, TString canvasname, TString strSystematics, std::vector<TH1F*>* histolist[2]){
  const TString ytitle = "Ratio to nominal";
  assert(histolist[0]->size()==histolist[1]->size());

  TDirectory* curdir = gDirectory;
  if (foutput) foutput->cd();

  HelperFunctions::replaceString(canvasname, "h1D_T_", "");
  {
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    const unsigned int stepsize = ncolors/histolist[0]->size();
    MELAout << "plotTH1Fs: Ncolors: " << ncolors << ", step size " << stepsize << endl;
    for (int is=0; is<2; is++){
      for (unsigned int ih=0; ih<histolist[is]->size(); ih++){
        int colorToUse = gStyle->GetColorPalette(ih*stepsize);
        histolist[is]->at(ih)->SetMarkerColor(colorToUse);
        histolist[is]->at(ih)->SetLineColor(colorToUse);
        histolist[is]->at(ih)->SetLineWidth(2);
      }
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

  TLegend* legend = new TLegend(0.38, 0.90-0.10/3.*float(histolist[0]->size()), 0.75, 0.90);
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
  unordered_map<TH1F*, TString> htitles;
  for (unsigned int is=0; is<2; is++){
    for (TH1F*& hh:*(histolist[is])){
      if (hh->GetName()[0]=='T') continue;

      const float minX=70;
      const float maxX=hh->GetXaxis()->GetBinUpEdge(hh->GetNbinsX());
      int binXlow = hh->GetXaxis()->FindBin(minX);
      int binXhigh = hh->GetXaxis()->FindBin(maxX);

      hh->GetXaxis()->SetRangeUser(minX, maxX);
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
      hh->GetYaxis()->SetTitle(ytitle);
      for (int ix=binXlow; ix<binXhigh; ix++){
        float bc = hh->GetBinContent(ix);
        if (bc!=0.){
          minY = std::min(bc, minY);
          maxY = std::max(bc, maxY);
        }
      }
      cout << "Min = " << minY << ", max = " << maxY << " after " << hh->GetName() << " (" << hh->GetTitle() << ")" << endl;
      if (is==0) legend->AddEntry(hh, hh->GetTitle(), "l");
      htitles[hh]=hh->GetTitle();
      hh->SetTitle("");
    }
  }
  for (unsigned int is=0; is<2; is++){
    for (TH1F*& hh:*(histolist[is])){
      if (hh->GetName()[0]=='T') continue;
      hh->GetYaxis()->SetRangeUser(0.5, 2.);
      if (first){
        hh->Draw("hist");
        first=false;
      }
      else hh->Draw("histsame");
    }
  }
  legend->Draw("same");
  pt->Draw();
  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".png"));
  canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".pdf"));
  foutput->WriteTObject(canvas);
  delete legend;
  canvas->Close();

  // Restore titles
  for (unsigned int is=0; is<2; is++){ for (TH1F*& hh:*(histolist[is])) hh->SetTitle(htitles[hh]); }
  curdir->cd();
}
void plotTH2Fs(TFile* foutput, TString coutdir, std::vector<TH2F*>& histolist){
  TDirectory* curdir = gDirectory;
  if (foutput) foutput->cd();

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

    float minZ=9e9, maxZ=-9e9;
    for (int ix=1; ix<=hh->GetNbinsX(); ix++){
      for (int iy=1; iy<=hh->GetNbinsY(); iy++){
        float bc = hh->GetBinContent(ix, iy);
        if (bc!=0.){
          minZ = std::min(minZ, bc);
          maxZ = std::max(maxZ, bc);
        }
      }
    }
    hh->GetZaxis()->SetRangeUser(minZ, maxZ);

    hh->Draw("colz");
    pt->Draw();

    canvas->RedrawAxis();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".png"));
    canvas->SaveAs(Form("%s/%s%s", coutdir.Data(), canvasname.Data(), ".pdf"));
    foutput->WriteTObject(canvas);
    delete pt;
    canvas->Close();
  }

  curdir->cd();
}

void plotProcessCheckStage_SystPairs(
  const Channel channel, const Category category, const ACHypothesis hypo, const SystematicVariationTypes syst,
  const unsigned int istage,
  const TString fixedDate,
  ProcessHandler::ProcessType proctype,
  const TString strGenerator
){
  if (channel==NChannels) return;
  if (!CheckSetTemplatesCategoryScheme(category)) return;
  ProcessHandler const* thePerProcessHandle=getOffshellProcessHandler(proctype);
  if (!thePerProcessHandle) return;
  if (!systematicAllowed(category, channel, thePerProcessHandle->getProcessType(), syst)) return;

  const int isyst = convertSystematicVariationTypeToInt(syst);
  if (syst==sNominal || isyst%2==1) return; // This means the code runs only when Up systematics are passed.

  const TString strChannel = getChannelName(channel);
  const TString strCategory = getCategoryName(category);
  const TString strACHypo = getACHypothesisName(hypo);
  const TString strStage = Form("Stage%i", istage);
  TString strSystematics = getSystematicsName(syst);
  replaceString(strSystematics, "Up", "Syst");

  const SystematicVariationTypes systDn = (const SystematicVariationTypes) (isyst-1);
  const SystematicVariationTypes systUp = syst;
  const TString strSystematics_Dn = getSystematicsName(systDn);
  const TString strSystematics_Up = getSystematicsName(systUp);

  // Setup the output directories
  TString sqrtsDir = Form("LHC_%iTeV/", theSqrts);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;
  TString cinput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Check_" + strStage + "/";
  TString coutput_common = user_output_dir + sqrtsDir + "Templates/" + strdate + "/Plots/" + strStage + "/"
    + strACHypo + "/" + strChannel + "/" + strCategory + "/" + strGenerator + "/" + strSystematics + "/";

  gSystem->Exec("mkdir -p " + coutput_common);
  TString INPUT_NOMINAL_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants.root",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    getSystematicsName(sNominal).Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString INPUT_DN_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants.root",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics_Dn.Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString INPUT_UP_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants.root",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics_Up.Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString OUTPUT_NAME = Form(
    "HtoZZ%s_%s_FinalTemplates_%s_%s_%s_Check%sDiscriminants",
    strChannel.Data(), strCategory.Data(),
    thePerProcessHandle->getProcessName().Data(),
    strSystematics.Data(),
    strGenerator.Data(),
    getACHypothesisName(hypo).Data()
  );
  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += "_plot.root";
  OUTPUT_LOG_NAME += "_plot.log";

  // Open the input files
  TString cinput[3] ={
    cinput_common + INPUT_NOMINAL_NAME,
    cinput_common + INPUT_DN_NAME,
    cinput_common + INPUT_UP_NAME,
  };
  TFile* finput[3] ={
    TFile::Open(cinput[0], "read"),
    TFile::Open(cinput[1], "read"),
    TFile::Open(cinput[2], "read")
  };
  if (!finput[0] || !finput[1] || !finput[2]) return;
  vector<TH1F*> h1dlist[3];
  MELAout << "Extracting 1D histograms..." << endl;
  for (unsigned int is=0; is<3; is++){
    HelperFunctions::extractHistogramsFromDirectory(finput[is], h1dlist[is]);
    MELAout << "N(1D histograms for syst " << (is==0 ? "Nominal" : (is==1 ? strSystematics_Dn : strSystematics_Up)) << ") = " << h1dlist[is].size() << endl;
    if (is>0){
      assert(h1dlist[is].size()==h1dlist[0].size());
      for (unsigned int ih=0; ih<h1dlist[0].size(); ih++) h1dlist[is][ih]->Divide(h1dlist[0][ih]);
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

  std::vector<TH1F*>* historatios[2]={
    &(h1dlist[1]),
    &(h1dlist[2])
  };
  plotSystRatioTH1Fs(foutput, coutput_common, Form("c_M4lDistribution_%s_RatioToNominal", thePerProcessHandle->getProcessName().Data()), strSystematics, historatios);

  for (unsigned int is=0; is<3; is++) finput[is]->Close();
  foutput->Close();
  MELAout.close();
}


#endif
