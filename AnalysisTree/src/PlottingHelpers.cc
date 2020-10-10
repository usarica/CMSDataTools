#include "HelperFunctionsCore.h"
#include "PlottingHelpers.h"
#include "TStyle.h"
#include "TText.h"


namespace PlottingHelpers{
  TCanvas* makeSquareCanvas(TString const& canvasname, bool is2D){
    TDirectory* curdir = gDirectory;

    TCanvas* canvas = nullptr;
    if (!is2D){
      canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
      canvas->cd();
      gStyle->SetOptStat(0);
      canvas->SetRightMargin(0.05);
    }
    else{
      canvas = new TCanvas(canvasname, "", 1000, 800);
      canvas->cd();
      gStyle->SetOptStat(0);
      canvas->SetRightMargin(0.12);
    }
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(2);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLeftMargin(0.17);
    canvas->SetTopMargin(0.07);
    canvas->SetBottomMargin(0.13);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameFillStyle(0);
    canvas->SetFrameBorderMode(0);

    curdir->cd();
    return canvas;
  }
  TLegend* makeLegend(float xlow, float ylow, float xhigh, float yhigh){
    TLegend* legend = new TLegend(xlow, ylow, xhigh, yhigh);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColor(1);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    return legend;
  }
  TPaveText* makeCMSLogo(bool hasData, PlottingHelpers::CMSLogoStep logoStep, float lumi, float sqrts){
    TText* text = nullptr;
    TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);

    text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
    text->SetTextSize(0.044);
    TString strLogoStep;
    if (!hasData) strLogoStep = "Simulation";
    if (logoStep == kPreliminary){
      if (strLogoStep!="") strLogoStep += " preliminary";
      else strLogoStep += "Preliminary";
    }
    else if (logoStep == kWIP){
      if (strLogoStep!="") strLogoStep += " (work in progress)";
      else strLogoStep += "Work in progress";
    }
    text = pt->AddText(0.165, 0.42, Form("#font[52]{%s}", strLogoStep.Data()));
    text->SetTextSize(0.0315);
    if (lumi>0.f && sqrts>0.f){
      TString cErgTev = Form("#font[42]{%s fb^{-1} %s TeV}", HelperFunctions::castValueToString(lumi, 1).data(), HelperFunctions::castValueToString(sqrts, 1).data());
      text = pt->AddText(0.82, 0.45, cErgTev);
      text->SetTextSize(0.0315);
    }

    return pt;
  }

}
