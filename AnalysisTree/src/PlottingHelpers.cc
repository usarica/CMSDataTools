#include <cassert>
#include <cmath>
#include "HelperFunctions.h"
#include "PlottingHelpers.h"
#include "MELAStreamHelpers.hh"
#include "TStyle.h"
#include "TROOT.h"
#include "TText.h"


using namespace std;
using namespace MELAStreamHelpers;


namespace PlottingHelpers{

  PlotCanvas::PlotCanvas(
    TString canvasname_,
    int npixels_per_pad_x, int npixels_per_pad_y,
    int npads_x, int npads_y,
    double margin_left_rel, double margin_right_rel, double margin_bottom_rel, double margin_top_rel,
    double space_x_rel, double space_y_rel,
    double vfrac_bottom_
  ) :
    canvasname(canvasname_),
    Nx(npads_x),
    Ny(npads_y),
    vfrac_bottom(vfrac_bottom_),

    canvas(nullptr)
  {
    TDirectory* curdir = gDirectory;

    // Disable disproportionate plotting if there is only one panel along the y axia.
    if (Ny==1) vfrac_bottom = 1;

    int lMargin_pixels = npixels_per_pad_x * margin_left_rel + 0.5; if (lMargin_pixels%2==1) lMargin_pixels++;
    int rMargin_pixels = npixels_per_pad_x * margin_right_rel + 0.5; if (rMargin_pixels%2==1) rMargin_pixels++;
    int space_x_pixels = npixels_per_pad_x * space_x_rel + 0.5; if (space_x_pixels%2==1) space_x_pixels++;
    int bMargin_pixels = npixels_per_pad_y * margin_bottom_rel + 0.5; if (bMargin_pixels%2==1) bMargin_pixels++;
    int tMargin_pixels = npixels_per_pad_y * margin_top_rel + 0.5; if (tMargin_pixels%2==1) tMargin_pixels++;
    int space_y_pixels = npixels_per_pad_y * space_y_rel + 0.5; if (space_y_pixels%2==1) space_y_pixels++;
    int npixels_bottompad_y = npixels_per_pad_y * vfrac_bottom + 0.5; if (npixels_bottompad_y%2==1) npixels_bottompad_y++;

    canvas_size_x = lMargin_pixels + rMargin_pixels + Nx*npixels_per_pad_x + (Nx-1)*space_x_pixels;
    canvas_size_y = bMargin_pixels + tMargin_pixels + (Ny-1)*(npixels_per_pad_y + space_y_pixels) + npixels_bottompad_y;

    // Convert quantities to absolute canvas size
    lMargin = double(lMargin_pixels)/double(canvas_size_x);
    rMargin = double(rMargin_pixels)/double(canvas_size_x);
    space_x = double(space_x_pixels)/double(canvas_size_x);
    hStep = double(npixels_per_pad_x)/double(canvas_size_x);
    hHalfOffset = (Nx==1 ? 0. : space_x/2.);
    //MELAout << "lMargin, rMargin, space_x, hStep, hHalfOffset = " << lMargin << ", " << rMargin << ", " << space_x << ", " << hStep << ", " << hHalfOffset << endl;

    bMargin = double(bMargin_pixels)/double(canvas_size_y);
    tMargin = double(tMargin_pixels)/double(canvas_size_y);
    space_y = double(space_y_pixels)/double(canvas_size_y);
    vStep = double(npixels_per_pad_y)/double(canvas_size_y);
    vHalfOffset = (Ny==1 ? 0. : space_y/2.);
    // Reset vfrac_bottom for pixel differences as well
    vfrac_bottom = double(npixels_bottompad_y)/double(npixels_per_pad_y);
    //MELAout << "bMargin, tMargin, space_y, vStep, vHalfOffset = " << bMargin << ", " << tMargin << ", " << space_y << ", " << vStep << ", " << vHalfOffset << endl;

    canvas = new TCanvas(canvasname, "", canvas_size_x, canvas_size_y);

    gStyle->SetOptStat(0);

    partitionCanvas();

    curdir->cd();
  }
  PlotCanvas::~PlotCanvas(){
    for (auto& text:texts) delete text;
    for (auto& legend:legends) delete legend;
    for (auto& pad:borderpanels) pad->Close();
    for (auto& vpad:insidepanels){ for (auto& pad:vpad) pad->Close(); }
    canvas->Close();
  }

  void PlotCanvas::partitionCanvas(){
    if (Nx<0 || Ny<0){
      MELAerr << "PlotCanvas::partitionCanvas: Invalid Nx, Ny = " << Nx << ", " << Ny << endl;
      assert(0);
    }

    double vposd=0, vposu=0, vmard=0, vmaru=0, vfactor=0;
    double hposl=0, hposr=0, hmarl=0, hmarr=0, hfactor=0;

    insidepanels.assign(Nx, std::vector<TPad*>(Ny, nullptr));

    canvas->cd();
    for (int i=0; i<Nx; i++){
      if (i==0){
        hposl = 0.;
        hposr = (i==Nx-1 ? 1. : lMargin + hStep + hHalfOffset);
        hfactor = hposr - hposl;
        hmarl = lMargin / hfactor;
        hmarr = (i==Nx-1 ? rMargin : hHalfOffset) / hfactor;
      }
      else if (i==Nx-1){
        hposl = hposr;
        hposr = 1.;
        hfactor = hposr - hposl;
        hmarl = hHalfOffset / hfactor;
        hmarr = rMargin / hfactor;
      }
      else{
        hposl = hposr;
        hposr = hposl + hStep + space_x;
        hfactor = hposr - hposl;
        hmarl = hHalfOffset/hfactor;
        hmarr = hHalfOffset/hfactor;
      }

      for (int j=0; j<Ny; j++){
        if (j==0){
          vposd = 0.;
          vposu = (j==Ny-1 ? 1. : bMargin + vStep*vfrac_bottom + vHalfOffset);
          vfactor = vposu - vposd;
          vmard = bMargin / vfactor;
          vmaru = (j==Ny-1 ? tMargin : vHalfOffset) / vfactor;
        }
        else if (j==Ny-1){
          vposd = vposu;
          vposu = 1.;
          vfactor = vposu - vposd;
          vmard = vHalfOffset / vfactor;
          vmaru = tMargin / vfactor;
        }
        else{
          vposd = vposu;
          vposu = vposd + vStep + space_y;
          vfactor = vposu - vposd;
          vmard = vHalfOffset / vfactor;
          vmaru = vHalfOffset / vfactor;
        }

        TPad* pad = nullptr;

        TString padname = canvasname + Form("_insidepanel_%i_%i", i, j);
        if (gROOT->FindObject(padname)) MELAerr << "PlotCanvas::partitionCanvas: There is already a pad named " << padname << "!" << endl;

        pad = new TPad(padname, "", hposl, vposd, hposr, vposu);

        pad->SetLeftMargin(hmarl);
        pad->SetRightMargin(hmarr);
        pad->SetBottomMargin(vmard);
        pad->SetTopMargin(vmaru);

        pad->SetFillColor(0);
        pad->SetBorderMode(0);
        pad->SetBorderSize(0);
        pad->SetFrameFillStyle(0);
        pad->SetFrameBorderMode(0);
        pad->SetFrameBorderSize(0);
        pad->SetTickx(1);
        pad->SetTicky(1);

        pad->Draw();

        insidepanels.at(i).at(j) = pad;
        canvas->cd();
      }
    }

    borderpanels.assign(4, nullptr);
    for (int iborder=0; iborder<4; iborder++){
      TPad* pad = nullptr;

      TString padname = canvasname + Form("_borderpanel_%i", iborder);
      if (gROOT->FindObject(padname)) MELAerr << "PlotCanvas::partitionCanvas: There is already a pad named " << padname << "!" << endl;

      switch (iborder){
      case 0:
      {
        // Bottom pad
        hposl = lMargin;
        hposr = 1.-rMargin;
        vposd = 0;
        double stdoffset = getStdPixelSize_XYTitle()/canvas_size_y*1.5;
        if (stdoffset>=bMargin){
          MELAerr << "PlotCanvas::partitionCanvas: Requested bottom offset " << stdoffset << ">= " << bMargin << ". Scaling to " << bMargin*0.75 << "." << endl;
          vposu = bMargin*0.75;
        }
        else vposu = stdoffset;
        break;
      }
      case 1:
      {
        // Left pad
        hposl = 0.;
        double stdoffset = getStdPixelSize_XYTitle()/canvas_size_x*1.5;
        if (stdoffset>=lMargin){
          MELAerr << "PlotCanvas::partitionCanvas: Requested left offset " << stdoffset << ">= " << lMargin << ". Scaling to " << lMargin*0.75 << "." << endl;
          hposr = lMargin*0.75;
        }
        else hposr = stdoffset;
        vposd = bMargin;
        vposu = 1.-tMargin;
        break;
      }
      case 2:
      {
        // Top pad
        hposl = lMargin;
        hposr = 1.-rMargin;
        vposd = 1.-tMargin*0.98;
        vposu = 1.;
        break;
      }
      default:
      {
        // Right pad
        hposl = 1.-rMargin*0.9;
        hposr = 1.;
        vposd = bMargin;
        vposu = 1.-tMargin;
        break;
      }
      }

      pad = new TPad(padname, "", hposl, vposd, hposr, vposu);

      pad->SetLeftMargin(hmarl);
      pad->SetRightMargin(hmarr);
      pad->SetBottomMargin(vmard);
      pad->SetTopMargin(vmaru);

      pad->SetFillColor(0);
      pad->SetBorderMode(0);
      pad->SetBorderSize(0);
      pad->SetFrameFillStyle(0);
      pad->SetFrameBorderMode(0);
      pad->SetFrameBorderSize(0);
      pad->SetTickx(1);
      pad->SetTicky(1);

      pad->Draw();

      borderpanels.at(iborder) = pad;
      canvas->cd();
    }
  }

  double PlotCanvas::getStdPixelSize_CMSLogo() const{
    return canvas_size_y*tMargin*0.8;
  }
  double PlotCanvas::getStdPixelSize_CMSLogoExtras() const{
    return getStdPixelSize_CMSLogo()*0.72;
  }
  double PlotCanvas::getStdPixelSize_XYTitle() const{
    return getStdPixelSize_CMSLogo()*0.91;
  }
  double PlotCanvas::getStdPixelSize_XYLabel() const{
    return getStdPixelSize_CMSLogo()*0.72;
  }

  int PlotCanvas::getStdFont_XYTitle(){ return 43; }
  int PlotCanvas::getStdFont_XYLabel(){ return 43; }

  double PlotCanvas::getStdOffset_XLabel() const{ return (0.007/0.8)*vStep/(Ny==1 ? 1. : bMargin + vStep + vHalfOffset); }
  double PlotCanvas::getStdOffset_YLabel() const{ return (0.007/0.78)*hStep/(Nx==1 ? 1. : lMargin + hStep + hHalfOffset); }

  void PlotCanvas::addCMSLogo(CMSLogoStep type, double sqrts, double lumi){
    TDirectory* curdir = gDirectory;

    borderpanels.at(2)->cd();

    // Draw the CMS label
    double cms_pixel_ysize = getStdPixelSize_CMSLogo();
    double cms_pixel_xsize = cms_pixel_ysize*2.2;
    TLatex* cmstxt = new TLatex(); addText(cmstxt);
    cmstxt->SetTextAlign(22);
    cmstxt->SetTextFont(63);
    cmstxt->SetTextSize(cms_pixel_ysize);
    cmstxt->DrawLatexNDC(0. + cms_pixel_xsize*0.5/((1.-rMargin-lMargin)*canvas_size_x), 0.55, "CMS");
    //MELAout << "CMS label relative position: " << 0. + cms_pixel_xsize*0.5/((1.-rMargin-lMargin)*canvas_size_x) << ", " << 0.55 << endl;
    //MELAout << "CMS label size: " << cms_pixel_xsize << ", " << cms_pixel_ysize << endl;

    // Draw the CMS extra label
    double cmsprelim_pixel_ysize = getStdPixelSize_CMSLogoExtras();
    double cmsprelim_relloffset_x = 0, cmsprelim_relloffset_y = 0;
    TString strprelim;
    if (type == kSimulation){
      double cmsprelim_pixel_xsize = cmsprelim_pixel_ysize*4.75;
      cmsprelim_relloffset_y = 0.55;
      cmsprelim_relloffset_x = (cms_pixel_xsize+cmsprelim_pixel_xsize/2.)/((1.-rMargin-lMargin)*canvas_size_x);
      strprelim = "Simulation";
    }
    else if (type == kPreliminary){
      double cmsprelim_pixel_xsize = cmsprelim_pixel_ysize*5.2;
      cmsprelim_relloffset_y = 0.5;
      cmsprelim_relloffset_x = (cms_pixel_xsize+cmsprelim_pixel_xsize/2.)/((1.-rMargin-lMargin)*canvas_size_x);
      strprelim = "Preliminary";
    }
    else if (type == kWIP){
      double cmsprelim_pixel_xsize = cmsprelim_pixel_ysize*7.5;
      cmsprelim_relloffset_y = 0.5;
      cmsprelim_relloffset_x = (cms_pixel_xsize+cmsprelim_pixel_xsize/2.)/((1.-rMargin-lMargin)*canvas_size_x);
      strprelim = "Work in progress";
    }
    if (cmsprelim_relloffset_x!=0. && cmsprelim_relloffset_y!=0.){
      TLatex* cmstxtmore = new TLatex(); addText(cmstxtmore);
      cmstxtmore->SetTextAlign(22);
      cmstxtmore->SetTextFont(53);
      cmstxtmore->SetTextSize(cmsprelim_pixel_ysize);
      cmstxtmore->DrawLatexNDC(cmsprelim_relloffset_x, cmsprelim_relloffset_y, strprelim);
    }

    // Draw the sqrts indicator
    if (lumi>0. && sqrts>0.){
      TString strSqrts = Form("%s fb^{-1} %s TeV", HelperFunctions::castValueToString(lumi, 1).data(), HelperFunctions::castValueToString(sqrts, 1).data());
      double cmssqrts_pixel_ysize = cmsprelim_pixel_ysize;
      double cmssqrts_pixel_xsize = cmssqrts_pixel_ysize*0.38*strSqrts.Length();
      TLatex* cmssqrts = new TLatex(); addText(cmssqrts);
      cmssqrts->SetTextAlign(22);
      cmssqrts->SetTextFont(PlotCanvas::getStdFont_XYTitle());
      cmssqrts->SetTextSize(cmsprelim_pixel_ysize);
      cmssqrts->DrawLatexNDC((1.-rMargin-lMargin-cmssqrts_pixel_xsize*0.5/canvas_size_x)/(1.-rMargin-lMargin), 0.55, strSqrts);
    }

    curdir->cd();
  }
  void PlotCanvas::addLegend(TLegend* const& legend){
    if (!HelperFunctions::checkListVariable(legends, legend)) legends.push_back(legend);
  }
  void PlotCanvas::addText(TLatex* const& text){
    if (!HelperFunctions::checkListVariable(texts, text)) texts.push_back(text);
  }
  void PlotCanvas::update(){
    for (auto& vpad:insidepanels){
      for (auto& pad:vpad){
        pad->RedrawAxis();
        pad->Modified();
        pad->Update();
      }
    }
    for (auto& pad:borderpanels){
      pad->Modified();
      pad->Update();
    }
    canvas->Modified();
    canvas->Update();
  }

  void PlotCanvas::save(TString outdir, TString strformat, TString newname){
    TString stroutput;
    stroutput = outdir + "/" + (newname=="" ? canvasname : newname) + "." + strformat;
    canvas->SaveAs(stroutput);
  }

  void get1DPlotYRange(std::vector<TH1F*> const& hlist, double const& factorYHigh, bool adjustYLow, double& ymin, double& ymax){
    ymax = -9e14;
    if (adjustYLow) ymin=9e14;
    else ymin = 0;

    for (auto const& hist:hlist){
      for (int ix=1; ix<=hist->GetNbinsX(); ix++){
        double bc = hist->GetBinContent(ix);
        double be = hist->GetBinError(ix);
        if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
        ymax = std::max(ymax, bc+be);
        double bclow = bc; if (be<=bclow) bclow -= be;
        if (adjustYLow && !(bc==0. && be==0.)) ymin = std::min(ymin, bclow);
      }
    }

    if (ymax>=0.) ymax *= (factorYHigh>0. ? factorYHigh : 1.5);
    else ymax /= (factorYHigh>0. ? factorYHigh : 1.5);
    ymin *= (ymin>=0. ? 0.95 : 1.05);
  }

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
