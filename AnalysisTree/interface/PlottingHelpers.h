#ifndef PLOTTINGHELPERS_H
#define PLOTTINGHELPERS_H

#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"


namespace PlottingHelpers{
  enum CMSLogoStep{
    kPaper,
    kSimulation,
    kPreliminary,
    kWIP
  };

  class PlotCanvas{
  protected:
    TString canvasname;
    int Nx;
    int Ny;
    double vfrac_bottom;
    double canvas_size_x;
    double canvas_size_y;
    double lMargin;
    double rMargin;
    double space_x;
    double hStep;
    double hHalfOffset;
    double bMargin;
    double tMargin;
    double space_y;
    double vStep;
    double vHalfOffset;

    TCanvas* canvas;
    std::vector<std::vector<TPad*>> insidepanels; // Order is from bottom-left to top-right
    std::vector<TPad*> borderpanels; // Order is bottom, left, top, right
    std::vector<TLatex*> texts;

    // These are objects created outside but specified as owned later on
    std::vector<TLegend*> legends;
    
    void partitionCanvas();

  public:
    PlotCanvas(
      TString canvasname,
      int npixels_per_pad_x, int npixels_per_pad_y,
      int npads_x, int npads_y,
      double margin_left_rel, double margin_right_rel, double margin_bottom_rel, double margin_top_rel, // These are relative to the pad size, not the canvas size!
      double space_x_rel, double space_y_rel,
      double vfrac_bottom_=1. // Size scale of the bottommost pad (e.g., when ratios are plotted)
    );
    virtual ~PlotCanvas();

    TCanvas* const& getCanvas() const{ return canvas; }
    std::vector<std::vector<TPad*>> const& getInsidePanels() const{ return insidepanels; }
    std::vector<TPad*> const& getBorderPanels() const{ return borderpanels; }

    double getStdPixelSize_CMSLogo() const;
    double getStdPixelSize_CMSLogoExtras() const;
    double getStdPixelSize_XYTitle() const;
    double getStdPixelSize_XYLabel() const;

    static int getStdFont_XYTitle();
    static int getStdFont_XYLabel();

    double getStdOffset_XLabel() const;
    double getStdOffset_YLabel() const;

    void addCMSLogo(CMSLogoStep type, double sqrts, double lumi); // sqrts or lumi = -1 disables their addition

    void addLegend(TLegend* const& legend);
    void addText(TLatex* const& text);

    void update();
  };

  TCanvas* makeSquareCanvas(TString const& canvasname, bool is2D);
  TLegend* makeLegend(float xlow, float ylow, float xhigh, float yhigh);
  TPaveText* makeCMSLogo();


}

#endif
