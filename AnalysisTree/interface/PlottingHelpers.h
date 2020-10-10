#ifndef PLOTTINGHELPERS_H
#define PLOTTINGHELPERS_H

#include "TDirectory.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"


namespace PlottingHelpers{
  enum CMSLogoStep{
    kPaper,
    kPreliminary,
    kWIP
  };

  TCanvas* makeSquareCanvas(TString const& canvasname, bool is2D);
  TLegend* makeLegend(float xlow, float ylow, float xhigh, float yhigh);
  TPaveText* makeCMSLogo();


}

#endif
