#ifndef COMMON_INCLUDES_H
#define COMMON_INCLUDES_H

#include <iostream>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "HelperFunctions.h"
#include "SampleHelpers.h"
#include "CJLSTSet.h"
#include "BaseTreeLooper.h"
#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ZXFakeRateHandler.h"
#include "ReweightingBuilder.h"
#include "ReweightingFunctions.h"
#include "SystematicsHelpers.h"
#include "ExtendedHistogram_1D.h"
#include "ExtendedHistogram_2D.h"
#include "ACHypothesisHelpers.h"
#include "TemplateHelpers.h"
#include "QuantFuncMathCore.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"

using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace CategorizationHelpers;
using namespace DiscriminantClasses;
using namespace ReweightingFunctions;
using namespace SystematicsHelpers;
using namespace ACHypothesisHelpers;
using namespace TemplateHelpers;
using namespace MELAStreamHelpers;

#endif
