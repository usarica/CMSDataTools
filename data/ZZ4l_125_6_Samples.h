#include <iostream>
#include <string>
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"

const float PI_VAL = TMath::Pi();

enum sample {
	kSig,
	kBkg,
	kBSI,
	kBSI25,

	kM4LSig,
	kM4LBkg,
	kM4LBSI,
	kM4LOn4Sig,
	kM4LOn4Bkg,
	kM4LOn4BSI,

	kNumSamples
};
TString sampleName[kNumSamples+7] = {
	"Sig",
	"Bkg",
	"BSI",
	"BSI25",
	"QQZZ",
	"M4LSig",
	"M4LBkg",
	"M4LBSI",
	"M4LOn4Sig",
	"M4LOn4Bkg",
	"M4LOn4BSI"
};
TString sample_suffix[kNumSamples+7]={
	"Sig",
	"Bkg",
	"BSI",
	"BSI25",
	"QQZZ",
	"M4LSig",
	"M4LBkg",
	"M4LBSI",
	"M4LOn4Sig",
	"M4LOn4Bkg",
	"M4LOn4BSI"
};
TString sample_suffix_gg2VV[kNumSamples-1]={
	"H125.6",
	"Continuum",
	"ContinuumInterfH125.6"
};
TString sample_suffix_MCFM[kNumSamples-1]={
	"SMH-MCFM67",
	"Contin-MCFM67",
	"SMHContinInterf-MCFM67",
	"BSMHContinInterf-MCFM67"
};

const float gi_phi2_phi4[kNumSamples][9]={ // g1-4; phia2,3; fa2, 3; fepspr
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0},
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0},
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0},
	{	1.0,	0,	0,	0,		0,	0,	0,	0,	0}
};

float luminosity[2]={5.051,19.712};
float BR_Table[2][3][4]={
	{
		{ 2.76E-04,	1.25E-04,	3.27E-05,	5.93E-05 }, // e+mu+tau, e+mu, 4e, 2e2mu: 125 GeV
		{ 2.91E-04,	1.32E-04,	3.45E-05,	6.27E-05 }, // e+mu+tau, e+mu, 4e, 2e2mu: 125.6 GeV
		{ 3.02E-04,	1.36E-04,	3.56E-05,	6.49E-05 } // e+mu+tau, e+mu, 4e, 2e2mu: 126 GeV
	},
	{
		{ 2.76E-04,	1.25E-04,	3.27E-05,	5.93E-05 }, // e+mu+tau, e+mu, 4e, 2e2mu: 125 GeV
		{ 2.91E-04,	1.32E-04,	3.45E-05,	6.27E-05 }, // e+mu+tau, e+mu, 4e, 2e2mu: 125.6 GeV
		{ 3.02E-04,	1.36E-04,	3.56E-05,	6.49E-05 } // e+mu+tau, e+mu, 4e, 2e2mu: 126 GeV
	}
};

float XSEC_Table[2][3]={
	{15.1326,14.987,14.8909},
	{19.2681,19.089,18.9708}
};

/*
float xsec_ggHZZ_BSI25_MCFM[2] = {0.0119633933965,0.01612206167475}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e, Gamma = 24.4106*Gamma_SM
float xsec_ggHZZ_BSI_MCFM[2] = {0.00900945471525,0.0118336229725}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_MCFM[2] = {0.00256261087575,0.00326849098575}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_MCFM[2] = {0.00687586094425,0.009192140617}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
*/
float xsec_ggHZZ_BSI25_MCFM[2] = {0.011959430747,0.016163693396}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e, Gamma = 24.4106*Gamma_SM
float xsec_ggHZZ_BSI_MCFM[2] = {0.009007274418,0.011840779163}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggHZZ_MCFM[2] = {0.002549983848,0.003255176587}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_MCFM[2] = {0.006885033938,0.009188219964}; // 4e+4mu+2e2mu, 2e2mu=4e+4mu=2x4e
float xsec_ggZZ_CERN[2][2] = {
	{0.00174,0.00348},
	{0.0048,0.01203}
}; // 4e or 4mu, mutliply by 2 for 2l2lp
float xsec_qqZZ_CERN[2][2] = {
	{0.06609,0.152},
	{0.07691,0.1767}
}; // 4e or 4mu


const int kNumBkg=6;
char hzz4lprefix[]="HZZ4lTree_";
TString sample_BackgroundFile[kNumBkg]={
	"ZZTo4mu",
	"ZZTo4e",
	"ZZTo2e2mu",
	"ZZTo2mu2tau",
	"ZZTo2e2tau",
	"ZZTo4tau"
};
//Change to relevant
TString user_dir="./UsingLastProduction/";
//string user_dir="/afs/cern.ch/work/u/usarica/WidthAnalysis/";
TString user_gg2VV_location="/scratch0/hep/ianderso/CJLST/ReprocessedTrees/lastProduction/";
//string user_gg2VV_location="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/";
TString user_folder[5]={
	"4mu",
	"4e",
	"2mu2e",
	"CR",
	"data"
};
