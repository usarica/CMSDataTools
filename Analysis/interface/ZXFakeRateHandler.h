#ifndef ZXFAKERATEHANDLER_H
#define ZXFAKERATEHANDLER_H

#include "HelperFunctions.h"


class ZXFakeRateHandler{
protected:
  struct FakeRateInput{
    TGraphErrors* tge;
    TGraphAsymmErrors* tgae;
    TH1F* h1f;

    FakeRateInput(TGraphErrors* tge_) : tge(tge_), tgae(nullptr), h1f(nullptr){}
    FakeRateInput(TGraphAsymmErrors* tgae_) : tge(nullptr), tgae(tgae_), h1f(nullptr){}
    FakeRateInput(TH1F* h1f_) : tge(nullptr), tgae(nullptr), h1f(h1f_){}
    ~FakeRateInput(){}
  };

  signed char useUpDn;

  TSpline3* FakeRateInterpolator_ZeeE_barrel;
  TSpline3* FakeRateInterpolator_ZeeE_endcap;
  TSpline3* FakeRateInterpolator_ZmmE_barrel;
  TSpline3* FakeRateInterpolator_ZmmE_endcap;

  TSpline3* FakeRateInterpolator_ZeeM_barrel;
  TSpline3* FakeRateInterpolator_ZeeM_endcap;
  TSpline3* FakeRateInterpolator_ZmmM_barrel;
  TSpline3* FakeRateInterpolator_ZmmM_endcap;

  TSpline3* convertInputToSpline(ZXFakeRateHandler::FakeRateInput const& frinput);

public:
  ZXFakeRateHandler(TString cinput, signed char useUpDn_=0);
  ~ZXFakeRateHandler();

  float eval(short const& Z1Flav, short const& LepId, float const& LepPt, float const& LepEta) const;

};


#endif
