#ifndef ZXFAKERATEHANDLER_H
#define ZXFAKERATEHANDLER_H

#include "HelperFunctions.h"
#include "CJLSTTree.h"


class ZXFakeRateHandler{
protected:
  // Struct to convert input format to something more sensible
  struct FakeRateInput{
    TGraphErrors* tge;
    TGraphAsymmErrors* tgae;
    TH1F* h1f;

    FakeRateInput(TGraphErrors* tge_) : tge(tge_), tgae(nullptr), h1f(nullptr){}
    FakeRateInput(TGraphAsymmErrors* tgae_) : tge(nullptr), tgae(tgae_), h1f(nullptr){}
    FakeRateInput(TH1F* h1f_) : tge(nullptr), tgae(nullptr), h1f(h1f_){}
    ~FakeRateInput(){}
  };

  // Data members
  signed char useUpDn;

  TSpline3* FakeRateInterpolator_ZeeE_barrel;
  TSpline3* FakeRateInterpolator_ZeeE_endcap;
  TSpline3* FakeRateInterpolator_ZmmE_barrel;
  TSpline3* FakeRateInterpolator_ZmmE_endcap;

  TSpline3* FakeRateInterpolator_ZeeM_barrel;
  TSpline3* FakeRateInterpolator_ZeeM_endcap;
  TSpline3* FakeRateInterpolator_ZmmM_barrel;
  TSpline3* FakeRateInterpolator_ZmmM_endcap;

  std::unordered_map<CJLSTTree*, short*> Z1FlavRegistry;
  std::unordered_map<CJLSTTree*, std::vector<short>* const*> LepIdRegistry;
  std::unordered_map<CJLSTTree*, std::vector<float>* const*> LepPtRegistry;
  std::unordered_map<CJLSTTree*, std::vector<float>* const*> LepEtaRegistry;

  // Functions
  TSpline3* convertInputToSpline(ZXFakeRateHandler::FakeRateInput const& frinput);

  float eval(short const& Z1Flav, short const& LepId, float const& LepPt, float const& LepEta) const;

public:
  ZXFakeRateHandler(TString cinput, signed char useUpDn_=0);
  ~ZXFakeRateHandler();

  void registerTree(CJLSTTree* tree);

  float getFakeRateWeight(CJLSTTree* tree) const;

};


#endif
