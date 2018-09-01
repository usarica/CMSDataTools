#include "ZXFakeRateHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace HelperFunctions;
using namespace MELAStreamHelpers;


ZXFakeRateHandler::ZXFakeRateHandler(TString cinput, ZXFakeRateHandler::FakeRateMethod FRMethod_, signed char useUpDn_):
  FRMethod(FRMethod_),
  useUpDn(useUpDn_),
  FakeRateInterpolator_ZeeE_barrel(nullptr),
  FakeRateInterpolator_ZeeE_endcap(nullptr),
  FakeRateInterpolator_ZmmE_barrel(nullptr),
  FakeRateInterpolator_ZmmE_endcap(nullptr),
  FakeRateInterpolator_ZeeM_barrel(nullptr),
  FakeRateInterpolator_ZeeM_endcap(nullptr),
  FakeRateInterpolator_ZmmM_barrel(nullptr),
  FakeRateInterpolator_ZmmM_endcap(nullptr)
{
  TDirectory* curdir = gDirectory; // Save current directory to return back to it later
  if (!gSystem->AccessPathName(cinput)){
    TFile* finput = TFile::Open(cinput, "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        finput->cd();
        vector<TGraphErrors*> tglist;
        vector<TGraphAsymmErrors*> tgalist;
        vector<TH1F*> h1flist;
        extractHistogramsFromDirectory(finput, tglist);
        extractHistogramsFromDirectory(finput, tgalist);
        extractHistogramsFromDirectory(finput, h1flist);
        if (!tglist.empty() || !tgalist.empty() || !h1flist.empty()){
          vector<FakeRateInput> frinputs;
          for (TGraphErrors* tg:tglist) frinputs.emplace_back(tg);
          for (TGraphAsymmErrors* tg:tgalist) frinputs.emplace_back(tg);
          for (TH1F* hh:h1flist) frinputs.emplace_back(hh);

          curdir->cd();
          if (frinputs.size()==4){ // Fake rates do not depend on real Z flavor
            for (FakeRateInput& frinput:frinputs){
              TSpline3* sptmp = convertInputToSpline(frinput);
              TString const spname = sptmp->GetName();
              if (spname.Contains("electron") && spname.Contains("EB")){
                FakeRateInterpolator_ZeeE_barrel = sptmp;
                FakeRateInterpolator_ZmmE_barrel = new TSpline3(*sptmp); FakeRateInterpolator_ZmmE_barrel->SetName(Form("%s_copy", spname.Data()));
              }
              else if (spname.Contains("electron") && spname.Contains("EE")){
                FakeRateInterpolator_ZeeE_endcap = sptmp;
                FakeRateInterpolator_ZmmE_endcap = new TSpline3(*sptmp); FakeRateInterpolator_ZmmE_endcap->SetName(Form("%s_copy", spname.Data()));
              }
              else if (spname.Contains("muon") && spname.Contains("EB")){
                FakeRateInterpolator_ZeeM_barrel = sptmp;
                FakeRateInterpolator_ZmmM_barrel = new TSpline3(*sptmp); FakeRateInterpolator_ZmmM_barrel->SetName(Form("%s_copy", spname.Data()));
              }
              else if (spname.Contains("muon") && spname.Contains("EE")){
                FakeRateInterpolator_ZeeM_endcap = sptmp;
                FakeRateInterpolator_ZmmM_endcap = new TSpline3(*sptmp); FakeRateInterpolator_ZmmM_endcap->SetName(Form("%s_copy", spname.Data()));
              }
              else{
                MELAerr << "ZXFakeRateHandler::ZXFakeRateHandler: Could not classify the fake rate type from the spline name " << spname << ". Aborting..." << endl;
                assert(0);
              }
            }
          }
          else if (frinputs.size()==8){ // Fake rates do depend on real Z flavor
            for (FakeRateInput& frinput:frinputs){
              TSpline3* sptmp = convertInputToSpline(frinput);
              TString const spname = sptmp->GetName();
              if (spname.Contains("Zee") && spname.Contains("electron") && spname.Contains("EB")) FakeRateInterpolator_ZeeE_barrel = sptmp;
              else if (spname.Contains("Zee") && spname.Contains("electron") && spname.Contains("EE")) FakeRateInterpolator_ZeeE_endcap = sptmp;
              else if (spname.Contains("Zee") && spname.Contains("muon") && spname.Contains("EB")) FakeRateInterpolator_ZeeM_barrel = sptmp;
              else if (spname.Contains("Zee") && spname.Contains("muon") && spname.Contains("EE")) FakeRateInterpolator_ZeeM_endcap = sptmp;
              else if (spname.Contains("Zmm") && spname.Contains("electron") && spname.Contains("EB")) FakeRateInterpolator_ZmmE_barrel = sptmp;
              else if (spname.Contains("Zmm") && spname.Contains("electron") && spname.Contains("EE")) FakeRateInterpolator_ZmmE_endcap = sptmp;
              else if (spname.Contains("Zmm") && spname.Contains("muon") && spname.Contains("EB")) FakeRateInterpolator_ZmmM_barrel = sptmp;
              else if (spname.Contains("Zmm") && spname.Contains("muon") && spname.Contains("EE")) FakeRateInterpolator_ZmmM_endcap = sptmp;
              else{
                MELAerr << "ZXFakeRateHandler::ZXFakeRateHandler: Could not classify the fake rate type from the spline name " << spname << ". Aborting..." << endl;
                assert(0);
              }
            }
          }
          else{
            MELAerr << "ZXFakeRateHandler::ZXFakeRateHandler: Size of fake rate inputs = " << frinputs.size() << " != 4 or 8. Aborting..." << endl;
            finput->Close();
            assert(0);
          }
        }
        else{
          MELAerr << "ZXFakeRateHandler::ZXFakeRateHandler: Could not extract input fake rates from file " << cinput << endl;
          finput->Close();
          assert(0);
        }
      }
      if (finput->IsOpen()) finput->Close();
    }
  }
  curdir->cd(); // Return back to the directory before opening the input file
}
ZXFakeRateHandler::~ZXFakeRateHandler(){
  delete FakeRateInterpolator_ZeeE_barrel;
  delete FakeRateInterpolator_ZeeE_endcap;
  delete FakeRateInterpolator_ZmmE_barrel;
  delete FakeRateInterpolator_ZmmE_endcap;

  delete FakeRateInterpolator_ZeeM_barrel;
  delete FakeRateInterpolator_ZeeM_endcap;
  delete FakeRateInterpolator_ZmmM_barrel;
  delete FakeRateInterpolator_ZmmM_endcap;
}

TSpline3* ZXFakeRateHandler::convertInputToSpline(ZXFakeRateHandler::FakeRateInput const& frinput){
  TSpline3* res = nullptr;
  if (frinput.tge){
    if (useUpDn!=0){ for (int ix=0; ix<frinput.tge->GetN(); ix++) frinput.tge->GetY()[ix] = std::max(frinput.tge->GetY()[ix] + frinput.tge->GetEY()[ix]*double(useUpDn), 0.); }
    res=convertGraphToSpline3(frinput.tge);
  }
  else if (frinput.tgae){
    if (useUpDn!=0){ for (int ix=0; ix<frinput.tgae->GetN(); ix++) frinput.tgae->GetY()[ix] = std::max(frinput.tgae->GetY()[ix] + (useUpDn>0. ? frinput.tgae->GetEYhigh()[ix] : frinput.tgae->GetEYlow()[ix])*double(useUpDn), 0.); }
    res=convertGraphToSpline3(frinput.tgae);
  }
  else if (frinput.h1f){
    TGraphErrors* tgtmp = makeGraphFromTH1(nullptr, frinput.h1f, Form("tg_%s", frinput.h1f->GetName()));
    if (useUpDn!=0){ for (int ix=0; ix<tgtmp->GetN(); ix++) tgtmp->GetY()[ix] = std::max(tgtmp->GetY()[ix] + tgtmp->GetEY()[ix]*double(useUpDn), 0.); }
    res=convertGraphToSpline3(tgtmp);
    delete tgtmp;
  }
  return res;
}

float ZXFakeRateHandler::eval(int const& CRFlag, short const& Z1Flav, short const& LepId, float const& LepPt, float const& LepEta) const{
  enum{
    CRZLLss=21,
    CRZLLos_2P2F=22,
    CRZLLos_3P1F=23
  };

  float res=0;
  if (FRMethod==mSS){
    if (test_bit(CRFlag, (unsigned int) CRZLLss)){
      TSpline3* spChosen=nullptr;
      unsigned short const absZ1Flav=std::abs(Z1Flav);
      unsigned short const absLepId=std::abs(LepId);
      if (absZ1Flav==121){ // Zee
        if (absLepId==11){
          if (LepEta<1.479) spChosen = FakeRateInterpolator_ZeeE_barrel;
          else spChosen = FakeRateInterpolator_ZeeE_endcap;
        }
        else if (absLepId==13){
          if (LepEta<1.2) spChosen = FakeRateInterpolator_ZeeM_barrel;
          else spChosen = FakeRateInterpolator_ZeeM_endcap;
        }
      }
      else{ // Zmm
        if (absLepId==11){
          if (LepEta<1.479) spChosen = FakeRateInterpolator_ZmmE_barrel;
          else spChosen = FakeRateInterpolator_ZmmE_endcap;
        }
        else if (absLepId==13){
          if (LepEta<1.2) spChosen = FakeRateInterpolator_ZmmM_barrel;
          else spChosen = FakeRateInterpolator_ZmmM_endcap;
        }
      }
      if (spChosen){
        const float xmin = spChosen->GetXmin();
        const float xmax = spChosen->GetXmax();
        float x = LepPt;
        if (x<xmin) x=xmin;
        else if (x>xmax) x=xmax;
        res = spChosen->Eval(x);
      }
    }
  }
  else{
    MELAout << "ZXFakeRateHandler::eval: Only SS method is implemented at this moment!" << endl;
    assert(0);
  }
  return res;
}

void ZXFakeRateHandler::registerTree(CJLSTTree* tree){
  if (!tree) return;

  int* CRFlagRef;
  short* Z1FlavRef;
  vector<short>* const* LepIdRef;
  vector<float>* const* LepPtRef;
  vector<float>* const* LepEtaRef;
  
  tree->bookBranch<BaseTree::BranchType_int_t>("CRflag");
  tree->bookBranch<BaseTree::BranchType_short_t>("Z1Flav");
  tree->bookBranch<BaseTree::BranchType_vshort_t>("LepLepId"); // Still have no idea why it is called LepLep...
  tree->bookBranch<BaseTree::BranchType_vfloat_t>("LepPt"); // See? This is not called LepLep!
  tree->bookBranch<BaseTree::BranchType_vfloat_t>("LepEta"); // See? This is not called LepLep either!

  tree->getValRef("CRflag", CRFlagRef);
  tree->getValRef("Z1Flav", Z1FlavRef);
  tree->getValRef("LepLepId", LepIdRef);
  tree->getValRef("LepPt", LepPtRef);
  tree->getValRef("LepEta", LepEtaRef);

  CRFlagRegistry[tree]=CRFlagRef;
  Z1FlavRegistry[tree]=Z1FlavRef;
  LepIdRegistry[tree]=LepIdRef;
  LepPtRegistry[tree]=LepPtRef;
  LepEtaRegistry[tree]=LepEtaRef;
}

float ZXFakeRateHandler::getFakeRateWeight(CJLSTTree* tree) const{
  float res=0;
  unordered_map<CJLSTTree*, int*>::const_iterator it_CRFlag=CRFlagRegistry.find(tree);
  if (it_CRFlag!=CRFlagRegistry.cend()){
    unordered_map<CJLSTTree*, short*>::const_iterator it_Z1Flav=Z1FlavRegistry.find(tree);
    unordered_map<CJLSTTree*, vector<short>* const*>::const_iterator it_LepId=LepIdRegistry.find(tree);
    unordered_map<CJLSTTree*, vector<float>* const*>::const_iterator it_LepPt=LepPtRegistry.find(tree);
    unordered_map<CJLSTTree*, vector<float>* const*>::const_iterator it_LepEta=LepEtaRegistry.find(tree);

    int const& CRFlag = *(it_CRFlag->second);
    short const& Z1Flav = *(it_Z1Flav->second);
    vector<short>* const& LepId = *(it_LepId->second);
    vector<float>* const& LepPt = *(it_LepPt->second);
    vector<float>* const& LepEta = *(it_LepEta->second);
    short Z2Flav = 1;
    if (LepId && LepPt && LepEta){
      if (LepId->size()>=3) res=1;
      for (unsigned int ilep=2; ilep<LepId->size(); ilep++){
        short const& lepid = LepId->at(ilep);
        res *= eval(CRFlag, Z1Flav, lepid, LepPt->at(ilep), LepEta->at(ilep));
        Z2Flav *= lepid;
      }
    }
    else{
      MELAerr << "ZXFakeRateHandler::getFakeRateWeight: Something went wrong! Lepton vector references are null." << endl;
      assert(0);
    }
    short const absZ1Flav=std::abs(Z1Flav);
    short const absZ2Flav=std::abs(Z2Flav);
    float scale=1;
    if (FRMethod==mSS){
      if (theDataPeriod=="2016"){
        /*
        // FROM RUN 1
        if (absZ1Flav==121 && absZ2Flav==121) scale=0.97;
        else if (absZ1Flav==169 && absZ2Flav==121) scale=0.98;
        else if (absZ1Flav==121 && absZ2Flav==169) scale=1.30;
        else if (absZ1Flav==169 && absZ2Flav==169) scale=1.22;
        */
        if (absZ1Flav==121 && absZ2Flav==121) scale=1.65436;
        else if (absZ1Flav==169 && absZ2Flav==121) scale=1.56623;
        else if (absZ1Flav==121 && absZ2Flav==169) scale=1.07378;
        else if (absZ1Flav==169 && absZ2Flav==169) scale=0.99124;
      }
      else if (theDataPeriod=="2017"){
        if (absZ1Flav==121 && absZ2Flav==121) scale=1.39427;
        else if (absZ1Flav==169 && absZ2Flav==121) scale=1.38892;
        else if (absZ1Flav==121 && absZ2Flav==169) scale=0.98636;
        else if (absZ1Flav==169 && absZ2Flav==169) scale=1.02401;
      }
      else{
        MELAerr << "ZXFakeRateHandler::getFakeRateWeight: Data period " << theDataPeriod << " has no OS/SS scales implemented!" << endl;
        assert(0);
      }
    }
    res *= scale;
  }
  return res;
}

TString ZXFakeRateHandler::TranslateFakeRateMethodToString(ZXFakeRateHandler::FakeRateMethod FRMethod_){
  if (FRMethod_==ZXFakeRateHandler::mSS) return "SS";
  else if (FRMethod_==ZXFakeRateHandler::mOS) return "OS";
  else return "";
}
ZXFakeRateHandler::FakeRateMethod ZXFakeRateHandler::TranslateFakeRateMethodToEnum(TString FRMethodName_){
  if (FRMethodName_.Contains("SS")) return ZXFakeRateHandler::mSS;
  else if (FRMethodName_.Contains("OS")) return ZXFakeRateHandler::mOS;
  else return ZXFakeRateHandler::NFakeRateMethods;
}
