#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "HelperFunctions.h"
#include "SampleHelpers.h"
#include "CJLSTSet.h"
#include "BaseTreeLooper.h"
#include "DiscriminantClasses.h"
#include "CategorizationHelpers.h"
#include "ReweightingBuilder.h"
#include "ReweightingFunctions.h"
#include "ExtendedBinning.h"
#include "MELAStreamHelpers.hh"
#include "Mela.h"


#ifndef doDebugKD
#define doDebugKD true
#endif
#ifndef doDebugKDExt
#define doDebugKDExt false
#endif

using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace CategorizationHelpers;
using namespace DiscriminantClasses;
using namespace ReweightingFunctions;
using namespace MELAStreamHelpers;


class EventAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool runEvent(CJLSTTree* tree, SimpleEntry& product);

public:
  float infTrackingVal;
  float supTrackingVal;

  EventAnalyzer(Channel channel_, Category category_) : BaseTreeLooper(), channel(channel_), category(category_){}
  EventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) : BaseTreeLooper(inTree), channel(channel_), category(category_){}
  EventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) : BaseTreeLooper(inTreeList), channel(channel_), category(category_){}
  EventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) : BaseTreeLooper(inTreeSet), channel(channel_), category(category_){}

};
bool EventAnalyzer::runEvent(CJLSTTree* tree, SimpleEntry& product){
  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    //float& GenHMass = *(valfloats["GenHMass"]);
    product.trackingval = ZZMass; product.setNamedVal("ZZMass", ZZMass);

    // Check if trackingval is between the range requested
    validProducts &= (product.trackingval>=this->infTrackingVal && product.trackingval<this->supTrackingVal);

    // Construct the weights
    float pure_reco_wgt = (*(valfloats["dataMCWeight"]))*(*(valfloats["trigEffWeight"]));
    float wgt = pure_reco_wgt;
    wgt *= tree->getAssociatedSet()->getPermanentWeight(tree);

    const unsigned int nCheckWeights=3;
    const TString strCheckWeights[nCheckWeights]={
      "KFactor_QCD_ggZZ_Nominal","KFactor_EW_qqZZ","KFactor_QCD_qqZZ_M"
    };
    for (unsigned int icw=0; icw<nCheckWeights; icw++){
      auto cwit=valfloats.find(strCheckWeights[icw]);
      if (cwit!=valfloats.cend()) wgt *= *(cwit->second);
    }

    auto pugenhep_it = Rewgtbuilders.find("PUGenHEPRewgt");
    if (pugenhep_it!=Rewgtbuilders.cend()){
      auto& pugenheprewgtBuilder = pugenhep_it->second;
      float pugenhep_wgt_sum = pugenheprewgtBuilder->getSumPostThresholdWeights(tree);
      float pugenhep_wgt = (pugenhep_wgt_sum!=0. ? pugenheprewgtBuilder->getPostThresholdWeight(tree)/pugenhep_wgt_sum : 0.); // Normalized to unit
      wgt *= pugenhep_wgt;
    }

    auto melarewgt_it = Rewgtbuilders.find("MELARewgt");
    if (melarewgt_it!=Rewgtbuilders.cend()){
      auto& melarewgtBuilder = melarewgt_it->second;
      float mela_wgt_sum = melarewgtBuilder->getSumPostThresholdWeights(tree);
      float mela_wgt = (mela_wgt_sum!=0. ? melarewgtBuilder->getPostThresholdWeight(tree)/mela_wgt_sum : 0.); // Normalized to unit
      mela_wgt *= melarewgtBuilder->getNormComponent(tree);
      wgt *= mela_wgt;
    }

    product.weight = wgt; product.setNamedVal("weight", wgt);
    if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.){
      if (wgt!=0.){
        MELAerr << "EventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded at mass " << ZZMass << " for tree " << tree->sampleIdentifier << "." << endl;
        exit(1);
      }
      validProducts=false;
    }

    // Compute the KDs
    // Reserve the special DjjVBF, DjjZH and DjjWH discriminants
    float DjjVBF=-1;
    float DjjWH=-1;
    float DjjZH=-1;
    float KDreq=-999;
    for (auto it=KDbuilders.cbegin(); it!=KDbuilders.cend(); it++){
      auto& KDbuilderpair = it->second;
      auto& KDbuilder = KDbuilderpair.first;
      auto& strKDVarsList = KDbuilderpair.second;
      vector<float> KDBuildVals; KDBuildVals.reserve(strKDVarsList.size());
      for (auto const& s : strKDVarsList) KDBuildVals.push_back(*(valfloats[s]));
      float KD = KDbuilder->update(KDBuildVals, ZZMass);
      validProducts &= !(std::isnan(KD) || std::isinf(KD));

      if (it->first=="DjjVBF") DjjVBF=KD;
      else if (it->first=="DjjZH") DjjZH=KD;
      else if (it->first=="DjjWH") DjjWH=KD;
      else if (it->first=="KD"){
        product.setNamedVal(it->first, KD);
        KDreq=KD;
      }
    }
    validProducts &= (KDreq>=0.);

    // Category check
    Category catFound = CategorizationHelpers::getCategory(DjjVBF, DjjZH, DjjWH, false);
    validProducts &= (category==Inclusive || category==catFound);

    // Channel check
    validProducts &= SampleHelpers::testChannel(channel, *(valshorts["Z1Flav"]), *(valshorts["Z2Flav"]));
  }

  return validProducts;
}


void constructSigSamples(TString sampleType, float sqrts, const std::vector<TString>& KDvars, CJLSTSet*& theSampleSet){
  vector<TString> samples;
  if (sampleType=="ggHPowheg" || sampleType=="Sig") samples.push_back("gg_Sig_POWHEG");
  if (sampleType=="VBFPowheg" || sampleType=="Sig") samples.push_back("VBF_Sig_POWHEG");
  if (sampleType=="ZHPowheg" || sampleType=="Sig") samples.push_back("ZH_Sig_POWHEG");
  if (sampleType=="WHPowheg" || sampleType=="Sig") samples.push_back("WH_Sig_POWHEG");
  if (sampleType=="WplusHPowheg" || sampleType=="Sig") samples.push_back("WplusH_Sig_POWHEG");
  if (sampleType=="WminusHPowheg" || sampleType=="Sig") samples.push_back("WminusH_Sig_POWHEG");
  if (sampleType=="ggHMCFM" || sampleType=="Sig") samples.push_back("gg_Sig_SM_MCFM");

  vector<TString> samplesList;
  SampleHelpers::getSamplesList(sqrts, samples, samplesList);
  theSampleSet = new CJLSTSet(samplesList);

  theSampleSet->bookXS(); // "xsec"
  theSampleSet->bookOverallEventWgt(); // Gen weigts "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("GenHMass", 0);
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for KDs
    for (auto& v:KDvars) tree->bookBranch<float>(v, 0);
    // Variables for reweighting on a case-by-case basis
    tree->bookBranch<float>("KFactor_QCD_ggZZ_Nominal", 1);
  }
}
void constructBkgSamples(TString sampleType, float sqrts, const std::vector<TString>& KDvars, CJLSTSet*& theSampleSet){
  vector<TString> samples;
  if (sampleType=="qqBkg" || sampleType=="Bkg") samples.push_back("qq_Bkg_Combined");
  if (sampleType=="ggBkg" || sampleType=="Bkg") samples.push_back("gg_Bkg_MCFM");

  vector<TString> samplesList;
  SampleHelpers::getSamplesList(sqrts, samples, samplesList);
  theSampleSet = new CJLSTSet(samplesList);

  theSampleSet->bookXS(); // "xsec"
  theSampleSet->bookOverallEventWgt(); // Gen weigts "PUWeight", "genHEPMCweight" and reco weights "dataMCWeight", "trigEffWeight"
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    // Book common variables needed for analysis
    tree->bookBranch<float>("GenHMass", 0);
    tree->bookBranch<float>("ZZMass", -1);
    tree->bookBranch<short>("Z1Flav", 0);
    tree->bookBranch<short>("Z2Flav", 0);
    // Variables for KDs
    for (auto& v:KDvars) tree->bookBranch<float>(v, 0);
    tree->bookBranch<float>("KFactor_QCD_ggZZ_Nominal", 1);
    tree->bookBranch<float>("KFactor_EW_qqZZ", 1);
    tree->bookBranch<float>("KFactor_QCD_qqZZ_M", 1);
  }
}
void constructSamples(TString sampleType, float sqrts, const std::vector<TString>& KDvars, CJLSTSet*& theSampleSet){
  if (sampleType=="qqBkg" || sampleType=="ggBkg" || sampleType=="Bkg") return constructBkgSamples(sampleType, sqrts, KDvars, theSampleSet);
  else return constructSigSamples(sampleType, sqrts, KDvars, theSampleSet);
}

class KDConstantByMass{
protected:
  float sqrts;
  TString strKD;
  unsigned int nstepsiter;
  std::vector<CJLSTSet::NormScheme> NormSchemeA;
  std::vector<CJLSTSet::NormScheme> NormSchemeB;

  void LoopForConstant(
    vector<SimpleEntry>(&index)[2],
    vector<unsigned int>(&indexboundaries)[2],
    TProfile*& px,
    TH1F*& hrec
  );

public:
  KDConstantByMass(float sqrts_, TString strKD_);

  void setNStepsIter(unsigned int nsteps){ nstepsiter=nsteps; }
  void setNormSchemeA(std::vector<CJLSTSet::NormScheme> scheme){ NormSchemeA=scheme; }
  void setNormSchemeB(std::vector<CJLSTSet::NormScheme> scheme){ NormSchemeB=scheme; }

  void run(
    vector<TString> strSamples[2], vector<vector<TString>> strMelaWgts[2],
    SampleHelpers::Channel channel, CategorizationHelpers::Category category,
    unsigned int divisor, const bool writeFinalTree, vector<pair<vector<float>, pair<float, float>>>* manualboundary_validity_pairs=0
  );

};

KDConstantByMass::KDConstantByMass(float sqrts_, TString strKD_) :
  sqrts(sqrts_), strKD(strKD_),
  nstepsiter(100)
{}


///////////////////
// Event helpers //
///////////////////
void KDConstantByMass::LoopForConstant(
  vector<SimpleEntry>(&index)[2],
  vector<unsigned int>(&indexboundaries)[2],
  TProfile*& px,
  TH1F*& hrec
){
  cout << "Begin KDConstantByMass::LoopForConstant" << endl;

  int nbins = indexboundaries[0].size()-1;

  for (int bin=0; bin<nbins; bin++){
    cout << "Bin " << bin << " / " << nbins << " is now being scrutinized..." << endl;

    // First find the average KD in each sample
    float sumKD[2]={ 0 }; float sumWgt[2]={ 0 }; float avgKD[2]={ 0 };
    std::vector<SimpleEntry>::iterator it_begin[2], it_end[2];
    for (unsigned int ih=0; ih<2; ih++){
      unsigned int const& evlow = indexboundaries[ih].at(bin);
      unsigned int const& evhigh = indexboundaries[ih].at(bin+1);
      it_begin[ih] = index[ih].begin()+evlow;
      it_end[ih] = index[ih].begin()+evhigh;
      cout << " - Scanning events [ " << evlow << " , " << evhigh << " ) for sample set " << ih << endl;
    }
    for (unsigned int ih=0; ih<2; ih++){
      for (std::vector<SimpleEntry>::iterator it_inst=it_begin[ih]; it_inst!=it_end[ih]; it_inst++){
        float const& KD = it_inst->namedfloats["KD"];
        float const& varTrack = it_inst->trackingval;
        float const& weight = it_inst->weight;
        sumKD[ih] += KD*weight;
        sumWgt[ih] += weight;

        if (px->GetXaxis()->GetBinLowEdge(bin+1)>varTrack || px->GetXaxis()->GetBinLowEdge(bin+2)<=varTrack) cerr
          << "Something terrible happened! " << varTrack << " is outside ["
          << px->GetXaxis()->GetBinLowEdge(bin+1) << " , " << px->GetXaxis()->GetBinLowEdge(bin+2)
          << "]" << endl;

        px->Fill(varTrack, varTrack, weight);
      }
      cout << " - Sum of weights for sample set " << ih << " = " << sumWgt[ih] << endl;
      avgKD[ih]=sumKD[ih]/sumWgt[ih];
      cout << " - Average KD for sample set " << ih << " = " << avgKD[ih] << endl;
    }

    float marginlow=1;
    float marginhigh=20;
    float Cfound=0;
    float centralConstant = (avgKD[0]+avgKD[1])*0.5; centralConstant = 1./(1./centralConstant-1.);
    unsigned int it=0;
    while (true){
      cout << " - Iteration " << it << " with margins = " << marginlow << ", " << marginhigh << endl;
      cout << "   - Checking c = [ " << centralConstant*(1.-marginlow) << " , " << centralConstant*(1.+marginhigh) << " ]" << endl;
      float mindiff=1;
      short sgnMinDiff=0;
      float finalFraction[2]={ 2, 2 };
      unsigned int nsteps = nstepsiter;
      if (it==0) nsteps*=1000;
      else if (it==1) nsteps*=100;
      else if (it==2) nsteps*=10;

      for (unsigned int step=0; step<=nsteps; step++){
        HelperFunctions::progressbar(step, nsteps);
        float testC = centralConstant*((1.-marginlow) + (marginhigh+marginlow)*(float(step))/((float) nsteps));
        if (testC<=0.) continue;

        float sumWgtAll[2]={ 0 };
        float sumWgtHalf[2]={ 0 };
        for (unsigned int ih=0; ih<2; ih++){
          for (std::vector<SimpleEntry>::iterator it_inst=it_begin[ih]; it_inst!=it_end[ih]; it_inst++){
            float const& KDold = it_inst->namedfloats["KD"];
            float const& weight = it_inst->weight;
            if (KDold==-999.) continue;
            float KD = KDold/(KDold+(1.-KDold)*testC);
            if (std::isnan(KD) || std::isinf(KD)){
              cerr << "Something went terribly wrong! KD is " << KD << endl;
              continue;
            }
            else if (KD<0.){
              cerr << "KD is invalid (" << KD << ")" << endl;
              continue;
            }
            sumWgtAll[ih] += weight;
            if (
              (KD>=0.5 && ih==0)
              ||
              (KD<0.5 && ih==1)
              ) sumWgtHalf[ih] += weight;
          }
          sumWgtHalf[ih]=sumWgtHalf[ih]/sumWgtAll[ih];
        }

        float sumWgtHalfDiff=sumWgtHalf[0]-sumWgtHalf[1];
        float absSumWgtHalfDiff=fabs(sumWgtHalfDiff);
        short sgnSumWgtHalfDiff = TMath::Sign(float(1), sumWgtHalfDiff);
        if (mindiff>absSumWgtHalfDiff){
          finalFraction[0] = sumWgtHalf[0];
          finalFraction[1] = sumWgtHalf[1];
          mindiff=absSumWgtHalfDiff;
          if (sgnMinDiff==0 && sgnSumWgtHalfDiff!=0) sgnMinDiff=sgnSumWgtHalfDiff;
          Cfound=testC;
        }
        if (sgnMinDiff!=sgnSumWgtHalfDiff) break;
      }
      cout << "  - New c found = " << Cfound << " (old was " << centralConstant << ")" << endl;
      cout << "  - Final fractions were = " << finalFraction[0] << " , " << finalFraction[1] << " (sign of difference: " << sgnMinDiff << ")" << endl;
      if (fabs(Cfound/centralConstant-1)<1e-4 && it>0) break;
      centralConstant=Cfound;
      if (it>2){
        marginhigh /= float(nsteps)/10.;
        marginlow /= float(nsteps)/10.;
      }
      else{
        marginhigh /= 5.;
        marginlow /= 2.;
      }
      Cfound=0;
      it++;
    }

    hrec->SetBinContent(bin+1, centralConstant);
  }

  cout << "End KDConstantByMass::LoopForConstant" << endl;
}
void KDConstantByMass::run(
  vector<TString> strSamples[2], vector<vector<TString>> strMelaWgts[2],
  SampleHelpers::Channel channel, CategorizationHelpers::Category category,
  unsigned int divisor, const bool writeFinalTree, vector<pair<vector<float>, pair<float, float>>>* manualboundary_validity_pairs
){
  if (strKD=="") return;

  cout << "Begin KDConstantByMass::run" << endl;
  for (unsigned int ih=0; ih<2; ih++){
    assert(strSamples[ih].size()==strMelaWgts[ih].size());
  }

  DiscriminantClasses::Type KDtype = DiscriminantClasses::getKDType(strKD);
  vector<TString> KDvars = DiscriminantClasses::getKDVars(KDtype);
  Discriminant* KDbuilder = constructKDFromType(KDtype);
  if (!KDbuilder) return;

  vector<SimpleEntry> index[2];

  gSystem->Exec("mkdir -p ./output/KDConstants");
  TString coutput = Form("KDConstant_m4l_%s", strKD.Data());
  if (channel!=SampleHelpers::NChannels) coutput += Form("_%s", SampleHelpers::getChannelName(channel).Data());
  if (category!=CategorizationHelpers::Inclusive) coutput += Form("_%s", CategorizationHelpers::getCategoryName(category).Data());
  if (sqrts>0.) coutput += Form("%.0fTeV", sqrts);
  TFile* foutput = TFile::Open(Form("./output/KDConstants/%s%s", coutput.Data(), ".root"), "recreate");

  int nEntries[2]={ 0 };
  float infimum=0;
  float supremum=sqrts*1000.;

  for (unsigned int ih=0; ih<2; ih++){
    vector<CJLSTSet*> theSets;
    theSets.assign(strSamples[ih].size(), nullptr);
    for (unsigned int ihs=0; ihs<strSamples[ih].size(); ihs++) constructSamples(strSamples[ih].at(ihs), 13, KDvars, theSets.at(ihs));
    auto melawgtcollit=strMelaWgts[ih].begin();
    std::vector<CJLSTSet::NormScheme> const& schemeSet = (ih==0 ? NormSchemeA : NormSchemeB);
    unsigned int iset=0;
    for (auto& theSet:theSets){
      for (auto& tree:theSet->getCJLSTTreeList()){ for (auto& strWgt:(*melawgtcollit)) tree->bookBranch<float>(strWgt, 0); }
      melawgtcollit++;

      theSet->setPermanentWeights((schemeSet.empty() ? CJLSTSet::NormScheme_None : schemeSet.at(iset)), true, false);

      for (auto& tree:theSet->getCJLSTTreeList()) tree->silenceUnused(); // Will no longer book another branch
      iset++;
    }

    // Setup GenHMass binning
    // Binning for PUGenHepRewgt
    ExtendedBinning GenHMassInclusiveBinning("GenHMass");

    // Construct PUGenHepRewgt
    vector<TString> strReweightingWeigths;
    strReweightingWeigths.push_back("PUWeight");
    strReweightingWeigths.push_back("genHEPMCweight");
    ReweightingBuilder* pugenheprewgtBuilder = new ReweightingBuilder(strReweightingWeigths, getSimpleWeight);
    pugenheprewgtBuilder->setWeightBinning(GenHMassInclusiveBinning);
    for (auto& theSet:theSets){ for (auto& tree:theSet->getCJLSTTreeList()) pugenheprewgtBuilder->setupWeightVariables(tree, 1.); }

    melawgtcollit=strMelaWgts[ih].begin();
    for (auto& theSet:theSets){
      ReweightingBuilder* melarewgtBuilder = nullptr;
      // Construct MELARewgt for each set
      if (!melawgtcollit->empty()){
        // Binning for MELARewgt
        ExtendedBinning GenHMassBinning("GenHMass");
        for (unsigned int is=0; is<theSet->getCJLSTTreeList().size()-1; is++){
          if (theSet->getCJLSTTreeList().at(is)->MHVal>0. && theSet->getCJLSTTreeList().at(is+1)->MHVal>0.){
            GenHMassBinning.addBinBoundary(
              0.5*(theSet->getCJLSTTreeList().at(is)->MHVal + theSet->getCJLSTTreeList().at(is+1)->MHVal)
            );
          }
        }
        if (GenHMassBinning.isValid()){
          GenHMassBinning.addBinBoundary(0);
          GenHMassBinning.addBinBoundary(sqrts*1000.);
        }
        ReweightingBuilder* melarewgtBuilder = new ReweightingBuilder(*melawgtcollit, ReweightingFunctions::getSimpleWeight);
        melarewgtBuilder->setDivideByNSample(true);
        melarewgtBuilder->setWeightBinning(GenHMassBinning);
        for (auto& tree:theSet->getCJLSTTreeList()) melarewgtBuilder->setupWeightVariables(tree, 0.999, 10);
      }

      EventAnalyzer theAnalyzer(theSet, channel, category);
      theAnalyzer.infTrackingVal=infimum;
      theAnalyzer.supTrackingVal=supremum;
      // Book common variables needed for analysis
      theAnalyzer.addConsumed<float>("dataMCWeight");
      theAnalyzer.addConsumed<float>("trigEffWeight");
      theAnalyzer.addConsumed<float>("GenHMass");
      theAnalyzer.addConsumed<float>("ZZMass");
      theAnalyzer.addConsumed<short>("Z1Flav");
      theAnalyzer.addConsumed<short>("Z2Flav");
      // Add discriminant builders
      theAnalyzer.addDiscriminantBuilder("KD", KDbuilder, KDvars);
      // Add reweighting builders
      theAnalyzer.addReweightingBuilder("PUGenHEPRewgt", pugenheprewgtBuilder);
      theAnalyzer.addReweightingBuilder("MELARewgt", melarewgtBuilder);

      // Loop
      theAnalyzer.setExternalProductList(&(index[ih]));
      theAnalyzer.loop(true, false, true);

      delete melarewgtBuilder;
      melawgtcollit++;
    }
    delete pugenheprewgtBuilder;
    for (auto& theSet:theSets) delete theSet;

    float firstVal=1000.*sqrts;
    float lastVal=0;
    cout << "Determining min/max for set " << ih << " (size=" << index[ih].size() << ")" << endl;
    for (auto const& ev:index[ih]){
      firstVal=std::min(firstVal, ev.trackingval);
      lastVal=std::max(lastVal, ev.trackingval);
    }
    firstVal = (float) ((int) firstVal); firstVal -= (float) (((int) firstVal)%10);
    lastVal = (float) ((int) (lastVal+0.5)); lastVal += (float) (10-((int) lastVal)%10);
    infimum = std::max(firstVal, infimum);
    supremum = std::min(lastVal, supremum);
    cout << "Unmodified Nproducts: " << index[ih].size() << ", firstVal: " << firstVal << ", lastVal: " << lastVal << endl;
  }
  for (unsigned int ih=0; ih<2; ih++){
    cout << "Cropping events for set " << ih << " based on min/max = " << infimum << " / " << supremum << endl;
    SimpleEntry::cropByTrueVal(index[ih], infimum, supremum);
    cout << "Sorting..." << endl;
    std::sort(index[ih].begin(), index[ih].end());
    nEntries[ih]=index[ih].size();
    cout << "Nentries remaining = " << nEntries[ih] << " | var = [ " << infimum << " , " << supremum << " ]" << endl;
  }

  vector<unsigned int> indexboundaries[2];
  for (unsigned int ih=0; ih<2; ih++) indexboundaries[ih].push_back(0);
  {
    unsigned int iit[2]={ divisor, divisor };
    while (iit[0]<index[0].size() && iit[1]<index[1].size()){
      if (
        (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
        <
        (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
        ){
        while (
          (iit[0]+1)<index[0].size() &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          <
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          &&
          (
          (index[0].at(iit[0]+1).trackingval+index[0].at(iit[0]).trackingval)*0.5
          <
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          ) iit[0]++;
      }
      else if (
        (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          >
        (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
        ){
        while (
          (iit[1]+1)<index[1].size() &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          >
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          >
          (index[1].at(iit[1]+1).trackingval+index[1].at(iit[1]).trackingval)*0.5
          )
          ) iit[1]++;
      }

      if (
        (index[0].size()-iit[0])<divisor
        ||
        (index[1].size()-iit[1])<divisor
        ) break;

      for (unsigned int ih=0; ih<2; ih++){ indexboundaries[ih].push_back(iit[ih]); iit[ih] += divisor; }
    }
  }
  for (unsigned int ih=0; ih<2; ih++) indexboundaries[ih].push_back(index[ih].size());
  for (unsigned int ih=0; ih<2; ih++) cout << "Final size of indexboundaries[" << ih << "] = " << indexboundaries[ih].size() << endl;

  unsigned int nbins=indexboundaries[0].size()-1;
  vector<float> binboundarylist;
  binboundarylist.push_back(infimum);
  for (unsigned int ix=1; ix<nbins; ix++){
    float binboundary = max(
      (index[0].at(indexboundaries[0].at(ix)-1).trackingval+index[0].at(indexboundaries[0].at(ix)).trackingval)*0.5
      ,
      (index[1].at(indexboundaries[1].at(ix)-1).trackingval+index[1].at(indexboundaries[1].at(ix)).trackingval)*0.5
    );

    cout << "Initial bin boundary for bin " << ix << ": " << binboundary << endl;

    bool skip=false;
    if (manualboundary_validity_pairs!=0){
      for (auto const& mbvpair : (*manualboundary_validity_pairs)){
        const pair<float, float>& valrange = mbvpair.second;
        skip = (
          (valrange.first<0 || valrange.first<=binboundary)
          &&
          (valrange.second<0 || valrange.second>binboundary)
          );
        if (skip) break;
      }
    }
    if (!skip) binboundarylist.push_back(binboundary);
  }
  binboundarylist.push_back(supremum);

  if (manualboundary_validity_pairs!=0){
    for (auto& mbvpair : (*manualboundary_validity_pairs)){
      vector<float>& bvals = mbvpair.first;
      for (unsigned int ib=0; ib<bvals.size(); ib++){
        float& bval = bvals.at(ib);
        cout << "Adding manual boundary " << bval << endl;
        addByLowest<float>(binboundarylist, bval, true);
      }
    }
  }

  nbins = binboundarylist.size()-1;
  float* binning = new float[nbins+1];
  for (unsigned int ix=0; ix<=nbins; ix++){
    binning[ix] = binboundarylist[ix];
    cout << "Boundary (" << ix << ")= " << binning[ix] << endl;
  }

  // Recalibrate index boundaries
  if (manualboundary_validity_pairs!=0){
    for (unsigned int ih=0; ih<2; ih++){
      indexboundaries[ih].clear();
      indexboundaries[ih].push_back(0);
      unsigned int ix=1;
      for (unsigned int ev=1; ev<index[ih].size(); ev++){
        if (ix==nbins) break;
        if (index[ih].at(ev).trackingval>=binning[ix]){
          indexboundaries[ih].push_back(ev);
          ix++;
        }
      }
      indexboundaries[ih].push_back(index[ih].size());
    }
  }

  foutput->cd();

  TH1F* h_varTrack_Constant = new TH1F("varReco_Constant", "", nbins, binning); h_varTrack_Constant->Sumw2();
  TProfile* p_varTrack = new TProfile("avg_varReco", "", nbins, binning); p_varTrack->Sumw2();
  delete[] binning;

  if (writeFinalTree){
    for (unsigned int bin=0; bin<nbins; bin++){
      for (unsigned int ih=0; ih<2; ih++){
        TTree* theFinalTree = new TTree(Form("Sample%i_Bin%i", ih, bin), "");
        SimpleEntry::writeToTree(index[ih].cbegin()+indexboundaries[ih].at(bin), index[ih].cbegin()+indexboundaries[ih].at(bin+1), theFinalTree);
        foutput->WriteTObject(theFinalTree);
        delete theFinalTree;
      }
    }
  }

  LoopForConstant(
    index, indexboundaries,
    p_varTrack,
    h_varTrack_Constant
  );

  TGraphErrors* gr = makeGraphFromTH1(p_varTrack, h_varTrack_Constant, "gr_varReco_Constant");
  foutput->WriteTObject(p_varTrack);
  foutput->WriteTObject(h_varTrack_Constant);
  foutput->WriteTObject(gr);
  delete gr;
  delete h_varTrack_Constant;
  delete p_varTrack;
  foutput->Close();

  delete KDbuilder;

  NormSchemeA.clear();
  NormSchemeB.clear();
  cout << "End KDConstantByMass::run" << endl;
}

/*
SPECIFIC COMMENT:
- Multiplies by Pmjj
*/
void getKDConstant_DjjZH(float sqrts=13, const bool writeFinalTree=false){
  float divisor=20000;
  TString strKD="DjjZH";

  vector<TString> strSamples[2];
  strSamples[0].push_back("ZHPowheg");
  strSamples[1].push_back("ggHPowheg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }
  {
    pair<float, float> valrange(230, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(245);
    //manualboundaries.push_back(300);
    //manualboundaries.push_back(500);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::NChannels, CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}
void getKDConstant_DjjWH(float sqrts=13, const bool writeFinalTree=false){
  float divisor=20000;
  TString strKD="DjjWH";

  vector<TString> strSamples[2];
  strSamples[0].push_back("WHPowheg");
  strSamples[1].push_back("ggHPowheg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }
  {
    pair<float, float> valrange(195, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(196);
    manualboundaries.push_back(209);
    manualboundaries.push_back(300);
    manualboundaries.push_back(500);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  {
    std::vector<CJLSTSet::NormScheme> NormSchemeA;
    NormSchemeA.assign(strSamples[0].size(), CJLSTSet::NormScheme_None);
    NormSchemeA.at(0)=CJLSTSet::NormScheme_XsecOverNgen_RelRenormToSumNgen;
    constProducer.setNormSchemeA(NormSchemeA);
  }
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::NChannels, CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_DjjVBF(float sqrts=13, const bool writeFinalTree=false){
  float divisor=40000;
  TString strKD="DjjVBF";

  vector<TString> strSamples[2];
  strSamples[0].push_back("VBFPowheg");
  strSamples[1].push_back("ggHPowheg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }
  {
    pair<float, float> valrange(750, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(770);
    manualboundaries.push_back(850);
    manualboundaries.push_back(1100);
    manualboundaries.push_back(1400);
    manualboundaries.push_back(2100);
    manualboundaries.push_back(2650);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::NChannels, CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_DjVBF(float sqrts=13, const bool writeFinalTree=false){
  float divisor=40000;
  TString strKD="DjVBF";

  vector<TString> strSamples[2];
  strSamples[0].push_back("VBFPowheg");
  strSamples[1].push_back("ggHPowheg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }
  {
    pair<float, float> valrange(1000, 2000);
    vector<float> manualboundaries;
    manualboundaries.push_back(1400);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
    ));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::NChannels, CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_Dbkgkin(TString const strchannel, float sqrts=13, const bool writeFinalTree=false){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;

  float divisor=21000;
  if (strchannel=="2l2l" || strchannel=="2e2mu") divisor = 50000;
  TString strKD="Dbkgkin";

  vector<TString> strSamples[2];
  strSamples[0].push_back("ggHPowheg");
  strSamples[0].push_back("ggHMCFM");
  strSamples[0].push_back("VBFPowheg");
  strSamples[1].push_back("qqBkg");
  strSamples[1].push_back("ggBkg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());
  // Reweight ggBkg decay kinematics to qqBkg
  strMelaWgts[1].at(0).push_back("xsec");
  strMelaWgts[1].at(0).push_back("p_Gen_QQB_BKG_MCFM");

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 142);
    vector<float> manualboundaries;
    manualboundaries.push_back(75); manualboundaries.push_back(85);
    manualboundaries.push_back(89); manualboundaries.push_back(93); manualboundaries.push_back(96);
    manualboundaries.push_back(100); manualboundaries.push_back(105); manualboundaries.push_back(110); manualboundaries.push_back(115);
    manualboundaries.push_back(120); manualboundaries.push_back(123); manualboundaries.push_back(135);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }
  {
    pair<float, float> valrange(600, 2500);
    vector<float> manualboundaries;
    manualboundaries.push_back(700); manualboundaries.push_back(900); manualboundaries.push_back(1100); manualboundaries.push_back(1400); manualboundaries.push_back(1900);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  {
    std::vector<CJLSTSet::NormScheme> NormSchemeB;
    NormSchemeB.assign(strSamples[1].size(), CJLSTSet::NormScheme_None);
    NormSchemeB.at(0)=CJLSTSet::NormScheme_XsecOverNgen;
    constProducer.setNormSchemeB(NormSchemeB);
  }
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::getChannelFromName(strchannel), CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_Dggbkgkin(TString const strchannel, float sqrts=13, const bool writeFinalTree=false){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;

  float divisor=21000;
  if (strchannel=="2l2l" || strchannel=="2e2mu") divisor = 50000;
  TString strKD="Dggbkgkin";

  vector<TString> strSamples[2];
  strSamples[0].push_back("ggHPowheg");
  strSamples[0].push_back("ggHMCFM");
  strSamples[1].push_back("ggBkg");
  vector<vector<TString>> strMelaWgts[2]; for (unsigned int ih=0; ih<2; ih++) strMelaWgts[ih].assign(strSamples[ih].size(), vector<TString>());

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 142);
    vector<float> manualboundaries;
    manualboundaries.push_back(75); manualboundaries.push_back(85);
    manualboundaries.push_back(89); manualboundaries.push_back(93); manualboundaries.push_back(96);
    manualboundaries.push_back(100); manualboundaries.push_back(105); manualboundaries.push_back(110); manualboundaries.push_back(115);
    manualboundaries.push_back(120); manualboundaries.push_back(123); manualboundaries.push_back(135);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }
  {
    pair<float, float> valrange(600, 2500);
    vector<float> manualboundaries;
    manualboundaries.push_back(700); manualboundaries.push_back(900); manualboundaries.push_back(1100); manualboundaries.push_back(1400); manualboundaries.push_back(1900);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }

  KDConstantByMass constProducer(sqrts, strKD);
  constProducer.run(
    strSamples, strMelaWgts,
    SampleHelpers::getChannelFromName(strchannel), CategorizationHelpers::Inclusive,
    divisor, writeFinalTree, &manualboundary_validity_pairs
  );
}


void generic_SmoothKDConstantProducer(
  float sqrts, TString strname, TString strcustomselection,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  bool useFaithfulSlopeFirst, bool useFaithfulSlopeSecond,
  vector<pair<pair<double, double>, unsigned int>>* addpoints=0
){
  const double xmin=0;
  const double xmax=(sqrts>0 ? (double) sqrts*1000. : 15000.);

  gSystem->Exec("mkdir -p ./output/KDConstants");
  TString cinput = Form("KDConstant_m4l_%s", strname.Data());
  if (strcustomselection!="") cinput += Form("_%s", strcustomselection.Data());
  if (sqrts>0.) cinput += Form("%.0fTeV", sqrts);

  TFile* finput = TFile::Open(Form("./output/KDConstants/%s%s", cinput.Data(), ".root"), "read");
  TFile* foutput = TFile::Open(Form("./output/KDConstants/Smooth%s%s", cinput.Data(), ".root"), "recreate");
  foutput->cd();

  TGraphErrors* tg = (TGraphErrors*) finput->Get("gr_varReco_Constant");
  foutput->WriteTObject(tg);

  if (addpoints!=0){ for (auto& prange : *addpoints) addPointsBetween(tg, prange.first.first, prange.first.second, prange.second); }

  int n = tg->GetN();
  double* xx = tg->GetX();
  double* ex = tg->GetEX();
  double* yy = tg->GetY();
  double* ey = tg->GetEY();

  TSpline3* sp;
  sp = convertGraphToSpline3(tg, useFaithfulSlopeFirst, useFaithfulSlopeSecond);
  double tglow = xx[0];
  double tghigh = xx[tg->GetN()-1];
  TF1* lowFcn = lowf(sp, xmin, tglow, true);
  TF1* highFcn = highf(sp, tghigh, xmax, false);
  lowFcn->SetNpx((int) (tglow-xmin)*5);
  highFcn->SetNpx((int) (xmax-tghigh)*5);

  vector<pair<double, double>> points;
  for (double xval=xmin; xval<tglow; xval+=1){
    double yval = lowFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }
  for (int ix=0; ix<n; ix++){
    addByLowest<double, double>(points, xx[ix], yy[ix]);
  }
  int tghigh_int = ((int) ((tghigh+1.)/100.+0.5))*100;
  if (tghigh>=(double) tghigh_int) tghigh_int+=100;
  for (double xval=tghigh_int; xval<=xmax; xval+=100){
    double yval = highFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }

  int nn_new = points.size();
  cout << "Number of new points: " << nn_new-n << endl;
  double* xy_new[2];
  for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
  for (int ix=0; ix<nn_new; ix++){
    xy_new[0][ix] = points.at(ix).first;
    xy_new[1][ix] = points.at(ix).second;
  }

  delete highFcn;
  delete lowFcn;
  delete sp;

  TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
  tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
  foutput->WriteTObject(tg_updated);

  sp = convertGraphToSpline3(tg_updated);
  foutput->WriteTObject(sp);
  delete sp;
  delete tg_updated;
  for (unsigned int i=0; i<2; i++) delete[] xy_new[i];

  foutput->Close();
  finput->Close();
}


void SmoothKDConstantProducer_DjjZH(){
  generic_SmoothKDConstantProducer(
    13, "DjjZH", "",
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    false, false
  );
}

void SmoothKDConstantProducer_DjjWH(){
  generic_SmoothKDConstantProducer(
    13, "DjjWH", "",
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    false, false
  );
}

void SmoothKDConstantProducer_DjjVBF(){
  TString strprod="VBF";

  vector<pair<pair<double, double>, unsigned int>> addpoints;
  {
    pair<double, double> xminmax(1000, 3000); unsigned int nadd=5;
    pair<pair<double, double>, unsigned int> addsingle(xminmax, nadd); addpoints.push_back(addsingle);
  }

  generic_SmoothKDConstantProducer(
    13, Form("Djj%s", strprod.Data()), "",
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    //&getFcn_a0plusa1overX,
    true, true,
    &addpoints
  );
}

void SmoothKDConstantProducer_DjVBF(){
  TString strprod="VBF";
  generic_SmoothKDConstantProducer(
    13, Form("Dj%s", strprod.Data()), "",
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    true, true
  );
}

void SmoothKDConstantProducer_Dbkgkin(TString strchannel){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;
  generic_SmoothKDConstantProducer(
    13, Form("Dbkgkin_%s", strchannel.Data()), "",
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    //&getFcn_a0plusa1overXN<6>,
    true, false
  );
}


void testConstant(TString fname, TString cfname, float massmin, float massmax, float KDcutval){
  TFile* cFile = TFile::Open(cfname, "read");
  TSpline3* spC = (TSpline3*) cFile->Get("sp_gr_varReco_Constant");

  TFile* treeFile = TFile::Open(fname, "read");
  TH1F* htmp = (TH1F*)treeFile->Get("varReco_Constant");
  const unsigned int nbins = htmp->GetNbinsX();

  TFile testFile("test.root", "recreate");
  TH1F hA("hA", "", 50, 0, 1); hA.Sumw2();
  TH1F hB("hB", "", 50, 0, 1); hB.Sumw2();
  float KD, ZZMass, weight;
  for (unsigned int ih=0; ih<2; ih++){
    treeFile->cd();
    float sumWgt=0;
    float sumWgtWithCut=0;
    vector<TTree*> treeList;
    for (unsigned int bin=0; bin<nbins; bin++){ treeList.push_back((TTree*)treeFile->Get(Form("Sample%i_Bin%i", ih, bin))); }
    for (auto& tree:treeList){
      tree->SetBranchAddress("KD", &KD);
      tree->SetBranchAddress("ZZMass", &ZZMass);
      tree->SetBranchAddress("weight", &weight);
      for (int ev=0; ev<tree->GetEntries(); ev++){
        tree->GetEntry(ev);
        if (ZZMass>=massmin && ZZMass<massmax){
          sumWgt += weight;
          if (ih==0){
            hA.Fill(KD, weight);
          }
          else{
            hB.Fill(KD, weight);
          }
          if (KD>=KDcutval) sumWgtWithCut += weight;
        }
      }
    }
    cout << "Sample " << ih << " cut efficiency: " << sumWgtWithCut << " / " << sumWgt << " = " << sumWgtWithCut/sumWgt << endl;
  }
  hA.Scale(1./hA.Integral());
  hB.Scale(1./hB.Integral());
  testFile.WriteTObject(&hA);
  testFile.WriteTObject(&hB);

  // Plot ROC curve
  TGraph* tgROC=HelperFunctions::createROCFromDistributions(&hA, &hB, "ROC");
  tgROC->GetYaxis()->SetTitle("Sample A efficiency");
  tgROC->GetXaxis()->SetTitle("Sample B efficiency");
  testFile.WriteTObject(tgROC);
  delete tgROC;

  testFile.Close();
  treeFile->Close();
  cFile->Close();
}


/*
g-constants for AC
*/

TGraph* getSingleTGraph(TString fname){
  TDirectory* curdir = gDirectory;

  cout << "Opening file " << fname << endl;

  TFile* finput = TFile::Open(Form("JHUGenXsec/%s%s", fname.Data(), ".root"), "read");
  TGraph* tgold = (TGraph*) finput->Get("Graph");
  double* xx = tgold->GetX();
  double* yy = tgold->GetY();

  vector<pair<double, double>> xyinterm;
  for (int ip=0; ip<tgold->GetN(); ip++){
    if (std::isnan(xx[ip]) || std::isinf(xx[ip])) continue;
    if (std::isnan(yy[ip]) || std::isinf(yy[ip])) continue;
    xyinterm.push_back(pair<double, double>(xx[ip], yy[ip]));
  }
  TGraph* tginterm = makeGraphFromPair(xyinterm, "tginterm");
  TSpline3* spinterm = convertGraphToSpline3(tginterm);
  delete tginterm;

  vector<pair<double, double>> xyfinal;
  for (int ip=0; ip<tgold->GetN(); ip++){
    if (std::isnan(xx[ip]) || std::isinf(xx[ip])) continue;
    if (std::isnan(yy[ip]) || std::isinf(yy[ip])) yy[ip] = spinterm->Eval(xx[ip]);
    xyfinal.push_back(pair<double, double>(xx[ip], yy[ip]));
  }

  curdir->cd();
  TGraph* result = makeGraphFromPair(xyfinal, Form("tgfinal_%s", fname.Data()));

  finput->Close();
  return result;
}

void generic_gConstantProducer(TString strprod, TString strhypo, bool useproddec=false){
  TString strSM = "SM";
  if (strhypo=="L1Zgs") strSM += "_photoncut";

  TFile* foutput = TFile::Open(Form("gConstant_%s%s_%s%s", strprod.Data(), (useproddec ? "_ProdDec" : ""), strhypo.Data(), ".root"), "recreate");
  TGraph* tgSMhypo = getSingleTGraph(Form("%s_%s", strprod.Data(), strSM.Data()));
  TGraph* tgBSMhypo = getSingleTGraph(Form("%s_%s", strprod.Data(), strhypo.Data()));

  if (useproddec){
    TGraph* tgSMhypodec = getSingleTGraph(Form("HZZ2e2mu_%s", strSM.Data()));
    TGraph* tgBSMhypodec = getSingleTGraph(Form("HZZ2e2mu_%s", strhypo.Data()));

    TGraph* tgSMhyponew = multiplyTGraphs(tgSMhypo, tgSMhypodec);
    TGraph* tgBSMhyponew = multiplyTGraphs(tgBSMhypo, tgBSMhypodec);

    delete tgSMhypo; tgSMhypo=tgSMhyponew;
    delete tgBSMhypo; tgBSMhypo=tgBSMhyponew;
  }

  TGraph* result = divideTGraphs(tgSMhypo, tgBSMhypo, 0.5, 0.5);
  foutput->WriteTObject(result);
  TSpline3* spresult = convertGraphToSpline3(result);
  foutput->WriteTObject(spresult);
  cout << "Result of spline at 125 GeV = " << spresult->Eval(125) << endl;
  delete spresult; delete result;
  foutput->Close();
}

void gConstantProducer(){
  generic_gConstantProducer("WH_Sig_POWHEG", "g2");
  generic_gConstantProducer("WH_Sig_POWHEG", "g4");
  generic_gConstantProducer("WH_Sig_POWHEG", "L1");
  generic_gConstantProducer("ZH_Sig_POWHEG", "g2");
  generic_gConstantProducer("ZH_Sig_POWHEG", "g4");
  generic_gConstantProducer("ZH_Sig_POWHEG", "L1");
  generic_gConstantProducer("ZH_Sig_POWHEG", "L1Zgs");
  generic_gConstantProducer("VBF", "g2");
  generic_gConstantProducer("VBF", "g4");
  generic_gConstantProducer("VBF", "L1");
  generic_gConstantProducer("VBF", "L1Zgs");
  generic_gConstantProducer("HZZ2e2mu", "g2");
  generic_gConstantProducer("HZZ2e2mu", "g4");
  generic_gConstantProducer("HZZ2e2mu", "L1");
  generic_gConstantProducer("HZZ2e2mu", "L1Zgs");
}

