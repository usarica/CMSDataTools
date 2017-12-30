#ifndef FIXTREEWEIGHTS_H
#define FIXTREEWEIGHTS_H

#include "common_includes.h"


TTree* fixTreeWeights(TTree* tree){
  const unsigned int nMarginalMax = 100;
  const unsigned int nMarginalMaxMult = 1000;
  const float nMarginalMaxFrac = 1./static_cast<float const>(nMarginalMaxMult);
  const unsigned int countThreshold=nMarginalMaxMult*nMarginalMax;

  const TString treename=tree->GetName();
  float ZZMass, weight;
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("weight", &weight);

  int const nbinsraw = (1000*theSqrts-70)/5;
  TH1F* hmass = new TH1F("hmass", "", nbinsraw, 70, 13000);
  TTree* newtree = tree->CloneTree(0);
  // Initial loop over the tree to count the events in each bin
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    hmass->Fill(ZZMass);
  }

  // Determine the final binning to set the weight thresholds
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Determining the final binning to set the weight thresholds"
    << endl;
  ExtendedBinning binning;
  binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(hmass->GetNbinsX()+1));
  vector<unsigned int> counts;
  unsigned int count=0;
  for (int bin=hmass->GetNbinsX(); bin>=0; bin--){
    count += hmass->GetBinContent(bin);
    if (count>countThreshold || bin==0){
      counts.push_back(count);
      binning.addBinBoundary(hmass->GetXaxis()->GetBinLowEdge(bin));
      count=0;
    }
  }
  delete hmass;
  std::reverse(counts.begin(), counts.end());
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;
  // These lines guarantee count>countThreshold in every bin
  if (counts.at(0)<countThreshold){
    counts.at(1) += counts.at(0);
    counts.erase(counts.begin());
    binning.removeBinLowEdge(1);
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "counts.size()=" << counts.size() << "=?" << "nbins=" << binning.getNbins()
    << endl;

  // Collect the count*nMarginalMaxFrac events with highest weights
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Collecting the count*" << nMarginalMaxFrac << " events with highest weights in " << binning.getNbins() << " bins"
    << endl;
  vector<vector<float>> wgtcollList;
  wgtcollList.assign(binning.getNbins(), vector<float>());
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int)binning.getNbins()){
      vector<float>& wgtcoll=wgtcollList.at(bin);
      const unsigned int maxPrunedSize = std::ceil(float(counts.at(bin))*nMarginalMaxFrac);
      if (wgtcoll.size()<maxPrunedSize) addByHighest(wgtcoll, fabs(weight), false);
      else if (wgtcoll.back()<fabs(weight)){
        addByHighest(wgtcoll, fabs(weight), false);
        wgtcoll.pop_back();
      }
    }
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Determining the weight thresholds"
    << endl;
  vector<float> wgtThresholds; wgtThresholds.reserve(binning.getNbins());
  for (auto const& wgtcoll:wgtcollList){
    unsigned int ns=wgtcoll.size();
    float threshold=0.5*(wgtcoll.at(ns-1)+wgtcoll.at(ns-2));
    if (wgtcoll.front()*5.<threshold) threshold=wgtcoll.front();
    else MELAout
      << "fixTreeWeights(" << treename << "): "
      << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
      << endl;
    wgtThresholds.push_back(threshold);
  }
  
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int)binning.getNbins() && fabs(weight)>wgtThresholds.at(bin)) weight = pow(wgtThresholds.at(bin), 2)/weight;
    newtree->Fill();
  }
  return newtree;
}


#endif
