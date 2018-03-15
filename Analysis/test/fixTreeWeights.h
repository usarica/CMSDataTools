#ifndef FIXTREEWEIGHTS_H
#define FIXTREEWEIGHTS_H

#include "common_includes.h"


TTree* fixTreeWeights(TTree* tree){
  if (!tree) return nullptr;
  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;

  float ZZMass, weight;
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("weight", &weight);
  const int nEntries = tree->GetEntries();
  TTree* newtree = tree->CloneTree(0);

  unsigned int nMarginalMax;
  if (nEntries>100000) nMarginalMax=100;
  else nMarginalMax=50;
  unsigned int nMarginalMaxMult;
  if (nEntries>100000) nMarginalMaxMult=1000;
  else if (nEntries>50000) nMarginalMaxMult=500;
  else nMarginalMaxMult=100;
  const float nMarginalMaxFrac = 1./static_cast<float>(nMarginalMaxMult);
  const unsigned int countThreshold=nMarginalMaxMult*nMarginalMax;

  const float massLow = 70.;
  const float massHigh = theSqrts*1000.;
  const float massWidth = 5.;
  int const nbinsraw = (massHigh-massLow)/massWidth+0.5;
  TH1F* hmass = new TH1F("hmass", "", nbinsraw, massLow, massHigh);
  // Initial loop over the tree to count the events in each bin
  for (int ev=0; ev<nEntries; ev++){
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
  if (counts.at(0)<countThreshold && counts.size()>1){
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
  for (int ev=0; ev<nEntries; ev++){
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
    if (threshold*5.>wgtcoll.front()) threshold=wgtcoll.front();
    else MELAout
      << "fixTreeWeights(" << treename << "): "
      << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
      << endl;
    wgtThresholds.push_back(threshold);
  }
  
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    int bin = binning.getBin(ZZMass);
    if (bin>=0 && bin<(int)binning.getNbins() && fabs(weight)>wgtThresholds.at(bin)) weight = pow(wgtThresholds.at(bin), 2)/weight;
    newtree->Fill();
  }
  return newtree;
}

// Assumes trackvar and weight are booked
// trimEdges==0 keeps underflow and overflow bins with no range restriction
// trimEdges==1 keeps underflow and overflow bins but restricts their range
// trimEdges==2 discards underflow and overflow bins
TTree* fixTreeWeights(TTree* tree, const ExtendedBinning& binning, float& trackvar, float& weight, int trimEdges){
  if (!tree) return nullptr;

  if (!tree) return nullptr;
  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Requested binning in " << binning.getNbins() << " bins: [ " << binning.getBinningVector() << " ]"
    << endl;

  const int nEntries = tree->GetEntries();
  TTree* newtree = tree->CloneTree(0);

  unsigned int nMarginalMax;
  if (nEntries>100000) nMarginalMax=100;
  else nMarginalMax=50;
  unsigned int nMarginalMaxMult;
  if (nEntries>100000) nMarginalMaxMult=1000;
  else if (nEntries>50000) nMarginalMaxMult=500;
  else nMarginalMaxMult=100;
  const float nMarginalMaxFrac = 1./static_cast<float>(nMarginalMaxMult);
  const unsigned int countThreshold=nMarginalMaxMult*nMarginalMax;

  int nbins=binning.getNbins()+2;
  vector<unsigned int> counts; counts.assign(nbins, 0);
  // Initial loop over the tree to count the events in each bin
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    int bin = 1+binning.getBin(trackvar);
    // Trim the tree
    bool doPass=false;
    switch (trimEdges){
    case 1:
    {
      double thrlow=2.*binning.getBinLowEdge(0)-binning.getBinLowEdge(1);
      double thrhigh=2.*binning.getBinLowEdge(nbins-2)-binning.getBinLowEdge(nbins-3);
      doPass=(trackvar<thrlow || trackvar>thrhigh);
      break;
    }
    case 2:
      doPass=(bin==0 || bin==nbins-1);
      break;
    }
    if (doPass) continue;
    counts.at(bin)=counts.at(bin)+1;
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Counts: " << counts
    << endl;

  // Collect the count*nMarginalMaxFrac events with highest weights
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Collecting the count*" << nMarginalMaxFrac << " events with highest weights in " << nbins << " bins"
    << endl;
  vector<vector<float>> wgtcollList;
  wgtcollList.assign(counts.size(), vector<float>());
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    int bin = 1+binning.getBin(trackvar);
    // Trim the tree
    bool doPass=false;
    switch (trimEdges){
    case 1:
    {
      double thrlow=2.*binning.getBinLowEdge(0)-binning.getBinLowEdge(1);
      double thrhigh=2.*binning.getBinLowEdge(nbins-2)-binning.getBinLowEdge(nbins-3);
      doPass=(trackvar<thrlow || trackvar>thrhigh);
      break;
    }
    case 2:
      doPass=(bin==0 || bin==nbins-1);
      break;
    }
    if (doPass) continue;
    vector<float>& wgtcoll=wgtcollList.at(bin);
    const unsigned int maxPrunedSize = std::ceil(float(counts.at(bin))*nMarginalMaxFrac);
    if (wgtcoll.size()<maxPrunedSize) addByHighest(wgtcoll, fabs(weight), false);
    else if (wgtcoll.back()<fabs(weight)){
      addByHighest(wgtcoll, fabs(weight), false);
      wgtcoll.pop_back();
    }
  }
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Determining the weight thresholds"
    << endl;
  vector<float> wgtThresholds; wgtThresholds.reserve(wgtcollList.size());
  for (auto const& wgtcoll:wgtcollList){
    unsigned int ns=wgtcoll.size();
    float threshold=0;
    if (ns>=3){
      threshold=0.5*(wgtcoll.at(ns-1)+wgtcoll.at(ns-2));
      if (threshold*5.>wgtcoll.front()) threshold=wgtcoll.front();
      else MELAout
        << "fixTreeWeights(" << treename << "): "
        << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
        << endl;
    }
    else if (ns>0) threshold=wgtcoll.front();
    wgtThresholds.push_back(threshold);
  }

  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    int bin = 1+binning.getBin(trackvar);
    // Trim the tree
    bool doPass=false;
    switch (trimEdges){
    case 1:
    {
      double thrlow=2.*binning.getBinLowEdge(0)-binning.getBinLowEdge(1);
      double thrhigh=2.*binning.getBinLowEdge(nbins-2)-binning.getBinLowEdge(nbins-3);
      doPass=(trackvar<thrlow || trackvar>thrhigh);
      break;
    }
    case 2:
      doPass=(bin==0 || bin==nbins-1);
      break;
    }
    if (doPass) continue;
    if (fabs(weight)>wgtThresholds.at(bin)) weight = pow(wgtThresholds.at(bin), 2)/weight;
    newtree->Fill();
  }
  return newtree;
}


#endif
