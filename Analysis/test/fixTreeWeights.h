#ifndef FIXTREEWEIGHTS_H
#define FIXTREEWEIGHTS_H

#include <cassert>
#include "common_includes.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"


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

  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;
  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Requested binning in " << binning.getLabel() << " with " << binning.getNbins() << " bins: [ " << binning.getBinningVector() << " ]"
    << endl;

  const bool allowOverUnderflows=(!binning.getLabel().Contains("m4l"));
  if (!allowOverUnderflows){
    MELAout
      << "fixTreeWeights(" << treename << "): "
      << "Will not allow under or overflow events."
      << endl;
    trimEdges=2;
  }

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

// Assumes trackvar and weight are booked
// trimEdges==0 keeps underflow and overflow bins with no range restriction
// trimEdges==1 keeps underflow and overflow bins but restricts their range
// trimEdges==2 discards underflow and overflow bins
void wipeLargeWeights(std::vector<TTree*>& treeList, const ExtendedBinning& binning, float& trackvar, float& weight, int trimEdges){
  std::vector<TTree*> res;
  if (treeList.empty()) return;
  unsigned int ntrees=treeList.size();
  res.reserve(ntrees);

  //const TString treename=tree->GetName();
  MELAout << "Begin wipeLargeWeights over " << ntrees << " trees." << endl;
  MELAout
    << "wipeLargeWeights: "
    << "Requested binning in " << binning.getLabel() << " with " << binning.getNbins() << " bins: [ " << binning.getBinningVector() << " ]"
    << endl;

  const bool allowOverUnderflows=(!binning.getLabel().Contains("m4l"));
  if (!allowOverUnderflows){
    MELAout
      << "wipeLargeWeights: "
      << "Will not allow under or overflow events."
      << endl;
    trimEdges=2;
  }

  // Store counts and sum of weights
  int const nbins=binning.getNbins()+2;
  vector<float> sumWList(treeList.size(), 0);
  vector<vector<unsigned int>> countsList;
  countsList.assign(treeList.size(), vector<unsigned int>(nbins, 0));

  // Determine if events can be wiped out simultaneously
  bool canWipeSimultaneous=true;
  for (unsigned int it=0; it<treeList.size(); it++){
    TTree* const& tree=treeList.at(it);
    vector<unsigned int>& counts=countsList.at(it);
    float& sumW=sumWList.at(it);

    // Initial loop over the tree to count the events in each bin
    const int nEntries = tree->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      int bin = 1+binning.getBin(trackvar);
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
      sumW+=weight;
    }
    MELAout
      << "wipeLargeWeights(" << tree->GetName() << "): "
      << "Counts: " << counts
      << endl;
    if (nEntries!=treeList.front()->GetEntries()) MELAout
      << "wipeLargeWeights(" << tree->GetName() << "): "
      << "Number of events do not match the front of the list."
      << endl;
    for (int bin=0; bin<nbins; bin++){
      canWipeSimultaneous &= (nEntries==treeList.front()->GetEntries() && counts.at(bin)==countsList.front().at(bin));
      if (counts.at(bin)!=countsList.front().at(bin)) MELAout
        << "wipeLargeWeights(" << tree->GetName() << "): "
        << "Number of counts in bin " << bin << " does not match: " << counts.at(bin) << " != " << countsList.front().at(bin)
        << endl;
    }
  }
  if (canWipeSimultaneous){
    MELAout
      << "wipeLargeWeights: "
      << "Trees satisfy simultaneous event removal."
      << endl;

    vector<int> vetoedEvents;
    for (unsigned int it=0; it<treeList.size(); it++){
      TTree* const& tree=treeList.at(it);
      vector<unsigned int>& counts=countsList.at(it);
      const int nEntries = tree->GetEntries();

      unsigned int nMarginalMax;
      if (nEntries>100000) nMarginalMax=100;
      else nMarginalMax=50;
      unsigned int nMarginalMaxMult;
      if (nEntries>100000) nMarginalMaxMult=1000;
      else if (nEntries>50000) nMarginalMaxMult=500;
      else nMarginalMaxMult=100;
      const float nMarginalMaxFrac = 1./static_cast<float>(nMarginalMaxMult);

      // Collect the count*nMarginalMaxFrac events with highest weights
      MELAout
        << "wipeLargeWeights(" << tree->GetName() << "): "
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
        << "wipeLargeWeights(" << tree->GetName() << "): "
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
            << "wipeLargeWeights(" << tree->GetName() << "): "
            << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
            << endl;
        }
        else if (ns>0) threshold=wgtcoll.front();
        wgtThresholds.push_back(threshold);
      }

      // Determine the vetoed events list
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
        if (fabs(weight)>wgtThresholds.at(bin)) addByLowest(vetoedEvents, ev, true);
      }
    }
    MELAout
      << "wipeLargeWeights: "
      << "All of these events will be vetoed: " << vetoedEvents
      << endl;
    for (unsigned int it=0; it<treeList.size(); it++){
      TTree* const& tree=treeList.at(it);
      vector<unsigned int>& counts=countsList.at(it);
      const int nEntries = tree->GetEntries();

      TTree* newtree = tree->CloneTree(0);
      float const& sumWold=sumWList.at(it);
      float sumWnew=0;
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
        if (TUtilHelpers::checkElementExists(ev, vetoedEvents)) continue;
        newtree->Fill();
        sumWnew += weight;
      }
      MELAout
        << "wipeLargeWeights: "
        << "Tree " << newtree->GetName() << " sum of weights changed from " << sumWold << " to " << sumWnew
        << endl;
      res.push_back(newtree);
    }
  }
  else{
    MELAout
      << "wipeLargeWeights: "
      << "Trees do not satisfy simultaneous event removal."
      << endl;

    for (unsigned int it=0; it<treeList.size(); it++){
      TTree* const& tree=treeList.at(it);
      vector<unsigned int>& counts=countsList.at(it);
      float const& sumWold=sumWList.at(it);

      const int nEntries = tree->GetEntries();
      TTree* newtree = tree->CloneTree(0);
      float sumWnew=0;

      unsigned int nMarginalMax;
      if (nEntries>100000) nMarginalMax=100;
      else nMarginalMax=50;
      unsigned int nMarginalMaxMult;
      if (nEntries>100000) nMarginalMaxMult=1000;
      else if (nEntries>50000) nMarginalMaxMult=500;
      else nMarginalMaxMult=100;
      const float nMarginalMaxFrac = 1./static_cast<float>(nMarginalMaxMult);

      // Collect the count*nMarginalMaxFrac events with highest weights
      MELAout
        << "wipeLargeWeights(" << tree->GetName() << "): "
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
        << "wipeLargeWeights(" << tree->GetName() << "): "
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
            << "wipeLargeWeights(" << tree->GetName() << "): "
            << "Threshold " << threshold << " is different from max. weight " << wgtcoll.front()
            << endl;
        }
        else if (ns>0) threshold=wgtcoll.front();
        wgtThresholds.push_back(threshold);
      }

      // Determine the vetoed events list
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
        if (fabs(weight)>wgtThresholds.at(bin)) continue;
        newtree->Fill();
        sumWnew += weight;
      }
      MELAout
        << "wipeLargeWeights: "
        << "Tree " << newtree->GetName() << " sum of weights changed from " << sumWold << " to " << sumWnew
        << endl;
      res.push_back(newtree);
    }
  }

  swap(res, treeList);
}


// Fix tree weights with the histogram given
TTree* fixTreeWeights(
  TH2F* hRatio,
  TTree* tree,
  float& xvar, float& yvar, float& weight, bool& flag
){
  if (!tree || !hRatio) return nullptr;

  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;

  const int nEntries = tree->GetEntries();
  TTree* newtree = tree->CloneTree(0);

  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Fixing tree weights with the histogram " << hRatio->GetName() << "."
    << endl;

  flag=true;
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    if (!flag) continue;

    int ix=hRatio->GetXaxis()->FindBin(xvar);
    int iy=hRatio->GetYaxis()->FindBin(yvar);
    weight *= hRatio->GetBinContent(ix, iy);

    if (weight>0.) newtree->Fill();
  }
  return newtree;
}
TTree* fixTreeWeights(
  TH3F* hRatio,
  TTree* tree,
  float& xvar, float& yvar, float& zvar, float& weight, bool& flag
){
  if (!tree || !hRatio) return nullptr;

  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;

  const int nEntries = tree->GetEntries();
  TTree* newtree = tree->CloneTree(0);

  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Fixing tree weights with the histogram " << hRatio->GetName() << "."
    << endl;

  flag=true;
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    if (!flag) continue;

    int ix=hRatio->GetXaxis()->FindBin(xvar);
    int iy=hRatio->GetYaxis()->FindBin(yvar);
    int iz=hRatio->GetZaxis()->FindBin(zvar);
    weight *= hRatio->GetBinContent(ix, iy, iz);

    if (weight>0.) newtree->Fill();
  }
  return newtree;
}

// Fix tree weights with the resolution
TTree* fixTreeWeights(
  CategorizationHelpers::Category const category,
  SampleHelpers::Channel const channel,
  ProcessHandler::ProcessType const proctype,
  SystematicsHelpers::SystematicVariationTypes const syst,
  TTree* tree,
  RooWorkspace* w,
  const ExtendedBinning& binning_mass,
  float& ZZMass, float& GenHMass, float& weight
){
  if (
    !tree || !w
    || !(
      syst==SystematicsHelpers::eLepScaleEleDn || syst==SystematicsHelpers::eLepScaleEleUp
      ||
      syst==SystematicsHelpers::eLepScaleMuDn || syst==SystematicsHelpers::eLepScaleMuUp
      ||
      syst==SystematicsHelpers::eLepResEleDn || syst==SystematicsHelpers::eLepResEleUp
      ||
      syst==SystematicsHelpers::eLepResMuDn || syst==SystematicsHelpers::eLepResMuUp
      )
    ) return nullptr;

  const TString treename=tree->GetName();
  MELAout << "Begin fixTreeWeights(" << treename << ")" << endl;

  TString strSystVarName = SystematicsHelpers::getSystematicsCombineName(category, channel, proctype, syst);
  HelperFunctions::replaceString(strSystVarName, "Down", "");
  HelperFunctions::replaceString(strSystVarName, "Dn", ""); // Just in case naming convention changes
  HelperFunctions::replaceString(strSystVarName, "Up", "");
  bool const doUp = (
    syst==SystematicsHelpers::eLepScaleEleUp
    ||
    syst==SystematicsHelpers::eLepScaleMuUp
    ||
    syst==SystematicsHelpers::eLepResEleUp
    ||
    syst==SystematicsHelpers::eLepResMuUp
    );

  RooAbsPdf* pdf = w->pdf("ResolutionModel");
  RooRealVar* ZZMassVar = w->var("mass");
  RooRealVar* GenHMassVar = w->var("MH"); // FIXME: Name to be revised later
  RooRealVar* systVar = w->var(strSystVarName);
  if (!pdf) MELAout
    << "fixTreeWeights(" << treename << "): "
    << "ResolutionModel could not be found!"
    << endl;
  if (!ZZMassVar) MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Reco. mass variable could not be found!"
    << endl;
  if (!GenHMassVar) MELAout
    << "fixTreeWeights(" << treename << "): "
    << "True mass variable could not be found!"
    << endl;
  if (!systVar) MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Systematic variation variable " << strSystVarName << " could not be found!"
    << endl;

  int const nbins=binning_mass.getNbins();
  double const thrlow=2.*binning_mass.getBinLowEdge(0)-binning_mass.getBinLowEdge(1);
  double const thrhigh=2.*binning_mass.getBinLowEdge(nbins)-binning_mass.getBinLowEdge(nbins-1);
  double const curthrlow=ZZMassVar->getMin();
  double const curthrhigh=ZZMassVar->getMax();
  ZZMassVar->setRange(thrlow, thrhigh); // Re-adjust the range of reco. mass

  bool const GenHMassWasConstant = GenHMassVar->getAttribute("Constant");
  GenHMassVar->setConstant(false);
  double const curthrlow_true=GenHMassVar->getMin();
  double const curthrhigh_true=GenHMassVar->getMax();
  GenHMassVar->setRange(0, theSqrts*1000.);

  const int nEntries = tree->GetEntries();
  TTree* newtree = tree->CloneTree(0);

  MELAout
    << "fixTreeWeights(" << treename << "): "
    << "Fixing tree weights with the pdf " << pdf->GetName() << "."
    << endl;

  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    if (ZZMass<thrlow || ZZMass>thrhigh) continue;

    ZZMassVar->setVal(ZZMass);
    GenHMassVar->setVal(GenHMass);
    systVar->setVal(0);
    double vnom = pdf->getVal(*ZZMassVar);
    systVar->setVal((doUp ? 1. : -1.));
    double vvar = pdf->getVal(*ZZMassVar);
    systVar->setVal(0);

    //if (ev%100==0) cout << "Event " << ev << " ZZMass / GenHMass = " << ZZMass << " / " << GenHMass << " | addiitonal wgt = " << vvar/vnom << endl;

    weight *= vvar/vnom;

    if (weight>0.) newtree->Fill();
  }

  ZZMassVar->setRange(curthrlow, curthrhigh);
  GenHMassVar->setRange(curthrlow_true, curthrhigh_true);
  GenHMassVar->setConstant(GenHMassWasConstant);
  return newtree;
}

#endif
