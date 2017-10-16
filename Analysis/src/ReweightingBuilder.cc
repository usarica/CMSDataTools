#include "ReweightingBuilder.h"
#include "SimpleEntry.h"

using namespace std;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<TString>&)) :
rule(infcn)
{ strWeights.push_back(inStrWeight); }

ReweightingBuilder::ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<TString>&)) :
rule(infcn),
strWeights(inStrWeights)
{}

float ReweightingBuilder::eval(CJLSTTree* theTree) const{ if (rule) return rule(theTree, strWeights); else return 0; }

void ReweightingBuilder::setupWeightVariables(CJLSTTree* theTree, const ExtendedBinning& binning, float fractionRequirement){
  std::vector<float> res;

  if (!theTree) return;
  weightBinning[theTree]=binning;
  TString strOrderingVal = binning.getLabel();

  const bool noBoundaries = !binning.isValid();
  const unsigned int ns = (!noBoundaries ? binning.getNbins() : 1);

  // Get sum of weights
  sumPostThrWeights[theTree]=std::vector<float>();
  sumEvents[theTree]=std::vector<unsigned int>();
  sumNonZeroWgtEvents[theTree]=std::vector<unsigned int>();
  vector<vector<SimpleEntry>> indexList;
  sumPostThrWeights[theTree].assign(ns, 0);
  sumEvents[theTree].assign(ns, 0);
  sumNonZeroWgtEvents[theTree].assign(ns, 0);
  indexList.assign(ns, vector<SimpleEntry>());

  int ev=0;
  while (theTree->getEvent(ev)){
    ev++;

    float wgt = this->eval(theTree);
    float orderingVal=-1;
    theTree->getVal(strOrderingVal, orderingVal);
    SimpleEntry theEntry(0, fabs(wgt), wgt);

    if (noBoundaries) addByLowest(indexList.front(), theEntry, false);
    else{
      int bin = binning.getBin(orderingVal);
      if (bin>=0 && bin<(int)ns) addByLowest(indexList.at(bin), theEntry, false);
    }
  }

  for (unsigned int ibin=0; ibin<ns; ibin++){
    const vector<SimpleEntry>& index = indexList.at(ibin);

    // First count the events
    sumEvents[theTree].at(ibin) = index.size();
    for (auto const& sentry:index){ if (sentry.trackingval!=0.) sumNonZeroWgtEvents[theTree].at(ibin)++; }

    // Find the threshold
    unsigned int index_entry = static_cast<int>((float(index.size())*fractionRequirement));
    if (index_entry==index.size() && index_entry>0) index_entry--;
    unsigned int index_entry_prev=index_entry; if (index_entry_prev>0) index_entry_prev--;
    float threshold = 0;
    if (!index.empty()) threshold = (index.at(index_entry_prev).trackingval + index.at(index_entry).trackingval)*0.5;
    if (threshold>0. && index.back().trackingval<threshold*5.) threshold = index.back().trackingval; // Prevent false-positives
    res.push_back(threshold);

    for (auto const& sentry:index){
      float weight = sentry.weight;
      if (fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
      sumPostThrWeights[theTree].at(ibin) += weight;
    }
  }

  weightThresholds[theTree] = res;
}
std::vector<float> ReweightingBuilder::getWeightThresholds(CJLSTTree* theTree){
  if (theTree && weightThresholds.find(theTree)!=weightThresholds.end()) return weightThresholds[theTree];
  else return vector<float>();
}

float ReweightingBuilder::getPostThresholdWeight(CJLSTTree* theTree) const{
  float weight = this->eval(theTree);
  auto binningIt = weightBinning.find(theTree);
  if (binningIt!=weightBinning.cend()){
    const ExtendedBinning& binning = binningIt->second;
    int bin=0;
    if (binning.isValid()){
      float orderingVal; theTree->getVal(binning.getLabel(), orderingVal);
      int bin = binning.getBin(orderingVal);
      if (bin<0 || bin>(int) binning.getNbins()) bin=-1;
    }
    if (bin>=0){
      const float& threshold=weightThresholds.find(theTree)->second.at(bin);
      if (fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
    }
  }
  return weight;
}

