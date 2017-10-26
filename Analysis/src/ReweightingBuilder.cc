#include "ReweightingBuilder.h"
#include "SimpleEntry.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<TString>&)) :
rule(infcn)
{ strWeights.push_back(inStrWeight); }

ReweightingBuilder::ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<TString>&)) :
rule(infcn),
strWeights(inStrWeights)
{}

float ReweightingBuilder::eval(CJLSTTree* theTree) const{ if (rule) return rule(theTree, strWeights); else return 0; }

void ReweightingBuilder::setupWeightVariables(CJLSTTree* theTree, float fractionRequirement){
  MELAout << "ReweightingBuilder[" << strWeights << "]::setupWeightVariables is called for tree " << theTree->sampleIdentifier << "." << endl;

  std::vector<float> res;

  if (!theTree) return;
  const ExtendedBinning& binning = this->weightBinning;
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
  const int nevents = theTree->getNEvents();
  MELAout << "\t- Ordering the " << nevents << " events";
  if (!noBoundaries) MELAout << " based on the " << binning.getNbins() << " bins: [ " << binning.getBinningVector() << " ]";
  MELAout << '.' << endl;
  while (theTree->getEvent(ev)){
    HelperFunctions::progressbar(ev, nevents);

    float wgt = this->eval(theTree);
    float orderingVal=-1;
    theTree->getVal(strOrderingVal, orderingVal);
    SimpleEntry theEntry(0, fabs(wgt), wgt);

    if (noBoundaries) addByLowest(indexList.front(), theEntry, false);
    else{
      int bin = binning.getBin(orderingVal);
      if (bin>=0 && bin<(int)ns) addByLowest(indexList.at(bin), theEntry, false);
    }

    ev++;
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
std::vector<float> ReweightingBuilder::getWeightThresholds(CJLSTTree* theTree) const{
  if (!theTree) return vector<float>();
  auto it = weightThresholds.find(theTree);

  if (it!=weightThresholds.cend()) return it->second;
  else return vector<float>();
}

void ReweightingBuilder::setWeightBinning(const ExtendedBinning& binning){ weightBinning=binning; }
int ReweightingBuilder::findBin(CJLSTTree* theTree) const{
  const ExtendedBinning& binning = weightBinning;
  int bin=0;
  if (binning.isValid()){
    float orderingVal; theTree->getVal(binning.getLabel(), orderingVal);
    int bin = binning.getBin(orderingVal);
    if (bin<0 || bin>(int) binning.getNbins()) bin=-1;
  }
  return bin;
}
float ReweightingBuilder::getPostThresholdWeight(CJLSTTree* theTree) const{
  float weight = this->eval(theTree);
  int bin=this->findBin(theTree);
  if (bin>=0){
    const float& threshold=weightThresholds.find(theTree)->second.at(bin);
    if (fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
  }
  return weight;
}
float ReweightingBuilder::getPostThresholdSumWeights(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin>=0) return sumPostThrWeights.find(theTree)->second.at(bin);
  else return 0;
}
unsigned int ReweightingBuilder::getSumEvents(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin>=0) return sumEvents.find(theTree)->second.at(bin);
  else return 0;
}
unsigned int ReweightingBuilder::getSumNonZeroWgtEvents(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin>=0) return sumNonZeroWgtEvents.find(theTree)->second.at(bin);
  else return 0;
}
