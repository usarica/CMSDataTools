#include "ReweightingBuilder.h"
#include "ReweightingFunctions.h"
#include "SimpleEntry.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(TString inStrWeight, float(*infcn)(CJLSTTree*, const std::vector<float*>&)) :
  allowNegativeWeights(true),
  divideByNSample(false),
  rule(infcn)
{
  strWeights.push_back(inStrWeight);
}
ReweightingBuilder::ReweightingBuilder(std::vector<TString> inStrWeights, float(*infcn)(CJLSTTree*, const std::vector<float*>&)) :
  allowNegativeWeights(true),
  divideByNSample(false),
  rule(infcn),
  strWeights(inStrWeights)
{}

std::vector<TString> const& ReweightingBuilder::getWeightVariables() const{ return strWeights; }

float ReweightingBuilder::eval(CJLSTTree* theTree) const{
  std::unordered_map<CJLSTTree*, std::vector<float*>>::const_iterator it = componentRefs.find(theTree);
  if (it==componentRefs.cend()){
    MELAerr << "ReweightingBuilder::eval: Could not find the weights " << strWeights << ". Call ReweightingBuilder::setupWeightVariables first!" << endl;
    return 0;
  }
  float divisor=1; if (divideByNSample) divisor=theTree->getNGenNoPU();
  float weight=0;
  if (rule && divisor!=0.) weight=rule(theTree, it->second)/divisor;
  if (!allowNegativeWeights && weight<0.) weight=0;
  return weight;
}

int ReweightingBuilder::findBin(CJLSTTree* theTree) const{
  const ExtendedBinning& binning = weightBinning;
  int bin=0;
  if (binning.isValid()){
    auto binningVarRefsIt = binningVarRefs.find(theTree);
    if (binningVarRefsIt!=binningVarRefs.cend()){
      float const& orderingVal=*(binningVarRefsIt->second);
      bin = binning.getBin(orderingVal);
      if (bin<0 || bin>=(int) binning.getNbins()) bin=-1;
    }
  }
  return bin;
}
void ReweightingBuilder::rejectNegativeWeights(const bool flag){ allowNegativeWeights = !flag; }
void ReweightingBuilder::setDivideByNSample(const bool flag){ divideByNSample = flag; }
void ReweightingBuilder::setWeightBinning(const ExtendedBinning& binning){ weightBinning=binning; }

void ReweightingBuilder::setupWeightVariables(CJLSTTree* theTree, float fractionRequirement, unsigned int minimumNevents){
  MELAout << "ReweightingBuilder[" << strWeights << "]::setupWeightVariables is called for tree " << theTree->sampleIdentifier << "." << endl;
  if (divideByNSample) MELAout << "ReweightingBuilder[" << strWeights << "]::setupWeightVariables will divide the weights by "
    << theTree->getNGenNoPU() << "." << endl;

  std::vector<float> res;

  if (!theTree) return;

  const ExtendedBinning& binning = this->weightBinning;
  TString strOrderingVal = binning.getLabel();
  if (strOrderingVal=="") return;
  float* orderingValRef = ReweightingFunctions::getWeightRef(theTree, strOrderingVal);
  if (!orderingValRef) return;

  // First link components
  componentRefs[theTree] = ReweightingFunctions::getWeightRefs(theTree, strWeights);
  binningVarRefs[theTree] = orderingValRef;

  // Now compute other weight-related variables
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
  for (auto& index:indexList) index.reserve(theTree->getNEvents()/ns+1);

  int ev=0;
  const int nevents = theTree->getNEvents();
  MELAout << "\t- Ordering the " << nevents << " events";
  if (!noBoundaries) MELAout << " over the " << ns << " bins: [ " << binning.getBinningVector() << " ]";
  MELAout << '.' << endl;
  while (theTree->getEvent(ev)){
    HelperFunctions::progressbar(ev, nevents);

    float wgt = this->eval(theTree);
    float const& orderingVal=*orderingValRef;
    SimpleEntry theEntry(0, fabs(wgt), wgt);

    int bin=0;
    if (!noBoundaries) bin = binning.getBin(orderingVal);
    if (bin>=0 && bin<(int) ns){
      // Accumulate the events
      indexList.at(bin).push_back(theEntry);
      sumEvents[theTree].at(bin)++;
      if (theEntry.trackingval!=0.) sumNonZeroWgtEvents[theTree].at(bin)++;
    }
    ev++;
  }

  for (unsigned int ibin=0; ibin<ns; ibin++){
    vector<SimpleEntry>& index = indexList.at(ibin);
    if (minimumNevents>sumNonZeroWgtEvents[theTree].at(ibin)){
      MELAout << "\t- Bin " << ibin << " has less number of events with non-zero weight than the requested number " << minimumNevents << ". Resetting the bin..." << endl;
      index.clear();
      sumNonZeroWgtEvents[theTree].at(ibin)=0;
      //sumEvents[theTree].at(ibin)=0;
    }
    const unsigned int nTotalPerBin = index.size();
    //MELAout << "\t- Looping over the " << nTotalPerBin << " events to find the threshold in bin " << ibin << endl;
    float threshold=0;
    if (nTotalPerBin>2){
      const unsigned int maxPrunedSize = std::max((unsigned int) (fractionRequirement>=0. ? std::ceil(float(nTotalPerBin)*(1.-fractionRequirement)) : nTotalPerBin), (unsigned int) 3);
      vector<SimpleEntry> indexPruned; indexPruned.reserve(maxPrunedSize);

      unsigned int itrk=0;
      for (auto& theEntry:index){
        //HelperFunctions::progressbar(itrk, nTotalPerBin);
        if (indexPruned.size()<maxPrunedSize) addByHighest(indexPruned, theEntry, false);
        else if (indexPruned.at(indexPruned.size()-1).trackingval<theEntry.trackingval){
          addByHighest(indexPruned, theEntry, false);
          indexPruned.pop_back();
        }
        itrk++;
      }

      // Find the threshold
      unsigned int index_entry = maxPrunedSize-2;
      unsigned int index_entry_prev=index_entry+1;
      threshold = (indexPruned.at(index_entry_prev).trackingval + indexPruned.at(index_entry).trackingval)*0.5;
      cout << "Threshold raw: " << threshold << endl;
      if (threshold>0. && indexPruned.front().trackingval<threshold*5.) threshold = indexPruned.front().trackingval; // Prevent false-positives
    }
    else if (nTotalPerBin==2) threshold = std::max(index.at(0).trackingval, index.at(1).trackingval);
    else if (nTotalPerBin==1) threshold = index.at(0).trackingval;
    res.push_back(threshold);

    //MELAout << "\t- Looping over the " << nTotalPerBin << " events to find the sum of weights after threshold " << threshold << " in bin " << ibin << endl;
    // Do a precise summation with the Kahan method
    KahanAccumulator<float> sum;
    for (auto& theEntry:index){
      float weight = theEntry.weight;
      if (fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
      sum += weight;
    }
    // Assign the sum
    sumPostThrWeights[theTree].at(ibin)=sum;
    MELAout << "\t- Threshold at bin " << ibin << ": " << threshold
      << ", sum of post-threshold weights: " << sumPostThrWeights[theTree].at(ibin)
      << ", Nevents: " << sumNonZeroWgtEvents[theTree].at(ibin) << " / " << sumEvents[theTree].at(ibin)
      << endl;
  }
  weightThresholds[theTree] = res;
}

std::vector<CJLSTTree*> ReweightingBuilder::getRegisteredTrees() const{
  std::vector<CJLSTTree*> res;
  for (auto it=componentRefs.cbegin(); it!=componentRefs.cend(); it++) res.push_back(it->first);
  return res;
}


std::vector<float> ReweightingBuilder::getWeightThresholds(CJLSTTree* theTree) const{
  if (!theTree) return vector<float>();
  auto it = weightThresholds.find(theTree);

  if (it!=weightThresholds.cend()) return it->second;
  else return vector<float>();
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
float ReweightingBuilder::getSumPostThresholdWeights(CJLSTTree* theTree) const{
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


float ReweightingBuilder::getSumAllPostThresholdWeights(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllPostThresholdWeights(bin);
}
float ReweightingBuilder::getSumAllPostThresholdWeights(int bin) const{
  if (bin<0) return 0;
  KahanAccumulator<float> sum;
  auto const& theMap=this->sumPostThrWeights;
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++) sum += it->second.at(bin);
  return sum;
}
unsigned int ReweightingBuilder::getSumAllEvents(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllEvents(bin);
}
unsigned int ReweightingBuilder::getSumAllEvents(int bin) const{
  if (bin<0) return 0;
  unsigned int sum(0);
  auto const& theMap=this->sumEvents;
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++) sum += it->second.at(bin);
  return sum;
}
unsigned int ReweightingBuilder::getSumAllNonZeroWgtEvents(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllNonZeroWgtEvents(bin);
}
unsigned int ReweightingBuilder::getSumAllNonZeroWgtEvents(int bin) const{
  if (bin<0) return 0;
  unsigned int sum(0);
  auto const& theMap=this->sumNonZeroWgtEvents;
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++) sum += it->second.at(bin);
  return sum;
}


float ReweightingBuilder::getNormComponent(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->ReweightingBuilder::getNormComponent(bin);
}
float ReweightingBuilder::getNormComponent(int bin) const{
  if (bin<0) return 0;

  vector<CJLSTTree*> const trees = this->getRegisteredTrees();
  KahanAccumulator<float> numerator;
  KahanAccumulator<float> denominator;

  // Numerator: sum{t} (sum_wgts_tj * N_tj / N_t)
  // Denominator: sum{t} (N_tj)
  for (auto& tree:trees){
    auto itSumWeights = this->sumPostThrWeights.find(tree); if (itSumWeights==this->sumPostThrWeights.cend()) continue;
    auto itNNonZeroWgtEvents = this->sumNonZeroWgtEvents.find(tree); if (itNNonZeroWgtEvents==this->sumNonZeroWgtEvents.cend()) continue;
    //auto itNEvents = this->sumEvents.find(tree); if (itNEvents==this->sumEvents.cend()) continue;

    auto const& sumwgts = itSumWeights->second.at(bin);
    auto const& nevtsnonzerowgt = itNNonZeroWgtEvents->second.at(bin); // N_tj
    //auto const& nevts = itNEvents->second.at(bin);

    // N_t = sum{j} (N_tj)
    unsigned int nSumNonZeroWgt = 0;
    for (auto const& nn:itNNonZeroWgtEvents->second) nSumNonZeroWgt += nn;
    unsigned int nSample = 1; if (divideByNSample) nSample=tree->getNGenNoPU();

    // sum_wgts_tj * N_tj / N_t
    float numerator_pertree = sumwgts * float(nevtsnonzerowgt);
    if (nSumNonZeroWgt!=0){
      if (!divideByNSample) numerator_pertree /= float(nSumNonZeroWgt);
      else numerator_pertree *= float(nSample)/float(nSumNonZeroWgt);
    }

    // Numerator += sum_wgts_tj * N_tj / N_t
    numerator += numerator_pertree;
    // Denominator += N_tj
    denominator += nevtsnonzerowgt;
  }

  float result=0;
  if (denominator!=0.) result = numerator/denominator;
  return result;
}
float ReweightingBuilder::getNorm() const{
  const ExtendedBinning& binning = this->weightBinning;
  const bool noBoundaries = !binning.isValid();
  const int ns = (!noBoundaries ? binning.getNbins() : 1);

  KahanAccumulator<float> sum;
  unsigned int sumN=0;
  for (int bin=0; bin<ns; bin++){
    sum += this->getNormComponent(bin);
    sumN += this->getSumAllNonZeroWgtEvents(bin);
  }
  if (sumN!=0) sum /= float(sumN);
  return sum;
}
