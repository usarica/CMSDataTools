#include "ReweightingBuilder.h"
#include "SimpleEntry.h"
#include "SampleHelpers.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(TString inStrWeight, ReweightingFunctions::ReweightingFunction_t infcn) :
  allowNegativeWeights(true),
  divideByNSample(false),
  rule(infcn)
{
  if (inStrWeight!="") strWeights.push_back(inStrWeight);
}
ReweightingBuilder::ReweightingBuilder(std::vector<TString> inStrWeights, ReweightingFunctions::ReweightingFunction_t infcn) :
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
  if (!allowNegativeWeights && weight<0.){
    weight=0;
    MELAerr << "ReweightingBuilder::eval: Negative weight encountered: "; for (auto& v:it->second) MELAerr << *v << ", "; MELAerr << endl;
    assert(0);
  }
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
  sumPostThrSqWeights[theTree]=std::vector<float>();
  weightThresholds[theTree]=std::vector<float>();
  sumEvents[theTree]=std::vector<unsigned int>();
  sumNonZeroWgtEvents[theTree]=std::vector<unsigned int>();
  sumPostThrWeights[theTree].assign(ns, 0);
  sumPostThrSqWeights[theTree].assign(ns, 0);
  if (fractionRequirement>=0.){
    weightThresholds[theTree].assign(ns, 0);
    sumEvents[theTree].assign(ns, 0);
    sumNonZeroWgtEvents[theTree].assign(ns, 0);
  }
  else{
    weightThresholds[theTree].assign(ns, -1);
    sumEvents[theTree].assign(ns, -1);
    sumNonZeroWgtEvents[theTree].assign(ns, -1);
    return;
  }

  vector<vector<SimpleEntry>> indexList;
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
    if (minimumNevents>sumNonZeroWgtEvents[theTree].at(ibin) && sumEvents[theTree].at(ibin)>0){
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

      //unsigned int itrk=0;
      for (auto& theEntry:index){
        //HelperFunctions::progressbar(itrk, nTotalPerBin);
        if (indexPruned.size()<maxPrunedSize) addByHighest(indexPruned, theEntry, false);
        else if (indexPruned.at(indexPruned.size()-1).trackingval<theEntry.trackingval){
          addByHighest(indexPruned, theEntry, false);
          indexPruned.pop_back();
        }
        //itrk++;
      }

      // Find the threshold
      unsigned int index_entry = maxPrunedSize-2;
      unsigned int index_entry_prev=index_entry+1;
      threshold = (indexPruned.at(index_entry_prev).trackingval + indexPruned.at(index_entry).trackingval)*0.5;
      cout << "Threshold raw: " << threshold << endl;
      if ((threshold>0. && indexPruned.front().trackingval<threshold*5.) || fractionRequirement>=1.) threshold = indexPruned.front().trackingval; // Prevent false-positives
    }
    else if (nTotalPerBin==2) threshold = std::max(index.at(0).trackingval, index.at(1).trackingval);
    else if (nTotalPerBin==1) threshold = index.at(0).trackingval;
    weightThresholds[theTree].at(ibin)=threshold;

    //MELAout << "\t- Looping over the " << nTotalPerBin << " events to find the sum of weights after threshold " << threshold << " in bin " << ibin << endl;
    // Do a precise summation with the Kahan method
    KahanAccumulator<float> sum;
    KahanAccumulator<float> sumsq;
    for (auto& theEntry:index){
      float weight = theEntry.weight;
      if (fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
      sum += weight;
      sumsq += pow(weight, 2);
    }
    // Assign the sum
    sumPostThrWeights[theTree].at(ibin)=sum;
    sumPostThrSqWeights[theTree].at(ibin)=sumsq;
    if (sumEvents[theTree].at(ibin)>0) MELAout << "\t- Threshold at bin " << ibin << ": " << threshold
      << ", sum of post-threshold weights: " << sumPostThrWeights[theTree].at(ibin) << " +- " << sqrt(sumPostThrSqWeights[theTree].at(ibin))
      << ", Nevents: " << sumNonZeroWgtEvents[theTree].at(ibin) << " / " << sumEvents[theTree].at(ibin)
      << endl;
  }
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
    if (threshold>=0. && fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
  }
  return weight;
}
float ReweightingBuilder::getPostThresholdSqWeight(CJLSTTree* theTree) const{
  float weight = this->eval(theTree);
  int bin=this->findBin(theTree);
  if (bin>=0){
    const float& threshold=weightThresholds.find(theTree)->second.at(bin);
    if (threshold>=0. && fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
  }
  return pow(weight, 2);
}
float ReweightingBuilder::getPostThresholdSqWeightInv(CJLSTTree* theTree) const{
  float weight = this->eval(theTree);
  int bin=this->findBin(theTree);
  if (bin>=0){
    const float& threshold=weightThresholds.find(theTree)->second.at(bin);
    if (threshold>=0. && fabs(weight)>threshold) weight = pow(threshold, 2)/weight;
  }
  return (weight!=0. ? pow(weight, -2) : 0.);
}
float ReweightingBuilder::getSumPostThresholdWeights(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin>=0) return sumPostThrWeights.find(theTree)->second.at(bin);
  else return 0;
}
float ReweightingBuilder::getSumPostThresholdSqWeights(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  if (bin>=0) return sumPostThrSqWeights.find(theTree)->second.at(bin);
  else return 0;
}
float ReweightingBuilder::getSumPostThresholdSqWeightInvs(CJLSTTree* theTree) const{
  float res = this->getSumPostThresholdSqWeights(theTree);
  return (res!=0. ? 1./res : 0.);
}
float ReweightingBuilder::getSumPostThresholdNeffs(CJLSTTree* theTree) const{
  float res = this->getSumPostThresholdSqWeightInvs(theTree);
  res *= pow(this->getSumPostThresholdWeights(theTree), 2);
  return res;
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
float ReweightingBuilder::getSumAllPostThresholdSqWeights(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllPostThresholdSqWeights(bin);
}
float ReweightingBuilder::getSumAllPostThresholdSqWeights(int bin) const{
  if (bin<0) return 0;
  KahanAccumulator<float> sum;
  auto const& theMap=this->sumPostThrSqWeights;
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++) sum += it->second.at(bin);
  return sum;
}
float ReweightingBuilder::getSumAllPostThresholdSqWeightInvs(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllPostThresholdSqWeightInvs(bin);
}
float ReweightingBuilder::getSumAllPostThresholdSqWeightInvs(int bin) const{
  if (bin<0) return 0;
  KahanAccumulator<float> sum;
  auto const& theMap=this->sumPostThrSqWeights;
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++){
    float const& component=it->second.at(bin);
    if (component!=0.) sum += 1./component;
  }
  return sum;
}
float ReweightingBuilder::getSumAllPostThresholdNeffs(CJLSTTree* theTree) const{
  int bin=this->findBin(theTree);
  return this->getSumAllPostThresholdNeffs(bin);
}
float ReweightingBuilder::getSumAllPostThresholdNeffs(int bin) const{
  if (bin<0) return 0;
  KahanAccumulator<float> sum;
  auto const& theMap=this->sumPostThrWeights;
  auto const& theMapSq=this->sumPostThrSqWeights;
  auto itSq=theMapSq.cbegin();
  for (auto it=theMap.cbegin(); it!=theMap.cend(); it++){
    float const& component=it->second.at(bin);
    float const& component_sq=itSq->second.at(bin);
    if (component_sq!=0.) sum += pow(component, 2)/component_sq;
    itSq++;
  }
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

  // Numerator: sum{t} ((S_tj/N_t) / (V_tj/N_t^2))
  // Denominator: sum{t} (1/(V_tj/N_t^2))
  for (auto& tree:trees){
    auto itSumWeights = this->sumPostThrWeights.find(tree); if (itSumWeights==this->sumPostThrWeights.cend()) continue;
    auto itSumSqWeights = this->sumPostThrSqWeights.find(tree); if (itSumSqWeights==this->sumPostThrSqWeights.cend()) continue;
    auto itNNonZeroWgtEvents = this->sumNonZeroWgtEvents.find(tree); if (itNNonZeroWgtEvents==this->sumNonZeroWgtEvents.cend()) continue;
    //auto itNEvents = this->sumEvents.find(tree); if (itNEvents==this->sumEvents.cend()) continue;

    float const& sumwgts = itSumWeights->second.at(bin); // S_tj
    float const& sumsqwgts = itSumSqWeights->second.at(bin); // V_tj
    unsigned int const& nevtsnonzerowgt = itNNonZeroWgtEvents->second.at(bin); // N_tj
    //auto const& nevts = itNEvents->second.at(bin);

    // N_t = sum{j} (N_tj)
    unsigned int nSumNonZeroWgt = 0;
    for (unsigned int const& nn:itNNonZeroWgtEvents->second) nSumNonZeroWgt += nn;

    float numerator_pertree = 0;
    float denominator_pertree = 0;
    if (nevtsnonzerowgt>0){
      // S_tj * N_tj / N_t
      //numerator_pertree = sumwgts * static_cast<float>(nevtsnonzerowgt*nSample)/static_cast<float>(nSumNonZeroWgt);

      unsigned int nSample = 1;
      if (divideByNSample) nSample = tree->getNGenNoPU();
      float scalefactor = static_cast<float>(nSample)/static_cast<float>(nSumNonZeroWgt);
      float sampleRawWeight = (ReweightingBuilder::useNeffInNormComponent ? pow(sumwgts, 2)/sumsqwgts : 1./sumsqwgts);
      denominator_pertree = sampleRawWeight/pow(scalefactor, 2);
      numerator_pertree = sumwgts*scalefactor*denominator_pertree;
    }
    /*
    if (divideByNSample && fabs(static_cast<float>(nSample)/static_cast<float>(nSumNonZeroWgt)-1.)>0.01)
    MELAerr
    << "ReweightingBuilder::getNormComponent: WARNING! N events with non-zero weights (" << nSumNonZeroWgt
    << ") is more than 1% different from N of sample (" << nSample
    << ") in reweighting bin " << bin << " of tree " << tree->sampleIdentifier
    << endl;
    */
    /*
    MELAout << "- Tree " << tree->sampleIdentifier << " bin " << bin << " norm: \n"
    << "- sum_wgts_tj: " << sumwgts << '\n'
    << "- sum_wgtsq_tj: " << sumsqwgts << '\n'
    << "- N_tj: " << nevtsnonzerowgt << '\n'
    << "- N_t: " << nSumNonZeroWgt << '\n'
    << "- Num per tree: " << numerator_pertree << '\n'
    << "- Den per tree: " << denominator_pertree
    << endl;
    */

    // Numerator += (S_tj/N_t) / (V_tj/N_t^2)
    numerator += numerator_pertree;
    // Denominator += N_tj
    //denominator += nevtsnonzerowgt;
    // Denominator += 1/(V_tj/N_t^2)
    denominator += denominator_pertree;
  }

  float result=0;
  if (denominator!=0.) result = numerator/denominator;
  //MELAout << "Final norm: " << numerator << " / " << denominator << " = " << result << endl;
  //exit(1);
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
    if (divideByNSample) sumN += this->getSumAllNonZeroWgtEvents(bin);
  }
  if (sumN!=0) sum /= float(sumN);
  return sum;
}

ExtendedBinning ReweightingBuilder::getTrueMassBinning(std::vector<CJLSTTree*> const& trees, float forceMass){
  ExtendedBinning GenHMassBinning("GenHMass");
  GenHMassBinning.addBinBoundary(0);
  GenHMassBinning.addBinBoundary(theSqrts*1000.);

  if (trees.empty()) return GenHMassBinning;
  SampleHelpers::makeGlobalMELA(theSqrts);

  if (forceMass<0.){ // Determine from the masses in the list
    vector<double> masslist;
    for (CJLSTTree* const& tree:trees){ if (tree->MHVal>0.) HelperFunctions::addByLowest<double>(masslist, (double) tree->MHVal, true); }
    if (masslist.size()==1){
      const double MHValfirst = masslist.back();
      GenHMassBinning.addBinBoundary(100);
      GenHMassBinning.addBinBoundary(MHValfirst-1);
      GenHMassBinning.addBinBoundary(MHValfirst+1);
      GenHMassBinning.addBinBoundary(160);
      GenHMassBinning.addBinBoundary(220);
      GenHMassBinning.addBinBoundary(450);
      GenHMassBinning.addBinBoundary(1300);
    }
    else if (!masslist.empty()){
      for (unsigned int is=0; is<masslist.size()-1; is++){
        double const& MHfirst = masslist.at(is);
        double const& MHsecond = masslist.at(is+1);
        double GHfirst = SampleHelpers::GlobalMELA->getHiggsWidthAtPoleMass(MHfirst);
        double GHsecond = SampleHelpers::GlobalMELA->getHiggsWidthAtPoleMass(MHsecond);
        double GHsum = GHfirst + GHsecond;
        double boundary;
        if (GHsum*10.<fabs(MHsecond - MHfirst)){
          boundary = (MHfirst*GHsecond + MHsecond*GHfirst)/GHsum;
          if (is==0 && MHfirst-4.*GHfirst>0.) GenHMassBinning.addBinBoundary(MHfirst-4.*GHfirst);
          if (MHfirst+4.*GHfirst<boundary && MHsecond-4.*GHsecond>boundary){
            GenHMassBinning.addBinBoundary(MHfirst+4.*GHfirst);
            GenHMassBinning.addBinBoundary(MHsecond-4.*GHsecond);
          }
          else GenHMassBinning.addBinBoundary(boundary);
        }
        else{
          boundary = (MHfirst + MHsecond)/2.;
          GenHMassBinning.addBinBoundary(boundary);
        }
      }
    }

  }
  else if (forceMass==125.){
    GenHMassBinning.addBinBoundary(124);
    GenHMassBinning.addBinBoundary(126);
    GenHMassBinning.addBinBoundary(160);
    GenHMassBinning.addBinBoundary(220);
    GenHMassBinning.addBinBoundary(450);
    GenHMassBinning.addBinBoundary(1300);
  }
  else{
    MELAerr << "ReweightingBuilder::getTrueMassBinning: Force mass " << forceMass << " is not implemented!" << endl;
    assert(0);
  }

  return GenHMassBinning;
}
