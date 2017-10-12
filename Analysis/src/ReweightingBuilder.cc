#include "ReweightingBuilder.h"
#include "SimpleEntry.h"

using namespace std;
using namespace HelperFunctions;


ReweightingBuilder::ReweightingBuilder(TString inStrWeight, ReweightingBuilder::BWScheme inscheme) :
theScheme(inscheme)
{
  strWeights.push_back(inStrWeight);
}

ReweightingBuilder::ReweightingBuilder(std::vector<TString> inStrWeights, ReweightingBuilder::BWScheme inscheme) :
theScheme(inscheme),
strWeights(inStrWeights)
{}

float ReweightingBuilder::eval(CJLSTTree* theTree) const{
  if (!theTree) return 0;

  float res=1;
  for (TString const& s : strWeights){
    float w=1;
    theTree->getVal(s, w);
    res *= w;
  }

  if (theScheme!=kDoNotConsider){
    float prop = theTree->getTrueBW();
    if (theScheme==kDivideSampleBW) res /= prop;
    else if (theScheme==kMultiplySampleBW) res *= prop;
  }

  return res;
}

std::vector<float> ReweightingBuilder::determineWeightThresholds(CJLSTTree* theTree, std::vector<std::pair<float, float>>& boundaries, float fractionRequirement, TString strOrderingVal){
  std::vector<float> res;

  if (!theTree) return res;
  const bool noBoundaries = boundaries.empty();

  vector<vector<SimpleEntry>> indexList;
  if (noBoundaries) indexList.push_back(vector<SimpleEntry>());
  else{
    indexList.reserve(boundaries.size());
    for (unsigned int ib=0; ib<boundaries.size(); ib++) indexList.push_back(vector<SimpleEntry>());
  }

  int ev=0;
  while (theTree->getEvent(ev)){
    ev++;

    float wgt = fabs(this->eval(theTree));
    float orderingVal=-1;
    theTree->getVal(strOrderingVal, orderingVal);
    SimpleEntry theEntry(0, wgt);

    if (noBoundaries) addByLowest(indexList.front(), theEntry, false);
    else{
      vector<vector<SimpleEntry>>::iterator indexIt=indexList.begin();
      vector<pair<float, float>>::iterator boundaryIt=boundaries.begin();
      while (boundaryIt!=boundaries.end() && indexIt!=indexList.end()){
        if (boundaryIt->first<=orderingVal && orderingVal<boundaryIt->second){
          addByLowest(*indexIt, theEntry, false);
          break;
        }
        indexIt++;
        boundaryIt++;
      }
    }
  }

  for (auto const& index:indexList){
    unsigned int index_entry = static_cast<int>((float(index.size())*fractionRequirement));
    if (index_entry==index.size() && index_entry>0) index_entry--;
    unsigned int index_entry_prev=index_entry; if (index_entry_prev>0) index_entry_prev--;
    float threshold = -1;
    if (!index.empty()) threshold = (index.at(index_entry_prev).trackingval + index.at(index_entry).trackingval)*0.5;
    if (threshold>0. && index.back().trackingval<threshold*5.) threshold = index.back().trackingval; // Prevent false-positives
    res.push_back(threshold);
  }

  weightThresholds[theTree] = res;
  return res;
}
std::vector<float> ReweightingBuilder::getWeightThresholds(CJLSTTree* theTree){
  if (theTree && weightThresholds.find(theTree)!=weightThresholds.end()) return weightThresholds[theTree];
  else return vector<float>();
}


