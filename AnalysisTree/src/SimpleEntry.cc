#include <iostream>
#include "SimpleEntry.h"
#include "SimpleEntry.hpp"
#include "SampleHelpersCore.h"


using namespace std;


SimpleEntry::SimpleEntry() : id(0), trackingval(0), weight(0) {}
SimpleEntry::SimpleEntry(int id_, float trackingval_, float weight_) : id(id_), trackingval(trackingval_), weight(weight_) {}
SimpleEntry::SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_) : id(id_), trackingval(trackingval_), weight(weight_), recoval(recoval_) {}
SimpleEntry::SimpleEntry(SimpleEntry const& other) :
  id(other.id),
  trackingval(other.trackingval),
  weight(other.weight),
  recoval(other.recoval),
  namedbools(other.namedbools),
  namedshorts(other.namedshorts),
  nameduints(other.nameduints),
  namedints(other.namedints),
  namedulongs(other.namedulongs),
  namedlongs(other.namedlongs),
  namedlonglongs(other.namedlonglongs),
  namedfloats(other.namedfloats),
  nameddoubles(other.nameddoubles),
  namedCMSLorentzVectors(other.namedCMSLorentzVectors),
  namedVbools(other.namedVbools),
  namedVshorts(other.namedVshorts),
  namedVuints(other.namedVuints),
  namedVints(other.namedVints),
  namedVulongs(other.namedVulongs),
  namedVlongs(other.namedVlongs),
  namedVlonglongs(other.namedVlonglongs),
  namedVfloats(other.namedVfloats),
  namedVdoubles(other.namedVdoubles),
  namedVCMSLorentzVectors(other.namedVCMSLorentzVectors)
{}
SimpleEntry::SimpleEntry(SimpleEntry&& other) : id(0), trackingval(0), weight(0)
{ this->swap(other); }

void SimpleEntry::swap(SimpleEntry& other){
  std::swap(id, other.id);
  std::swap(trackingval, other.trackingval);
  std::swap(weight, other.weight);
  std::swap(recoval, other.recoval);
  std::swap(namedbools, other.namedbools);
  std::swap(namedshorts, other.namedshorts);
  std::swap(nameduints, other.nameduints);
  std::swap(namedints, other.namedints);
  std::swap(namedulongs, other.namedulongs);
  std::swap(namedlongs, other.namedlongs);
  std::swap(namedlonglongs, other.namedlonglongs);
  std::swap(namedfloats, other.namedfloats);
  std::swap(nameddoubles, other.nameddoubles);
  std::swap(namedCMSLorentzVectors, other.namedCMSLorentzVectors);
  std::swap(namedVbools, other.namedVbools);
  std::swap(namedVshorts, other.namedVshorts);
  std::swap(namedVuints, other.namedVuints);
  std::swap(namedVints, other.namedVints);
  std::swap(namedVulongs, other.namedVulongs);
  std::swap(namedVlongs, other.namedVlongs);
  std::swap(namedVlonglongs, other.namedVlonglongs);
  std::swap(namedVfloats, other.namedVfloats);
  std::swap(namedVdoubles, other.namedVdoubles);
  std::swap(namedVCMSLorentzVectors, other.namedVCMSLorentzVectors);
}
SimpleEntry& SimpleEntry::operator=(const SimpleEntry& other){
  SimpleEntry tmp(other);
  this->swap(tmp);
  return *this;
}

bool SimpleEntry::operator != (const SimpleEntry& other)const{ return trackingval!=other.trackingval; }
bool SimpleEntry::operator == (const SimpleEntry& other)const{ return trackingval==other.trackingval; }
bool SimpleEntry::operator > (const SimpleEntry& other)const{ return trackingval>other.trackingval; }
bool SimpleEntry::operator >= (const SimpleEntry& other)const{ return trackingval>=other.trackingval; }
bool SimpleEntry::operator < (const SimpleEntry& other)const{ return trackingval<other.trackingval; }
bool SimpleEntry::operator <= (const SimpleEntry& other)const{ return trackingval<=other.trackingval; }

void SimpleEntry::cropByTrueVal(std::vector<SimpleEntry>& vec, float minval, float maxval){
  std::vector<unsigned int> erasepos;
  unsigned int pos=0;
  for (std::vector<SimpleEntry>::const_iterator it=vec.cbegin(); it!=vec.cend(); it++){
    if (it->trackingval<minval || it->trackingval>maxval) erasepos.push_back(pos);
    pos++;
  }
  for (std::vector<unsigned int>::reverse_iterator ipos=erasepos.rbegin(); ipos!=erasepos.rend(); ipos++) vec.erase(vec.begin()+(*ipos));
}
void SimpleEntry::writeToTree(std::vector<SimpleEntry>::const_iterator const& vecBegin, std::vector<SimpleEntry>::const_iterator const& vecEnd, TTree* const& tree){
  if (!tree) return;
  SimpleEntry commonEntry;
  for (std::vector<SimpleEntry>::const_iterator it=vecBegin; it!=vecEnd; it++){
    SimpleEntry const& entry = *it;
    for (auto itb=entry.namedbools.begin(); itb!=entry.namedbools.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedshorts.begin(); itb!=entry.namedshorts.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.nameduints.begin(); itb!=entry.nameduints.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedints.begin(); itb!=entry.namedints.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedulongs.begin(); itb!=entry.namedulongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedlongs.begin(); itb!=entry.namedlongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedlonglongs.begin(); itb!=entry.namedlonglongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedfloats.begin(); itb!=entry.namedfloats.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.nameddoubles.begin(); itb!=entry.nameddoubles.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedCMSLorentzVectors.begin(); itb!=entry.namedCMSLorentzVectors.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVbools.begin(); itb!=entry.namedVbools.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVshorts.begin(); itb!=entry.namedVshorts.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVuints.begin(); itb!=entry.namedVuints.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVints.begin(); itb!=entry.namedVints.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVulongs.begin(); itb!=entry.namedVulongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVlongs.begin(); itb!=entry.namedVlongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVlonglongs.begin(); itb!=entry.namedVlonglongs.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVfloats.begin(); itb!=entry.namedVfloats.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVdoubles.begin(); itb!=entry.namedVdoubles.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    for (auto itb=entry.namedVCMSLorentzVectors.begin(); itb!=entry.namedVCMSLorentzVectors.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    if (it==vecBegin){
      for (auto itb=commonEntry.namedbools.begin(); itb!=commonEntry.namedbools.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedbools[itb->first]);
      for (auto itb=commonEntry.namedshorts.begin(); itb!=commonEntry.namedshorts.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedshorts[itb->first]);
      for (auto itb=commonEntry.nameduints.begin(); itb!=commonEntry.nameduints.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.nameduints[itb->first]);
      for (auto itb=commonEntry.namedints.begin(); itb!=commonEntry.namedints.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedints[itb->first]);
      for (auto itb=commonEntry.namedulongs.begin(); itb!=commonEntry.namedulongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedulongs[itb->first]);
      for (auto itb=commonEntry.namedlongs.begin(); itb!=commonEntry.namedlongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedlongs[itb->first]);
      for (auto itb=commonEntry.namedlonglongs.begin(); itb!=commonEntry.namedlonglongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedlonglongs[itb->first]);
      for (auto itb=commonEntry.namedfloats.begin(); itb!=commonEntry.namedfloats.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedfloats[itb->first]);
      for (auto itb=commonEntry.nameddoubles.begin(); itb!=commonEntry.nameddoubles.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.nameddoubles[itb->first]);
      for (auto itb=commonEntry.namedCMSLorentzVectors.begin(); itb!=commonEntry.namedCMSLorentzVectors.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedCMSLorentzVectors[itb->first]);
      for (auto itb=commonEntry.namedVbools.begin(); itb!=commonEntry.namedVbools.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVbools[itb->first]);
      for (auto itb=commonEntry.namedVshorts.begin(); itb!=commonEntry.namedVshorts.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVshorts[itb->first]);
      for (auto itb=commonEntry.namedVuints.begin(); itb!=commonEntry.namedVuints.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVuints[itb->first]);
      for (auto itb=commonEntry.namedVints.begin(); itb!=commonEntry.namedVints.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVints[itb->first]);
      for (auto itb=commonEntry.namedVulongs.begin(); itb!=commonEntry.namedVulongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVulongs[itb->first]);
      for (auto itb=commonEntry.namedVlongs.begin(); itb!=commonEntry.namedVlongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVlongs[itb->first]);
      for (auto itb=commonEntry.namedVlonglongs.begin(); itb!=commonEntry.namedVlonglongs.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVlonglongs[itb->first]);
      for (auto itb=commonEntry.namedVfloats.begin(); itb!=commonEntry.namedVfloats.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVfloats[itb->first]);
      for (auto itb=commonEntry.namedVdoubles.begin(); itb!=commonEntry.namedVdoubles.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVdoubles[itb->first]);
      for (auto itb=commonEntry.namedVCMSLorentzVectors.begin(); itb!=commonEntry.namedVCMSLorentzVectors.end(); itb++) SampleHelpers::putBranch(tree, itb->first, commonEntry.namedVCMSLorentzVectors[itb->first]);
    }
    tree->Fill();
  }
}

void SimpleEntry::print(){
  cout << "Simple entry:" << endl;
  cout << " - Id = " << id << endl;
  cout << " - Weight = " << weight << endl;
  cout << " - Trueval: " << trackingval << endl;
  cout << " - Recoval: ";
  for (auto& v : recoval) cout << v << " ";
  cout << endl;
}

void ExtBin::addEvent(SimpleEntry evt){ collection.push_back(evt); }

