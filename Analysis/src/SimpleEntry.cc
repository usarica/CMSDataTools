#include <iostream>
#include "SimpleEntry.h"
#include "SimpleEntry.hpp"


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
  namedfloats(other.namedfloats),
  nameddoubles(other.nameddoubles),
  namedVbools(other.namedVbools),
  namedVshorts(other.namedVshorts),
  namedVuints(other.namedVuints),
  namedVints(other.namedVints),
  namedVfloats(other.namedVfloats),
  namedVdoubles(other.namedVdoubles)
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
  std::swap(namedfloats, other.namedfloats);
  std::swap(nameddoubles, other.nameddoubles);
  std::swap(namedVbools, other.namedVbools);
  std::swap(namedVshorts, other.namedVshorts);
  std::swap(namedVuints, other.namedVuints);
  std::swap(namedVints, other.namedVints);
  std::swap(namedVfloats, other.namedVfloats);
  std::swap(namedVdoubles, other.namedVdoubles);
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
  vector<unsigned int> erasepos;
  unsigned int pos=0;
  for (std::vector<SimpleEntry>::iterator it=vec.begin(); it!=vec.end(); it++){
    if (it->trackingval<minval || it->trackingval>maxval) erasepos.push_back(pos);
    pos++;
  }
  for (int ipos=(int) erasepos.size()-1; ipos>=0; ipos--) vec.erase(vec.begin()+erasepos.at(ipos));
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

