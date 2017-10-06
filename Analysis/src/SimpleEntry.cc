#include "SimpleEntry.h"
#include <iostream>

using namespace std;


SimpleEntry::SimpleEntry() : id(0), trackingval(0), weight(0) {}
SimpleEntry::SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_) : id(id_), trackingval(trackingval_), recoval(recoval_), weight(weight_) {}

bool SimpleEntry::operator != (const SimpleEntry& other)const{ return trackingval!=other.trackingval; }
bool SimpleEntry::operator == (const SimpleEntry& other)const{ return trackingval==other.trackingval; }
bool SimpleEntry::operator > (const SimpleEntry& other)const{ return trackingval>other.trackingval; }
bool SimpleEntry::operator >= (const SimpleEntry& other)const{ return trackingval>=other.trackingval; }
bool SimpleEntry::operator < (const SimpleEntry& other)const{ return trackingval<other.trackingval; }
bool SimpleEntry::operator <= (const SimpleEntry& other)const{ return trackingval<=other.trackingval; }

void SimpleEntry::setNamedVal(TString strname, unsigned int& val){ nameduints[strname]=val; }
void SimpleEntry::setNamedVal(TString strname, short& val){ namedshorts[strname]=val; }
void SimpleEntry::setNamedVal(TString strname, int& val){ namedints[strname]=val; }
void SimpleEntry::setNamedVal(TString strname, float& val){ namedfloats[strname]=val; }
void SimpleEntry::setNamedVal(TString strname, double& val){ nameddoubles[strname]=val; }

void SimpleEntry::getNamedVal(TString strname, unsigned int& val){ val = nameduints[strname]; }
void SimpleEntry::getNamedVal(TString strname, short& val){ val = namedshorts[strname]; }
void SimpleEntry::getNamedVal(TString strname, int& val){ val = namedints[strname]; }
void SimpleEntry::getNamedVal(TString strname, float& val){ val = namedfloats[strname]; }
void SimpleEntry::getNamedVal(TString strname, double& val){ val = nameddoubles[strname]; }

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

