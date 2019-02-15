#include <iostream>
#include "SimpleEntry.h"
#include "SimpleEntry.hpp"
#include "SampleHelpersCore.h"


using namespace std;


SimpleEntry::SimpleEntry() : id(0), trackingval(0), weight(0) {}
SimpleEntry::SimpleEntry(int id_, float trackingval_, float weight_) : id(id_), trackingval(trackingval_), weight(weight_) {}
SimpleEntry::SimpleEntry(int id_, float trackingval_, std::vector<float> recoval_, float weight_) : id(id_), trackingval(trackingval_), weight(weight_), recoval(recoval_) {}
SimpleEntry::SimpleEntry(SimpleEntry const& other) :
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) named##name_t##s(other.named##name_t##s),
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) namedV##name_t##s(other.namedV##name_t##s),
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) namedVV##name_t##s(other.namedVV##name_t##s),
  SIMPLE_DATA_OUTPUT_DIRECTIVES
  VECTOR_DATA_OUTPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
  id(other.id),
  trackingval(other.trackingval),
  weight(other.weight),
  recoval(other.recoval)
{}
SimpleEntry::SimpleEntry(SimpleEntry&& other) : id(0), trackingval(0), weight(0)
{ this->swap(other); }

void SimpleEntry::swap(SimpleEntry& other){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) std::swap(named##name_t##s, other.named##name_t##s);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) std::swap(namedV##name_t##s, other.namedV##name_t##s);
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) std::swap(namedVV##name_t##s, other.namedVV##name_t##s);
  SIMPLE_DATA_OUTPUT_DIRECTIVES
  VECTOR_DATA_OUTPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
  std::swap(id, other.id);
  std::swap(trackingval, other.trackingval);
  std::swap(weight, other.weight);
  std::swap(recoval, other.recoval);
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
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.named##name_t##s.begin(); itb!=entry.named##name_t##s.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedV##name_t##s.begin(); itb!=entry.namedV##name_t##s.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=entry.namedVV##name_t##s.begin(); itb!=entry.namedVV##name_t##s.end(); itb++) commonEntry.setNamedVal(itb->first, itb->second);
    SIMPLE_DATA_OUTPUT_DIRECTIVES
    VECTOR_DATA_OUTPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
    if (it==vecBegin){
#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.named##name_t##s.begin(); itb!=commonEntry.named##name_t##s.end(); itb++) SampleHelpers::putBranch(tree, itb->first, itb->second);
#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedV##name_t##s.begin(); itb!=commonEntry.namedV##name_t##s.end(); itb++) SampleHelpers::putBranch(tree, itb->first, itb->second);
#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) for (auto itb=commonEntry.namedVV##name_t##s.begin(); itb!=commonEntry.namedVV##name_t##s.end(); itb++) SampleHelpers::putBranch(tree, itb->first, itb->second);
      SIMPLE_DATA_OUTPUT_DIRECTIVES
      VECTOR_DATA_OUTPUT_DIRECTIVES
      DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES
#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE
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

