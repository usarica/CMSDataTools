#ifndef _COMPILE_STANDALONE_

#include "TSystem.h"
#include "TDirectory.h"
#include "BaseEDMInputTree.h"
#include "BaseEDMInputTree.hpp"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


BaseEDMInputTree::BaseEDMInputTree() :
  BaseTree()
{}
BaseEDMInputTree::BaseEDMInputTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname) :
  BaseTree(cinput, treename, failedtreename, countersname)
{}
BaseEDMInputTree::~BaseEDMInputTree(){
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) HelperFunctions::cleanUnorderedMap(bridge##name##s);
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) HelperFunctions::cleanUnorderedMap(bridgeV##name##s);
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) HelperFunctions::cleanUnorderedMap(bridgeVV##name##s);
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}

void BaseEDMInputTree::synchronizeEDMBranches(){
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) for (auto& it:bridge##name##s){ if (it.second) it.second->synchronize(); }
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto& it:bridgeV##name##s){ if (it.second) it.second->synchronize(); }
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto& it:bridgeVV##name##s){ if (it.second) it.second->synchronize(); }
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}

void BaseEDMInputTree::print() const{
  BaseTree::print();

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) for (auto const& it:bridge##name##s){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto const& it:bridgeV##name##s){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) for (auto const& it:bridgeVV##name##s){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}


bool BaseEDMInputTree::getSelectedEvent(int ev){
  bool result = BaseTree::getSelectedEvent(ev);
  synchronizeEDMBranches();
  return result;
}
bool BaseEDMInputTree::getFailedEvent(int ev){
  bool result = BaseTree::getFailedEvent(ev);
  synchronizeEDMBranches();
  return result;
}

void BaseEDMInputTree::resetBranches(){
  BaseTree::resetBranches();

#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) this->resetEDMBranch<BaseTree::BranchType_##name##_t>();
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) this->resetEDMBranch<BaseTree::BranchType_v##name##_t>();
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) this->resetEDMBranch<BaseTree::BranchType_vv##name##_t>();
  SIMPLE_DATA_INPUT_DIRECTIVES
  VECTOR_DATA_INPUT_DIRECTIVES
  DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
}


#endif
