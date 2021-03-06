#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "IvyBase.h"
#include "IvyBase.hpp"


using namespace std;


IvyBase::IvyBase() :
  verbosity(TVar::ERROR),
  currentTree(nullptr),
  eventCache_currentTTree(nullptr),
  eventCache_currentEvent(-1)
{}
IvyBase::~IvyBase(){}


void IvyBase::defineConsumedSloppy(TString name){
  if (std::find(this->sloppyConsumes.begin(), this->sloppyConsumes.end(), name)==this->sloppyConsumes.end()){
    this->sloppyConsumes.push_back(name);
    if (this->verbosity>=TVar::INFO) MELAout << "IvyBase::defineConsumedSloppy: Consumed " << name << " will be treated sloppily." << endl;
  }
}

bool IvyBase::linkConsumes(BaseTree* tree){
  bool process = tree->isValid();
  if (!process && verbosity>=TVar::ERROR) MELAerr << "IvyBase::linkConsumes: Tree " << tree->sampleIdentifier << " is already invalid." << endl;
  if (process){
#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) process &= this->linkConsumed<type>(tree);
#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) process &= this->linkConsumed<type*>(tree);
#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) process &= this->linkConsumed<type*>(tree);
    SIMPLE_DATA_INPUT_DIRECTIVES
    VECTOR_DATA_INPUT_DIRECTIVES
    DOUBLEVECTOR_DATA_INPUT_DIRECTIVES
#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE
    // Silence unused branches
    tree->silenceUnused();
  }
  if (!process && verbosity>=TVar::ERROR) MELAerr << "IvyBase::linkConsumes: Linking failed for some reason for tree " << tree->sampleIdentifier << endl;
  return process;
}

bool IvyBase::cacheEvent(){
  return (
    this->currentTree
    &&
    this->currentTree->getCurrentEventInfo(this->eventCache_currentTTree, this->eventCache_currentEvent)
    );
}
bool IvyBase::isAlreadyCached() const{
  return (
    this->currentTree
    &&
    this->currentTree->isSameEvent(this->eventCache_currentTTree, this->eventCache_currentEvent)
    );
}
void IvyBase::resetCache(){
  this->eventCache_currentTTree = nullptr;
  this->eventCache_currentEvent = -1;
}

bool IvyBase::wrapTree(BaseTree* tree){
  this->currentTree = tree;
  this->eventCache_currentTTree = nullptr;
  this->eventCache_currentEvent = -1;
  if (!(this->currentTree)){
    if (this->verbosity>=TVar::ERROR) MELAerr << "IvyBase::wrapTree: The input tree is null!" << endl;
    return false;
  }
  return this->linkConsumes(this->currentTree);
}
