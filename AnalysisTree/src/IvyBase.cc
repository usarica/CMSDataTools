#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "IvyBase.h"
#include "IvyBase.hpp"


using namespace std;


IvyBase::IvyBase() : verbosity(TVar::ERROR), currentTree(nullptr) {}
IvyBase::~IvyBase(){}


void IvyBase::defineConsumedSloppy(TString name){
  if (std::find(this->sloppyConsumes.begin(), this->sloppyConsumes.end(), name)==this->sloppyConsumes.end()){
    this->sloppyConsumes.push_back(name);
    if (this->verbosity>=TVar::INFO) MELAout << "IvyBase::defineConsumedSloppy: Consumed " << name << " will be treated sloppily." << endl;
  }
}

bool IvyBase::linkConsumes(BaseTree* tree){
  bool process = tree->isValid();
  if (process){
    process &= this->linkConsumed<bool>(tree);
    process &= this->linkConsumed<short>(tree);
    process &= this->linkConsumed<unsigned int>(tree);
    process &= this->linkConsumed<int>(tree);
    process &= this->linkConsumed<unsigned long>(tree);
    process &= this->linkConsumed<long>(tree);
    process &= this->linkConsumed<unsigned long long>(tree);
    process &= this->linkConsumed<long long>(tree);
    process &= this->linkConsumed<float>(tree);
    process &= this->linkConsumed<double>(tree);
    process &= this->linkConsumed<std::string>(tree);
    process &= this->linkConsumed<TString>(tree);
    process &= this->linkConsumed<CMSLorentzVector>(tree);
    process &= this->linkConsumed<std::vector<bool>*>(tree);
    process &= this->linkConsumed<std::vector<short>*>(tree);
    process &= this->linkConsumed<std::vector<unsigned int>*>(tree);
    process &= this->linkConsumed<std::vector<int>*>(tree);
    process &= this->linkConsumed<std::vector<unsigned long>*>(tree);
    process &= this->linkConsumed<std::vector<long>*>(tree);
    process &= this->linkConsumed<std::vector<unsigned long long>*>(tree);
    process &= this->linkConsumed<std::vector<long long>*>(tree);
    process &= this->linkConsumed<std::vector<float>*>(tree);
    process &= this->linkConsumed<std::vector<double>*>(tree);
    process &= this->linkConsumed<std::vector<std::string>*>(tree);
    process &= this->linkConsumed<std::vector<TString>*>(tree);
    process &= this->linkConsumed<std::vector<CMSLorentzVector>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<bool>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<short>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<unsigned int>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<int>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<unsigned long>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<long>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<unsigned long long>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<long long>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<float>>*>(tree);
    process &= this->linkConsumed<std::vector<std::vector<double>>*>(tree);
    // Silence unused branches
    tree->silenceUnused();
  }
  if (!process && verbosity>=TVar::ERROR) MELAerr << "IvyBase::linkConsumes: Linking failed for some reason for tree " << tree->sampleIdentifier << endl;
  return process;
}

bool IvyBase::wrapTree(BaseTree* tree){
  this->currentTree = tree;
  if (!(this->currentTree)){
    if (this->verbosity>=TVar::ERROR) MELAerr << "IvyBase::wrapTree: The input tree is null!" << endl;
    return false;
  }
  return this->linkConsumes(this->currentTree);
}
