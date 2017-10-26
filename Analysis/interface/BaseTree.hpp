#ifndef BASETREE_HPP
#define BASETREE_HPP

#include "BaseTree.h"


template<> bool BaseTree::getBranchCIterator<std::pair<short, short>*>(TString branchname, std::unordered_map<TString, std::pair<short, short>*>::iterator& it){
  auto& theMap = valshorts;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<short>*>(TString branchname, std::unordered_map<TString, std::vector<short>*>::iterator& it){
  auto& theMap = valVshorts;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned int, unsigned int>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::iterator& it){
  auto& theMap = valuints;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned int>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned int>*>::iterator& it){
  auto& theMap = valVuints;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<int, int>*>(TString branchname, std::unordered_map<TString, std::pair<int, int>*>::iterator& it){
  auto& theMap = valints;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<int>*>(TString branchname, std::unordered_map<TString, std::vector<int>*>::iterator& it){
  auto& theMap = valVints;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<float, float>*>(TString branchname, std::unordered_map<TString, std::pair<float, float>*>::iterator& it){
  auto& theMap = valfloats;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<float>*>(TString branchname, std::unordered_map<TString, std::vector<float>*>::iterator& it){
  auto& theMap = valVfloats;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<double, double>*>(TString branchname, std::unordered_map<TString, std::pair<double, double>*>::iterator& it){
  auto& theMap = valdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<double>*>(TString branchname, std::unordered_map<TString, std::vector<double>*>::iterator& it){
  auto& theMap = valVdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}

template<> bool BaseTree::getBranchCIterator<std::pair<short, short>*>(TString branchname, std::unordered_map<TString, std::pair<short, short>*>::const_iterator& it) const{
  auto const& theMap = valshorts;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<short>*>(TString branchname, std::unordered_map<TString, std::vector<short>*>::const_iterator& it) const{
  auto const& theMap = valVshorts;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned int, unsigned int>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned int, unsigned int>*>::const_iterator& it) const{
  auto const& theMap = valuints;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned int>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned int>*>::const_iterator& it) const{
  auto const& theMap = valVuints;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<int, int>*>(TString branchname, std::unordered_map<TString, std::pair<int, int>*>::const_iterator& it) const{
  auto const& theMap = valints;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<int>*>(TString branchname, std::unordered_map<TString, std::vector<int>*>::const_iterator& it) const{
  auto const& theMap = valVints;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<float, float>*>(TString branchname, std::unordered_map<TString, std::pair<float, float>*>::const_iterator& it) const{
  auto const& theMap = valfloats;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<float>*>(TString branchname, std::unordered_map<TString, std::vector<float>*>::const_iterator& it) const{
  auto const& theMap = valVfloats;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<double, double>*>(TString branchname, std::unordered_map<TString, std::pair<double, double>*>::const_iterator& it) const{
  auto const& theMap = valdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<double>*>(TString branchname, std::unordered_map<TString, std::vector<double>*>::const_iterator& it) const{
  auto const& theMap = valVdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}


template<> void BaseTree::resetBranch<BaseTree::BranchType_short_t>(){ for (auto& it:valshorts){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vshort_t>(){ for (auto& it:valVshorts){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_uint_t>(){ for (auto& it:valuints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vuint_t>(){ for (auto& it:valVuints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_int_t>(){ for (auto& it:valints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vint_t>(){ for (auto& it:valVints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_float_t>(){ for (auto& it:valfloats){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vfloat_t>(){ for (auto& it:valVfloats){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_double_t>(){ for (auto& it:valdoubles){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vdouble_t>(){ for (auto& it:valVdoubles){ if (it.second) it.second->clear(); } }

template<> void BaseTree::removeBranch<BaseTree::BranchType_short_t>(TString branchname){ for (auto& it:valshorts){ if (it.first){ delete it.second; it.second=0; } } valshorts.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vshort_t>(TString branchname){ for (auto& it:valVshorts){ if (it.first){ delete it.second; it.second=0; } } valVshorts.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_uint_t>(TString branchname){ for (auto& it:valuints){ if (it.first){ delete it.second; it.second=0; } } valuints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vuint_t>(TString branchname){ for (auto& it:valVuints){ if (it.first){ delete it.second; it.second=0; } } valVuints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_int_t>(TString branchname){ for (auto& it:valints){ if (it.first){ delete it.second; it.second=0; } } valints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vint_t>(TString branchname){ for (auto& it:valVints){ if (it.first){ delete it.second; it.second=0; } } valVints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_float_t>(TString branchname){ for (auto& it:valfloats){ if (it.first){ delete it.second; it.second=0; } } valfloats.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ for (auto& it:valVfloats){ if (it.first){ delete it.second; it.second=0; } } valVfloats.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_double_t>(TString branchname){ for (auto& it:valdoubles){ if (it.first){ delete it.second; it.second=0; } } valdoubles.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ for (auto& it:valVdoubles){ if (it.first){ delete it.second; it.second=0; } } valVdoubles.erase(branchname); }

template<> bool BaseTree::bookBranch<short>(TString branchname, short valdef){
  if (valshorts.find(branchname)==valshorts.end()) valshorts[branchname] = new std::pair<short, short>(valdef, valdef);
  else{ valshorts[branchname]->first=valdef; valshorts[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valshorts[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valshorts[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<unsigned int>(TString branchname, unsigned int valdef){
  if (valuints.find(branchname)==valuints.end()) valuints[branchname] = new std::pair<unsigned int, unsigned int>(valdef, valdef);
  else{ valuints[branchname]->first=valdef; valuints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valuints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valuints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<int>(TString branchname, int valdef){
  if (valints.find(branchname)==valints.end()) valints[branchname] = new std::pair<int, int>(valdef, valdef);
  else{ valints[branchname]->first=valdef; valints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<float>(TString branchname, float valdef){
  if (valfloats.find(branchname)==valfloats.end()) valfloats[branchname] = new std::pair<float, float>(valdef, valdef);
  else{ valfloats[branchname]->first=valdef; valfloats[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valfloats[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valfloats[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<double>(TString branchname, double valdef){
  if (valdoubles.find(branchname)==valdoubles.end()) valdoubles[branchname] = new std::pair<double, double>(valdef, valdef);
  else{ valdoubles[branchname]->first=valdef; valdoubles[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valdoubles[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valdoubles[branchname]->first));
  return true;
}

template<> bool BaseTree::bookBranch<std::vector<short>*>(TString branchname, std::vector<short>*/* valdef*/){
  valVshorts[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVshorts[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVshorts[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>*/* valdef*/){
  valVuints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVuints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVuints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<int>*>(TString branchname, std::vector<int>*/* valdef*/){
  valVints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<float>*>(TString branchname, std::vector<float>*/* valdef*/){
  valVfloats[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVfloats[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVfloats[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<double>*>(TString branchname, std::vector<double>*/* valdef*/){
  valVdoubles[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVdoubles[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVdoubles[branchname]));
  return true;
}

template<> bool BaseTree::bookBranch<BaseTree::BranchType_short_t>(TString branchname){ return this->bookBranch<short>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_uint_t>(TString branchname){ return this->bookBranch<unsigned int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_int_t>(TString branchname){ return this->bookBranch<int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_float_t>(TString branchname){ return this->bookBranch<float>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_double_t>(TString branchname){ return this->bookBranch<double>(branchname, 0); }

template<> bool BaseTree::bookBranch<BaseTree::BranchType_vshort_t>(TString branchname){ return this->bookBranch<std::vector<short>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vuint_t>(TString branchname){ return this->bookBranch<std::vector<unsigned int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vint_t>(TString branchname){ return this->bookBranch<std::vector<int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ return this->bookBranch<std::vector<float>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ return this->bookBranch<std::vector<double>*>(branchname, 0); }

template<> void BaseTree::getVal<short>(TString branchname, short& val) const{
  typedef short itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<short> const*>(TString branchname, std::vector<short> const*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<unsigned int>(TString branchname, unsigned int& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<unsigned int> const*>(TString branchname, std::vector<unsigned int> const*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<int>(TString branchname, int& val) const{
  typedef int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<int> const*>(TString branchname, std::vector<int> const*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<float>(TString branchname, float& val) const{
  typedef float itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<float> const*>(TString branchname, std::vector<float> const*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<double>(TString branchname, double& val) const{
  typedef double itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<double> const*>(TString branchname, std::vector<double> const*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}

template<> void BaseTree::setVal<short>(TString branchname, short const& val){
  typedef short itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<short>*>(TString branchname, std::vector<short>* const& val){
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<unsigned int>(TString branchname, unsigned int const& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>* const& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<int>(TString branchname, int const& val){
  typedef int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<int>*>(TString branchname, std::vector<int>* const& val){
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<float>(TString branchname, float const& val){
  typedef float itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<float>*>(TString branchname, std::vector<float>* const& val){
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<double>(TString branchname, double const& val){
  typedef double itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<double>*>(TString branchname, std::vector<double>* const& val){
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}

template<> void BaseTree::getValRef<short>(TString branchname, short*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::vector<short>>(TString branchname, std::vector<short>*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<unsigned int>(TString branchname, unsigned int*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::vector<unsigned int>>(TString branchname, std::vector<unsigned int>*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<int>(TString branchname, int*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::vector<int>>(TString branchname, std::vector<int>*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<float>(TString branchname, float*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::vector<float>>(TString branchname, std::vector<float>*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<double>(TString branchname, double*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::vector<double>>(TString branchname, std::vector<double>*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}


#endif
