#ifndef BASETREE_HPP
#define BASETREE_HPP

#include "SampleHelpersCore.h"
#include "BaseTree.h"


#define SIMPLE_DATA_INPUT_DIRECTIVE(name, type, default_value) \
template<> bool BaseTree::getBranchCIterator<std::pair<type, type>*>(TString branchname, std::unordered_map<TString, std::pair<type, type>*>::iterator& it){ \
  auto& theMap = val##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.end()); \
} \
template<> bool BaseTree::getBranchCIterator<std::pair<type, type>*>(TString branchname, std::unordered_map<TString, std::pair<type, type>*>::const_iterator& it) const{ \
  auto const& theMap = val##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.cend()); \
} \
template<> void BaseTree::resetBranch<BaseTree::BranchType_##name##_t>(){ for (auto& it:val##name##s){ if (it.second){ it.second->first=it.second->second; } } } \
template<> void BaseTree::removeBranch<BaseTree::BranchType_##name##_t>(TString branchname){ for (auto& it:val##name##s){ if (it.first==branchname){ delete it.second; it.second=0; } } val##name##s.erase(branchname); } \
template<> bool BaseTree::bookBranch<type>(TString branchname, type valdef){ \
  if (val##name##s.find(branchname)==val##name##s.end()) val##name##s[branchname] = new std::pair<type, type>(valdef, valdef); \
  else{ val##name##s[branchname]->first=valdef; val##name##s[branchname]->second=valdef; } \
  SampleHelpers::bookBranch(tree, branchname, &(val##name##s[branchname]->first)); \
  SampleHelpers::bookBranch(failedtree, branchname, &(val##name##s[branchname]->first)); \
  return true; \
} \
template<> bool BaseTree::bookBranch<BaseTree::BranchType_##name##_t>(TString branchname){ return this->bookBranch<type>(branchname, default_value); } \
template<> bool BaseTree::putBranch<type>(TString branchname, type valdef){ \
  if (val##name##s.find(branchname)==val##name##s.end()) val##name##s[branchname] = new std::pair<type, type>(valdef, valdef); \
  else{ val##name##s[branchname]->first=valdef; val##name##s[branchname]->second=valdef; } \
  SampleHelpers::putBranch(tree, branchname, val##name##s[branchname]->first); \
  return true; \
} \
template<> bool BaseTree::putBranch<BaseTree::BranchType_##name##_t>(TString branchname){ return this->putBranch<type>(branchname, default_value); } \
template<> void BaseTree::getVal<type>(TString branchname, type& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it; \
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; } \
} \
template<> void BaseTree::setVal<type>(TString branchname, type const& val){ \
  typedef type itType; \
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it; \
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; } \
} \
template<> void BaseTree::getValRef<type>(TString branchname, type*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it; \
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); } \
} \
template<> void BaseTree::getValRef<type>(TString branchname, type*& val){ \
  typedef type itType; \
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it; \
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); } \
} \


#define VECTOR_DATA_INPUT_DIRECTIVE(name, type) \
template<> bool BaseTree::getBranchCIterator<type*>(TString branchname, std::unordered_map<TString, type*>::iterator& it){ \
  auto& theMap = valV##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.end()); \
} \
template<> bool BaseTree::getBranchCIterator<type*>(TString branchname, std::unordered_map<TString, type*>::const_iterator& it) const{ \
  auto const& theMap = valV##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.cend()); \
} \
template<> void BaseTree::resetBranch<BaseTree::BranchType_v##name##_t>(){ for (auto& it:valV##name##s){ if (it.second) it.second->clear(); } } \
template<> void BaseTree::removeBranch<BaseTree::BranchType_v##name##_t>(TString branchname){ for (auto& it:valV##name##s){ if (it.first==branchname){ delete it.second; it.second=0; } } valV##name##s.erase(branchname); } \
template<> bool BaseTree::bookBranch<type*>(TString branchname, type*/* valdef*/){ \
  valV##name##s[branchname] = nullptr; \
  SampleHelpers::bookBranch(tree, branchname, &(valV##name##s[branchname])); \
  SampleHelpers::bookBranch(failedtree, branchname, &(valV##name##s[branchname])); \
  return true; \
} \
template<> bool BaseTree::bookBranch<BaseTree::BranchType_v##name##_t>(TString branchname){ return this->bookBranch<type*>(branchname, nullptr); } \
template<> bool BaseTree::putBranch<type*>(TString branchname, type*/* valdef*/){ \
  if (valV##name##s.find(branchname)==valV##name##s.end()) valV##name##s[branchname] = new type(); \
  else valV##name##s[branchname]->clear(); \
  SampleHelpers::putBranch(tree, branchname, *(valV##name##s[branchname])); \
  return true; \
} \
template<> bool BaseTree::putBranch<BaseTree::BranchType_v##name##_t>(TString branchname){ return this->putBranch<type*>(branchname, nullptr); } \
template<> void BaseTree::getVal<type const*>(TString branchname, type const*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = it->second; \
} \
template<> void BaseTree::setVal<type*>(TString branchname, type* const& val){ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end()); \
} \
template<> void BaseTree::getValRef<type>(TString branchname, type*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = it->second; \
} \
template<> void BaseTree::getValRef<type* const>(TString branchname, type* const*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = &(it->second); \
} \
template<> void BaseTree::getValRef<type*>(TString branchname, type**& val){ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = &(it->second); \
} \


#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(name, type) \
template<> bool BaseTree::getBranchCIterator<type*>(TString branchname, std::unordered_map<TString, type*>::iterator& it){ \
  auto& theMap = valVV##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.end()); \
} \
template<> bool BaseTree::getBranchCIterator<type*>(TString branchname, std::unordered_map<TString, type*>::const_iterator& it) const{ \
  auto const& theMap = valVV##name##s; \
  it = theMap.find(branchname); \
  return (it!=theMap.cend()); \
} \
template<> void BaseTree::resetBranch<BaseTree::BranchType_vv##name##_t>(){ for (auto& it:valVV##name##s){ if (it.second) it.second->clear(); } } \
template<> void BaseTree::removeBranch<BaseTree::BranchType_vv##name##_t>(TString branchname){ for (auto& it:valVV##name##s){ if (it.first==branchname){ delete it.second; it.second=0; } } valVV##name##s.erase(branchname); } \
template<> bool BaseTree::bookBranch<type*>(TString branchname, type*/* valdef*/){ \
  valVV##name##s[branchname] = nullptr; \
  SampleHelpers::bookBranch(tree, branchname, &(valVV##name##s[branchname])); \
  SampleHelpers::bookBranch(failedtree, branchname, &(valVV##name##s[branchname])); \
  return true; \
} \
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vv##name##_t>(TString branchname){ return this->bookBranch<type*>(branchname, nullptr); } \
template<> bool BaseTree::putBranch<type*>(TString branchname, type*/* valdef*/){ \
  if (valVV##name##s.find(branchname)==valVV##name##s.end()) valVV##name##s[branchname] = new type(); \
  else valVV##name##s[branchname]->clear(); \
  SampleHelpers::putBranch(tree, branchname, *(valVV##name##s[branchname])); \
  return true; \
} \
template<> bool BaseTree::putBranch<BaseTree::BranchType_vv##name##_t>(TString branchname){ return this->putBranch<type*>(branchname, nullptr); } \
template<> void BaseTree::getVal<type const*>(TString branchname, type const*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = it->second; \
} \
template<> void BaseTree::setVal<type*>(TString branchname, type* const& val){ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end()); \
} \
template<> void BaseTree::getValRef<type>(TString branchname, type*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = it->second; \
} \
template<> void BaseTree::getValRef<type* const>(TString branchname, type* const*& val) const{ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::const_iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = &(it->second); \
} \
template<> void BaseTree::getValRef<type*>(TString branchname, type**& val){ \
  typedef type itType; \
  std::unordered_map<TString, itType*>::iterator it; \
  if (this->getBranchCIterator<itType*>(branchname, it)) val = &(it->second); \
} \


SIMPLE_DATA_INPUT_DIRECTIVES
VECTOR_DATA_INPUT_DIRECTIVES
DOUBLEVECTOR_DATA_INPUT_DIRECTIVES


#undef SIMPLE_DATA_INPUT_DIRECTIVE
#undef VECTOR_DATA_INPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_INPUT_DIRECTIVE


#endif
