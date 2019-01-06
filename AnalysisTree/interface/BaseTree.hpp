#ifndef BASETREE_HPP
#define BASETREE_HPP

#include "SampleHelpersCore.h"
#include "BaseTree.h"


template<> bool BaseTree::getBranchCIterator<std::pair<bool, bool>*>(TString branchname, std::unordered_map<TString, std::pair<bool, bool>*>::iterator& it){
  auto& theMap = valbools;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<bool>*>(TString branchname, std::unordered_map<TString, std::vector<bool>*>::iterator& it){
  auto& theMap = valVbools;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<bool>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<bool>>*>::iterator& it){
  auto& theMap = valVVbools;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<short>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<short>>*>::iterator& it){
  auto& theMap = valVVshorts;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned int>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned int>>*>::iterator& it){
  auto& theMap = valVVuints;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<int>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<int>>*>::iterator& it){
  auto& theMap = valVVints;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned long, unsigned long>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>::iterator& it){
  auto& theMap = valulongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned long>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned long>*>::iterator& it){
  auto& theMap = valVulongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned long>>*>::iterator& it){
  auto& theMap = valVVulongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<long, long>*>(TString branchname, std::unordered_map<TString, std::pair<long, long>*>::iterator& it){
  auto& theMap = vallongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<long>*>(TString branchname, std::unordered_map<TString, std::vector<long>*>::iterator& it){
  auto& theMap = valVlongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<long>>*>::iterator& it){
  auto& theMap = valVVlongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned long long, unsigned long long>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>::iterator& it){
  auto& theMap = valulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned long long>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned long long>*>::iterator& it){
  auto& theMap = valVulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned long long>>*>::iterator& it){
  auto& theMap = valVVulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<long long, long long>*>(TString branchname, std::unordered_map<TString, std::pair<long long, long long>*>::iterator& it){
  auto& theMap = vallonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<long long>*>(TString branchname, std::unordered_map<TString, std::vector<long long>*>::iterator& it){
  auto& theMap = valVlonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<long long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<long long>>*>::iterator& it){
  auto& theMap = valVVlonglongs;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<float>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<float>>*>::iterator& it){
  auto& theMap = valVVfloats;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<double>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<double>>*>::iterator& it){
  auto& theMap = valVVdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<std::string, std::string>*>(TString branchname, std::unordered_map<TString, std::pair<std::string, std::string>*>::iterator& it){
  auto& theMap = valstrings;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::string>*>(TString branchname, std::unordered_map<TString, std::vector<std::string>*>::iterator& it){
  auto& theMap = valVstrings;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<TString, TString>*>(TString branchname, std::unordered_map<TString, std::pair<TString, TString>*>::iterator& it){
  auto& theMap = valTStrings;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<TString>*>(TString branchname, std::unordered_map<TString, std::vector<TString>*>::iterator& it){
  auto& theMap = valVTStrings;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::pair<CMSLorentzVector, CMSLorentzVector>*>(TString branchname, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>::iterator& it){
  auto& theMap = valCMSLorentzVectors;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}
template<> bool BaseTree::getBranchCIterator<std::vector<CMSLorentzVector>*>(TString branchname, std::unordered_map<TString, std::vector<CMSLorentzVector>*>::iterator& it){
  auto& theMap = valVCMSLorentzVectors;
  it = theMap.find(branchname);
  return (it!=theMap.end());
}

template<> bool BaseTree::getBranchCIterator<std::pair<bool, bool>*>(TString branchname, std::unordered_map<TString, std::pair<bool, bool>*>::const_iterator& it) const{
  auto const& theMap = valbools;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<bool>*>(TString branchname, std::unordered_map<TString, std::vector<bool>*>::const_iterator& it) const{
  auto const& theMap = valVbools;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<bool>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<bool>>*>::const_iterator& it) const{
  auto const& theMap = valVVbools;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<short>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<short>>*>::const_iterator& it) const{
  auto const& theMap = valVVshorts;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned int>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned int>>*>::const_iterator& it) const{
  auto const& theMap = valVVuints;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<int>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<int>>*>::const_iterator& it) const{
  auto const& theMap = valVVints;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned long, unsigned long>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned long, unsigned long>*>::const_iterator& it) const{
  auto const& theMap = valulongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned long>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned long>*>::const_iterator& it) const{
  auto const& theMap = valVulongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned long>>*>::const_iterator& it) const{
  auto const& theMap = valVVulongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<long, long>*>(TString branchname, std::unordered_map<TString, std::pair<long, long>*>::const_iterator& it) const{
  auto const& theMap = vallongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<long>*>(TString branchname, std::unordered_map<TString, std::vector<long>*>::const_iterator& it) const{
  auto const& theMap = valVlongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<long>>*>::const_iterator& it) const{
  auto const& theMap = valVVlongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<unsigned long long, unsigned long long>*>(TString branchname, std::unordered_map<TString, std::pair<unsigned long long, unsigned long long>*>::const_iterator& it) const{
  auto const& theMap = valulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<unsigned long long>*>(TString branchname, std::unordered_map<TString, std::vector<unsigned long long>*>::const_iterator& it) const{
  auto const& theMap = valVulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<unsigned long long>>*>::const_iterator& it) const{
  auto const& theMap = valVVulonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<long long, long long>*>(TString branchname, std::unordered_map<TString, std::pair<long long, long long>*>::const_iterator& it) const{
  auto const& theMap = vallonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<long long>*>(TString branchname, std::unordered_map<TString, std::vector<long long>*>::const_iterator& it) const{
  auto const& theMap = valVlonglongs;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<long long>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<long long>>*>::const_iterator& it) const{
  auto const& theMap = valVVlonglongs;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<float>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<float>>*>::const_iterator& it) const{
  auto const& theMap = valVVfloats;
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
template<> bool BaseTree::getBranchCIterator<std::vector<std::vector<double>>*>(TString branchname, std::unordered_map<TString, std::vector<std::vector<double>>*>::const_iterator& it) const{
  auto const& theMap = valVVdoubles;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<std::string, std::string>*>(TString branchname, std::unordered_map<TString, std::pair<std::string, std::string>*>::const_iterator& it) const{
  auto const& theMap = valstrings;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<std::string>*>(TString branchname, std::unordered_map<TString, std::vector<std::string>*>::const_iterator& it) const{
  auto const& theMap = valVstrings;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<TString, TString>*>(TString branchname, std::unordered_map<TString, std::pair<TString, TString>*>::const_iterator& it) const{
  auto const& theMap = valTStrings;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<TString>*>(TString branchname, std::unordered_map<TString, std::vector<TString>*>::const_iterator& it) const{
  auto const& theMap = valVTStrings;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::pair<CMSLorentzVector, CMSLorentzVector>*>(TString branchname, std::unordered_map<TString, std::pair<CMSLorentzVector, CMSLorentzVector>*>::const_iterator& it) const{
  auto const& theMap = valCMSLorentzVectors;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}
template<> bool BaseTree::getBranchCIterator<std::vector<CMSLorentzVector>*>(TString branchname, std::unordered_map<TString, std::vector<CMSLorentzVector>*>::const_iterator& it) const{
  auto const& theMap = valVCMSLorentzVectors;
  it = theMap.find(branchname);
  return (it!=theMap.cend());
}


template<> void BaseTree::resetBranch<BaseTree::BranchType_bool_t>(){ for (auto& it:valbools){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vbool_t>(){ for (auto& it:valVbools){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvbool_t>(){ for (auto& it:valVVbools){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_short_t>(){ for (auto& it:valshorts){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vshort_t>(){ for (auto& it:valVshorts){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvshort_t>(){ for (auto& it:valVVshorts){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_uint_t>(){ for (auto& it:valuints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vuint_t>(){ for (auto& it:valVuints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvuint_t>(){ for (auto& it:valVVuints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_int_t>(){ for (auto& it:valints){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vint_t>(){ for (auto& it:valVints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvint_t>(){ for (auto& it:valVVints){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_ulong_t>(){ for (auto& it:valulongs){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vulong_t>(){ for (auto& it:valVulongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvulong_t>(){ for (auto& it:valVVulongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_long_t>(){ for (auto& it:vallongs){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vlong_t>(){ for (auto& it:valVlongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvlong_t>(){ for (auto& it:valVVlongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_ulonglong_t>(){ for (auto& it:valulonglongs){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vulonglong_t>(){ for (auto& it:valVulonglongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvulonglong_t>(){ for (auto& it:valVVulonglongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_longlong_t>(){ for (auto& it:vallonglongs){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vlonglong_t>(){ for (auto& it:valVlonglongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvlonglong_t>(){ for (auto& it:valVVlonglongs){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_float_t>(){ for (auto& it:valfloats){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vfloat_t>(){ for (auto& it:valVfloats){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvfloat_t>(){ for (auto& it:valVVfloats){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_double_t>(){ for (auto& it:valdoubles){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vdouble_t>(){ for (auto& it:valVdoubles){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vvdouble_t>(){ for (auto& it:valVVdoubles){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_string_t>(){ for (auto& it:valstrings){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vstring_t>(){ for (auto& it:valVstrings){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_TString_t>(){ for (auto& it:valTStrings){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vTString_t>(){ for (auto& it:valVTStrings){ if (it.second) it.second->clear(); } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_CMSLorentzVector_t>(){ for (auto& it:valCMSLorentzVectors){ if (it.second){ it.second->first=it.second->second; } } }
template<> void BaseTree::resetBranch<BaseTree::BranchType_vCMSLorentzVector_t>(){ for (auto& it:valVCMSLorentzVectors){ if (it.second) it.second->clear(); } }


template<> void BaseTree::removeBranch<BaseTree::BranchType_bool_t>(TString branchname){ for (auto& it:valbools){ if (it.first==branchname){ delete it.second; it.second=0; } } valbools.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vbool_t>(TString branchname){ for (auto& it:valVbools){ if (it.first==branchname){ delete it.second; it.second=0; } } valVbools.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvbool_t>(TString branchname){ for (auto& it:valVVbools){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVbools.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_short_t>(TString branchname){ for (auto& it:valshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } valshorts.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vshort_t>(TString branchname){ for (auto& it:valVshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } valVshorts.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvshort_t>(TString branchname){ for (auto& it:valVVshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVshorts.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_uint_t>(TString branchname){ for (auto& it:valuints){ if (it.first==branchname){ delete it.second; it.second=0; } } valuints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vuint_t>(TString branchname){ for (auto& it:valVuints){ if (it.first==branchname){ delete it.second; it.second=0; } } valVuints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvuint_t>(TString branchname){ for (auto& it:valVVuints){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVuints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_int_t>(TString branchname){ for (auto& it:valints){ if (it.first==branchname){ delete it.second; it.second=0; } } valints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vint_t>(TString branchname){ for (auto& it:valVints){ if (it.first==branchname){ delete it.second; it.second=0; } } valVints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvint_t>(TString branchname){ for (auto& it:valVVints){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVints.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_ulong_t>(TString branchname){ for (auto& it:valulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valulongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vulong_t>(TString branchname){ for (auto& it:valVulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVulongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvulong_t>(TString branchname){ for (auto& it:valVVulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVulongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_long_t>(TString branchname){ for (auto& it:vallongs){ if (it.first==branchname){ delete it.second; it.second=0; } } vallongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vlong_t>(TString branchname){ for (auto& it:valVlongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVlongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvlong_t>(TString branchname){ for (auto& it:valVVlongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVlongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_ulonglong_t>(TString branchname){ for (auto& it:valulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valulonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vulonglong_t>(TString branchname){ for (auto& it:valVulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVulonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvulonglong_t>(TString branchname){ for (auto& it:valVVulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVulonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_longlong_t>(TString branchname){ for (auto& it:vallonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } vallonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vlonglong_t>(TString branchname){ for (auto& it:valVlonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVlonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvlonglong_t>(TString branchname){ for (auto& it:valVVlonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVlonglongs.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_float_t>(TString branchname){ for (auto& it:valfloats){ if (it.first==branchname){ delete it.second; it.second=0; } } valfloats.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ for (auto& it:valVfloats){ if (it.first==branchname){ delete it.second; it.second=0; } } valVfloats.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvfloat_t>(TString branchname){ for (auto& it:valVVfloats){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVfloats.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_double_t>(TString branchname){ for (auto& it:valdoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } valdoubles.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ for (auto& it:valVdoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } valVdoubles.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vvdouble_t>(TString branchname){ for (auto& it:valVVdoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } valVVdoubles.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_string_t>(TString branchname){ for (auto& it:valstrings){ if (it.first==branchname){ delete it.second; it.second=0; } } valstrings.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vstring_t>(TString branchname){ for (auto& it:valVstrings){ if (it.first==branchname){ delete it.second; it.second=0; } } valVstrings.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_TString_t>(TString branchname){ for (auto& it:valTStrings){ if (it.first==branchname){ delete it.second; it.second=0; } } valTStrings.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vTString_t>(TString branchname){ for (auto& it:valVTStrings){ if (it.first==branchname){ delete it.second; it.second=0; } } valVTStrings.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_CMSLorentzVector_t>(TString branchname){ for (auto& it:valCMSLorentzVectors){ if (it.first==branchname){ delete it.second; it.second=0; } } valCMSLorentzVectors.erase(branchname); }
template<> void BaseTree::removeBranch<BaseTree::BranchType_vCMSLorentzVector_t>(TString branchname){ for (auto& it:valVCMSLorentzVectors){ if (it.first==branchname){ delete it.second; it.second=0; } } valVCMSLorentzVectors.erase(branchname); }


template<> bool BaseTree::bookBranch<bool>(TString branchname, bool valdef){
  if (valbools.find(branchname)==valbools.end()) valbools[branchname] = new std::pair<bool, bool>(valdef, valdef);
  else{ valbools[branchname]->first=valdef; valbools[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valbools[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valbools[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<bool>*>(TString branchname, std::vector<bool>*/* valdef*/){
  valVbools[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVbools[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVbools[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<bool>>*>(TString branchname, std::vector<std::vector<bool>>*/* valdef*/){
  valVVbools[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVbools[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVbools[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<short>(TString branchname, short valdef){
  if (valshorts.find(branchname)==valshorts.end()) valshorts[branchname] = new std::pair<short, short>(valdef, valdef);
  else{ valshorts[branchname]->first=valdef; valshorts[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valshorts[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valshorts[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<short>*>(TString branchname, std::vector<short>*/* valdef*/){
  valVshorts[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVshorts[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVshorts[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<short>>*>(TString branchname, std::vector<std::vector<short>>*/* valdef*/){
  valVVshorts[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVshorts[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVshorts[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<unsigned int>(TString branchname, unsigned int valdef){
  if (valuints.find(branchname)==valuints.end()) valuints[branchname] = new std::pair<unsigned int, unsigned int>(valdef, valdef);
  else{ valuints[branchname]->first=valdef; valuints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valuints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valuints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>*/* valdef*/){
  valVuints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVuints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVuints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<unsigned int>>*>(TString branchname, std::vector<std::vector<unsigned int>>*/* valdef*/){
  valVVuints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVuints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVuints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<int>(TString branchname, int valdef){
  if (valints.find(branchname)==valints.end()) valints[branchname] = new std::pair<int, int>(valdef, valdef);
  else{ valints[branchname]->first=valdef; valints[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valints[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valints[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<int>*>(TString branchname, std::vector<int>*/* valdef*/){
  valVints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<int>>*>(TString branchname, std::vector<std::vector<int>>*/* valdef*/){
  valVVints[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVints[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVints[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<unsigned long>(TString branchname, unsigned long valdef){
  if (valulongs.find(branchname)==valulongs.end()) valulongs[branchname] = new std::pair<unsigned long, unsigned long>(valdef, valdef);
  else{ valulongs[branchname]->first=valdef; valulongs[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valulongs[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valulongs[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<unsigned long>*>(TString branchname, std::vector<unsigned long>*/* valdef*/){
  valVulongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVulongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVulongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<unsigned long>>*>(TString branchname, std::vector<std::vector<unsigned long>>*/* valdef*/){
  valVVulongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVulongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVulongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<long>(TString branchname, long valdef){
  if (vallongs.find(branchname)==vallongs.end()) vallongs[branchname] = new std::pair<long, long>(valdef, valdef);
  else{ vallongs[branchname]->first=valdef; vallongs[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(vallongs[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(vallongs[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<long>*>(TString branchname, std::vector<long>*/* valdef*/){
  valVlongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVlongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVlongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<long>>*>(TString branchname, std::vector<std::vector<long>>*/* valdef*/){
  valVVlongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVlongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVlongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<unsigned long long>(TString branchname, unsigned long long valdef){
  if (valulonglongs.find(branchname)==valulonglongs.end()) valulonglongs[branchname] = new std::pair<unsigned long long, unsigned long long>(valdef, valdef);
  else{ valulonglongs[branchname]->first=valdef; valulonglongs[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valulonglongs[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valulonglongs[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<unsigned long long>*>(TString branchname, std::vector<unsigned long long>*/* valdef*/){
  valVulonglongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVulonglongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVulonglongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::vector<std::vector<unsigned long long>>*/* valdef*/){
  valVVulonglongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVulonglongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVulonglongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<long long>(TString branchname, long long valdef){
  if (vallonglongs.find(branchname)==vallonglongs.end()) vallonglongs[branchname] = new std::pair<long long, long long>(valdef, valdef);
  else{ vallonglongs[branchname]->first=valdef; vallonglongs[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(vallonglongs[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(vallonglongs[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<long long>*>(TString branchname, std::vector<long long>*/* valdef*/){
  valVlonglongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVlonglongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVlonglongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<long long>>*>(TString branchname, std::vector<std::vector<long long>>*/* valdef*/){
  valVVlonglongs[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVlonglongs[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVlonglongs[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<float>(TString branchname, float valdef){
  if (valfloats.find(branchname)==valfloats.end()) valfloats[branchname] = new std::pair<float, float>(valdef, valdef);
  else{ valfloats[branchname]->first=valdef; valfloats[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valfloats[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valfloats[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<float>*>(TString branchname, std::vector<float>*/* valdef*/){
  valVfloats[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVfloats[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVfloats[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<float>>*>(TString branchname, std::vector<std::vector<float>>*/* valdef*/){
  valVVfloats[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVfloats[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVfloats[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<double>(TString branchname, double valdef){
  if (valdoubles.find(branchname)==valdoubles.end()) valdoubles[branchname] = new std::pair<double, double>(valdef, valdef);
  else{ valdoubles[branchname]->first=valdef; valdoubles[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valdoubles[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valdoubles[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<double>*>(TString branchname, std::vector<double>*/* valdef*/){
  valVdoubles[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVdoubles[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVdoubles[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::vector<double>>*>(TString branchname, std::vector<std::vector<double>>*/* valdef*/){
  valVVdoubles[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVVdoubles[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVVdoubles[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<std::string>(TString branchname, std::string valdef){
  if (valstrings.find(branchname)==valstrings.end()) valstrings[branchname] = new std::pair<std::string, std::string>(valdef, valdef);
  else{ valstrings[branchname]->first=valdef; valstrings[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valstrings[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valstrings[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<std::string>*>(TString branchname, std::vector<std::string>*/* valdef*/){
  valVstrings[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVstrings[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVstrings[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<TString>(TString branchname, TString valdef){
  if (valTStrings.find(branchname)==valTStrings.end()) valTStrings[branchname] = new std::pair<TString, TString>(valdef, valdef);
  else{ valTStrings[branchname]->first=valdef; valTStrings[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valTStrings[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valTStrings[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<TString>*>(TString branchname, std::vector<TString>*/* valdef*/){
  valVTStrings[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVTStrings[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVTStrings[branchname]));
  return true;
}
template<> bool BaseTree::bookBranch<CMSLorentzVector>(TString branchname, CMSLorentzVector valdef){
  if (valCMSLorentzVectors.find(branchname)==valCMSLorentzVectors.end()) valCMSLorentzVectors[branchname] = new std::pair<CMSLorentzVector, CMSLorentzVector>(valdef, valdef);
  else{ valCMSLorentzVectors[branchname]->first=valdef; valCMSLorentzVectors[branchname]->second=valdef; }
  SampleHelpers::bookBranch(tree, branchname, &(valCMSLorentzVectors[branchname]->first));
  SampleHelpers::bookBranch(failedtree, branchname, &(valCMSLorentzVectors[branchname]->first));
  return true;
}
template<> bool BaseTree::bookBranch<std::vector<CMSLorentzVector>*>(TString branchname, std::vector<CMSLorentzVector>*/* valdef*/){
  valVCMSLorentzVectors[branchname] = nullptr;
  SampleHelpers::bookBranch(tree, branchname, &(valVCMSLorentzVectors[branchname]));
  SampleHelpers::bookBranch(failedtree, branchname, &(valVCMSLorentzVectors[branchname]));
  return true;
}


template<> bool BaseTree::bookBranch<BaseTree::BranchType_bool_t>(TString branchname){ return this->bookBranch<bool>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vbool_t>(TString branchname){ return this->bookBranch<std::vector<bool>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvbool_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<bool>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_short_t>(TString branchname){ return this->bookBranch<short>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vshort_t>(TString branchname){ return this->bookBranch<std::vector<short>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvshort_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<short>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_uint_t>(TString branchname){ return this->bookBranch<unsigned int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vuint_t>(TString branchname){ return this->bookBranch<std::vector<unsigned int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvuint_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<unsigned int>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_int_t>(TString branchname){ return this->bookBranch<int>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vint_t>(TString branchname){ return this->bookBranch<std::vector<int>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvint_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<int>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_ulong_t>(TString branchname){ return this->bookBranch<unsigned long>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vulong_t>(TString branchname){ return this->bookBranch<std::vector<unsigned long>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvulong_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<unsigned long>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_long_t>(TString branchname){ return this->bookBranch<long>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vlong_t>(TString branchname){ return this->bookBranch<std::vector<long>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvlong_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<long>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_ulonglong_t>(TString branchname){ return this->bookBranch<unsigned long long>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vulonglong_t>(TString branchname){ return this->bookBranch<std::vector<unsigned long long>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvulonglong_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<unsigned long long>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_longlong_t>(TString branchname){ return this->bookBranch<long long>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vlonglong_t>(TString branchname){ return this->bookBranch<std::vector<long long>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvlonglong_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<long long>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_float_t>(TString branchname){ return this->bookBranch<float>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ return this->bookBranch<std::vector<float>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvfloat_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<float>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_double_t>(TString branchname){ return this->bookBranch<double>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ return this->bookBranch<std::vector<double>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vvdouble_t>(TString branchname){ return this->bookBranch<std::vector<std::vector<double>>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_string_t>(TString branchname){ return this->bookBranch<std::string>(branchname, ""); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vstring_t>(TString branchname){ return this->bookBranch<std::vector<std::string>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_TString_t>(TString branchname){ return this->bookBranch<TString>(branchname, ""); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vTString_t>(TString branchname){ return this->bookBranch<std::vector<TString>*>(branchname, 0); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_CMSLorentzVector_t>(TString branchname){ return this->bookBranch<CMSLorentzVector>(branchname, CMSLorentzVector(0, 0, 0, 0)); }
template<> bool BaseTree::bookBranch<BaseTree::BranchType_vCMSLorentzVector_t>(TString branchname){ return this->bookBranch<std::vector<CMSLorentzVector>*>(branchname, 0); }


template<> bool BaseTree::putBranch<bool>(TString branchname, bool valdef){
  if (valbools.find(branchname)==valbools.end()) valbools[branchname] = new std::pair<bool, bool>(valdef, valdef);
  else{ valbools[branchname]->first=valdef; valbools[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valbools[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<bool>*>(TString branchname, std::vector<bool>*/* valdef*/){
  if (valVbools.find(branchname)==valVbools.end()) valVbools[branchname] = new std::vector<bool>();
  else valVbools[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVbools[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<bool>>*>(TString branchname, std::vector<std::vector<bool>>*/* valdef*/){
  if (valVVbools.find(branchname)==valVVbools.end()) valVVbools[branchname] = new std::vector<std::vector<bool>>();
  else valVVbools[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVbools[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<short>(TString branchname, short valdef){
  if (valshorts.find(branchname)==valshorts.end()) valshorts[branchname] = new std::pair<short, short>(valdef, valdef);
  else{ valshorts[branchname]->first=valdef; valshorts[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valshorts[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<short>*>(TString branchname, std::vector<short>*/* valdef*/){
  if (valVshorts.find(branchname)==valVshorts.end()) valVshorts[branchname] = new std::vector<short>();
  else valVshorts[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVshorts[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<short>>*>(TString branchname, std::vector<std::vector<short>>*/* valdef*/){
  if (valVVshorts.find(branchname)==valVVshorts.end()) valVVshorts[branchname] = new std::vector<std::vector<short>>();
  else valVVshorts[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVshorts[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<unsigned int>(TString branchname, unsigned int valdef){
  if (valuints.find(branchname)==valuints.end()) valuints[branchname] = new std::pair<unsigned int, unsigned int>(valdef, valdef);
  else{ valuints[branchname]->first=valdef; valuints[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valuints[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>*/* valdef*/){
  if (valVuints.find(branchname)==valVuints.end()) valVuints[branchname] = new std::vector<unsigned int>();
  else valVuints[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVuints[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<unsigned int>>*>(TString branchname, std::vector<std::vector<unsigned int>>*/* valdef*/){
  if (valVVuints.find(branchname)==valVVuints.end()) valVVuints[branchname] = new std::vector<std::vector<unsigned int>>();
  else valVVuints[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVuints[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<int>(TString branchname, int valdef){
  if (valints.find(branchname)==valints.end()) valints[branchname] = new std::pair<int, int>(valdef, valdef);
  else{ valints[branchname]->first=valdef; valints[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valints[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<int>*>(TString branchname, std::vector<int>*/* valdef*/){
  if (valVints.find(branchname)==valVints.end()) valVints[branchname] = new std::vector<int>();
  else valVints[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVints[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<int>>*>(TString branchname, std::vector<std::vector<int>>*/* valdef*/){
  if (valVVints.find(branchname)==valVVints.end()) valVVints[branchname] = new std::vector<std::vector<int>>();
  else valVVints[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVints[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<unsigned long>(TString branchname, unsigned long valdef){
  if (valulongs.find(branchname)==valulongs.end()) valulongs[branchname] = new std::pair<unsigned long, unsigned long>(valdef, valdef);
  else{ valulongs[branchname]->first=valdef; valulongs[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valulongs[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<unsigned long>*>(TString branchname, std::vector<unsigned long>*/* valdef*/){
  if (valVulongs.find(branchname)==valVulongs.end()) valVulongs[branchname] = new std::vector<unsigned long>();
  else valVulongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVulongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<unsigned long>>*>(TString branchname, std::vector<std::vector<unsigned long>>*/* valdef*/){
  if (valVVulongs.find(branchname)==valVVulongs.end()) valVVulongs[branchname] = new std::vector<std::vector<unsigned long>>();
  else valVVulongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVulongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<long>(TString branchname, long valdef){
  if (vallongs.find(branchname)==vallongs.end()) vallongs[branchname] = new std::pair<long, long>(valdef, valdef);
  else{ vallongs[branchname]->first=valdef; vallongs[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, vallongs[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<long>*>(TString branchname, std::vector<long>*/* valdef*/){
  if (valVlongs.find(branchname)==valVlongs.end()) valVlongs[branchname] = new std::vector<long>();
  else valVlongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVlongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<long>>*>(TString branchname, std::vector<std::vector<long>>*/* valdef*/){
  if (valVVlongs.find(branchname)==valVVlongs.end()) valVVlongs[branchname] = new std::vector<std::vector<long>>();
  else valVVlongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVlongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<unsigned long long>(TString branchname, unsigned long long valdef){
  if (valulonglongs.find(branchname)==valulonglongs.end()) valulonglongs[branchname] = new std::pair<unsigned long long, unsigned long long>(valdef, valdef);
  else{ valulonglongs[branchname]->first=valdef; valulonglongs[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valulonglongs[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<unsigned long long>*>(TString branchname, std::vector<unsigned long long>*/* valdef*/){
  if (valVulonglongs.find(branchname)==valVulonglongs.end()) valVulonglongs[branchname] = new std::vector<unsigned long long>();
  else valVulonglongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVulonglongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::vector<std::vector<unsigned long long>>*/* valdef*/){
  if (valVVulonglongs.find(branchname)==valVVulonglongs.end()) valVVulonglongs[branchname] = new std::vector<std::vector<unsigned long long>>();
  else valVVulonglongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVulonglongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<long long>(TString branchname, long long valdef){
  if (vallonglongs.find(branchname)==vallonglongs.end()) vallonglongs[branchname] = new std::pair<long long, long long>(valdef, valdef);
  else{ vallonglongs[branchname]->first=valdef; vallonglongs[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, vallonglongs[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<long long>*>(TString branchname, std::vector<long long>*/* valdef*/){
  if (valVlonglongs.find(branchname)==valVlonglongs.end()) valVlonglongs[branchname] = new std::vector<long long>();
  else valVlonglongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVlonglongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<long long>>*>(TString branchname, std::vector<std::vector<long long>>*/* valdef*/){
  if (valVVlonglongs.find(branchname)==valVVlonglongs.end()) valVVlonglongs[branchname] = new std::vector<std::vector<long long>>();
  else valVVlonglongs[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVlonglongs[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<float>(TString branchname, float valdef){
  if (valfloats.find(branchname)==valfloats.end()) valfloats[branchname] = new std::pair<float, float>(valdef, valdef);
  else{ valfloats[branchname]->first=valdef; valfloats[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valfloats[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<float>*>(TString branchname, std::vector<float>*/* valdef*/){
  if (valVfloats.find(branchname)==valVfloats.end()) valVfloats[branchname] = new std::vector<float>();
  else valVfloats[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVfloats[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<float>>*>(TString branchname, std::vector<std::vector<float>>*/* valdef*/){
  if (valVVfloats.find(branchname)==valVVfloats.end()) valVVfloats[branchname] = new std::vector<std::vector<float>>();
  else valVVfloats[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVfloats[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<double>(TString branchname, double valdef){
  if (valdoubles.find(branchname)==valdoubles.end()) valdoubles[branchname] = new std::pair<double, double>(valdef, valdef);
  else{ valdoubles[branchname]->first=valdef; valdoubles[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valdoubles[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<double>*>(TString branchname, std::vector<double>*/* valdef*/){
  if (valVdoubles.find(branchname)==valVdoubles.end()) valVdoubles[branchname] = new std::vector<double>();
  else valVdoubles[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVdoubles[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::vector<double>>*>(TString branchname, std::vector<std::vector<double>>*/* valdef*/){
  if (valVVdoubles.find(branchname)==valVVdoubles.end()) valVVdoubles[branchname] = new std::vector<std::vector<double>>();
  else valVVdoubles[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVVdoubles[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<std::string>(TString branchname, std::string valdef){
  if (valstrings.find(branchname)==valstrings.end()) valstrings[branchname] = new std::pair<std::string, std::string>(valdef, valdef);
  else{ valstrings[branchname]->first=valdef; valstrings[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valstrings[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<std::string>*>(TString branchname, std::vector<std::string>*/* valdef*/){
  if (valVstrings.find(branchname)==valVstrings.end()) valVstrings[branchname] = new std::vector<std::string>();
  else valVstrings[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVstrings[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<TString>(TString branchname, TString valdef){
  if (valTStrings.find(branchname)==valTStrings.end()) valTStrings[branchname] = new std::pair<TString, TString>(valdef, valdef);
  else{ valTStrings[branchname]->first=valdef; valTStrings[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valTStrings[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<TString>*>(TString branchname, std::vector<TString>*/* valdef*/){
  if (valVTStrings.find(branchname)==valVTStrings.end()) valVTStrings[branchname] = new std::vector<TString>();
  else valVTStrings[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVTStrings[branchname]));
  return true;
}
template<> bool BaseTree::putBranch<CMSLorentzVector>(TString branchname, CMSLorentzVector valdef){
  if (valCMSLorentzVectors.find(branchname)==valCMSLorentzVectors.end()) valCMSLorentzVectors[branchname] = new std::pair<CMSLorentzVector, CMSLorentzVector>(valdef, valdef);
  else{ valCMSLorentzVectors[branchname]->first=valdef; valCMSLorentzVectors[branchname]->second=valdef; }
  SampleHelpers::putBranch(tree, branchname, valCMSLorentzVectors[branchname]->first);
  return true;
}
template<> bool BaseTree::putBranch<std::vector<CMSLorentzVector>*>(TString branchname, std::vector<CMSLorentzVector>*/* valdef*/){
  if (valVCMSLorentzVectors.find(branchname)==valVCMSLorentzVectors.end()) valVCMSLorentzVectors[branchname] = new std::vector<CMSLorentzVector>();
  else valVCMSLorentzVectors[branchname]->clear();
  SampleHelpers::putBranch(tree, branchname, *(valVCMSLorentzVectors[branchname]));
  return true;
}


template<> bool BaseTree::putBranch<BaseTree::BranchType_bool_t>(TString branchname){ return this->putBranch<bool>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vbool_t>(TString branchname){ return this->putBranch<std::vector<bool>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvbool_t>(TString branchname){ return this->putBranch<std::vector<std::vector<bool>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_short_t>(TString branchname){ return this->putBranch<short>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vshort_t>(TString branchname){ return this->putBranch<std::vector<short>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvshort_t>(TString branchname){ return this->putBranch<std::vector<std::vector<short>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_uint_t>(TString branchname){ return this->putBranch<unsigned int>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vuint_t>(TString branchname){ return this->putBranch<std::vector<unsigned int>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvuint_t>(TString branchname){ return this->putBranch<std::vector<std::vector<unsigned int>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_int_t>(TString branchname){ return this->putBranch<int>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vint_t>(TString branchname){ return this->putBranch<std::vector<int>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvint_t>(TString branchname){ return this->putBranch<std::vector<std::vector<int>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_ulong_t>(TString branchname){ return this->putBranch<unsigned long>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vulong_t>(TString branchname){ return this->putBranch<std::vector<unsigned long>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvulong_t>(TString branchname){ return this->putBranch<std::vector<std::vector<unsigned long>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_long_t>(TString branchname){ return this->putBranch<long>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vlong_t>(TString branchname){ return this->putBranch<std::vector<long>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvlong_t>(TString branchname){ return this->putBranch<std::vector<std::vector<long>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_ulonglong_t>(TString branchname){ return this->putBranch<unsigned long long>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vulonglong_t>(TString branchname){ return this->putBranch<std::vector<unsigned long long>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvulonglong_t>(TString branchname){ return this->putBranch<std::vector<std::vector<unsigned long long>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_longlong_t>(TString branchname){ return this->putBranch<long long>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vlonglong_t>(TString branchname){ return this->putBranch<std::vector<long long>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvlonglong_t>(TString branchname){ return this->putBranch<std::vector<std::vector<long long>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_float_t>(TString branchname){ return this->putBranch<float>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ return this->putBranch<std::vector<float>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvfloat_t>(TString branchname){ return this->putBranch<std::vector<std::vector<float>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_double_t>(TString branchname){ return this->putBranch<double>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ return this->putBranch<std::vector<double>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vvdouble_t>(TString branchname){ return this->putBranch<std::vector<std::vector<double>>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_string_t>(TString branchname){ return this->putBranch<std::string>(branchname, ""); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vstring_t>(TString branchname){ return this->putBranch<std::vector<std::string>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_TString_t>(TString branchname){ return this->putBranch<TString>(branchname, ""); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vTString_t>(TString branchname){ return this->putBranch<std::vector<TString>*>(branchname, 0); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_CMSLorentzVector_t>(TString branchname){ return this->putBranch<CMSLorentzVector>(branchname, CMSLorentzVector(0, 0, 0, 0)); }
template<> bool BaseTree::putBranch<BaseTree::BranchType_vCMSLorentzVector_t>(TString branchname){ return this->putBranch<std::vector<CMSLorentzVector>*>(branchname, 0); }


template<> void BaseTree::getVal<bool>(TString branchname, bool& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<bool> const*>(TString branchname, std::vector<bool> const*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::vector<std::vector<bool>> const*>(TString branchname, std::vector<std::vector<bool>> const*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
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
template<> void BaseTree::getVal<std::vector<std::vector<short>> const*>(TString branchname, std::vector<std::vector<short>> const*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
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
template<> void BaseTree::getVal<std::vector<std::vector<unsigned int>> const*>(TString branchname, std::vector<std::vector<unsigned int>> const*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
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
template<> void BaseTree::getVal<std::vector<std::vector<int>> const*>(TString branchname, std::vector<std::vector<int>> const*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<unsigned long>(TString branchname, unsigned long& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<unsigned long> const*>(TString branchname, std::vector<unsigned long> const*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::vector<std::vector<unsigned long>> const*>(TString branchname, std::vector<std::vector<unsigned long>> const*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<long>(TString branchname, long& val) const{
  typedef long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<long> const*>(TString branchname, std::vector<long> const*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::vector<std::vector<long>> const*>(TString branchname, std::vector<std::vector<long>> const*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<unsigned long long>(TString branchname, unsigned long long& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<unsigned long long> const*>(TString branchname, std::vector<unsigned long long> const*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::vector<std::vector<unsigned long long>> const*>(TString branchname, std::vector<std::vector<unsigned long long>> const*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<long long>(TString branchname, long long& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<long long> const*>(TString branchname, std::vector<long long> const*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::vector<std::vector<long long>> const*>(TString branchname, std::vector<std::vector<long long>> const*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
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
template<> void BaseTree::getVal<std::vector<std::vector<float>> const*>(TString branchname, std::vector<std::vector<float>> const*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
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
template<> void BaseTree::getVal<std::vector<std::vector<double>> const*>(TString branchname, std::vector<std::vector<double>> const*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<std::string>(TString branchname, std::string& val) const{
  typedef std::string itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<std::string> const*>(TString branchname, std::vector<std::string> const*& val) const{
  typedef std::string itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<TString>(TString branchname, TString& val) const{
  typedef TString itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<TString> const*>(TString branchname, std::vector<TString> const*& val) const{
  typedef TString itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getVal<CMSLorentzVector>(TString branchname, CMSLorentzVector& val) const{
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=tmp->first; }
}
template<> void BaseTree::getVal<std::vector<CMSLorentzVector> const*>(TString branchname, std::vector<CMSLorentzVector> const*& val) const{
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}


template<> void BaseTree::setVal<bool>(TString branchname, bool const& val){
  typedef bool itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<bool>*>(TString branchname, std::vector<bool>* const& val){
  typedef bool itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::vector<std::vector<bool>>*>(TString branchname, std::vector<std::vector<bool>>* const& val){
  typedef bool itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
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
template<> void BaseTree::setVal<std::vector<std::vector<short>>*>(TString branchname, std::vector<std::vector<short>>* const& val){
  typedef short itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
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
template<> void BaseTree::setVal<std::vector<std::vector<unsigned int>>*>(TString branchname, std::vector<std::vector<unsigned int>>* const& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
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
template<> void BaseTree::setVal<std::vector<std::vector<int>>*>(TString branchname, std::vector<std::vector<int>>* const& val){
  typedef int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<unsigned long>(TString branchname, unsigned long const& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<unsigned long>*>(TString branchname, std::vector<unsigned long>* const& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::vector<std::vector<unsigned long>>*>(TString branchname, std::vector<std::vector<unsigned long>>* const& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<long>(TString branchname, long const& val){
  typedef long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<long>*>(TString branchname, std::vector<long>* const& val){
  typedef long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::vector<std::vector<long>>*>(TString branchname, std::vector<std::vector<long>>* const& val){
  typedef long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<unsigned long long>(TString branchname, unsigned long long const& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<unsigned long long>*>(TString branchname, std::vector<unsigned long long>* const& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::vector<std::vector<unsigned long long>>* const& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<long long>(TString branchname, long long const& val){
  typedef long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<long long>*>(TString branchname, std::vector<long long>* const& val){
  typedef long long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::vector<std::vector<long long>>*>(TString branchname, std::vector<std::vector<long long>>* const& val){
  typedef long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
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
template<> void BaseTree::setVal<std::vector<std::vector<float>>*>(TString branchname, std::vector<std::vector<float>>* const& val){
  typedef float itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
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
template<> void BaseTree::setVal<std::vector<std::vector<double>>*>(TString branchname, std::vector<std::vector<double>>* const& val){
  typedef double itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<std::string>(TString branchname, std::string const& val){
  typedef std::string itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<std::string>*>(TString branchname, std::vector<std::string>* const& val){
  typedef std::string itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<TString>(TString branchname, TString const& val){
  typedef TString itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<TString>*>(TString branchname, std::vector<TString>* const& val){
  typedef TString itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}
template<> void BaseTree::setVal<CMSLorentzVector>(TString branchname, CMSLorentzVector const& val){
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) tmp->first=val; }
}
template<> void BaseTree::setVal<std::vector<CMSLorentzVector>*>(TString branchname, std::vector<CMSLorentzVector>* const& val){
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it) && it->second && val) it->second->assign(val->begin(), val->end());
}


template<> void BaseTree::getValRef<bool>(TString branchname, bool*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<short>(TString branchname, short*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned int>(TString branchname, unsigned int*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<int>(TString branchname, int*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned long>(TString branchname, unsigned long*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<long>(TString branchname, long*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned long long>(TString branchname, unsigned long long*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<long long>(TString branchname, long long*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<float>(TString branchname, float*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<double>(TString branchname, double*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::string>(TString branchname, std::string*& val) const{
  typedef std::string itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<TString>(TString branchname, TString*& val) const{
  typedef TString itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<CMSLorentzVector>(TString branchname, CMSLorentzVector*& val) const{
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}

template<> void BaseTree::getValRef<bool>(TString branchname, bool*& val){
  typedef bool itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<short>(TString branchname, short*& val){
  typedef short itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned int>(TString branchname, unsigned int*& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<int>(TString branchname, int*& val){
  typedef int itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned long>(TString branchname, unsigned long*& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<long>(TString branchname, long*& val){
  typedef long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<unsigned long long>(TString branchname, unsigned long long*& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<long long>(TString branchname, long long*& val){
  typedef long long itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<float>(TString branchname, float*& val){
  typedef float itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<double>(TString branchname, double*& val){
  typedef double itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<std::string>(TString branchname, std::string*& val){
  typedef std::string itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<TString>(TString branchname, TString*& val){
  typedef TString itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}
template<> void BaseTree::getValRef<CMSLorentzVector>(TString branchname, CMSLorentzVector*& val){
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::pair<itType, itType>*>::iterator it;
  if (this->getBranchCIterator<std::pair<itType, itType>*>(branchname, it)){ auto& tmp = it->second; if (tmp) val=&(tmp->first); }
}

template<> void BaseTree::getValRef<std::vector<bool>>(TString branchname, std::vector<bool>*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<short>>(TString branchname, std::vector<short>*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<unsigned int>>(TString branchname, std::vector<unsigned int>*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<int>>(TString branchname, std::vector<int>*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<unsigned long>>(TString branchname, std::vector<unsigned long>*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<long>>(TString branchname, std::vector<long>*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<unsigned long long>>(TString branchname, std::vector<unsigned long long>*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<long long>>(TString branchname, std::vector<long long>*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<float>>(TString branchname, std::vector<float>*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<double>>(TString branchname, std::vector<double>*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::string>>(TString branchname, std::vector<std::string>*& val) const{
  typedef std::string itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<TString>>(TString branchname, std::vector<TString>*& val) const{
  typedef TString itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<CMSLorentzVector>>(TString branchname, std::vector<CMSLorentzVector>*& val) const{
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<bool>>>(TString branchname, std::vector<std::vector<bool>>*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<short>>>(TString branchname, std::vector<std::vector<short>>*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned int>>>(TString branchname, std::vector<std::vector<unsigned int>>*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<int>>>(TString branchname, std::vector<std::vector<int>>*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long>>>(TString branchname, std::vector<std::vector<unsigned long>>*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<long>>>(TString branchname, std::vector<std::vector<long>>*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long long>>>(TString branchname, std::vector<std::vector<unsigned long long>>*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<long long>>>(TString branchname, std::vector<std::vector<long long>>*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<float>>>(TString branchname, std::vector<std::vector<float>>*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}
template<> void BaseTree::getValRef<std::vector<std::vector<double>>>(TString branchname, std::vector<std::vector<double>>*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = it->second;
}

template<> void BaseTree::getValRef<std::vector<bool>* const>(TString branchname, std::vector<bool>* const*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<short>* const>(TString branchname, std::vector<short>* const*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned int>* const>(TString branchname, std::vector<unsigned int>* const*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<int>* const>(TString branchname, std::vector<int>* const*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned long>* const>(TString branchname, std::vector<unsigned long>* const*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<long>* const>(TString branchname, std::vector<long>* const*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned long long>* const>(TString branchname, std::vector<unsigned long long>* const*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<long long>* const>(TString branchname, std::vector<long long>* const*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<float>* const>(TString branchname, std::vector<float>* const*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<double>* const>(TString branchname, std::vector<double>* const*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::string>* const>(TString branchname, std::vector<std::string>* const*& val) const{
  typedef std::string itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<TString>* const>(TString branchname, std::vector<TString>* const*& val) const{
  typedef TString itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<CMSLorentzVector>* const>(TString branchname, std::vector<CMSLorentzVector>* const*& val) const{
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::vector<itType>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<bool>>* const>(TString branchname, std::vector<std::vector<bool>>* const*& val) const{
  typedef bool itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<short>>* const>(TString branchname, std::vector<std::vector<short>>* const*& val) const{
  typedef short itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned int>>* const>(TString branchname, std::vector<std::vector<unsigned int>>* const*& val) const{
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<int>>* const>(TString branchname, std::vector<std::vector<int>>* const*& val) const{
  typedef int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long>>* const>(TString branchname, std::vector<std::vector<unsigned long>>* const*& val) const{
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<long>>* const>(TString branchname, std::vector<std::vector<long>>* const*& val) const{
  typedef long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long long>>* const>(TString branchname, std::vector<std::vector<unsigned long long>>* const*& val) const{
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<long long>>* const>(TString branchname, std::vector<std::vector<long long>>* const*& val) const{
  typedef long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<float>>* const>(TString branchname, std::vector<std::vector<float>>* const*& val) const{
  typedef float itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<double>>* const>(TString branchname, std::vector<std::vector<double>>* const*& val) const{
  typedef double itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::const_iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}

template<> void BaseTree::getValRef<std::vector<bool>*>(TString branchname, std::vector<bool>**& val){
  typedef bool itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<short>*>(TString branchname, std::vector<short>**& val){
  typedef short itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>**& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<int>*>(TString branchname, std::vector<int>**& val){
  typedef int itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned long>*>(TString branchname, std::vector<unsigned long>**& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<long>*>(TString branchname, std::vector<long>**& val){
  typedef long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<unsigned long long>*>(TString branchname, std::vector<unsigned long long>**& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<long long>*>(TString branchname, std::vector<long long>**& val){
  typedef long long itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<float>*>(TString branchname, std::vector<float>**& val){
  typedef float itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<double>*>(TString branchname, std::vector<double>**& val){
  typedef double itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::string>*>(TString branchname, std::vector<std::string>**& val){
  typedef std::string itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<TString>*>(TString branchname, std::vector<TString>**& val){
  typedef TString itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<CMSLorentzVector>*>(TString branchname, std::vector<CMSLorentzVector>**& val){
  typedef CMSLorentzVector itType;
  std::unordered_map<TString, std::vector<itType>*>::iterator it;
  if (this->getBranchCIterator<std::vector<itType>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<bool>>*>(TString branchname, std::vector<std::vector<bool>>**& val){
  typedef bool itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<short>>*>(TString branchname, std::vector<std::vector<short>>**& val){
  typedef short itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned int>>*>(TString branchname, std::vector<std::vector<unsigned int>>**& val){
  typedef unsigned int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<int>>*>(TString branchname, std::vector<std::vector<int>>**& val){
  typedef int itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long>>*>(TString branchname, std::vector<std::vector<unsigned long>>**& val){
  typedef unsigned long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<long>>*>(TString branchname, std::vector<std::vector<long>>**& val){
  typedef long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::vector<std::vector<unsigned long long>>**& val){
  typedef unsigned long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<long long>>*>(TString branchname, std::vector<std::vector<long long>>**& val){
  typedef long long itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<float>>*>(TString branchname, std::vector<std::vector<float>>**& val){
  typedef float itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}
template<> void BaseTree::getValRef<std::vector<std::vector<double>>*>(TString branchname, std::vector<std::vector<double>>**& val){
  typedef double itType;
  std::unordered_map<TString, std::vector<std::vector<itType>>*>::iterator it;
  if (this->getBranchCIterator<std::vector<std::vector<itType>>*>(branchname, it)) val = &(it->second);
}


#endif
