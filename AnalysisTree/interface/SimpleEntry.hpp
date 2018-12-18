#ifndef SIMPLEENTRY_HPP
#define SIMPLEENTRY_HPP

#include "SimpleEntry.h"


template<> std::unordered_map<TString, bool> const& SimpleEntry::getNamedMap<bool>() const{ return namedbools; }
template<> std::unordered_map<TString, short> const& SimpleEntry::getNamedMap<short>() const{ return namedshorts; }
template<> std::unordered_map<TString, unsigned int> const& SimpleEntry::getNamedMap<unsigned int>() const{ return nameduints; }
template<> std::unordered_map<TString, int> const& SimpleEntry::getNamedMap<int>() const{ return namedints; }
template<> std::unordered_map<TString, unsigned long> const& SimpleEntry::getNamedMap<unsigned long>() const{ return namedulongs; }
template<> std::unordered_map<TString, long> const& SimpleEntry::getNamedMap<long>() const{ return namedlongs; }
template<> std::unordered_map<TString, long long> const& SimpleEntry::getNamedMap<long long>() const{ return namedlonglongs; }
template<> std::unordered_map<TString, float> const& SimpleEntry::getNamedMap<float>() const{ return namedfloats; }
template<> std::unordered_map<TString, double> const& SimpleEntry::getNamedMap<double>() const{ return nameddoubles; }
template<> std::unordered_map<TString, CMSLorentzVector> const& SimpleEntry::getNamedMap<CMSLorentzVector>() const{ return namedCMSLorentzVectors; }
template<> std::unordered_map<TString, std::vector<bool>> const& SimpleEntry::getNamedMap<std::vector<bool>>() const{ return namedVbools; }
template<> std::unordered_map<TString, std::vector<short>> const& SimpleEntry::getNamedMap<std::vector<short>>() const{ return namedVshorts; }
template<> std::unordered_map<TString, std::vector<unsigned int>> const& SimpleEntry::getNamedMap<std::vector<unsigned int>>() const{ return namedVuints; }
template<> std::unordered_map<TString, std::vector<int>> const& SimpleEntry::getNamedMap<std::vector<int>>() const{ return namedVints; }
template<> std::unordered_map<TString, std::vector<unsigned long>> const& SimpleEntry::getNamedMap<std::vector<unsigned long>>() const{ return namedVulongs; }
template<> std::unordered_map<TString, std::vector<long>> const& SimpleEntry::getNamedMap<std::vector<long>>() const{ return namedVlongs; }
template<> std::unordered_map<TString, std::vector<long long>> const& SimpleEntry::getNamedMap<std::vector<long long>>() const{ return namedVlonglongs; }
template<> std::unordered_map<TString, std::vector<float>> const& SimpleEntry::getNamedMap<std::vector<float>>() const{ return namedVfloats; }
template<> std::unordered_map<TString, std::vector<double>> const& SimpleEntry::getNamedMap<std::vector<double>>() const{ return namedVdoubles; }
template<> std::unordered_map<TString, std::vector<CMSLorentzVector>> const& SimpleEntry::getNamedMap<std::vector<CMSLorentzVector>>() const{ return namedVCMSLorentzVectors; }

template<> std::unordered_map<TString, bool>& SimpleEntry::getNamedMap<bool>(){ return namedbools; }
template<> std::unordered_map<TString, short>& SimpleEntry::getNamedMap<short>(){ return namedshorts; }
template<> std::unordered_map<TString, unsigned int>& SimpleEntry::getNamedMap<unsigned int>(){ return nameduints; }
template<> std::unordered_map<TString, int>& SimpleEntry::getNamedMap<int>(){ return namedints; }
template<> std::unordered_map<TString, unsigned long>& SimpleEntry::getNamedMap<unsigned long>(){ return namedulongs; }
template<> std::unordered_map<TString, long>& SimpleEntry::getNamedMap<long>(){ return namedlongs; }
template<> std::unordered_map<TString, long long>& SimpleEntry::getNamedMap<long long>(){ return namedlonglongs; }
template<> std::unordered_map<TString, float>& SimpleEntry::getNamedMap<float>(){ return namedfloats; }
template<> std::unordered_map<TString, double>& SimpleEntry::getNamedMap<double>(){ return nameddoubles; }
template<> std::unordered_map<TString, CMSLorentzVector>& SimpleEntry::getNamedMap<CMSLorentzVector>(){ return namedCMSLorentzVectors; }
template<> std::unordered_map<TString, std::vector<bool>>& SimpleEntry::getNamedMap<std::vector<bool>>(){ return namedVbools; }
template<> std::unordered_map<TString, std::vector<short>>& SimpleEntry::getNamedMap<std::vector<short>>(){ return namedVshorts; }
template<> std::unordered_map<TString, std::vector<unsigned int>>& SimpleEntry::getNamedMap<std::vector<unsigned int>>(){ return namedVuints; }
template<> std::unordered_map<TString, std::vector<int>>& SimpleEntry::getNamedMap<std::vector<int>>(){ return namedVints; }
template<> std::unordered_map<TString, std::vector<unsigned long>>& SimpleEntry::getNamedMap<std::vector<unsigned long>>(){ return namedVulongs; }
template<> std::unordered_map<TString, std::vector<long>>& SimpleEntry::getNamedMap<std::vector<long>>(){ return namedVlongs; }
template<> std::unordered_map<TString, std::vector<long long>>& SimpleEntry::getNamedMap<std::vector<long long>>(){ return namedVlonglongs; }
template<> std::unordered_map<TString, std::vector<float>>& SimpleEntry::getNamedMap<std::vector<float>>(){ return namedVfloats; }
template<> std::unordered_map<TString, std::vector<double>>& SimpleEntry::getNamedMap<std::vector<double>>(){ return namedVdoubles; }
template<> std::unordered_map<TString, std::vector<CMSLorentzVector>>& SimpleEntry::getNamedMap<std::vector<CMSLorentzVector>>(){ return namedVCMSLorentzVectors; }

template<typename T> void SimpleEntry::setNamedVal(TString strname, T const& val){
  std::unordered_map<TString, T>& theMap = this->getNamedMap<T>();
  theMap[strname] = val;
}
template void SimpleEntry::setNamedVal<bool>(TString strname, bool const& val);
template void SimpleEntry::setNamedVal<short>(TString strname, short const& val);
template void SimpleEntry::setNamedVal<unsigned int>(TString strname, unsigned int const& val);
template void SimpleEntry::setNamedVal<int>(TString strname, int const& val);
template void SimpleEntry::setNamedVal<unsigned long>(TString strname, unsigned long const& val);
template void SimpleEntry::setNamedVal<long>(TString strname, long const& val);
template void SimpleEntry::setNamedVal<long long>(TString strname, long long const& val);
template void SimpleEntry::setNamedVal<float>(TString strname, float const& val);
template void SimpleEntry::setNamedVal<double>(TString strname, double const& val);
template void SimpleEntry::setNamedVal<CMSLorentzVector>(TString strname, CMSLorentzVector const& val);
template void SimpleEntry::setNamedVal<std::vector<bool>>(TString strname, std::vector<bool> const& val);
template void SimpleEntry::setNamedVal<std::vector<short>>(TString strname, std::vector<short> const& val);
template void SimpleEntry::setNamedVal<std::vector<unsigned int>>(TString strname, std::vector<unsigned int> const& val);
template void SimpleEntry::setNamedVal<std::vector<int>>(TString strname, std::vector<int> const& val);
template void SimpleEntry::setNamedVal<std::vector<unsigned long>>(TString strname, std::vector<unsigned long> const& val);
template void SimpleEntry::setNamedVal<std::vector<long>>(TString strname, std::vector<long> const& val);
template void SimpleEntry::setNamedVal<std::vector<long long>>(TString strname, std::vector<long long> const& val);
template void SimpleEntry::setNamedVal<std::vector<float>>(TString strname, std::vector<float> const& val);
template void SimpleEntry::setNamedVal<std::vector<double>>(TString strname, std::vector<double> const& val);
template void SimpleEntry::setNamedVal<std::vector<CMSLorentzVector>>(TString strname, std::vector<CMSLorentzVector> const& val);

template<typename T> void SimpleEntry::getNamedVal(TString strname, T& val) const{
  typename std::unordered_map<TString, T>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<T>(strname, this->getNamedMap<T>(), it)) val = it->second;
}
template void SimpleEntry::getNamedVal<bool>(TString strname, bool& val) const;
template void SimpleEntry::getNamedVal<short>(TString strname, short& val) const;
template void SimpleEntry::getNamedVal<unsigned int>(TString strname, unsigned int& val) const;
template void SimpleEntry::getNamedVal<int>(TString strname, int& val) const;
template void SimpleEntry::getNamedVal<unsigned long>(TString strname, unsigned long& val) const;
template void SimpleEntry::getNamedVal<long>(TString strname, long& val) const;
template void SimpleEntry::getNamedVal<long long>(TString strname, long long& val) const;
template void SimpleEntry::getNamedVal<float>(TString strname, float& val) const;
template void SimpleEntry::getNamedVal<double>(TString strname, double& val) const;
template void SimpleEntry::getNamedVal<CMSLorentzVector>(TString strname, CMSLorentzVector& val) const;
template void SimpleEntry::getNamedVal<std::vector<bool>>(TString strname, std::vector<bool>& val) const;
template void SimpleEntry::getNamedVal<std::vector<short>>(TString strname, std::vector<short>& val) const;
template void SimpleEntry::getNamedVal<std::vector<unsigned int>>(TString strname, std::vector<unsigned int>& val) const;
template void SimpleEntry::getNamedVal<std::vector<int>>(TString strname, std::vector<int>& val) const;
template void SimpleEntry::getNamedVal<std::vector<unsigned long>>(TString strname, std::vector<unsigned long>& val) const;
template void SimpleEntry::getNamedVal<std::vector<long>>(TString strname, std::vector<long>& val) const;
template void SimpleEntry::getNamedVal<std::vector<long long>>(TString strname, std::vector<long long>& val) const;
template void SimpleEntry::getNamedVal<std::vector<float>>(TString strname, std::vector<float>& val) const;
template void SimpleEntry::getNamedVal<std::vector<double>>(TString strname, std::vector<double>& val) const;
template void SimpleEntry::getNamedVal<std::vector<CMSLorentzVector>>(TString strname, std::vector<CMSLorentzVector>& val) const;


#endif
