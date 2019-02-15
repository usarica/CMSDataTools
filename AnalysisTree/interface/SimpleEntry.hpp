#ifndef SIMPLEENTRY_HPP
#define SIMPLEENTRY_HPP

#include "SimpleEntry.h"


template<typename T> void SimpleEntry::setNamedVal(TString strname, T const& val){
  std::unordered_map<TString, T>& theMap = this->getNamedMap<T>();
  theMap[strname] = val;
}
template<typename T> void SimpleEntry::getNamedVal(TString strname, T& val) const{
  typename std::unordered_map<TString, T>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator<T>(strname, this->getNamedMap<T>(), it)) val = it->second;
}


#define SIMPLE_DATA_OUTPUT_DIRECTIVE(name_t, type) \
template<> std::unordered_map<TString, type> const& SimpleEntry::getNamedMap<type>() const{ return named##name_t##s; } \
template<> std::unordered_map<TString, type>& SimpleEntry::getNamedMap<type>(){ return named##name_t##s; } \
template void SimpleEntry::setNamedVal<type>(TString strname, type const& val); \
template void SimpleEntry::getNamedVal<type>(TString strname, type& val) const; \


#define VECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) \
template<> std::unordered_map<TString, type> const& SimpleEntry::getNamedMap<type>() const{ return namedV##name_t##s; } \
template<> std::unordered_map<TString, type>& SimpleEntry::getNamedMap<type>(){ return namedV##name_t##s; } \
template void SimpleEntry::setNamedVal<type>(TString strname, type const& val); \
template void SimpleEntry::getNamedVal<type>(TString strname, type& val) const; \


#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(name_t, type) \
template<> std::unordered_map<TString, type> const& SimpleEntry::getNamedMap<type>() const{ return namedVV##name_t##s; } \
template<> std::unordered_map<TString, type>& SimpleEntry::getNamedMap<type>(){ return namedVV##name_t##s; } \
template void SimpleEntry::setNamedVal<type>(TString strname, type const& val); \
template void SimpleEntry::getNamedVal<type>(TString strname, type& val) const; \


SIMPLE_DATA_OUTPUT_DIRECTIVES
VECTOR_DATA_OUTPUT_DIRECTIVES
DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES


#undef SIMPLE_DATA_OUTPUT_DIRECTIVE
#undef VECTOR_DATA_OUTPUT_DIRECTIVE
#undef DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE


#endif
