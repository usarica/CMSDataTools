#ifndef STRINGSTREAMER_H
#define STRINGSTREAMER_H

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TString.h"


class StringStreamer{
protected:
  std::ofstream theFile;

public:
  StringStreamer(const char* fname);
  ~StringStreamer();

  template<typename T> StringStreamer& operator<<(const T& val);
  StringStreamer& operator<<(std::ostream& (*fcn)(std::ostream&)); // For endl

};

template<typename T> StringStreamer& StringStreamer::operator<<(const T& val){
  theFile << val;
  std::cout << val;
  return *this;
}
template StringStreamer& StringStreamer::operator<< <bool>(const bool& val);
template StringStreamer& StringStreamer::operator<< <short>(const short& val);
template StringStreamer& StringStreamer::operator<< <unsigned int>(const unsigned int& val);
template StringStreamer& StringStreamer::operator<< <float>(const float& val);
template StringStreamer& StringStreamer::operator<< <double>(const double& val);
template StringStreamer& StringStreamer::operator<< <char*>(char* const& val);
template StringStreamer& StringStreamer::operator<< <char const*>(char const* const& val);
template StringStreamer& StringStreamer::operator<< <char>(const char& val);
template StringStreamer& StringStreamer::operator<< <std::string>(const std::string& val);
template StringStreamer& StringStreamer::operator<< <TString>(const TString& val);


#endif
