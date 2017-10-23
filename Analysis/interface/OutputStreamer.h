#ifndef OUTPUTSTREAMER_H
#define OUTPUTSTREAMER_H

#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <iomanip>
#include <string>
#include "TString.h"


class OutputStreamer{
protected:
  std::ofstream theFile;
  std::ostream* stdout_ptr;

public:
  OutputStreamer(const char* fname, std::ios_base::openmode fmode = std::ios_base::out, bool printError=false);
  ~OutputStreamer();

  template<typename T> OutputStreamer& operator<<(const T& val);
  OutputStreamer& operator<<(std::ostream& (*fcn)(std::ostream&));
  OutputStreamer& operator<<(std::ios& (*fcn)(std::ios&));
  OutputStreamer& operator<<(std::ios_base& (*fcn)(std::ios_base&));

  std::streamsize width() const;
  std::streamsize width(std::streamsize wide);

  char fill() const;
  char fill(char fillch);

  void open(const char* fname, std::ios_base::openmode fmode = std::ios_base::out);
  void close();

  template<typename T> void writeCentered(const T& val, char fillch=' ', std::streamsize gapsize=0);

};

template<typename T> OutputStreamer& OutputStreamer::operator<<(const T& val){
  theFile << val;
  if (stdout_ptr) *stdout_ptr << val;
  return *this;
}
template OutputStreamer& OutputStreamer::operator<< <bool>(const bool& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned short>(const unsigned short& val);
template OutputStreamer& OutputStreamer::operator<< <short>(const short& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned int>(const unsigned int& val);
template OutputStreamer& OutputStreamer::operator<< <int>(const int& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned long>(const unsigned long& val);
template OutputStreamer& OutputStreamer::operator<< <long>(const long& val);
template OutputStreamer& OutputStreamer::operator<< <float>(const float& val);
template OutputStreamer& OutputStreamer::operator<< <double>(const double& val);
template OutputStreamer& OutputStreamer::operator<< <long double>(const long double& val);

template OutputStreamer& OutputStreamer::operator<< <char*>(char* const& val);
template OutputStreamer& OutputStreamer::operator<< <char const*>(char const* const& val);
template OutputStreamer& OutputStreamer::operator<< <char>(const char& val);
template OutputStreamer& OutputStreamer::operator<< <signed char*>(signed char* const& val);
template OutputStreamer& OutputStreamer::operator<< <signed char const*>(signed char const* const& val);
template OutputStreamer& OutputStreamer::operator<< <signed char>(const signed char& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned char*>(unsigned char* const& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned char const*>(unsigned char const* const& val);
template OutputStreamer& OutputStreamer::operator<< <unsigned char>(const unsigned char& val);

template OutputStreamer& OutputStreamer::operator<< <std::string>(std::string const& val);
template OutputStreamer& OutputStreamer::operator<< <TString>(TString const& val);
template OutputStreamer& OutputStreamer::operator<< <std::streambuf*>(std::streambuf* const& val);
template OutputStreamer& OutputStreamer::operator<< <void*>(void* const& val);

template<typename T> void OutputStreamer::writeCentered(const T& val, char fillch, std::streamsize gapsize){
  char deffillch = this->fill(fillch);

  std::stringstream tmpss;
  tmpss << val;
  std::string tmpstr = tmpss.str();
  std::streamsize strlength = (std::streamsize) tmpstr.length();

  if (strlength>gapsize) *this << std::setw(gapsize) << "";
  else{
    std::streamsize leftgap = (gapsize+strlength)/2;
    std::streamsize rightgap = gapsize-leftgap;
    *this << std::setw(leftgap) << tmpstr << std::setw(rightgap) << "";
  }

  this->fill(deffillch);
}
template void OutputStreamer::writeCentered<bool>(const bool& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned short>(const unsigned short& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<short>(const short& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned int>(const unsigned int& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<int>(const int& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned long>(const unsigned long& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<float>(const float& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<double>(const double& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<long double>(const long double& val, char fillch, std::streamsize gapsize);

template void OutputStreamer::writeCentered<char*>(char* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<char const*>(char const* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<char>(char const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<signed char*>(signed char* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<signed char const*>(signed char const* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<signed char>(signed char const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned char*>(unsigned char* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned char const*>(unsigned char const* const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<unsigned char>(unsigned char const& val, char fillch, std::streamsize gapsize);

template void OutputStreamer::writeCentered<std::string>(std::string const& val, char fillch, std::streamsize gapsize);
template void OutputStreamer::writeCentered<TString>(TString const& val, char fillch, std::streamsize gapsize);



#endif
