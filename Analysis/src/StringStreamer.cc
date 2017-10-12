#include "StringStreamer.h"


using namespace std;


StringStreamer::StringStreamer(const char* fname) : theFile(fname, ios::out){}
StringStreamer::~StringStreamer(){ theFile.close(); }

StringStreamer& StringStreamer::operator<<(std::ostream& (*fcn)(std::ostream&)){
  fcn(theFile);
  fcn(cout);
  return *this;
}
