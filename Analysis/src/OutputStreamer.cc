#include "OutputStreamer.h"


using namespace std;


OutputStreamer::OutputStreamer(const char* fname, std::ios_base::openmode fmode, bool printError) :
theFile(fname, fmode)
{
  if (!printError) stdout_ptr = &std::cout;
  else stdout_ptr = &std::cerr;
}
OutputStreamer::~OutputStreamer(){ this->close(); }

OutputStreamer& OutputStreamer::operator<<(std::ostream& (*fcn)(std::ostream&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
OutputStreamer& OutputStreamer::operator<<(std::ios& (*fcn)(std::ios&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
OutputStreamer& OutputStreamer::operator<<(std::ios_base& (*fcn)(std::ios_base&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}

std::streamsize OutputStreamer::width() const{ return theFile.width(); }
std::streamsize OutputStreamer::width(std::streamsize wide){
  if (stdout_ptr) stdout_ptr->width(wide);
  return theFile.width(wide);
}

char OutputStreamer::fill() const{ return theFile.fill(); }
char OutputStreamer::fill(char fillch){
  if (stdout_ptr) stdout_ptr->fill(fillch);
  return theFile.fill(fillch);
}

void OutputStreamer::close(){
  theFile.close();
}
void OutputStreamer::open(const char* fname, std::ios_base::openmode fmode){
  this->close();
  theFile.open(fname, fmode);
}
