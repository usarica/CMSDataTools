#include "ExtendedBinning.h"
#include "HelperFunctions.h"


ExtendedBinning::ExtendedBinning(const TString label_) : label(label_){}
ExtendedBinning::ExtendedBinning(const unsigned int nbins, const double min, const double max, const TString label_) : label(label_){
  if (max>min && nbins>0){
    const double inc = (max-min)/((const double) nbins);
    vbinlow.reserve(nbins+1);
    for (unsigned int i=0; i<=nbins; i++) vbinlow.push_back(min+inc*(double(i)));
  }
}
ExtendedBinning::ExtendedBinning(const double* abinlow, const TString label_) : label(label_){
  if (abinlow!=nullptr){
    const int np = sizeof(abinlow)/sizeof(abinlow[0]);
    vbinlow = std::vector<double>(abinlow, abinlow+np);
  }
}
ExtendedBinning::ExtendedBinning(const std::vector<double>& vbinlow_, const TString label_) : vbinlow(vbinlow_), label(label_){}
ExtendedBinning::ExtendedBinning(ExtendedBinning const& other) : vbinlow(other.vbinlow), label(other.label){}

bool ExtendedBinning::isValid() const{
  return (vbinlow.size()>1);
}

void ExtendedBinning::setLabel(const TString label_){ label=label_; }
TString ExtendedBinning::getLabel() const{ return label; }

double* ExtendedBinning::getBinning(){ return vbinlow.data(); }
const double* ExtendedBinning::getBinning() const{ return vbinlow.data(); }
std::vector<double> ExtendedBinning::getBinningVector(){ return vbinlow; }
const std::vector<double>& ExtendedBinning::getBinningVector() const{ return vbinlow; }

unsigned int ExtendedBinning::getNbins() const{
  return (this->isValid() ? vbinlow.size()-1 : 0);
}

int ExtendedBinning::getBin(double val) const{
  if (!this->isValid()) return -1;
  for (int bin=0; bin<(int)(vbinlow.size()-1); bin++){
    if (vbinlow.at(bin)<=val && val<vbinlow.at(bin+1)) return bin;
  }
  if (val>=vbinlow.back()) return vbinlow.size();
  else return -1;
}
double ExtendedBinning::getBinLowEdge(const int bin) const{
  if (bin>=0 && bin<(int)vbinlow.size()) return vbinlow.at(bin);
  else if (bin<0 && vbinlow.size()>0) return vbinlow.at(0);
  else if (vbinlow.size()>0) return vbinlow.at(vbinlow.size()-1);
  else return -1;
}
double ExtendedBinning::getBinHighEdge(const int bin) const{
  if (bin>=-1 && bin<(int)vbinlow.size()-1) return vbinlow.at(bin+1);
  else if (bin<-1 && vbinlow.size()>0) return vbinlow.at(0);
  else if (vbinlow.size()>0) return vbinlow.at(vbinlow.size()-1);
  else return -1;
}

void ExtendedBinning::addBinBoundary(double boundary){
  HelperFunctions::addByLowest<double>(vbinlow, boundary, true);
}
void ExtendedBinning::removeBinLowEdge(const int bin){
  if (bin>=0 && bin<(int) vbinlow.size()) vbinlow.erase(vbinlow.begin()+bin);
}

