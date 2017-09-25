#include "ExtendedBinning.h"


ExtendedBinning::ExtendedBinning() : isvalid(false){}
ExtendedBinning::ExtendedBinning(const unsigned int nbins, const double min, const double max, const TString label_) : label(label_), isvalid(nbins>0){
  if (isvalid){
    const double inc = (max-min)/((const double)nbins);
    vbinlow.reserve(nbins);
    for (unsigned int i=0; i<=nbins; i++) vbinlow.push_back(min+inc*(double(i)));
  }
}
ExtendedBinning::ExtendedBinning(const double* abinlow, const TString label_) : label(label_), isvalid(abinlow!=nullptr){
  if (isvalid){
    const int nbins = sizeof(abinlow)/sizeof(abinlow[0]);
    vbinlow = std::vector<double>(abinlow, abinlow+nbins);
    isvalid=true;
  }
}
ExtendedBinning::ExtendedBinning(const std::vector<double>& vbinlow_, const TString label_) : vbinlow(vbinlow_), label(label_), isvalid(vbinlow.size()>0){}

void ExtendedBinning::setLabel(const TString label_){ label=label_; }

double* ExtendedBinning::getBinning(){ return vbinlow.data(); }
const double* ExtendedBinning::getBinning()const{ return vbinlow.data(); }
const unsigned int ExtendedBinning::getNbins()const{ return vbinlow.size(); }

int ExtendedBinning::getBin(double val)const{
  for (unsigned int bin=0; bin<vbinlow.size()-1; bin++){
    if (vbinlow.at(bin)<=val && val<vbinlow.at(bin+1)) return bin;
  }
  if (val>=vbinlow.back()) return vbinlow.size();
  else return -1;
}
double ExtendedBinning::getBinLowEdge(const int bin)const{
  if (bin>=0 && bin<(int)vbinlow.size()) return vbinlow.at(bin);
  else if (bin<0 && vbinlow.size()>0) return vbinlow.at(0);
  else if (vbinlow.size()>0) return vbinlow.at(vbinlow.size()-1);
  else return -1;
}
double ExtendedBinning::getBinHighEdge(const int bin)const{
  if (bin>=-1 && bin<(int)vbinlow.size()-1) return vbinlow.at(bin+1);
  else if (bin<-1 && vbinlow.size()>0) return vbinlow.at(0);
  else if (vbinlow.size()>0) return vbinlow.at(vbinlow.size()-1);
  else return -1;
}

