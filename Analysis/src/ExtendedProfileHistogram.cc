#include "ExtendedProfileHistogram.h"
#include "HelperFunctions.h"


using namespace std;
using namespace HelperFunctions;


ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, unsigned int fillAverageQ_) :
  fillAverageQ(fillAverageQ_),
  xbinning(bX),
  sumWQ(nullptr), sumWQsq(nullptr)
{
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(1, std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  if (fillAverageQ>0){
    sumWQ=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    sumWQsq=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    for (unsigned int iQ=0; iQ<fillAverageQ; iQ++){
      sumWQ[iQ]=sumW;
      sumWQsq[iQ]=sumW;
    }
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, unsigned int fillAverageQ_) :
  fillAverageQ(fillAverageQ_),
  xbinning(bX), ybinning(bY),
  sumWQ(nullptr), sumWQsq(nullptr)
{
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  if (fillAverageQ>0){
    sumWQ=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    sumWQsq=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    for (unsigned int iQ=0; iQ<fillAverageQ; iQ++){
      sumWQ[iQ]=sumW;
      sumWQsq[iQ]=sumW;
    }
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ, unsigned int fillAverageQ_) :
  fillAverageQ(fillAverageQ_),
  xbinning(bX), ybinning(bY), zbinning(bZ),
  sumWQ(nullptr), sumWQsq(nullptr)
{
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(zbinning.getNbins(), double(0.)); }
  sumWsq=sumW;
  if (fillAverageQ>0){
    sumWQ=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    sumWQsq=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    for (unsigned int iQ=0; iQ<fillAverageQ; iQ++){
      sumWQ[iQ]=sumW;
      sumWQsq[iQ]=sumW;
    }
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedProfileHistogram const& other) :
  fillAverageQ(other.fillAverageQ),
  xbinning(other.xbinning),
  ybinning(other.ybinning),
  zbinning(other.zbinning),
  sumW(other.sumW),
  sumWsq(other.sumWsq)
{
  if (fillAverageQ>0){
    sumWQ=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    sumWQsq=new std::vector<std::vector<std::vector<double>>>[fillAverageQ];
    for (unsigned int iQ=0; iQ<fillAverageQ; iQ++){
      sumWQ[iQ]=other.sumWQ[iQ];
      sumWQsq[iQ]=other.sumWQsq[iQ];
    }
  }
}
ExtendedProfileHistogram::~ExtendedProfileHistogram(){
  if (fillAverageQ>0){
    delete[] sumWQ;
    delete[] sumWQsq;
  }
}


double ExtendedProfileHistogram::getBinSumW(int ix, int iy, int iz) const{
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sumW.at(ix).at(iy).at(iz);
}
double ExtendedProfileHistogram::getBinSumWsq(int ix, int iy, int iz) const{
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sumWsq.at(ix).at(iy).at(iz);
}
double ExtendedProfileHistogram::getBinAvgQ(unsigned int iQ, int ix, int iy, int iz) const{
  if (fillAverageQ==0 || iQ>=fillAverageQ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  const std::vector<std::vector<std::vector<double>>>& sWQ = sumWQ[iQ];
  return (sumW.at(ix).at(iy).at(iz)!=0. ? sWQ.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) : 0.);
}
double ExtendedProfileHistogram::getBinSigmaQ(unsigned int iQ, int ix, int iy, int iz) const{
  if (fillAverageQ==0 || iQ>=fillAverageQ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  const std::vector<std::vector<std::vector<double>>>& sWQ = sumWQ[iQ];
  const std::vector<std::vector<std::vector<double>>>& sWQsq = sumWQsq[iQ];
  return sqrt(max(0., (sumW.at(ix).at(iy).at(iz)!=0. ? sWQsq.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) - pow(sWQ.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz), 2) : 0.)));
}

std::vector<std::vector<std::vector<double>>>& ExtendedProfileHistogram::getSumWQContainer(unsigned int iQ){
  assert(fillAverageQ!=0 && iQ<fillAverageQ);
  return sumWQ[iQ];
}
std::vector<std::vector<std::vector<double>>>& ExtendedProfileHistogram::getSumWQsqContainer(unsigned int iQ){
  assert(fillAverageQ!=0 && iQ<fillAverageQ);
  return sumWQsq[iQ];
}
std::vector<std::vector<std::vector<double>>> const& ExtendedProfileHistogram::getSumWQContainer(unsigned int iQ) const{
  assert(fillAverageQ!=0 && iQ<fillAverageQ);
  return sumWQ[iQ];
}
std::vector<std::vector<std::vector<double>>> const& ExtendedProfileHistogram::getSumWQsqContainer(unsigned int iQ) const{
  assert(fillAverageQ!=0 && iQ<fillAverageQ);
  return sumWQsq[iQ];
}


void ExtendedProfileHistogram::fill(double x, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (ybinning.isValid()) return;
  if (zbinning.isValid()) return;
  sumW.at(ix).at(iy).at(iz)+=w;
  sumWsq.at(ix).at(iy).at(iz)+=w*w;
}
void ExtendedProfileHistogram::fill(double x, double y, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (zbinning.isValid()) return;
  sumW.at(ix).at(iy).at(iz)+=w;
  sumWsq.at(ix).at(iy).at(iz)+=w*w;
}
void ExtendedProfileHistogram::fill(double x, double y, double z, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (!zbinning.isValid()) return;
  else{ iz=zbinning.getBin(z); if (iz<0 || iz>=(int) zbinning.getNbins()) return; }
  sumW.at(ix).at(iy).at(iz)+=w;
  sumWsq.at(ix).at(iy).at(iz)+=w*w;
}


void ExtendedProfileHistogram::fillQ(unsigned int iQ, double x, double Q, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (ybinning.isValid()) return;
  if (zbinning.isValid()) return;
  std::vector<std::vector<std::vector<double>>>& sWQ = sumWQ[iQ];
  std::vector<std::vector<std::vector<double>>>& sWQsq = sumWQsq[iQ];
  sWQ.at(ix).at(iy).at(iz)+=w*Q;
  sWQsq.at(ix).at(iy).at(iz)+=w*Q*Q;
}
void ExtendedProfileHistogram::fillQ(unsigned int iQ, double x, double y, double Q, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (zbinning.isValid()) return;
  std::vector<std::vector<std::vector<double>>>& sWQ = sumWQ[iQ];
  std::vector<std::vector<std::vector<double>>>& sWQsq = sumWQsq[iQ];
  sWQ.at(ix).at(iy).at(iz)+=w*Q;
  sWQsq.at(ix).at(iy).at(iz)+=w*Q*Q;
}
void ExtendedProfileHistogram::fillQ(unsigned int iQ, double x, double y, double z, double Q, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (!zbinning.isValid()) return;
  else{ iz=zbinning.getBin(z); if (iz<0 || iz>=(int) zbinning.getNbins()) return; }
  std::vector<std::vector<std::vector<double>>>& sWQ = sumWQ[iQ];
  std::vector<std::vector<std::vector<double>>>& sWQsq = sumWQsq[iQ];
  sWQ.at(ix).at(iy).at(iz)+=w*Q;
  sWQsq.at(ix).at(iy).at(iz)+=w*Q*Q;
}
