#include "ExtendedProfileHistogram.h"
#include "HelperFunctions.h"


using namespace std;
using namespace HelperFunctions;


ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, bool fillAverageXYZ_) : fillAverageXYZ(fillAverageXYZ_), xbinning(bX){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(1, std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  if (fillAverageXYZ){
    sumWX=sumW;
    sumWXsq=sumW;
    sumWY=sumW;
    sumWYsq=sumW;
    sumWZ=sumW;
    sumWZsq=sumW;
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, bool fillAverageXYZ_) : fillAverageXYZ(fillAverageXYZ_), xbinning(bX), ybinning(bY){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  if (fillAverageXYZ){
    sumWX=sumW;
    sumWXsq=sumW;
    sumWY=sumW;
    sumWYsq=sumW;
    sumWZ=sumW;
    sumWZsq=sumW;
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ, bool fillAverageXYZ_) : fillAverageXYZ(fillAverageXYZ_), xbinning(bX), ybinning(bY), zbinning(bZ){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(zbinning.getNbins(), double(0.)); }
  sumWsq=sumW;
  if (fillAverageXYZ){
    sumWX=sumW;
    sumWXsq=sumW;
    sumWY=sumW;
    sumWYsq=sumW;
    sumWZ=sumW;
    sumWZsq=sumW;
  }
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedProfileHistogram const& other) :
  fillAverageXYZ(other.fillAverageXYZ),
  xbinning(other.xbinning),
  ybinning(other.ybinning),
  zbinning(other.zbinning),
  sumW(other.sumW),
  sumWsq(other.sumWsq),
  sumWX(other.sumWX),
  sumWXsq(other.sumWXsq),
  sumWY(other.sumWY),
  sumWYsq(other.sumWYsq),
  sumWZ(other.sumWZ),
  sumWZsq(other.sumWZsq)
{}

double ExtendedProfileHistogram::getBinSumW(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sumW.at(ix).at(iy).at(iz);
}
double ExtendedProfileHistogram::getBinSumWsq(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sumWsq.at(ix).at(iy).at(iz);
}
double ExtendedProfileHistogram::getBinAvgX(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW.at(ix).at(iy).at(iz)!=0. ? sumWX.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) : 0.);
}
double ExtendedProfileHistogram::getBinSigmaX(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW.at(ix).at(iy).at(iz)!=0. ? sumWXsq.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) - pow(sumWX.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz), 2) : 0.)));
}
double ExtendedProfileHistogram::getBinAvgY(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW.at(ix).at(iy).at(iz)!=0. ? sumWY.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) : 0.);
}
double ExtendedProfileHistogram::getBinSigmaY(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW.at(ix).at(iy).at(iz)!=0. ? sumWYsq.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) - pow(sumWY.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz), 2) : 0.)));
}
double ExtendedProfileHistogram::getBinAvgZ(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW.at(ix).at(iy).at(iz)!=0. ? sumWZ.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) : 0.);
}
double ExtendedProfileHistogram::getBinSigmaZ(int ix, int iy, int iz){
  if (!fillAverageXYZ) return 0;
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW.at(ix).at(iy).at(iz)!=0. ? sumWZsq.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz) - pow(sumWZ.at(ix).at(iy).at(iz) / sumW.at(ix).at(iy).at(iz), 2) : 0.)));
}

void ExtendedProfileHistogram::fill(double x, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (ybinning.isValid()) return;
  if (zbinning.isValid()) return;
  sumW.at(ix).at(iy).at(iz)+=w;
  sumWsq.at(ix).at(iy).at(iz)+=w*w;
  if (!fillAverageXYZ) return;
  sumWX.at(ix).at(iy).at(iz)+=w*x;
  sumWXsq.at(ix).at(iy).at(iz)+=w*x*x;
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
  if (!fillAverageXYZ) return;
  sumWX.at(ix).at(iy).at(iz)+=w*x;
  sumWXsq.at(ix).at(iy).at(iz)+=w*x*x;
  sumWY.at(ix).at(iy).at(iz)+=w*y;
  sumWYsq.at(ix).at(iy).at(iz)+=w*y*y;
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
  if (!fillAverageXYZ) return;
  sumWX.at(ix).at(iy).at(iz)+=w*x;
  sumWXsq.at(ix).at(iy).at(iz)+=w*x*x;
  sumWY.at(ix).at(iy).at(iz)+=w*y;
  sumWYsq.at(ix).at(iy).at(iz)+=w*y*y;
  sumWZ.at(ix).at(iy).at(iz)+=w*z;
  sumWZsq.at(ix).at(iy).at(iz)+=w*z*z;
}
