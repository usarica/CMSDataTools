#include "ExtendedProfileHistogram.h"
#include "HelperFunctions.h"


using namespace std;
using namespace HelperFunctions;


ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX) : xbinning(bX){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(1, std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  sumWX=sumW;
  sumWXsq=sumW;
  sumWY=sumW;
  sumWYsq=sumW;
  sumWZ=sumW;
  sumWZsq=sumW;
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY) : xbinning(bX), ybinning(bY){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(1, double(0.)); }
  sumWsq=sumW;
  sumWX=sumW;
  sumWXsq=sumW;
  sumWY=sumW;
  sumWYsq=sumW;
  sumWZ=sumW;
  sumWZsq=sumW;
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedBinning bX, ExtendedBinning bY, ExtendedBinning bZ) : xbinning(bX), ybinning(bY), zbinning(bZ){
  sumW.assign(xbinning.getNbins(), std::vector<std::vector<double>>()); for (auto& v:sumW){ v.assign(ybinning.getNbins(), std::vector<double>()); for (auto& vv:v) vv.assign(zbinning.getNbins(), double(0.)); }
  sumWsq=sumW;
  sumWX=sumW;
  sumWXsq=sumW;
  sumWY=sumW;
  sumWYsq=sumW;
  sumWZ=sumW;
  sumWZsq=sumW;
}
ExtendedProfileHistogram::ExtendedProfileHistogram(ExtendedProfileHistogram const& other) :
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
  return sumW[ix][iy][iz];
}
double ExtendedProfileHistogram::getBinSumWsq(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sumWsq[ix][iy][iz];
}
double ExtendedProfileHistogram::getBinAvgX(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW[ix][iy][iz]!=0. ? sumWX[ix][iy][iz] / sumW[ix][iy][iz] : 0.);
}
double ExtendedProfileHistogram::getBinSigmaX(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW[ix][iy][iz]!=0. ? sumWXsq[ix][iy][iz] / sumW[ix][iy][iz] - pow(sumWX[ix][iy][iz] / sumW[ix][iy][iz], 2) : 0.)));
}
double ExtendedProfileHistogram::getBinAvgY(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW[ix][iy][iz]!=0. ? sumWY[ix][iy][iz] / sumW[ix][iy][iz] : 0.);
}
double ExtendedProfileHistogram::getBinSigmaY(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW[ix][iy][iz]!=0. ? sumWYsq[ix][iy][iz] / sumW[ix][iy][iz] - pow(sumWY[ix][iy][iz] / sumW[ix][iy][iz], 2) : 0.)));
}
double ExtendedProfileHistogram::getBinAvgZ(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return (sumW[ix][iy][iz]!=0. ? sumWZ[ix][iy][iz] / sumW[ix][iy][iz] : 0.);
}
double ExtendedProfileHistogram::getBinSigmaZ(int ix, int iy, int iz){
  if (ix<0 || (xbinning.isValid() && ix>=(int) xbinning.getNbins()) || (!xbinning.isValid() && ix>0)) return 0;
  if (iy<0 || (ybinning.isValid() && iy>=(int) ybinning.getNbins()) || (!ybinning.isValid() && iy>0)) return 0;
  if (iz<0 || (zbinning.isValid() && iz>=(int) zbinning.getNbins()) || (!zbinning.isValid() && iz>0)) return 0;
  return sqrt(max(0., (sumW[ix][iy][iz]!=0. ? sumWZsq[ix][iy][iz] / sumW[ix][iy][iz] - pow(sumWZ[ix][iy][iz] / sumW[ix][iy][iz], 2) : 0.)));
}

void ExtendedProfileHistogram::fill(double x, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (ybinning.isValid()) return;
  if (zbinning.isValid()) return;
  sumW[ix][iy][iz]+=w;
  sumWsq[ix][iy][iz]+=w*w;
  sumWX[ix][iy][iz]+=w*x;
  sumWXsq[ix][iy][iz]+=w*x*x;
}
void ExtendedProfileHistogram::fill(double x, double y, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (zbinning.isValid()) return;
  sumW[ix][iy][iz]+=w;
  sumWsq[ix][iy][iz]+=w*w;
  sumWX[ix][iy][iz]+=w*x;
  sumWXsq[ix][iy][iz]+=w*x*x;
  sumWY[ix][iy][iz]+=w*y;
  sumWYsq[ix][iy][iz]+=w*y*y;
}
void ExtendedProfileHistogram::fill(double x, double y, double z, double w){
  int ix=0, iy=0, iz=0;
  if (!xbinning.isValid()) return;
  else{ ix=xbinning.getBin(x); if (ix<0 || ix>=(int) xbinning.getNbins()) return; }
  if (!ybinning.isValid()) return;
  else{ iy=ybinning.getBin(y); if (iy<0 || iy>=(int) ybinning.getNbins()) return; }
  if (!zbinning.isValid()) return;
  else{ iz=zbinning.getBin(z); if (iz<0 || iz>=(int) zbinning.getNbins()) return; }
  sumW[ix][iy][iz]+=w;
  sumWsq[ix][iy][iz]+=w*w;
  sumWX[ix][iy][iz]+=w*x;
  sumWXsq[ix][iy][iz]+=w*x*x;
  sumWY[ix][iy][iz]+=w*y;
  sumWYsq[ix][iy][iz]+=w*y*y;
  sumWZ[ix][iy][iz]+=w*z;
  sumWZsq[ix][iy][iz]+=w*z*z;
}
