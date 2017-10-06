#ifndef CALCHELPERS_H
#define CALCHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TF1.h"
#include "TSpline.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "StdExtensions.h"

namespace CalcHelpers{
  using namespace std;

  // Return A, B, C for A + B*x + C*x**2
  vector<double> solveQuadratic(vector<pair<double, double>> xy){
    vector<double> res;
    if (xy.size()==3){
      double& x1 = xy[0].first;
      double& x2 = xy[1].first;
      double& x3 = xy[2].first;
      double det = pow(x1, 2)*(x3-x2) + pow(x2, 2)*(x1-x3) + pow(x3, 2)*(x2-x1);
      double Ainv[3][3]={
        { x2*x3*(x3-x2), x1*x3*(x1-x3), x1*x2*(x2-x1) },
        { (x2-x3)*(x2+x3), (x3-x1)*(x3+x1), (x1-x2)*(x1+x2) },
        { -(x2-x3), -(x3-x1), -(x1-x2) }
      };
      for (unsigned int i=0; i<3; i++){
        double sum=0;
        for (unsigned int j=0; j<3; j++) sum += Ainv[i][j]/det*xy[j].second;
        res.push_back(sum);
      }
    }
    return res;
  }
  // Return A, B, C for A*exp(-(x-B)**2/(2*C))
  vector<double> solveGaussian(vector<pair<double, double>> xy){
    vector<double> res;
    if (xy.size()==3){
      vector<pair<double, double>> xlny;
      for (auto ipair : xy) xlny.push_back(pair<double, double>(ipair.first, log(ipair.second)));
      vector<double> qc = solveQuadratic(xlny);

      double A = exp(qc[0]-pow(qc[1], 2)/qc[2]);
      double B = -qc[1]*0.5/qc[2];
      double C = -0.5/qc[2];
      res.push_back(A);
      res.push_back(B);
      res.push_back(C);
    }
    return res;
  }

  /* SPECIFIC COMMENT: Convert a TGraph to a TSpline5 */
  TSpline5* convertGraphToSpline5(TGraph* tg, double* dfirst=0, double* dlast=0){
    unsigned int nbins = tg->GetN();
    double* xy[2]={
      tg->GetX(),
      tg->GetY()
    };
    double derivative_first = (xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0]);
    double derivative_last = (xy[1][nbins-1]-xy[1][nbins-2])/(xy[0][nbins-1]-xy[0][nbins-2]);
    TSpline5* spline = new TSpline5("spline", tg, "b1e1", derivative_first, derivative_last);
    spline->SetName(Form("sp_%s", tg->GetName()));
    if (dfirst!=0) *dfirst=derivative_first;
    if (dlast!=0) *dlast=derivative_last;
    return spline;
  }
  TSpline3* convertTSpline5ToTspline3(TSpline5* sp){
    double xmin = sp->GetXmin();
    double xmax = sp->GetXmax();
    const int nbins=500;
    double xyval[2][nbins+1];
    double interval = (xmax-xmin)/nbins;
    for (int ix=0; ix<=nbins; ix++){
      xyval[0][ix] = xmin + ix*interval;
      xyval[1][ix] = sp->Eval(xyval[0][ix]);
    }
    double dfirst = sp->Derivative(xmin);
    double dlast = sp->Derivative(xmax);
    TSpline3* sp_new = new TSpline3("spline", xyval[0], xyval[1], nbins+1, "b1e1", dfirst, dlast);
    sp_new->SetName(sp->GetName());
    sp_new->SetTitle(sp->GetTitle());
    return sp_new;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1/x */
  TF1* getFcn_a0plusa1overX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    a0 = y+s*x;
    a1 = -s*pow(x, 2);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]/x", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
  }


  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0/(1+a1*x) */
  TF1* getFcn_a0overOneplusa1timesX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);
    double k = y/s;

    double a0, a1;
    a1 = -1./(x + k);
    a0 = y*(1.+a1*x);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]/(1+[1]*x)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0/x**2-a1/x */
  TF1* getFcn_a0overX2minusa1overX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    a0 = -y*pow(x, 2)-s*pow(x, 3);
    a1 = -2.*y*x-s*pow(x, 2);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]/x/x-[1]/x", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1/x-(a1/x)**2 */
  TF1* getFcn_a0plusXPinvminusXpsqinv(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    double disc = 1.+8.*x*s;
    if (disc>0.){
      a1 = x/4.*(1.+sqrt(disc));
      a0 = y-a1/x+pow(a1/x, 2);
      cout << x << '\t' << a0 << '\t' << a1 << '\t' << s << '\t' << y << endl;

      TString fcnName;
      if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
      else fcnName = Form("highFcn_%s", sp->GetName());
      TF1* fcn = new TF1(fcnName, "[0]+[1]/x-pow([1]/x, 2)", xmin, xmax);
      fcn->SetParameter(0, a0);
      fcn->SetParameter(1, a1);
      return fcn;
    }
    else return getFcn_a0plusa1overX(sp, xmin, xmax, useLowBound);
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x */
  TF1* getFcn_a0plusa1timesX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    a0 = y-s*x;
    a1 = s;

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]*x", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x+a2*x**2 */
  TF1* getFcn_a0plusa1timesXplusa2overX2(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s, d;
    double xh, sh;
    if (useLowBound){
      x = sp->GetXmin();
      xh = x+0.001;
    }
    else{
      x = sp->GetXmax();
      xh = x-0.001;
    }
    y = sp->Eval(x);
    s = sp->Derivative(x);
    sh = sp->Derivative(xh);
    d = (sh-s)/(xh-x);

    double a0, a1, a2;
    a2 = d*pow(x, 4)/6.;
    a1 = s+2.*a2/pow(x, 3);
    a0 = y-a1*x-a2/pow(x, 2);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]*x+[2]/pow(x, 2)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);
    fcn->SetParameter(2, a2);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x+a2/x**2 */
  TF1* getFcn_a0plusa1timesXplusa2timesX2(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s, d;
    double xh, sh;
    if (useLowBound){
      x = sp->GetXmin();
      xh = x+0.001;
    }
    else{
      x = sp->GetXmax();
      xh = x-0.001;
    }
    y = sp->Eval(x);
    s = sp->Derivative(x);
    sh = sp->Derivative(xh);
    d = (sh-s)/(xh-x);

    double a0, a1, a2;
    a2 = d/2.;
    a1 = s-2.*a2*x;
    a0 = y-a1*x-a2*pow(x, 2);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]*x+[2]*pow(x, 2)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);
    fcn->SetParameter(2, a2);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get external input added by a0 and multiplied by a1 */
  TGraph* getPatch_a0a1(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
    double x, y, s;
    double y_patch, s_patch;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);
    y_patch = sppatch->Eval(x);
    s_patch = sppatch->Derivative(x);

    double a0, a1;
    a1 = s/s_patch;
    a0 = y-a1*y_patch;

    int n = (xmax-xmin+0.5);
    cout << "Patcher n: " << n << endl;
    double xwidth = (xmax-xmin)/n;
    double* xy[2];
    for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
    for (int ip=0; ip<=n; ip++){
      double xval = xmin+xwidth*ip;
      double yval = sppatch->Eval(xval);
      yval = a0+a1*yval;
      xy[0][ip]=xval;
      xy[1][ip]=yval;
    }
    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
    fcn->SetName(fcnName);

    for (unsigned int i=0; i<2; i++) delete[] xy[i];
    return fcn;
  }

  /* SPECIFIC COMMENT: Get external input added by a0+a1/x */
  TGraph* getPatch_a0plusa1overX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);
    y -= sppatch->Eval(x);
    s -= sppatch->Derivative(x);

    double a0, a1;
    a0 = y+s*x;
    a1 = -s*pow(x, 2);

    cout << "getPatch_a0plusa1overX: Checking if patch will flip: " << endl;
    double xext;
    if (useLowBound) xext = x;
    else xext = x;
    if (1.-a1/pow(xext, 2)/sppatch->Derivative(xext)<0.){
      cout << "Sign flips at " << xext << endl;
      if (!forceOutput) return 0;
    }
    else cout << "Sign does not flip at " << xext << endl;

    int n = (xmax-xmin+0.5);
    cout << "Patcher n: " << n << endl;
    double xwidth = (xmax-xmin)/n;
    double* xy[2];
    for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
    for (int ip=0; ip<=n; ip++){
      double xval = xmin+xwidth*ip;
      if (xval==0.) xval += xwidth*0.5;
      double yval = sppatch->Eval(xval);
      yval += (a0+a1/xval);
      xy[0][ip]=xval;
      xy[1][ip]=yval;
    }
    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
    fcn->SetName(fcnName);

    for (unsigned int i=0; i<2; i++) delete[] xy[i];
    return fcn;
  }

  /* SPECIFIC COMMENT: Get external input added by a0+a1*x */
  TGraph* getPatch_a0plusa1timesX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);
    y -= sppatch->Eval(x);
    s -= sppatch->Derivative(x);

    double a0, a1;
    a1 = s;
    a0 = y-a1*x;

    cout << "getPatch_a0plusa1timesX: Checking if patch will flip: " << endl;
    double xext;
    if (useLowBound) xext = x;
    else xext = x;
    if (1.+a1/sppatch->Derivative(xext)<0.){
      cout << "Sign flips at " << xext << endl;
      if (!forceOutput) return 0;
    }
    else cout << "Sign does not flip at " << xext << endl;

    int n = (xmax-xmin+0.5);
    cout << "Patcher n: " << n << endl;
    double xwidth = (xmax-xmin)/n;
    double* xy[2];
    for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
    for (int ip=0; ip<=n; ip++){
      double xval = xmin+xwidth*ip;
      double yval = sppatch->Eval(xval);
      yval += (a0+a1*xval);
      xy[0][ip]=xval;
      xy[1][ip]=yval;
    }
    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
    fcn->SetName(fcnName);

    for (unsigned int i=0; i<2; i++) delete[] xy[i];
    return fcn;
  }

  /* SPECIFIC COMMENT: Get external input added by a0*exp(a1*x) */
  TGraph* getPatch_a0expa1timesX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);
    y -= sppatch->Eval(x);
    s -= sppatch->Derivative(x);

    double a0, a1;
    a1 = s/y;
    a0 = y/exp(a1*x);

    int n = (xmax-xmin+0.5);
    cout << "Patcher n: " << n << endl;
    double xwidth = (xmax-xmin)/n;
    double* xy[2];
    for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
    for (int ip=0; ip<=n; ip++){
      double xval = xmin+xwidth*ip;
      double yval = sppatch->Eval(xval);
      yval += (a0*exp(a1*xval));
      xy[0][ip]=xval;
      xy[1][ip]=yval;
    }
    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
    fcn->SetName(fcnName);

    for (unsigned int i=0; i<2; i++) delete[] xy[i];
    return fcn;
  }

  /* SPECIFIC COMMENT: THIS FUNCTION TAKES A TGRAPHERRORS AND MODIFIES IT DIRECTLY. THE OPTIONAL STD::VECTOR IS FOR FIXING <P> FOR CERTAIN X-VALUES. */
  void regularizeSlice(TGraphErrors* tgSlice, std::vector<double>* fixedX=0, double omitbelow=0., int nIter_=-1, double threshold_ = -1){
    unsigned int nbins_slice = tgSlice->GetN();
    double* xexyey_slice[4]={
      tgSlice->GetX(),
      tgSlice->GetEX(),
      tgSlice->GetY(),
      tgSlice->GetEY()
    };

    double* xexyey_linear[4];
    for (unsigned int ix=0; ix<4; ix++){
      xexyey_linear[ix] = new double[nbins_slice];
      for (unsigned int iy=0; iy<nbins_slice; iy++){
        if (ix<2) xexyey_linear[ix][iy] = xexyey_slice[ix][iy];
        else if (ix==3) xexyey_linear[ix][iy] = exp(xexyey_slice[ix][iy]);
        else xexyey_linear[ix][iy] = xexyey_slice[ix][iy]*xexyey_linear[ix-1][iy];
      }
    }
    TGraphErrors* tgSlice_linear = new TGraphErrors(nbins_slice, xexyey_linear[0], xexyey_linear[2], xexyey_linear[1], xexyey_linear[3]);
    double integral_in=tgSlice_linear->Integral();
    delete tgSlice_linear;
    for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_linear[ix];


    double* xexyey_mod[4];
    for (unsigned int ix=0; ix<4; ix++){
      xexyey_mod[ix] = new double[nbins_slice];
      for (unsigned int iy=0; iy<nbins_slice; iy++) xexyey_mod[ix][iy] = xexyey_slice[ix][iy];
    }
    unsigned int bin_first = 2, bin_last = nbins_slice-1;

    std::vector<int> fixedBins;
    if (fixedX!=0){
      for (unsigned int ifx=0; ifx<fixedX->size(); ifx++){
        double requestedVal = fixedX->at(ifx);
        double distance=1e15;
        int bin_to_fix=-1;
        for (unsigned int bin=0; bin<nbins_slice; bin++){ if (distance>fabs(xexyey_mod[0][bin]-requestedVal)){ bin_to_fix = bin; distance = fabs(xexyey_mod[0][bin]-requestedVal); } }
        if (bin_to_fix>=0) fixedBins.push_back(bin_to_fix);
        cout << "Requested to fix bin " << bin_to_fix << endl;
      }
    }
    if (omitbelow>0.){
      for (unsigned int bin=0; bin<nbins_slice; bin++){
        if (xexyey_mod[0][bin]<omitbelow) fixedBins.push_back(bin);
        cout << "Requested to fix bin " << bin << endl;
      }
    }

    double* xx_second;
    double* yy_second;

    int nIter = (nIter_<0 ? 1000 : nIter_);
    for (int it=0; it<nIter; it++){
      double threshold = (threshold_<0. ? 0.01 : threshold_);
      for (unsigned int binIt = bin_first; binIt<=bin_last; binIt++){
        bool doFix=false;
        for (unsigned int ifx=0; ifx<fixedBins.size(); ifx++){
          if ((int)(binIt-1)==fixedBins.at(ifx)){ doFix=true; /*cout << "Iteration " << it << " is fixing bin " << (binIt-1) << endl; */break; }
        }
        if (doFix) continue;

        int ctr = 0;
        int nbins_second = nbins_slice-1;
        xx_second = new double[nbins_second];
        yy_second = new double[nbins_second];
        for (unsigned int bin = 1; bin<=nbins_slice; bin++){
          if (bin==binIt) continue;
          xx_second[ctr] = xexyey_mod[0][bin-1];
          yy_second[ctr] = xexyey_mod[2][bin-1];
          ctr++;
        }

        TGraph* interpolator = new TGraph(nbins_second, xx_second, yy_second);
        double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
        double derivative_last = (yy_second[nbins_second-1]-yy_second[nbins_second-2])/(xx_second[nbins_second-1]-xx_second[nbins_second-2]);
        TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

        double center = xexyey_mod[0][binIt-1];
        double val = spline->Eval(center);
        if (fabs(xexyey_mod[2][binIt-1]-val)>threshold*val) xexyey_mod[2][binIt-1]=val;

        delete spline;
        delete interpolator;

        delete[] yy_second;
        delete[] xx_second;
      }
    }

    for (unsigned int ix=0; ix<4; ix++){
      xexyey_linear[ix] = new double[nbins_slice];
      for (unsigned int iy=0; iy<nbins_slice; iy++){
        if (ix<2) xexyey_linear[ix][iy] = xexyey_mod[ix][iy];
        else if (ix==3) xexyey_linear[ix][iy] = exp(xexyey_mod[ix][iy]);
        else xexyey_linear[ix][iy] = xexyey_mod[ix][iy]*xexyey_linear[ix-1][iy];
      }
    }
    tgSlice_linear = new TGraphErrors(nbins_slice, xexyey_linear[0], xexyey_linear[2], xexyey_linear[1], xexyey_linear[3]);
    double integral_out=tgSlice_linear->Integral();
    delete tgSlice_linear;
    for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_linear[ix];

    double scale = integral_out / integral_in;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      xexyey_slice[2][iy] = xexyey_mod[2][iy]*scale;
      xexyey_slice[3][iy] *= xexyey_mod[2][iy]/xexyey_slice[2][iy];
    }

    for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_mod[ix];
  }

  TGraphErrors* removePointsBetween(TGraphErrors* tgSlice, double xmin, double xmax){
    const unsigned int nbins_slice = tgSlice->GetN();
    double* xexyey_slice[4]={
      tgSlice->GetX(),
      tgSlice->GetEX(),
      tgSlice->GetY(),
      tgSlice->GetEY()
    };

    double xexyey[4][nbins_slice];
    unsigned int ctr=0;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (xexyey_slice[0][iy]<=xmax && xexyey_slice[0][iy]>=xmin) continue;
      for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
      cout << "Point " << ctr << " X: " << xexyey[0][ctr] << endl;
      ctr++;
    }
    cout << "removePointsBetween: " << "Starting number of points " << nbins_slice << " final number " << ctr << endl;
    TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tgSlice_new->SetName(tgSlice->GetName());
    tgSlice_new->SetTitle(tgSlice->GetTitle());
    return tgSlice_new;
  }

  TGraphErrors* replacePointsBetween(TGraphErrors* tgSlice, double xmin, double xmax){
    const unsigned int nbins_slice = tgSlice->GetN();
    double* xexyey_slice[4]={
      tgSlice->GetX(),
      tgSlice->GetEX(),
      tgSlice->GetY(),
      tgSlice->GetEY()
    };

    double xexyey[4][nbins_slice];
    unsigned int lowbin=0, highbin=0;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (xexyey_slice[0][iy]<xmin) lowbin=iy;
      if (xexyey_slice[0][iy]>xmax){ highbin=iy; break; }
    }
    cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
    cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (xexyey_slice[0][iy]<=xmax && xexyey_slice[0][iy]>=xmin){
        for (unsigned int ix=0; ix<4; ix++){
          double xlow = xexyey_slice[0][lowbin];
          double xhigh = xexyey_slice[0][highbin];
          double ylow = xexyey_slice[ix][lowbin];
          double yhigh = xexyey_slice[ix][highbin];
          double val;
          if (ix==0) val = xexyey_slice[0][iy];
          else if (ix==2) val = ylow + (yhigh-ylow)/(xhigh-xlow)*(xexyey_slice[0][iy]-xlow);
          else val = xexyey_slice[ix][iy]*(xexyey[ix-1][iy]/xexyey_slice[ix-1][iy]);
          xexyey[ix][iy] = val;
        }
      }
      else{
        for (unsigned int ix=0; ix<4; ix++) xexyey[ix][iy] = xexyey_slice[ix][iy];
      }
      cout << "Point " << iy << " X: " << xexyey[0][iy] << endl;
    }
    TGraphErrors* tgSlice_new = new TGraphErrors(nbins_slice, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tgSlice_new->SetName(tgSlice->GetName());
    tgSlice_new->SetTitle(tgSlice->GetTitle());
    return tgSlice_new;
  }

  TGraphErrors* addPoint(TGraphErrors* tgSlice, double x){
    const unsigned int nbins_slice = tgSlice->GetN();
    double* xexyey_slice[4]={
      tgSlice->GetX(),
      tgSlice->GetEX(),
      tgSlice->GetY(),
      tgSlice->GetEY()
    };

    double xexyey[4][nbins_slice+1];
    unsigned int lowbin=0, highbin=0;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (xexyey_slice[0][iy]<x) lowbin=iy;
      if (xexyey_slice[0][iy]>x){ highbin=iy; break; }
    }
    cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
    cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

    int ctr=0;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
      ctr++;
      if (iy==lowbin){
        for (unsigned int ix=0; ix<4; ix++){
          double ylow = xexyey_slice[ix][lowbin];
          double yhigh = xexyey_slice[ix][highbin];
          double val = (yhigh+ylow)*0.5;
          xexyey[ix][ctr] = val;
        }
        ctr++;
      }
    }
    TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tgSlice_new->SetName(tgSlice->GetName());
    tgSlice_new->SetTitle(tgSlice->GetTitle());
    return tgSlice_new;
  }
  void addPoint(TGraph*& tg, double x, double y){
    TString strname = tg->GetName();
    TString strtitle = tg->GetTitle();
    TString strxtitle = tg->GetXaxis()->GetTitle();
    TString strytitle = tg->GetYaxis()->GetTitle();

    bool hasErrors = false;

    vector<double> xarray;
    vector<double> yarray;
    xarray.push_back(x);
    yarray.push_back(y);
    for (int ip=0; ip<tg->GetN(); ip++){
      if (tg->GetX()[ip]!=x){
        xarray.push_back(tg->GetX()[ip]);
        yarray.push_back(tg->GetY()[ip]);
      }
    }
    vector<pair<double, int>> xorder;
    for (unsigned int ip=0; ip<xarray.size(); ip++) addByLowest<double, int>(xorder, xarray.at(ip), ip);

    double* xynew[2];
    for (unsigned int i=0; i<2; i++) xynew[i] = new double[xorder.size()];
    for (unsigned int ip=0; ip<xarray.size(); ip++){
      unsigned int pos = xorder[ip].second;
      xynew[0][ip] = xarray[pos];
      xynew[1][ip] = yarray[pos];
    }

    delete tg;

    tg = new TGraph(xorder.size(), xynew[0], xynew[1]);
    tg->SetName(strname);
    tg->SetTitle(strtitle);
    tg->GetXaxis()->SetTitle(strxtitle);
    tg->GetYaxis()->SetTitle(strytitle);
    for (unsigned int i=0; i<2; i++) delete[] xynew[i];
  }
  void addPoint(TGraphErrors*& tg, double x, double y, double ex, double ey){
    TString strname = tg->GetName();
    TString strtitle = tg->GetTitle();
    TString strxtitle = tg->GetXaxis()->GetTitle();
    TString strytitle = tg->GetYaxis()->GetTitle();

    vector<double> xarray;
    vector<double> yarray;
    vector<double> exarray;
    vector<double> eyarray;
    xarray.push_back(x);
    yarray.push_back(y);
    exarray.push_back(ex);
    eyarray.push_back(ey);
    for (int ip=0; ip<tg->GetN(); ip++){
      if (tg->GetX()[ip]!=x){
        xarray.push_back(tg->GetX()[ip]);
        yarray.push_back(tg->GetY()[ip]);
        exarray.push_back(tg->GetEX()[ip]);
        eyarray.push_back(tg->GetEY()[ip]);
      }
    }
    vector<pair<double, int>> xorder;
    for (unsigned int ip=0; ip<xarray.size(); ip++) addByLowest<double, int>(xorder, xarray.at(ip), ip);

    double* xynew[4];
    for (unsigned int i=0; i<4; i++) xynew[i] = new double[xorder.size()];
    for (unsigned int ip=0; ip<xarray.size(); ip++){
      unsigned int pos = xorder[ip].second;
      xynew[0][ip] = xarray[pos];
      xynew[1][ip] = yarray[pos];
      xynew[2][ip] = exarray[pos];
      xynew[3][ip] = eyarray[pos];
    }

    delete tg;

    tg = new TGraphErrors(xorder.size(), xynew[0], xynew[1], xynew[2], xynew[3]);
    tg->SetName(strname);
    tg->SetTitle(strtitle);
    tg->GetXaxis()->SetTitle(strxtitle);
    tg->GetYaxis()->SetTitle(strytitle);
    for (unsigned int i=0; i<4; i++) delete[] xynew[i];
  }

  TGraphErrors* addPointAfterBin(TGraphErrors* tgSlice, int abin){
    const unsigned int nbins_slice = tgSlice->GetN();
    double* xexyey_slice[4]={
      tgSlice->GetX(),
      tgSlice->GetEX(),
      tgSlice->GetY(),
      tgSlice->GetEY()
    };

    double xexyey[4][nbins_slice+1];
    unsigned int lowbin=abin, highbin=abin+1;
    cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
    cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

    int ctr=0;
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
      ctr++;
      if (iy==lowbin){
        for (unsigned int ix=0; ix<4; ix++){
          double ylow = xexyey_slice[ix][lowbin];
          double yhigh = xexyey_slice[ix][highbin];
          double val = (yhigh+ylow)*0.5;
          xexyey[ix][ctr] = val;
        }
        ctr++;
        cout << "Adding additional point at " << xexyey[0][ctr-1] << '\t' << xexyey[2][ctr-1] << endl;
      }
    }
    TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tgSlice_new->SetName(tgSlice->GetName());
    tgSlice_new->SetTitle(tgSlice->GetTitle());
    return tgSlice_new;
  }

}

#endif
