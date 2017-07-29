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

namespace CalcHelpers{
  using namespace std;

  struct SimpleEntry{
    int id;
    float trackingval;
    vector<float> recoval;
    float weight;

    SimpleEntry() : id(0), trackingval(0), weight(0) {}
    SimpleEntry(int id_, float trackingval_, vector<float> recoval_, float weight_=1) : id(id_), trackingval(trackingval_), recoval(recoval_), weight(weight_) {}

    bool operator != (const SimpleEntry& other)const{ return trackingval!=other.trackingval; }
    bool operator == (const SimpleEntry& other)const{ return trackingval==other.trackingval; }
    bool operator > (const SimpleEntry& other)const{ return trackingval>other.trackingval; }
    bool operator >= (const SimpleEntry& other)const{ return trackingval>=other.trackingval; }
    bool operator < (const SimpleEntry& other)const{ return trackingval<other.trackingval; }
    bool operator <= (const SimpleEntry& other)const{ return trackingval<=other.trackingval; }

    static void cropByTrueVal(vector<SimpleEntry>& vec, float minval, float maxval){
      vector<unsigned int> erasepos;
      unsigned int pos=0;
      for (vector<SimpleEntry>::iterator it=vec.begin(); it!=vec.end(); it++){
        if (it->trackingval<minval || it->trackingval>maxval) erasepos.push_back(pos);
        pos++;
      }
      for (int ipos=(int)erasepos.size()-1; ipos>=0; ipos--) vec.erase(vec.begin()+erasepos.at(ipos));
    }
    void print(){
      cout << "Simple entry:" << endl;
      cout << " - Id = " << id << endl;
      cout << " - Weight = " << weight << endl;
      cout << " - Trueval: " << trackingval << endl;
      cout << " - Recoval: ";
      for (auto& v : recoval) cout << v << " ";
      cout << endl;
    }
  };
  struct ExtBin{
    double binlow;
    double binhigh;

    vector<SimpleEntry> collection;

    void addEvent(SimpleEntry evt){
      collection.push_back(evt);
    }
  };

  template<typename T> void appendVector(std::vector<T>& a, std::vector<T>& b){ a.insert(a.end(), b.begin(), b.end()); }

  template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique){
    bool inserted = false;
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it>val || (!unique && *it==val)){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
    if (!inserted) valArray.push_back(val);
  }

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index){
    bool inserted = false;
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=val){
        inserted=true;
        if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
        break;
      }
    }
    if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
  }

  template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false){
    if (consecutive){
      bool inserted = false;
      typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
      typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
      for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
        if ((*it).first>=(*inbegin).first){
          inserted=true;
          if ((*it).second!=(*inbegin).second) valArray.insert(it, inbegin, inend);
          break;
        }
      }
      if (!inserted) appendVector<pair<T, U>>(valArray, inArray);
    }
    else if (!inputordered){
      for (typename std::vector<std::pair<T, U>>::iterator init = inArray.begin(); init<inArray.end(); init++){
        bool inserted = false;
        for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
          if ((*it).first>=(*init).first){
            inserted=true;
            if ((*it).second!=(*init).second) valArray.insert(it, *init);
            break;
          }
        }
        if (!inserted) valArray.push_back(*init);
      }
    }
    else if (inArray.size()>0){
      typename std::vector<std::pair<T, U>>::iterator infirst = inArray.begin();
      typename std::vector<std::pair<T, U>>::iterator inlast = inArray.end()-1;
      typename std::vector<std::pair<T, U>>::iterator valfirst = valArray.begin();
      typename std::vector<std::pair<T, U>>::iterator vallast = valArray.end()-1;
      while ((*valfirst).first<(*infirst).first) valfirst++;
      while ((*vallast).first>=(*inlast).first) vallast--;
      vallast++;
      inlast++;

      for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
        bool inserted = false;
        for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
          if ((*it).first>=(*init).first){
            inserted=true;
            if ((*it).second!=(*init).second) valArray.insert(it, *init);
            break;
          }
        }
        if (!inserted) valArray.insert(vallast, *init);
      }
    }
  }

  template<typename T> bool checkListVariable(const vector<T>& list, const T& var){
    for (unsigned int v=0; v<list.size(); v++){
      if (list.at(v)==var) return true; // Look for exact match
    }
    return false;
  }


  void splitOption(const string rawoption, string& wish, string& value, char delimiter){
    size_t posEq = rawoption.find(delimiter);
    if (posEq!=string::npos){
      wish=rawoption;
      value=rawoption.substr(posEq+1);
      wish.erase(wish.begin()+posEq, wish.end());
    }
    else{
      wish="";
      value=rawoption;
    }
  }
  void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
    string suboption=rawoption, result=rawoption;
    string remnant;
    while (result!=""){
      splitOption(suboption, result, remnant, delimiter);
      if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
      suboption = remnant;
    }
    if (remnant!="") splitoptions.push_back(remnant);
  }


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

  /* SPECIFIC COMMENT: Convert a TGraph to a TSpline3 */
  TSpline3* convertGraphToSpline3(TGraph* tg, double* dfirst=0, double* dlast=0){
    unsigned int nbins = tg->GetN();
    double* xy[2]={
      tg->GetX(),
      tg->GetY()
    };
    TSpline3* spline = new TSpline3("spline", tg, "b2e2", 0, 0);
    spline->SetName(Form("sp_%s", tg->GetName()));
    if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
    if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
    return spline;
  }

  /* SPECIFIC COMMENT: Convert a TGraph to a TSpline3 */
  TSpline3* convertGraphToSpline3_FaithfulSlopes(TGraph* tg, double* dfirst=0, double* dlast=0){
    unsigned int nbins = tg->GetN();
    double* xy[2]={
      tg->GetX(),
      tg->GetY()
    };
    double derivative_first = (xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0]);
    double derivative_last = (xy[1][nbins-1]-xy[1][nbins-2])/(xy[0][nbins-1]-xy[0][nbins-2]);
    cout << "End derivatives to be preserved: " << derivative_first << " , " << derivative_last << endl;
    TSpline3* spline = new TSpline3("spline", tg, "b1e1", derivative_first, derivative_last);
    spline->SetName(Form("sp_%s", tg->GetName()));
    if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
    if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
    return spline;
  }

  /* SPECIFIC COMMENT: Convert a TGraph to a TSpline3 */
  TSpline3* convertGraphToSpline3_MightFaithfulSlopes(TGraph* tg, bool faithfulFirst, bool faithfulSecond, double* dfirst=0, double* dlast=0){
    unsigned int nbins = tg->GetN();
    double* xy[2]={
      tg->GetX(),
      tg->GetY()
    };

    double derivative_first=0;
    double derivative_last=0;
    TString spopt="";
    if (faithfulFirst){
      spopt += "b1";
      derivative_first = (xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0]);
    }
    else spopt += "b2";
    if (faithfulFirst){
      spopt += "e1";
      derivative_last = (xy[1][nbins-1]-xy[1][nbins-2])/(xy[0][nbins-1]-xy[0][nbins-2]);
    }
    else spopt += "e2";
    cout << "End derivatives to be preserved (" << faithfulFirst << " , " << faithfulSecond << "): " << derivative_first << " , " << derivative_last << endl;
    TSpline3* spline = new TSpline3("spline", tg, spopt, derivative_first, derivative_last);
    spline->SetName(Form("sp_%s", tg->GetName()));
    if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
    if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
    return spline;
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

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*exp(x) */
  TF1* getFcn_a0plusa1expX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    a0 = y-s;
    a1 = s*exp(-x);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]*exp(x)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
  }

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0*exp(a1*x) */
  TF1* getFcn_a0timesexpa1X(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    double a0, a1;
    // a0*a1*exp(a1*x)=s a0*exp(a1*x)=y
    a1=s/y;
    a0=y*exp(-a1*x);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]*exp([1]*x)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);

    return fcn;
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

  /* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1/x**N */
  template<int N> TF1* getFcn_a0plusa1overXN(TSpline3* sp, double xmin, double xmax, bool useLowBound){
    double x, y, s;
    if (useLowBound) x = sp->GetXmin();
    else x = sp->GetXmax();
    y = sp->Eval(x);
    s = sp->Derivative(x);

    // y=a0+a1/x^N
    // s=-N*a1/x^(N+1)
    double a0, a1;
    a0 = y+s*x/(double(N));
    a1 = -s*pow(x, N+1)/(double(N));

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, Form("[0]+[1]/pow(x, %i)", N), xmin, xmax);
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

  TGraphErrors* makeGraphFromTH1(TH1* hx, TH1* hy, TString name){
    if (hx->GetNbinsX()!=hy->GetNbinsX()){
      cerr << "Number of bins for x coordinate != those for y" << endl;
      assert(0);
    }
    unsigned int nbins = hx->GetNbinsX();
    double* xexyey[4];
    for (unsigned int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (unsigned int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hx->GetBinContent(bin+1);
      xexyey[1][bin] = hx->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hy->GetBinContent(bin+1);
      xexyey[3][bin] = hy->GetBinError(bin+1);
    }
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(name);
    for (unsigned int ix=0; ix<4; ix++) delete[] xexyey[ix];
    return tg;
  }
  template<typename T> TGraph* makeGraphFromPair(vector<pair<T, T>> points, TString name){
    if (points.size()==0) return 0;
    unsigned int nbins = points.size();
    double* xy[2];
    for (unsigned int ix=0; ix<2; ix++) xy[ix] = new double[nbins];
    for (unsigned int bin=0; bin<nbins; bin++){
      xy[0][bin] = points[bin].first;
      xy[1][bin] = points[bin].second;
    }
    TGraph* tg = new TGraph(nbins, xy[0], xy[1]);
    tg->SetName(name);
    for (unsigned int ix=0; ix<2; ix++) delete[] xy[ix];
    return tg;
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

  void addPointsBetween(TGraph*& tgOriginal, double xmin, double xmax, unsigned int nadd){
    const unsigned int np = tgOriginal->GetN();
    double* xy[2]={
      tgOriginal->GetX(),
      tgOriginal->GetY()
    };

    unsigned int pfirst=0, plast=np-1;
    vector<pair<double, double>> tmpPoints;
    for (unsigned int ip=0; ip<np; ip++){
      if (xy[0][ip]<xmin) pfirst=ip;
      if (xy[0][ip]>xmax){ plast=ip; break; }
    }
    for (unsigned int ip=pfirst; ip<=plast; ip++) tmpPoints.push_back(pair<double, double>(xy[0][ip], xy[1][ip]));

    TGraph* tgtmp = makeGraphFromPair(tmpPoints, "tgtmp");
    TSpline3* sptmp = convertGraphToSpline3(tgtmp);

    double addWidth = (xmax-xmin)/double(nadd+1);
    for (unsigned int it=1; it<=nadd; it++){
      double xnew = xmin + addWidth*double(it);
      double ynew = sptmp->Eval(xnew);
      addPoint(tgOriginal, xnew, ynew);
    }

    delete sptmp;
    delete tgtmp;
  }
  void addPointsBetween(TGraphErrors*& tgOriginal, double xmin, double xmax, unsigned int nadd){
    const unsigned int np = tgOriginal->GetN();
    double* xy[2]={
      tgOriginal->GetX(),
      tgOriginal->GetY()
    };

    unsigned int pfirst=0, plast=np-1;
    vector<pair<double, double>> tmpPoints;
    for (unsigned int ip=0; ip<np; ip++){
      if (xy[0][ip]<xmin) pfirst=ip;
      if (xy[0][ip]>xmax){ plast=ip; break; }
    }
    for (unsigned int ip=pfirst; ip<=plast; ip++) tmpPoints.push_back(pair<double, double>(xy[0][ip], xy[1][ip]));

    TGraph* tgtmp = makeGraphFromPair(tmpPoints, "tgtmp");
    TSpline3* sptmp = convertGraphToSpline3(tgtmp);

    double addWidth = (xmax-xmin)/double(nadd+1);
    for (unsigned int it=1; it<=nadd; it++){
      double xnew = xmin + addWidth*double(it);
      double ynew = sptmp->Eval(xnew);
      addPoint(tgOriginal, xnew, ynew, 0, 0);
    }

    delete sptmp;
    delete tgtmp;
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

  TGraph* multiplyTGraphs(TGraph* tgfirst, TGraph* tgsecond){
    TSpline3* spfirst = convertGraphToSpline3(tgfirst);
    TSpline3* spsecond = convertGraphToSpline3(tgsecond);

    vector<pair<double, double>> xy;
    double* xx = tgfirst->GetX();
    for (int ip=0; ip<tgfirst->GetN(); ip++) addByLowest<double, double>(xy, xx[ip], 0);
    xx = tgsecond->GetX();
    for (int ip=0; ip<tgsecond->GetN(); ip++) addByLowest<double, double>(xy, xx[ip], 0);

    vector<pair<double, double>> xynew;
    for (unsigned int ip=0; ip<xy.size(); ip++){
      double xval = xy.at(ip).first;
      double yval = spfirst->Eval(xval);
      yval *= spsecond->Eval(xval);
      xynew.push_back(pair<double, double>(xval, yval));
    }

    return makeGraphFromPair(xynew, Form("%s_x_%s", tgfirst->GetName(), tgsecond->GetName()));
  }
  TGraph* divideTGraphs(TGraph* tgnum, TGraph* tgdenom, double powernum=1, double powerdenom=1){
    TSpline3* spnum = convertGraphToSpline3(tgnum);
    TSpline3* spdenom = convertGraphToSpline3(tgdenom);

    vector<pair<double, double>> xy;
    double* xxnum = tgnum->GetX();
    for (int ip=0; ip<tgnum->GetN(); ip++) addByLowest<double, double>(xy, xxnum[ip], 0);
    double* xxdenom = tgdenom->GetX();
    for (int ip=0; ip<tgdenom->GetN(); ip++) addByLowest<double, double>(xy, xxdenom[ip], 0);

    vector<pair<double, double>> xynew;
    for (unsigned int ip=0; ip<xy.size(); ip++){
      double xval = xy.at(ip).first;
      double numerator = spnum->Eval(xval);
      double denominator = spdenom->Eval(xval);
      if (denominator!=0.){
        double result = pow(numerator, powernum)/pow(denominator, powerdenom);
        if (std::isnan(result) || std::isinf(result)){
          cerr << "Division gave NaN! ";
          cerr << "Numerator = " << numerator << ", denominator = " << denominator << " at " << xval << endl;
          cerr << " - TGraph numerator = " << tgnum->Eval(xval) << ", denominator = " << tgdenom->Eval(xval) << " at " << xval << endl;
        }
        else xynew.push_back(pair<double, double>(xval, result));
      }
    }
    return makeGraphFromPair(xynew, Form("%s_over_%s", tgnum->GetName(), tgdenom->GetName()));
  }


  // Explicit instantiations
  template void addByLowest<SimpleEntry>(std::vector<SimpleEntry>& valArray, SimpleEntry val, bool unique);
  template void addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, double val, int index);
  template void addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, double val, double index);
  template void addByLowest<double, int>(std::vector<std::pair<double, int>>& valArray, std::vector<std::pair<double, int>>& inArray, bool consecutive, bool inputordered);
  template void addByLowest<double, double>(std::vector<std::pair<double, double>>& valArray, std::vector<std::pair<double, double>>& inArray, bool consecutive, bool inputordered);
  template bool checkListVariable<string>(const vector<string>& list, const string& var);
  template bool checkListVariable<double>(const vector<double>& list, const double& var);
  template TGraph* makeGraphFromPair<float>(vector<pair<float, float>> points, TString name);
  template TGraph* makeGraphFromPair<double>(vector<pair<double, double>> points, TString name);

  template TF1* getFcn_a0plusa1overXN<1>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<2>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<3>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<4>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<5>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<6>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<7>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<8>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<9>(TSpline3* sp, double xmin, double xmax, bool useLowBound);
  template TF1* getFcn_a0plusa1overXN<10>(TSpline3* sp, double xmin, double xmax, bool useLowBound);

}

#endif
