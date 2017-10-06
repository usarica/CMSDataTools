#include "HelperFunctions.h"


using namespace std;


void HelperFunctions::splitOption(const std::string rawoption, std::string& wish, std::string& value, char delimiter){
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
void HelperFunctions::splitOptionRecursive(const std::string rawoption, std::vector<std::string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*exp(x) */
TF1* HelperFunctions::getFcn_a0plusa1expX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
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
TF1* HelperFunctions::getFcn_a0timesexpa1X(TSpline3* sp, double xmin, double xmax, bool useLowBound){
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


TSpline3* HelperFunctions::convertGraphToSpline3(TGraph* tg, bool faithfulFirst, bool faithfulSecond, double* dfirst, double* dlast){
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
  //cout << "End derivatives to be preserved (" << faithfulFirst << " , " << faithfulSecond << "): " << derivative_first << " , " << derivative_last << endl;
  TSpline3* spline = new TSpline3("spline", tg, spopt, derivative_first, derivative_last);
  spline->SetName(Form("sp_%s", tg->GetName()));
  if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
  if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
  return spline;
}

TGraphErrors* HelperFunctions::makeGraphFromTH1(TH1* hx, TH1* hy, TString name){
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

TGraph* HelperFunctions::multiplyTGraphs(TGraph* tgfirst, TGraph* tgsecond){
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
TGraph* HelperFunctions::divideTGraphs(TGraph* tgnum, TGraph* tgdenom, double powernum, double powerdenom){
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

TGraphErrors* HelperFunctions::addPoint(TGraphErrors* tgSlice, double x){
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
void HelperFunctions::addPoint(TGraph*& tg, double x, double y){
  TString strname = tg->GetName();
  TString strtitle = tg->GetTitle();
  TString strxtitle = tg->GetXaxis()->GetTitle();
  TString strytitle = tg->GetYaxis()->GetTitle();

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
void HelperFunctions::addPoint(TGraphErrors*& tg, double x, double y, double ex, double ey){
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

template<> void HelperFunctions::addPointsBetween<TGraph>(TGraph*& tgOriginal, double xmin, double xmax, unsigned int nadd){
  const unsigned int np = tgOriginal->GetN();
  double* xy[2]={
    tgOriginal->GetX(),
    tgOriginal->GetY()
  };

  unsigned int pfirst=0, plast=np-1;
  std::vector<std::pair<double, double>> tmpPoints;
  for (unsigned int ip=0; ip<np; ip++){
    if (xy[0][ip]<xmin) pfirst=ip;
    if (xy[0][ip]>xmax){ plast=ip; break; }
  }
  for (unsigned int ip=pfirst; ip<=plast; ip++) tmpPoints.push_back(std::pair<double, double>(xy[0][ip], xy[1][ip]));

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
template<> void HelperFunctions::addPointsBetween<TGraphErrors>(TGraphErrors*& tgOriginal, double xmin, double xmax, unsigned int nadd){
  const unsigned int np = tgOriginal->GetN();
  double* xy[2]={
    tgOriginal->GetX(),
    tgOriginal->GetY()
  };

  unsigned int pfirst=0, plast=np-1;
  std::vector<std::pair<double, double>> tmpPoints;
  for (unsigned int ip=0; ip<np; ip++){
    if (xy[0][ip]<xmin) pfirst=ip;
    if (xy[0][ip]>xmax){ plast=ip; break; }
  }
  for (unsigned int ip=pfirst; ip<=plast; ip++) tmpPoints.push_back(std::pair<double, double>(xy[0][ip], xy[1][ip]));

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
