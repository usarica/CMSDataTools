#include <cstring>
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


template<> void HelperFunctions::castStringToValue(std::string const name, bool& val){
  std::string namelower=name;
  std::transform(namelower.begin(), namelower.end(), namelower.begin(), ::tolower);
  if (namelower=="true" || namelower=="t") val=1;
  else if (namelower=="false" || namelower=="f") val=0;
  else{ std::stringstream ss(name); ss >> val; }
}

TString HelperFunctions::todaysdate(){
  TString result;
  time_t now = time(0);
  tm* ltm = localtime(&now);
  int year = 1900 + ltm->tm_year-2000;
  int month = 1 + ltm->tm_mon;
  int day = ltm->tm_mday;
  TString strday = (day>=10 ? Form("%i", day) : Form("%i%i", 0, day));
  TString strmonth = (month>=10 ? Form("%i", month) : Form("%i%i", 0, month));
  TString stryear = Form("%i", year);
  result = stryear + strmonth + strday;
  return result;
}

void HelperFunctions::progressbar(unsigned int val, unsigned int tot){
  unsigned int percent=std::ceil(0.01*tot);
  tot=(tot==0 ? tot : tot-1);
  if (val%percent==0 || val==tot){
    unsigned int percent_done = val/percent;
    if (val==tot) percent_done=100;
    cout << "[ " << setw(3) << percent_done << "% | ";
    for (unsigned int k=0; k<percent_done; k++) cout << '=';
    if (percent_done<100) cout << '>';
    else cout << ' ';
    for (unsigned int k=percent_done; k<100; k++) cout << ' ';
    cout << "| ]";
    cout << flush;
    if (val==tot) cout << endl;
    else cout << '\r';
  }
}

bool HelperFunctions::test_bit(int mask, unsigned int iBit){ return (mask >> iBit) & 1; }

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

void HelperFunctions::splitOption(const TString rawoption, TString& wish, TString& value, char delimiter){
  std::string const s_rawoption = rawoption.Data();
  std::string s_wish, s_value;
  splitOption(s_rawoption, s_wish, s_value, delimiter);
  wish=s_wish.c_str();
  value=s_value.c_str();
}
void HelperFunctions::splitOptionRecursive(const TString rawoption, std::vector<TString>& splitoptions, char delimiter){
  std::string const s_rawoption = rawoption.Data();
  std::vector<std::string> s_splitoptions;
  splitOptionRecursive(s_rawoption, s_splitoptions, delimiter);
  for (std::string const& s:s_splitoptions) splitoptions.push_back(s.c_str());
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

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*exp(x/a2) */
TF1* HelperFunctions::getFcn_a0plusa1timesexpXovera2(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, xp, dx, y, b, c, d;
  double s, ss;
  if (useLowBound){
    x = sp->GetXmin();
    sp->GetCoeff(0, xp, y, b, c, d); // Hopefully xp=x; we don't care about y.
  }
  else{
    x = sp->GetXmax();
    sp->GetCoeff(sp->GetNp()-2, xp, y, b, c, d); // We don't care about y.
  }
  dx=x-xp;
  y = sp->Eval(x);
  s = sp->Derivative(x); // == b + 2.*c*dx + 3.*d*dx*dx
  ss = 2.*c + 6.*d*dx;

  cout << "s = " << s << " =? b + 2.*c*dx + 3.*d*dx*dx = " << b + 2.*c*dx + 3.*d*dx*dx << endl;
  cout << "b = " << b << ", c = " << c << ", d = " << d << endl;
  cout << "x = " << x << ", xp = " << xp << endl;

  double a0, a1, a2;
  // y=a0+a1*exp(x/a2) => s=a1/a2*exp(x/a2), ss=a1/pow(a2, 2)*exp(x/a2)
  TF1* fcn=nullptr;
  if (ss!=0.){
    a2=s/ss;
    a1=s*a2*exp(-x/a2);
    a0=y-a1*exp(-x/a2);

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    fcn = new TF1(fcnName, "[0]+[1]*exp(x/[2])", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);
    fcn->SetParameter(2, a2);
  }
  else{
    MELAout << "HelperFunctions::getFcn_a0plusa1timesexpXovera2: Second derivative was 0! Forcing constant interpolation." << endl;
    a0=y;

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    fcn = new TF1(fcnName, "[0]", xmin, xmax);
    fcn->SetParameter(0, a0);
  }
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
  TSpline3* spline = new TSpline3("spline", tg, spopt, derivative_first, derivative_last);
  spline->SetName(Form("sp_%s", tg->GetName()));
  if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
  if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
  return spline;
}

void HelperFunctions::convertTGraphErrorsToTH1F(TGraphErrors* tg, TH1F* histo){
  if (tg==0){
    MELAerr << "convertTGraphErrorsToTH1F: TGraph is 0!" << endl;
    histo=0;
    return;
  }

  double* xx = tg->GetX();
  double* yy = tg->GetY();
  double* ex = tg->GetEX();
  double* ey = tg->GetEY();
  // This is just stupid.
  /*
  double* ex_up = tg->GetEXhigh();
  double* ex_dn = tg->GetEXlow();
  double* ey_up = tg->GetEYhigh();
  double* ey_dn = tg->GetEYlow();
  */
  const int nbins = tg->GetN();
  double xarray[nbins+1];
  /*
  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex_dn[ix];
  xarray[nbins] = xx[nbins-1]+ex_up[nbins-1];
  */
  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex[ix];
  xarray[nbins] = xx[nbins-1]+ex[nbins-1];
  gROOT->cd();
  histo = new TH1F(Form("h_%s", tg->GetName()), "", nbins, xarray);
  for (int ix=1; ix<=nbins; ix++){
    //MELAout << "x, y = " << xarray[ix-1] << " " << yy[ix-1] << endl;
    histo->SetBinContent(ix, yy[ix-1]);
    //histo->SetBinError(ix, sqrt((pow(ey_dn[ix-1], 2)+pow(ey_up[ix-1], 2))*0.5));
    histo->SetBinError(ix, ey[ix-1]);
  }
}
void HelperFunctions::convertTGraphAsymmErrorsToTH1F(TGraphAsymmErrors* tg, TH1F* histo){
  if (tg==0){
    MELAerr << "convertTGraphAsymmErrorsToTH1F: TGraph is 0!" << endl;
    histo=0;
    return;
  }

  double* xx = tg->GetX();
  double* yy = tg->GetY();
  //double* ex = tg->GetEX();
  //double* ey = tg->GetEY();
  double* ex_up = tg->GetEXhigh();
  double* ex_dn = tg->GetEXlow();
  double* ey_up = tg->GetEYhigh();
  double* ey_dn = tg->GetEYlow();
  const int nbins = tg->GetN();
  double xarray[nbins+1];

  for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex_dn[ix];
  xarray[nbins] = xx[nbins-1]+ex_up[nbins-1];

  //for (int ix=0; ix<nbins; ix++) xarray[ix] = xx[ix]-ex[ix];
  //xarray[nbins] = xx[nbins-1]+ex[nbins-1];
  gROOT->cd();
  histo = new TH1F(Form("h_%s", tg->GetName()), "", nbins, xarray);
  for (int ix=1; ix<=nbins; ix++){
    //MELAout << "x, y = " << xarray[ix-1] << " " << yy[ix-1] << endl;
    histo->SetBinContent(ix, yy[ix-1]);
    histo->SetBinError(ix, sqrt((pow(ey_dn[ix-1], 2)+pow(ey_up[ix-1], 2))*0.5));
    //histo->SetBinError(ix, ey[ix-1]);
  }
}

TGraph* HelperFunctions::createROCFromDistributions(TH1* hA, TH1* hB, TString name){
  if (!hA || !hB) return nullptr;
  assert(hA->GetNbinsX()==hB->GetNbinsX());
  const int nbins=hA->GetNbinsX();
  double integral[2]={ hA->Integral(0, nbins), hB->Integral(0, nbins) };
  vector<pair<float, float>> sumWgtsPerBin; sumWgtsPerBin.assign(nbins, pair<float, float>(0, 0));
  for (int ix=1; ix<=nbins; ix++){
    for (int jx=ix; jx<=nbins; jx++){
      sumWgtsPerBin.at(ix-1).second += hA->GetBinContent(jx)/integral[0];
      sumWgtsPerBin.at(ix-1).first += hB->GetBinContent(jx)/integral[1];
    }
  }
  TGraph* tg=HelperFunctions::makeGraphFromPair(sumWgtsPerBin, name);
  HelperFunctions::addPoint(tg, 0, 0);
  tg->GetYaxis()->SetRangeUser(0, 1);
  tg->GetYaxis()->SetTitle(TString(hA->GetTitle())+" eff.");
  tg->GetXaxis()->SetRangeUser(0, 1);
  tg->GetXaxis()->SetTitle(TString(hB->GetTitle())+" eff.");
  return tg;
}

TGraphErrors* HelperFunctions::makeGraphFromTH1(TH1* hx, TH1* hy, TString name){
  if (!hy) return nullptr;
  if (hx && hx->GetNbinsX()!=hy->GetNbinsX()){
    MELAerr << "Number of bins for x coordinate != those for y" << endl;
    assert(0);
  }
  unsigned int nbins = hy->GetNbinsX();
  double* xexyey[4];
  for (unsigned int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    if (hx){
      xexyey[0][bin] = hx->GetBinContent(bin+1);
      xexyey[1][bin] = hx->GetBinError(bin+1);
    }
    else{
      xexyey[0][bin] = hy->GetBinCenter(bin+1);
      xexyey[1][bin] = 0;
    }

    //MELAout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
    xexyey[2][bin] = hy->GetBinContent(bin+1);
    xexyey[3][bin] = hy->GetBinError(bin+1);
  }
  TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey[ix];
  return tg;
}

TGraph* HelperFunctions::addTGraphs(TGraph* tgfirst, TGraph* tgsecond){
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
    yval += spsecond->Eval(xval);
    xynew.push_back(pair<double, double>(xval, yval));
  }

  return makeGraphFromPair(xynew, Form("%s_plus_%s", tgfirst->GetName(), tgsecond->GetName()));
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
        MELAerr << "Division gave NaN! ";
        MELAerr << "Numerator = " << numerator << ", denominator = " << denominator << " at " << xval << endl;
        MELAerr << " - TGraph numerator = " << tgnum->Eval(xval) << ", denominator = " << tgdenom->Eval(xval) << " at " << xval << endl;
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
  //MELAout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
  //MELAout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

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

TGraph* HelperFunctions::genericPatcher(
  TGraph* tg, const TString newname,
  double const xmin, double const xmax,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  bool useFaithfulSlopeFirst, bool useFaithfulSlopeSecond,
  vector<pair<pair<double, double>, unsigned int>>* addpoints
){
  if (addpoints!=0){ for (auto& prange : *addpoints) addPointsBetween(tg, prange.first.first, prange.first.second, prange.second); }

  int n = tg->GetN();
  double* xx = tg->GetX();
  double* yy = tg->GetY();

  TSpline3* sp = convertGraphToSpline3(tg, useFaithfulSlopeFirst, useFaithfulSlopeSecond);
  double tglow = xx[0];
  double tghigh = xx[tg->GetN()-1];
  TF1* lowFcn = lowf(sp, xmin, tglow, true);
  TF1* highFcn = highf(sp, tghigh, xmax, false);
  lowFcn->SetNpx((int) (tglow-xmin)*5);
  highFcn->SetNpx((int) (xmax-tghigh)*5);
  delete sp;

  vector<pair<double, double>> points;
  for (double xval=xmin; xval<tglow; xval+=1){
    double yval = lowFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }
  delete lowFcn;

  for (int ix=0; ix<n; ix++) addByLowest<double, double>(points, xx[ix], yy[ix]);

  int tghigh_int = ((int) ((tghigh+1.)/100.+0.5))*100;
  if (tghigh>=(double) tghigh_int) tghigh_int+=100;
  for (double xval=tghigh_int; xval<=xmax; xval+=100){
    double yval = highFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }
  delete highFcn;

  TGraph* res = makeGraphFromPair(points, newname);
  return res;
}


void HelperFunctions::regularizeSlice(TGraph* tgSlice, std::vector<double>* fixedX, double omitbelow, double omitabove, int nIter_, double threshold_){
  unsigned int nbins_slice = tgSlice->GetN();
  double* xy_slice[2]={
    tgSlice->GetX(),
    tgSlice->GetY()
  };

  double* xy_mod[2];
  for (unsigned int ix=0; ix<2; ix++){
    xy_mod[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++) xy_mod[ix][iy] = xy_slice[ix][iy];
  }
  unsigned int bin_first = 2, bin_last = nbins_slice-1;

  std::vector<int> fixedBins;
  if (fixedX!=0){
    for (unsigned int ifx=0; ifx<fixedX->size(); ifx++){
      double requestedVal = fixedX->at(ifx);
      double distance=1e15;
      int bin_to_fix=-1;
      for (unsigned int bin=0; bin<nbins_slice; bin++){ if (distance>fabs(xy_mod[0][bin]-requestedVal)){ bin_to_fix = bin; distance = fabs(xy_mod[0][bin]-requestedVal); } }
      if (bin_to_fix>=0) fixedBins.push_back(bin_to_fix);
      //MELAout << "Requested to fix bin " << bin_to_fix << endl;
    }
  }
  for (unsigned int bin=0; bin<nbins_slice; bin++){
    if (xy_mod[0][bin]<omitbelow) fixedBins.push_back(bin);
    //MELAout << "Requested to fix bin " << bin << endl;
  }
  for (unsigned int bin=0; bin<nbins_slice; bin++){
    if (xy_mod[0][bin]>omitabove) fixedBins.push_back(bin);
    //MELAout << "Requested to fix bin " << bin << endl;
  }

  double* xx_second;
  double* yy_second;

  int nIter = (nIter_<0 ? 1000 : nIter_);
  for (int it=0; it<nIter; it++){
    double threshold = (threshold_<0. ? 0.01 : threshold_);
    for (unsigned int binIt = bin_first; binIt<=bin_last; binIt++){
      bool doFix=false;
      for (unsigned int ifx=0; ifx<fixedBins.size(); ifx++){
        if ((int) (binIt-1)==fixedBins.at(ifx)){ doFix=true; /*MELAout << "Iteration " << it << " is fixing bin " << (binIt-1) << endl; */break; }
      }
      if (doFix) continue;

      int ctr = 0;
      int nbins_second = nbins_slice-1;
      xx_second = new double[nbins_second];
      yy_second = new double[nbins_second];
      for (unsigned int bin = 1; bin<=nbins_slice; bin++){
        if (bin==binIt) continue;
        xx_second[ctr] = xy_mod[0][bin-1];
        yy_second[ctr] = xy_mod[1][bin-1];
        ctr++;
      }

      TGraph* interpolator = new TGraph(nbins_second, xx_second, yy_second);
      //double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
      //double derivative_last = (yy_second[nbins_second-1]-yy_second[nbins_second-2])/(xx_second[nbins_second-1]-xx_second[nbins_second-2]);
      //TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);
      TSpline3* spline = new TSpline3("spline", interpolator, "b2e2", 0, 0);

      double center = xy_mod[0][binIt-1];
      double val = spline->Eval(center);
      if (val!=0. && fabs(xy_mod[1][binIt-1]/val-1.)>threshold) xy_mod[1][binIt-1]=val;

      delete spline;
      delete interpolator;

      delete[] yy_second;
      delete[] xx_second;
    }
  }

  for (unsigned int iy=0; iy<nbins_slice; iy++) xy_slice[1][iy] = xy_mod[1][iy];
  for (unsigned int ix=0; ix<2; ix++) delete[] xy_mod[ix];
}
void HelperFunctions::regularizeSlice(TGraphErrors* tgSlice, std::vector<double>* fixedX, double omitbelow, double omitabove, int nIter_, double threshold_, double acceleration_){
  const unsigned int ndim=4;
  unsigned int nbins_slice = tgSlice->GetN();
  double* xy_slice[ndim]={
    tgSlice->GetX(),
    tgSlice->GetY(),
    tgSlice->GetEX(),
    tgSlice->GetEY()
  };

  double* xy_mod[ndim];
  for (unsigned int ix=0; ix<ndim; ix++){
    xy_mod[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++) xy_mod[ix][iy] = xy_slice[ix][iy];
  }
  unsigned int bin_first = 0, bin_last = nbins_slice-1;

  std::vector<int> fixedBins;
  if (fixedX){
    for (unsigned int ifx=0; ifx<fixedX->size(); ifx++){
      double requestedVal = fixedX->at(ifx);
      double distance=1e15;
      int bin_to_fix=-1;
      for (unsigned int bin=0; bin<nbins_slice; bin++){ if (distance>fabs(xy_mod[0][bin]-requestedVal)){ bin_to_fix = bin; distance = fabs(xy_mod[0][bin]-requestedVal); } }
      if (bin_to_fix>=0) fixedBins.push_back(bin_to_fix);
    }
  }
  for (unsigned int bin=0; bin<nbins_slice; bin++){ if (xy_mod[0][bin]<omitbelow) fixedBins.push_back(bin); }
  for (unsigned int bin=0; bin<nbins_slice; bin++){ if (xy_mod[0][bin]>omitabove) fixedBins.push_back(bin); }

  double* xx_second;
  double* yy_second;

  int nIter = (nIter_<0 ? 1000 : nIter_);
  double threshold = (threshold_<0. ? 0.5 : threshold_);
  double acceleration = (acceleration_<0. ? 0.5 : acceleration_);
  for (int it=0; it<nIter; it++){
    vector<double> residuals;

    // First loop: Collect unbiased residuals
    unsigned int binIt = bin_first;
    while (binIt<=bin_last){
      double residual=0;
      bool doFix=false;
      for (unsigned int ifx=0; ifx<fixedBins.size(); ifx++){
        if ((int) binIt==fixedBins.at(ifx)){ doFix=true; break; }
      }
      if (!doFix){
        int ctr = 0;
        int nbins_second = nbins_slice-1;
        xx_second = new double[nbins_second];
        yy_second = new double[nbins_second];
        for (unsigned int bin = 0; bin<nbins_slice; bin++){
          if (bin==binIt) continue;
          xx_second[ctr] = xy_mod[0][bin];
          yy_second[ctr] = xy_mod[1][bin];
          ctr++;
        }

        TGraph* interpolator = new TGraph(nbins_second, xx_second, yy_second);
        TSpline3* spline = convertGraphToSpline3(interpolator, (binIt!=bin_first), (binIt!=bin_last)); spline->SetName("tmpspline");

        double center = xy_mod[0][binIt];
        double val;
        if (binIt!=bin_first && binIt!=bin_last) val = spline->Eval(center);
        else if (binIt==bin_first){
          double deriv = spline->Derivative(xx_second[0]);
          val = spline->Eval(xx_second[0]) + deriv*(center-xx_second[0]);
        }
        else{
          double deriv = spline->Derivative(xx_second[nbins_second-1]);
          val = spline->Eval(xx_second[nbins_second-1]) + deriv*(center-xx_second[nbins_second-1]);
        }
        residual = xy_mod[1][binIt]-val;

        delete spline;
        delete interpolator;
        delete[] yy_second;
        delete[] xx_second;
      }
      residuals.push_back(residual);
      binIt++;
    }

    // Second loop: Shift values according to the unbiased residuals
    binIt = bin_first;
    for (double const& residual:residuals){
      if (fabs(residual)>threshold*xy_mod[3][binIt]){
        double valMove = -residual*acceleration;
        if (xy_mod[3][binIt]>0. && fabs(valMove)>xy_mod[3][binIt]) valMove = xy_mod[3][binIt] * (valMove<0. ? double(-1) : double(1)); // Slow down the movement to allow better convergence
        double newval = xy_mod[1][binIt] + valMove;
        if (newval<=0.) newval = xy_mod[1][binIt];
        xy_mod[3][binIt] *= newval / xy_mod[1][binIt];
        xy_mod[1][binIt] = newval;
      }
      binIt++;
    }
  }

  for (unsigned int iy=0; iy<nbins_slice; iy++){
    xy_slice[1][iy] = xy_mod[1][iy];
    xy_slice[3][iy] = xy_mod[3][iy];
  }
  for (unsigned int ix=0; ix<ndim; ix++) delete[] xy_mod[ix];
}

float HelperFunctions::calculateEfficiencyError(
  float const sumW, float const sumWAll,
  float const sumWsq, float const sumWsqAll
  ){
  float const& sumWp=sumW;
  float const& sumWsqp=sumWsq;
  float const sumWm = sumWAll-sumWp;
  float const sumWsqm = sumWsqAll-sumWsqp;
  float numerator, denominator;
  float ratio=0;
  if (sumWAll!=0.){
    numerator = sqrt(sumWsqp*pow(sumWm, 2) + sumWsqm*pow(sumWp, 2));
    denominator = pow(sumWAll, 2);
    ratio = numerator/denominator;
  }
  return ratio;
}


template<> void HelperFunctions::replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1) strinput.Replace(ipos, strTakeOut.Length(), strPutIn);
}
template<> void HelperFunctions::replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1) strinput.Replace(ipos, strlen(strTakeOut), strPutIn);
}
template<> void HelperFunctions::replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos) strinput.replace(ipos, strTakeOut.length(), strPutIn);
}
template<> void HelperFunctions::replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos) strinput.replace(ipos, strlen(strTakeOut), strPutIn);
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

  MELAout << "HelperFunctions::addPointsBetween: Adding " << nadd << " points between " << xy[0][pfirst] << " and " << xy[0][plast] << endl;

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

  MELAout << "HelperFunctions::addPointsBetween: Adding " << nadd << " points between " << xy[0][pfirst] << " and " << xy[0][plast] << endl;

  double addWidth = (xmax-xmin)/double(nadd+1);
  for (unsigned int it=1; it<=nadd; it++){
    double xnew = xmin + addWidth*double(it);
    double ynew = sptmp->Eval(xnew);
    addPoint(tgOriginal, xnew, ynew, 0, 0);
  }

  delete sptmp;
  delete tgtmp;
}

template<> double HelperFunctions::evaluateTObject<TH1F>(TH1F* obj, float val){
  int bin = obj->GetXaxis()->FindBin(val);
  if (bin>obj->GetXaxis()->GetNbins()) bin=obj->GetXaxis()->GetNbins();
  else if (bin<1) bin=1;
  return obj->GetBinContent(bin);
}
template<> double HelperFunctions::evaluateTObject<TGraph>(TGraph* obj, float val){
  double* xx = obj->GetX();
  int n = obj->GetN();
  if (val>xx[n-1]) val=xx[n-1];
  else if (val<xx[0]) val=xx[0];
  return obj->Eval(val);
}

template<> void HelperFunctions::regularizeHistogram<TH1F>(TH1F* histo, int nIter_, double threshold_, double acceleration_, unsigned int /*iaxis_*/){
  const int nbinsx = histo->GetNbinsX();

  double xy[4][nbinsx]={ { 0 } };
  double abssumerr=0;
  for (int ix=1; ix<=nbinsx; ix++){
    xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
    xy[1][ix-1] = histo->GetBinContent(ix);
    xy[3][ix-1] = histo->GetBinError(ix);
    abssumerr += xy[3][ix-1];
  }
  if (abssumerr==0.){
    TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
    tg->SetName("tg_tmp");
    regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_);
    for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, tg->GetY()[ix-1]);
    delete tg;
  }
  else{
    TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
    tg->SetName("tg_tmp");
    regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_, acceleration_);
    for (int ix=1; ix<=nbinsx; ix++){
      histo->SetBinContent(ix, tg->GetY()[ix-1]);
      histo->SetBinError(ix, tg->GetEY()[ix-1]);
    }
    delete tg;
  }
}
template<> void HelperFunctions::regularizeHistogram<TH2F>(TH2F* histo, int nIter_, double threshold_, double acceleration_, unsigned int iaxis_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();

  switch (iaxis_){
  case 0:
  {
    for (int iy=1; iy<=nbinsy; iy++){
      double xy[4][nbinsx]={ { 0 } };
      double abssumerr=0;
      for (int ix=1; ix<=nbinsx; ix++){
        xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
        xy[1][ix-1] = histo->GetBinContent(ix, iy);
        xy[3][ix-1] = histo->GetBinError(ix, iy);
        abssumerr += xy[3][ix-1];
      }
      if (abssumerr==0.){
        TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
        tg->SetName("tg_tmp");
        regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_);
        for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, tg->GetY()[ix-1]);
        delete tg;
      }
      else{
        TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
        tg->SetName("tg_tmp");
        regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_, acceleration_);
        for (int ix=1; ix<=nbinsx; ix++){
          histo->SetBinContent(ix, iy, tg->GetY()[ix-1]);
          histo->SetBinError(ix, iy, tg->GetEY()[ix-1]);
        }
        delete tg;
      }
    }
    break;
  }
  case 1:
  {
    for (int ix=1; ix<=nbinsx; ix++){
      double xy[4][nbinsy]={ { 0 } };
      double abssumerr=0;
      for (int iy=1; iy<=nbinsy; iy++){
        xy[0][iy-1] = histo->GetYaxis()->GetBinCenter(iy);
        xy[1][iy-1] = histo->GetBinContent(ix, iy);
        xy[3][iy-1] = histo->GetBinError(ix, iy);
        abssumerr += xy[3][iy-1];
      }
      if (abssumerr==0.){
        TGraph* tg = new TGraph(nbinsy, xy[0], xy[1]);
        tg->SetName("tg_tmp");
        regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsy-2]+xy[0][nbinsy-1])*0.5, nIter_, threshold_);
        for (int iy=1; iy<=nbinsy; iy++) histo->SetBinContent(ix, iy, tg->GetY()[iy-1]);
        delete tg;
      }
      else{
        TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
        tg->SetName("tg_tmp");
        regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsy-2]+xy[0][nbinsy-1])*0.5, nIter_, threshold_, acceleration_);
        for (int iy=1; iy<=nbinsy; iy++){
          histo->SetBinContent(ix, iy, tg->GetY()[iy-1]);
          histo->SetBinError(ix, iy, tg->GetEY()[iy-1]);
        }
        delete tg;
      }
    }
    break;
  }
  }
}
template<> void HelperFunctions::regularizeHistogram<TH3F>(TH3F* histo, int nIter_, double threshold_, double acceleration_, unsigned int iaxis_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();
  const int nbinsz = histo->GetNbinsZ();

  switch (iaxis_){
  case 0:
  {
    for (int iy=1; iy<=nbinsy; iy++){
      for (int iz=1; iz<=nbinsz; iz++){
        double xy[4][nbinsx]={ { 0 } };
        double abssumerr=0;
        for (int ix=1; ix<=nbinsx; ix++){
          xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
          xy[1][ix-1] = histo->GetBinContent(ix, iy, iz);
          xy[3][ix-1] = histo->GetBinError(ix, iy, iz);
          abssumerr += xy[3][ix-1];
        }
        if (abssumerr==0.){
          TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_);
          for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, iz, tg->GetY()[ix-1]);
          delete tg;
        }
        else{
          TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsx-2]+xy[0][nbinsx-1])*0.5, nIter_, threshold_, acceleration_);
          for (int ix=1; ix<=nbinsx; ix++){
            histo->SetBinContent(ix, iy, iz, tg->GetY()[ix-1]);
            histo->SetBinError(ix, iy, iz, tg->GetEY()[ix-1]);
          }
          delete tg;
        }
      }
    }
    break;
  }
  case 1:
  {
    for (int iz=1; iz<=nbinsz; iz++){
      for (int ix=1; ix<=nbinsx; ix++){
        double xy[4][nbinsy]={ { 0 } };
        double abssumerr=0;
        for (int iy=1; iy<=nbinsy; iy++){
          xy[0][iy-1] = histo->GetYaxis()->GetBinCenter(iy);
          xy[1][iy-1] = histo->GetBinContent(ix, iy);
          xy[3][iy-1] = histo->GetBinError(ix, iy, iz);
          abssumerr += xy[3][iy-1];
        }
        if (abssumerr==0.){
          TGraph* tg = new TGraph(nbinsy, xy[0], xy[1]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsy-2]+xy[0][nbinsy-1])*0.5, nIter_, threshold_);
          for (int iy=1; iy<=nbinsy; iy++) histo->SetBinContent(ix, iy, tg->GetY()[iy-1]);
          delete tg;
        }
        else{
          TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsy-2]+xy[0][nbinsy-1])*0.5, nIter_, threshold_, acceleration_);
          for (int iy=1; iy<=nbinsy; iy++){
            histo->SetBinContent(ix, iy, iz, tg->GetY()[iy-1]);
            histo->SetBinError(ix, iy, iz, tg->GetEY()[iy-1]);
          }
          delete tg;
        }
      }
    }
    break;
  }
  case 2:
  {
    for (int ix=1; ix<=nbinsx; ix++){
      for (int iy=1; iy<=nbinsy; iy++){
        double xy[4][nbinsz]={ { 0 } };
        double abssumerr=0;
        for (int iz=1; iz<=nbinsz; iz++){
          xy[0][iz-1] = histo->GetZaxis()->GetBinCenter(iz);
          xy[1][iz-1] = histo->GetBinContent(ix, iy, iz);
          xy[3][iz-1] = histo->GetBinError(ix, iy, iz);
          abssumerr += xy[3][iz-1];
        }
        if (abssumerr==0.){
          TGraph* tg = new TGraph(nbinsz, xy[0], xy[1]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsz-2]+xy[0][nbinsz-1])*0.5, nIter_, threshold_);
          for (int iz=1; iz<=nbinsz; iz++) histo->SetBinContent(ix, iy, iz, tg->GetY()[iz-1]);
          delete tg;
        }
        else{
          TGraphErrors* tg = new TGraphErrors(nbinsx, xy[0], xy[1], xy[2], xy[3]);
          tg->SetName("tg_tmp");
          regularizeSlice(tg, 0, (xy[0][0]+xy[0][1])*0.5, (xy[0][nbinsz-2]+xy[0][nbinsz-1])*0.5, nIter_, threshold_, acceleration_);
          for (int iz=1; iz<=nbinsz; iz++){
            histo->SetBinContent(ix, iy, iz, tg->GetY()[iz-1]);
            histo->SetBinError(ix, iy, iz, tg->GetEY()[iz-1]);
          }
          delete tg;
        }
      }
    }
    break;
  }
  }
}

template<> void HelperFunctions::conditionalizeHistogram<TH2F>(TH2F* histo, unsigned int axis, std::vector<std::pair<TH2F*, float>> const* conditionalsReference){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral=1;
      double width = (ix>=1 && ix<=histo->GetNbinsX() ? histo->GetXaxis()->GetBinWidth(ix) : 1.);
      if (!conditionalsReference) integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1)/width;
      else{ for (std::pair<TH2F*, float> const& hh:(*conditionalsReference)) integral *= pow(hh.first->Integral(ix, ix, 0, hh.first->GetNbinsY()+1)/width, hh.second); }
      if (integral==0.) continue; // All bins across y are 0.
      for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
  else{
    for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
      double integral = 1;
      double width = (iy>=1 && iy<=histo->GetNbinsY() ? histo->GetYaxis()->GetBinWidth(iy) : 1.);
      if (!conditionalsReference) integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy)/width;
      else{ for (std::pair<TH2F*, float> const& hh:(*conditionalsReference)) integral *= pow(hh.first->Integral(0, hh.first->GetNbinsX()+1, iy, iy)/width, hh.second); }
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
}
template<> void HelperFunctions::conditionalizeHistogram<TH3F>(TH3F* histo, unsigned int axis, std::vector<std::pair<TH3F*, float>> const* conditionalsReference){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral=1;
      double width = (ix>=1 && ix<=histo->GetNbinsX() ? histo->GetXaxis()->GetBinWidth(ix) : 1.);
      if (!conditionalsReference) integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1, 0, histo->GetNbinsZ()+1)/width;
      else{ for (std::pair<TH3F*, float> const& hh:(*conditionalsReference)) integral *= pow(hh.first->Integral(ix, ix, 0, hh.first->GetNbinsY()+1, 0, hh.first->GetNbinsZ()+1)/width, hh.second); }
      if (integral==0.) continue; // All bins across y are 0.
      for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
        for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
  else if (axis==1){
    for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
      double integral=1;
      double width = (iy>=1 && iy<=histo->GetNbinsY() ? histo->GetYaxis()->GetBinWidth(iy) : 1.);
      if (!conditionalsReference) integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy, 0, histo->GetNbinsZ()+1)/width;
      else{ for (std::pair<TH3F*, float> const& hh:(*conditionalsReference)) integral *= pow(hh.first->Integral(0, hh.first->GetNbinsX()+1, iy, iy, 0, hh.first->GetNbinsZ()+1)/width, hh.second); }
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
  else{
    for (int iz=0; iz<=histo->GetNbinsZ()+1; iz++){
      double integral=1;
      double width = (iz>=1 && iz<=histo->GetNbinsZ() ? histo->GetZaxis()->GetBinWidth(iz) : 1.);
      if (!conditionalsReference) integral = histo->Integral(0, histo->GetNbinsX()+1, 0, histo->GetNbinsY()+1, iz, iz)/width;
      else{ for (std::pair<TH3F*, float> const& hh:(*conditionalsReference)) integral *= pow(hh.first->Integral(0, hh.first->GetNbinsX()+1, 0, hh.first->GetNbinsY()+1, iz, iz)/width, hh.second); }
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
          histo->SetBinContent(ix, iy, iz, histo->GetBinContent(ix, iy, iz)/integral);
          histo->SetBinError(ix, iy, iz, histo->GetBinError(ix, iy, iz)/integral);
        }
      }
    }
  }
}

template<> void HelperFunctions::wipeOverUnderFlows<TH1F>(TH1F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    hwipe->SetBinContent(binx, 0);
    hwipe->SetBinError(binx, 0);
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}
template<> void HelperFunctions::wipeOverUnderFlows<TH2F>(TH2F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      hwipe->SetBinContent(binx, biny, 0);
      hwipe->SetBinError(binx, biny, 0);
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}
template<> void HelperFunctions::wipeOverUnderFlows<TH3F>(TH3F* hwipe, bool rescale){
  double integral = hwipe->Integral(0, hwipe->GetNbinsX()+1, 0, hwipe->GetNbinsY()+1, 0, hwipe->GetNbinsZ()+1);
  for (int binx=0; binx<=hwipe->GetNbinsX()+1; binx++){
    if (binx>=1 && binx<=hwipe->GetNbinsX()) continue;
    for (int biny=0; biny<=hwipe->GetNbinsY()+1; biny++){
      if (biny>=1 && biny<=hwipe->GetNbinsY()) continue;
      for (int binz=0; binz<=hwipe->GetNbinsZ()+1; binz++){
        if (binz>=1 && binz<=hwipe->GetNbinsZ()) continue;
        hwipe->SetBinContent(binx, biny, binz, 0);
        hwipe->SetBinError(binx, biny, binz, 0);
      }
    }
  }
  double wipeScale = hwipe->Integral();
  wipeScale = integral / wipeScale;
  if (rescale) hwipe->Scale(wipeScale);
}

template<> void HelperFunctions::divideBinWidth<TH1F>(TH1F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    histo->SetBinContent(binx, histo->GetBinContent(binx)/binwidthX);
    histo->SetBinError(binx, histo->GetBinError(binx)/binwidthX);
  }
}
template<> void HelperFunctions::divideBinWidth<TH2F>(TH2F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      float binwidth=binwidthX*binwidthY;
      histo->SetBinContent(binx, biny, histo->GetBinContent(binx, biny)/binwidth);
      histo->SetBinError(binx, biny, histo->GetBinError(binx, biny)/binwidth);
    }
  }
}
template<> void HelperFunctions::divideBinWidth<TH3F>(TH3F* histo){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  TAxis const* zaxis = histo->GetZaxis();
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      for (int binz=1; binz<=histo->GetNbinsZ(); binz++){
        float binwidthZ = zaxis->GetBinWidth(binz);
        float binwidth=binwidthX*binwidthY*binwidthZ;
        histo->SetBinContent(binx, biny, binz, histo->GetBinContent(binx, biny, binz)/binwidth);
        histo->SetBinError(binx, biny, binz, histo->GetBinError(binx, biny, binz)/binwidth);
      }
    }
  }
}

template<> float HelperFunctions::computeIntegral<TH1F>(TH1F* histo, bool useWidth){
  TAxis const* xaxis = histo->GetXaxis();
  float sum=0;
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    float bincontent = histo->GetBinContent(binx);
    if (useWidth) bincontent *= binwidthX;
    sum += bincontent;
  }
  return sum;
}
template<> float HelperFunctions::computeIntegral<TH2F>(TH2F* histo, bool useWidth){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  float sum=0;
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      float binwidth=binwidthX*binwidthY;
      float bincontent = histo->GetBinContent(binx, biny);
      if (useWidth) bincontent *= binwidth;
      sum += bincontent;
    }
  }
  return sum;
}
template<> float HelperFunctions::computeIntegral<TH3F>(TH3F* histo, bool useWidth){
  TAxis const* xaxis = histo->GetXaxis();
  TAxis const* yaxis = histo->GetYaxis();
  TAxis const* zaxis = histo->GetZaxis();
  float sum=0;
  for (int binx=1; binx<=histo->GetNbinsX(); binx++){
    float binwidthX = xaxis->GetBinWidth(binx);
    for (int biny=1; biny<=histo->GetNbinsY(); biny++){
      float binwidthY = yaxis->GetBinWidth(biny);
      for (int binz=1; binz<=histo->GetNbinsZ(); binz++){
        float binwidthZ = zaxis->GetBinWidth(binz);
        float binwidth=binwidthX*binwidthY*binwidthZ;
        float bincontent = histo->GetBinContent(binx, biny, binz);
        if (useWidth) bincontent *= binwidth;
        sum += bincontent;
      }
    }
  }
  return sum;
}

template<> void HelperFunctions::divideHistograms<TH1F>(TH1F* hnum, TH1F* hden, TH1F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  for (int binx=0; binx<=nbinsx+1; binx++){
    float sumW = hnum->GetBinContent(binx);
    float sumWAll = hden->GetBinContent(binx);
    float sumWsq = pow(hnum->GetBinError(binx), 2);
    float sumWsqAll = pow(hden->GetBinError(binx), 2);
    float bincontent=0;
    float binerror=0;
    if (sumWAll!=0.) bincontent = sumW/sumWAll;
    if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
    else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
    if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
      bincontent=0;
      binerror=0;
    }
    hAssign->SetBinContent(binx, bincontent);
    hAssign->SetBinError(binx, binerror);
  }
}
template<> void HelperFunctions::divideHistograms<TH2F>(TH2F* hnum, TH2F* hden, TH2F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  if (hnum->GetNbinsY()!=hden->GetNbinsY()) return; const int nbinsy = hnum->GetNbinsY();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      float sumW = hnum->GetBinContent(binx, biny);
      float sumWAll = hden->GetBinContent(binx, biny);
      float sumWsq = pow(hnum->GetBinError(binx, biny), 2);
      float sumWsqAll = pow(hden->GetBinError(binx, biny), 2);
      float bincontent=0;
      float binerror=0;
      if (sumWAll!=0.) bincontent = sumW/sumWAll;
      if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
      else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
      if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
        bincontent=0;
        binerror=0;
      }
      hAssign->SetBinContent(binx, biny, bincontent);
      hAssign->SetBinError(binx, biny, binerror);
    }
  }
}
template<> void HelperFunctions::divideHistograms<TH3F>(TH3F* hnum, TH3F* hden, TH3F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  if (hnum->GetNbinsY()!=hden->GetNbinsY()) return; const int nbinsy = hnum->GetNbinsY();
  if (hnum->GetNbinsZ()!=hden->GetNbinsZ()) return; const int nbinsz = hnum->GetNbinsZ();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      for (int binz=0; binz<=nbinsz+1; binz++){
        float sumW = hnum->GetBinContent(binx, biny, binz);
        float sumWAll = hden->GetBinContent(binx, biny, binz);
        float sumWsq = pow(hnum->GetBinError(binx, biny, binz), 2);
        float sumWsqAll = pow(hden->GetBinError(binx, biny, binz), 2);
        float bincontent=0;
        float binerror=0;
        if (sumWAll!=0.) bincontent = sumW/sumWAll;
        if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
        else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
        if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
          bincontent=0;
          binerror=0;
        }
        hAssign->SetBinContent(binx, biny, binz, bincontent);
        hAssign->SetBinError(binx, biny, binz, binerror);
      }
    }
  }
}

template <> void HelperFunctions::symmetrizeHistogram<TH1F>(TH1F* histo, unsigned int const axis){
  if (!histo) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = nbinsx/2;
  for (int ix=1; ix<=ixlast; ix++){
    int jx = nbinsx-ix+1;
    float cfirst = histo->GetBinContent(ix);
    float csecond = histo->GetBinContent(jx);
    float c = (cfirst+csecond)/2.;
    histo->SetBinContent(ix, c);
    histo->SetBinContent(jx, c);
  }
}
template <> void HelperFunctions::symmetrizeHistogram<TH2F>(TH2F* histo, unsigned int const axis){
  if (!histo || axis>1) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = (axis==0 ? nbinsx/2 : nbinsx);
  const int nbinsy = histo->GetNbinsY();
  const int iylast = (axis==1 ? nbinsy/2 : nbinsy);
  for (int ix=1; ix<=ixlast; ix++){
    int jx = (axis==0 ? nbinsx-ix+1 : ix);
    for (int iy=1; iy<=iylast; iy++){
      int jy = (axis==1 ? nbinsy-iy+1 : iy);
      float cfirst = histo->GetBinContent(ix, iy);
      float csecond = histo->GetBinContent(jx, jy);
      float c = (cfirst+csecond)/2.;
      histo->SetBinContent(ix, iy, c);
      histo->SetBinContent(jx, jy, c);
    }
  }
}
template <> void HelperFunctions::symmetrizeHistogram<TH3F>(TH3F* histo, unsigned int const axis){
  if (!histo || axis>2) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = (axis==0 ? nbinsx/2 : nbinsx);
  const int nbinsy = histo->GetNbinsY();
  const int iylast = (axis==1 ? nbinsy/2 : nbinsy);
  const int nbinsz = histo->GetNbinsZ();
  const int izlast = (axis==2 ? nbinsz/2 : nbinsz);
  for (int ix=1; ix<=ixlast; ix++){
    int jx = (axis==0 ? nbinsx-ix+1 : ix);
    for (int iy=1; iy<=iylast; iy++){
      int jy = (axis==1 ? nbinsy-iy+1 : iy);
      for (int iz=1; iz<=izlast; iz++){
        int jz = (axis==2 ? nbinsz-iz+1 : iz);
        float cfirst = histo->GetBinContent(ix, iy, iz);
        float csecond = histo->GetBinContent(jx, jy, jz);
        float c = (cfirst+csecond)/2.;
        histo->SetBinContent(ix, iy, iz, c);
        histo->SetBinContent(jx, jy, jz, c);
      }
    }
  }
}

template <> void HelperFunctions::antisymmetrizeHistogram<TH1F>(TH1F* histo, unsigned int const axis){
  if (!histo) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = nbinsx/2;
  for (int ix=1; ix<=ixlast; ix++){
    int jx = nbinsx-ix+1;
    float cfirst = histo->GetBinContent(ix);
    float csecond = histo->GetBinContent(jx);
    float c = (cfirst-csecond)/2.;
    histo->SetBinContent(ix, c);
    histo->SetBinContent(jx, -c);
  }
}
template <> void HelperFunctions::antisymmetrizeHistogram<TH2F>(TH2F* histo, unsigned int const axis){
  if (!histo || axis>1) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = (axis==0 ? nbinsx/2 : nbinsx);
  const int nbinsy = histo->GetNbinsY();
  const int iylast = (axis==1 ? nbinsy/2 : nbinsy);
  for (int ix=1; ix<=ixlast; ix++){
    int jx = (axis==0 ? nbinsx-ix+1 : ix);
    for (int iy=1; iy<=iylast; iy++){
      int jy = (axis==1 ? nbinsy-iy+1 : iy);
      float cfirst = histo->GetBinContent(ix, iy);
      float csecond = histo->GetBinContent(jx, jy);
      float c = (cfirst-csecond)/2.;
      histo->SetBinContent(ix, iy, c);
      histo->SetBinContent(jx, jy, -c);
    }
  }
}
template <> void HelperFunctions::antisymmetrizeHistogram<TH3F>(TH3F* histo, unsigned int const axis){
  if (!histo || axis>2) return;
  const int nbinsx = histo->GetNbinsX();
  const int ixlast = (axis==0 ? nbinsx/2 : nbinsx);
  const int nbinsy = histo->GetNbinsY();
  const int iylast = (axis==1 ? nbinsy/2 : nbinsy);
  const int nbinsz = histo->GetNbinsZ();
  const int izlast = (axis==2 ? nbinsz/2 : nbinsz);
  for (int ix=1; ix<=ixlast; ix++){
    int jx = (axis==0 ? nbinsx-ix+1 : ix);
    for (int iy=1; iy<=iylast; iy++){
      int jy = (axis==1 ? nbinsy-iy+1 : iy);
      for (int iz=1; iz<=izlast; iz++){
        int jz = (axis==2 ? nbinsz-iz+1 : iz);
        float cfirst = histo->GetBinContent(ix, iy, iz);
        float csecond = histo->GetBinContent(jx, jy, jz);
        float c = (cfirst-csecond)/2.;
        histo->SetBinContent(ix, iy, iz, c);
        histo->SetBinContent(jx, jy, jz, -c);
      }
    }
  }
}

void HelperFunctions::CopyFile(TString fname, TTree*(*fcnTree)(TTree*), TDirectory*(*fcnDirectory)(TDirectory*)){
  // Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory* target = gDirectory;
  TFile* f = TFile::Open(fname, "read");
  if (!f || f->IsZombie()){
    MELAerr << "HelperFunctions::CopyFile: Cannot copy file " << fname << endl;
    target->cd();
    if (f && f->IsOpen()){ f->Close(); f=nullptr; }
    delete f;
    return;
  }
  target->cd();
  CopyDirectory(f, fcnTree, fcnDirectory);
  f->Close();
  target->cd();
}
void HelperFunctions::CopyDirectory(TDirectory* source, TTree*(*fcnTree)(TTree*), TDirectory*(*fcnDirectory)(TDirectory*)){
  // Copy all objects and subdirs of directory source as a subdir of the current directory
  source->ls();
  TDirectory* savdir = gDirectory;
  TDirectory* adir;
  if (dynamic_cast<TFile*>(source)==nullptr) adir = savdir->mkdir(source->GetName());
  else adir=savdir;
  adir->cd();
  // Loop on all entries of this directory
  TKey* key;
  TIter nextkey(source->GetListOfKeys());
  vector<TString> copiedTrees;
  while ((key = (TKey*)nextkey())){
    MELAout << "HelperFunctions::CopyDirectory: Copying key " << key->GetName() << endl;
    const char* classname = key->GetClassName();
    TClass* cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())){
      MELAout << "- HelperFunctions::CopyDirectory: Key " << key->GetName() << " is a directory." << endl;
      source->cd(key->GetName());
      TDirectory* subdir = gDirectory;
      if (fcnDirectory) subdir = fcnDirectory(subdir);
      adir->cd();
      CopyDirectory(subdir, fcnTree, fcnDirectory);
      adir->cd();
    }
    else if (cl->InheritsFrom(TTree::Class())){
      MELAout << "- HelperFunctions::CopyDirectory: Key " << key->GetName() << " is a tree." << endl;
      TTree* T = (TTree*)source->Get(key->GetName());
      bool alreadyCopied=false;
      for (auto& k:copiedTrees){
        if (k==key->GetName()){
          alreadyCopied=true;
          break;
        }
      }
      adir->cd();
      if (!alreadyCopied){
        TTree* newT=nullptr;
        if (fcnTree) newT = fcnTree(T);
        if (!newT) newT = T->CloneTree(-1, "fast");
        if (newT){
          MELAout << "- HelperFunctions::CopyDirectory: A new tree " << newT->GetName() << " is created." << endl;
          newT->Write(nullptr, TObject::kWriteDelete);
          copiedTrees.push_back(key->GetName());
        }
        //delete newT;
      }
      else MELAout << "- HelperFunctions::CopyDirectory: A copy of " << key->GetName() << " already exists." << endl;
    }
    else{
      source->cd();
      TObject* obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void HelperFunctions::extractTreesFromDirectory(TDirectory* source, std::vector<TTree*>& res, bool doClone){
  TDirectory* target = gDirectory;

  source->cd();
  // Loop on all entries of this directory
  TKey* key;
  TIter nextkey(source->GetListOfKeys());
  vector<TString> copiedTrees;
  while ((key = (TKey*) nextkey())){
    MELAout << "HelperFunctions::GetTreesInDirectory: Acquiring key " << key->GetName() << endl;
    const char* classname = key->GetClassName();
    TClass* cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())){
      MELAout << "- HelperFunctions::GetTreesInDirectory: Key " << key->GetName() << " is a directory." << endl;
      source->cd(key->GetName());
      TDirectory* subdir = gDirectory;
      source->cd();
      extractTreesFromDirectory(subdir, res, doClone);
    }
    else if (cl->InheritsFrom(TTree::Class())){
      MELAout << "- HelperFunctions::GetTreesInDirectory: Key " << key->GetName() << " is a tree." << endl;
      TTree* T = (TTree*) source->Get(key->GetName());
      bool alreadyCopied=false;
      for (auto& k:copiedTrees){
        if (k==key->GetName()){
          MELAout << "- HelperFunctions::GetTreesInDirectory: Tree " << key->GetName() << " already copied." << endl;
          alreadyCopied=true;
          break;
        }
      }
      if (!alreadyCopied){
        TTree* newT;
        if (doClone) newT = T->CloneTree(-1, "fast");
        else newT = T;
        if (newT){
          copiedTrees.push_back(key->GetName());
          res.push_back(newT);
        }
      }
    }
  }

  target->cd();
}

