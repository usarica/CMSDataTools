#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


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

  vector<pair<float, float>> sumWgtsPerBin; sumWgtsPerBin.assign(hA->GetNbinsX(), pair<float, float>(0, 0));
  for (int ix=1; ix<=hA->GetNbinsX(); ix++){
    for (int jx=ix; jx<=hA->GetNbinsX(); jx++){
      sumWgtsPerBin.at(ix-1).second += hA->GetBinContent(jx);
      sumWgtsPerBin.at(ix-1).first += hB->GetBinContent(jx);
    }
  }
  TGraph* tg=HelperFunctions::makeGraphFromPair(sumWgtsPerBin, name);
  HelperFunctions::addPoint(tg, 0, 0);
  //tg->GetYaxis()->SetTitle("Sample A efficiency");
  //tg->GetXaxis()->SetTitle("Sample B efficiency");
  return tg;
}

TGraphErrors* HelperFunctions::makeGraphFromTH1(TH1* hx, TH1* hy, TString name){
  if (hx->GetNbinsX()!=hy->GetNbinsX()){
    MELAerr << "Number of bins for x coordinate != those for y" << endl;
    assert(0);
  }
  unsigned int nbins = hx->GetNbinsX();
  double* xexyey[4];
  for (unsigned int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xexyey[0][bin] = hx->GetBinContent(bin+1);
    xexyey[1][bin] = hx->GetBinError(bin+1);

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

void HelperFunctions::regularizeSlice(TGraph* tgSlice, std::vector<double>* fixedX, double omitbelow, int nIter_, double threshold_){
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
  if (omitbelow>0.){
    for (unsigned int bin=0; bin<nbins_slice; bin++){
      if (xy_mod[0][bin]<omitbelow) fixedBins.push_back(bin);
      //MELAout << "Requested to fix bin " << bin << endl;
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
      double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
      double derivative_last = (yy_second[nbins_second-1]-yy_second[nbins_second-2])/(xx_second[nbins_second-1]-xx_second[nbins_second-2]);
      TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double center = xy_mod[0][binIt-1];
      double val = spline->Eval(center);
      if (fabs(xy_mod[1][binIt-1]-val)>threshold*val && val>0.) xy_mod[1][binIt-1]=val;

      delete spline;
      delete interpolator;

      delete[] yy_second;
      delete[] xx_second;
    }
  }

  for (unsigned int iy=0; iy<nbins_slice; iy++) xy_slice[1][iy] = xy_mod[1][iy];
  for (unsigned int ix=0; ix<2; ix++) delete[] xy_mod[ix];
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

template<> void HelperFunctions::regularizeHistogram<TH2F>(TH2F* histo, int nIter_, double threshold_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();

  for (int iy=1; iy<=nbinsy; iy++){
    //MELAout << "regularizeHistogram::Bin " << iy << " being regularized..." << endl;
    double xy[2][nbinsx];
    for (int ix=1; ix<=nbinsx; ix++){
      xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
      xy[1][ix-1] = histo->GetBinContent(ix, iy);
    }

    TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
    tg->SetName("tg_tmp");
    regularizeSlice(tg, 0, 0, nIter_, threshold_);
    for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, tg->Eval(xy[0][ix-1]));
  }
  conditionalizeHistogram(histo, 0);
}
template<> void HelperFunctions::regularizeHistogram<TH3F>(TH3F* histo, int nIter_, double threshold_){
  const int nbinsx = histo->GetNbinsX();
  const int nbinsy = histo->GetNbinsY();
  const int nbinsz = histo->GetNbinsZ();

  for (int iy=1; iy<=nbinsy; iy++){
    for (int iz=1; iz<=nbinsz; iz++){
      //MELAout << "regularizeHistogram::Bin " << iy << ", " << iz << " being regularized..." << endl;
      double xy[2][nbinsx];
      for (int ix=1; ix<=nbinsx; ix++){
        xy[0][ix-1] = histo->GetXaxis()->GetBinCenter(ix);
        xy[1][ix-1] = histo->GetBinContent(ix, iy, iz);
      }

      TGraph* tg = new TGraph(nbinsx, xy[0], xy[1]);
      tg->SetName("tg_tmp");
      regularizeSlice(tg, 0, 0, nIter_, threshold_);
      for (int ix=1; ix<=nbinsx; ix++) histo->SetBinContent(ix, iy, iz, tg->Eval(xy[0][ix-1]));
    }
  }
  conditionalizeHistogram(histo, 0);
}

template<> void HelperFunctions::conditionalizeHistogram<TH2F>(TH2F* histo, unsigned int axis){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1);
      if (integral==0.) continue; // All bins across y are 0.
      for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
  else{
    for (int iy=0; iy<=histo->GetNbinsY()+1; iy++){
      double integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy);
      if (integral==0.) continue; // All bins across y are 0.
      for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
        histo->SetBinContent(ix, iy, histo->GetBinContent(ix, iy)/integral);
        histo->SetBinError(ix, iy, histo->GetBinError(ix, iy)/integral);
      }
    }
  }
}
template<> void HelperFunctions::conditionalizeHistogram<TH3F>(TH3F* histo, unsigned int axis){
  if (axis==0){
    for (int ix=0; ix<=histo->GetNbinsX()+1; ix++){
      double integral = histo->Integral(ix, ix, 0, histo->GetNbinsY()+1, 0, histo->GetNbinsZ()+1);
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
      double integral = histo->Integral(0, histo->GetNbinsX()+1, iy, iy, 0, histo->GetNbinsZ()+1);
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
      double integral = histo->Integral(0, histo->GetNbinsX()+1, 0, histo->GetNbinsY()+1, iz, iz);
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
