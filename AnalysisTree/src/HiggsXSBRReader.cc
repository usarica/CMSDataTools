#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include "StdExtensions.h"
#include "HelperFunctions.h"
#include "HiggsXSBRReader.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


HiggsXSBRReader::HiggsXSBRReader(TString fname, TString partial_width_type){
  std::unordered_map<TString, std::vector<double>> type_vallist_map;
  std::vector<TString> column_list;

  ifstream file;
  file.open(fname.Data());
  bool firstLine=true;
  while (!file.eof()){
    std::string line;
    std::getline(file, line);
    if (line.empty()) continue;
    std::vector<std::string> splitline;
    HelperFunctions::splitOptionRecursive(line, splitline, ',', false);
    if (firstLine){
      for (auto const& str:splitline){
        TString tstr = str.data();
        column_list.push_back(tstr);
        type_vallist_map[tstr] = std::vector<double>();
      }
      if (std::find(column_list.begin(), column_list.end(), partial_width_type)==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type '" << partial_width_type << "' does not exist." << endl;
      }
      if (std::find(column_list.begin(), column_list.end(), "mass")==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type 'mass' does not exist." << endl;
      }
      if (std::find(column_list.begin(), column_list.end(), "total_width")==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type 'total_width' does not exist." << endl;
      }
      firstLine=false;
    }
    else{
      if (column_list.size()!=splitline.size()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: ERROR! column_list.size() (" << column_list.size() << ") != splitline.size() (" << splitline.size() << ")" << endl;
      }
      for (size_t ic=0; ic<column_list.size(); ic++){
        type_vallist_map[column_list.at(ic)].push_back(std::stod(splitline.at(ic)));
      }
    }
  }
  file.close();

  masses = type_vallist_map["mass"];
  total_widths = type_vallist_map["total_width"];
  partial_widths = type_vallist_map[partial_width_type];
  size_t nmasses = masses.size();
  if (partial_width_type!="total_width"){ for (size_t irow=0; irow<nmasses; irow++) partial_widths.at(irow) *= total_widths.at(irow); }
  {
    double dbegin = (partial_widths.at(1)-partial_widths.front())/(masses.at(1)-masses.front());
    double cB = (partial_widths.at(nmasses-1)-partial_widths.at(nmasses-2))/(pow(masses.at(nmasses-1), 3)-pow(masses.at(nmasses-2), 3));
    double dend = 3.*cB*pow(masses.at(nmasses-1), 2);
    sp_partial_width = TSpline3("sp", masses.data(), partial_widths.data(), nmasses, "b1e1", dbegin, dend);
  }
  {
    double dbegin = (total_widths.at(1)-total_widths.front())/(masses.at(1)-masses.front());
    double cB = (total_widths.at(nmasses-1)-total_widths.at(nmasses-2))/(pow(masses.at(nmasses-1), 3)-pow(masses.at(nmasses-2), 3));
    double dend = 3.*cB*pow(masses.at(nmasses-1), 2);
    sp_total_width = TSpline3("sp", masses.data(), total_widths.data(), nmasses, "b1e1", dbegin, dend);
  }
}

float HiggsXSBRReader::eval_partial_width(float const& mass) const{
  if (mass<=masses.back()) return sp_partial_width.Eval(mass);
  else{
    size_t npoints = masses.size();
    double cB = (partial_widths.at(npoints-1)-partial_widths.at(npoints-2))/(pow(masses.at(npoints-1), 3)-pow(masses.at(npoints-2), 3));
    double cA = partial_widths.at(npoints-1) - cB*pow(masses.at(npoints-1), 3);
    return cA + cB*pow(mass, 3);
  }
}
float HiggsXSBRReader::eval_total_width(float const& mass) const{
  if (mass<=masses.back()) return sp_total_width.Eval(mass);
  else{
    size_t npoints = masses.size();
    double cB = (total_widths.at(npoints-1)-total_widths.at(npoints-2))/(pow(masses.at(npoints-1), 3)-pow(masses.at(npoints-2), 3));
    double cA = total_widths.at(npoints-1) - cB*pow(masses.at(npoints-1), 3);
    return cA + cB*pow(mass, 3);
  }
}
