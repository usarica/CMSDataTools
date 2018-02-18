#include "ExtendedHistogram.h"


ExtendedHistogram::ExtendedHistogram() : name(""), title(""){}
ExtendedHistogram::ExtendedHistogram(const TString name_, const TString title_) : name(name_), title(title_){}
ExtendedHistogram::ExtendedHistogram(ExtendedHistogram const& other) : name(other.name), title(other.title){}
ExtendedHistogram::~ExtendedHistogram(){}

const TString& ExtendedHistogram::getName() const{ return name; }
TString ExtendedHistogram::getName(){ return name; }
const TString& ExtendedHistogram::getTitle() const{ return title; }
TString ExtendedHistogram::getTitle(){ return title; }

void ExtendedHistogram::setNameTitle(const TString name_, const TString title_){ name=name_; title=title_; }
