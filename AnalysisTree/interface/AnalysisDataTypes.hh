#ifndef ANALYSISDATATYPES_HH
#define ANALYSISDATATYPES_HH


// First argument is name, second argument is type, third argument in simple types is default value
#define SIMPLE_DATA_INPUT_DIRECTIVES \
SIMPLE_DATA_INPUT_DIRECTIVE(bool, bool, false) \
SIMPLE_DATA_INPUT_DIRECTIVE(uchar, unsigned char, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(char, char, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(ushort, unsigned short, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(short, short, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(uint, unsigned int, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(int, int, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(ulong, unsigned long, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(long, long, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(ulonglong, unsigned long long, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(longlong, long long, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(float, float, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(double, double, 0) \
SIMPLE_DATA_INPUT_DIRECTIVE(TBits, TBits, TBits()) \
SIMPLE_DATA_INPUT_DIRECTIVE(string, std::string, "") \
SIMPLE_DATA_INPUT_DIRECTIVE(TString, TString, "") \
SIMPLE_DATA_INPUT_DIRECTIVE(CMSLorentzVector, CMSLorentzVector, CMSLorentzVector(0, 0, 0, 0))

#define VECTOR_DATA_INPUT_DIRECTIVES \
VECTOR_DATA_INPUT_DIRECTIVE(bool, std::vector<bool>) \
VECTOR_DATA_INPUT_DIRECTIVE(uchar, std::vector<unsigned char>) \
VECTOR_DATA_INPUT_DIRECTIVE(char, std::vector<char>) \
VECTOR_DATA_INPUT_DIRECTIVE(ushort, std::vector<unsigned short>) \
VECTOR_DATA_INPUT_DIRECTIVE(short, std::vector<short>) \
VECTOR_DATA_INPUT_DIRECTIVE(uint, std::vector<unsigned int>) \
VECTOR_DATA_INPUT_DIRECTIVE(int, std::vector<int>) \
VECTOR_DATA_INPUT_DIRECTIVE(ulong, std::vector<unsigned long>) \
VECTOR_DATA_INPUT_DIRECTIVE(long, std::vector<long>) \
VECTOR_DATA_INPUT_DIRECTIVE(ulonglong, std::vector<unsigned long long>) \
VECTOR_DATA_INPUT_DIRECTIVE(longlong, std::vector<long long>) \
VECTOR_DATA_INPUT_DIRECTIVE(float, std::vector<float>) \
VECTOR_DATA_INPUT_DIRECTIVE(double, std::vector<double>) \
/*VECTOR_DATA_INPUT_DIRECTIVE(TBits, std::vector<TBits>)*/ \
VECTOR_DATA_INPUT_DIRECTIVE(string, std::vector<std::string>) \
VECTOR_DATA_INPUT_DIRECTIVE(TString, std::vector<TString>) \
VECTOR_DATA_INPUT_DIRECTIVE(CMSLorentzVector, std::vector<CMSLorentzVector>)

#define DOUBLEVECTOR_DATA_INPUT_DIRECTIVES \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(bool, std::vector<std::vector<bool>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(uchar, std::vector<std::vector<unsigned char>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(char, std::vector<std::vector<char>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(ushort, std::vector<std::vector<unsigned short>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(short, std::vector<std::vector<short>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(uint, std::vector<std::vector<unsigned int>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(int, std::vector<std::vector<int>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(ulong, std::vector<std::vector<unsigned long>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(long, std::vector<std::vector<long>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(ulonglong, std::vector<std::vector<unsigned long long>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(longlong, std::vector<std::vector<long long>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(float, std::vector<std::vector<float>>) \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(double, std::vector<std::vector<double>>) \
/*DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(TBits, std::vector<std::vector<TBits>>)*/ \
/*DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(string, std::vector<std::vector<std::string>>)*/ \
/*DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(TString, std::vector<std::vector<TString>>)*/ \
DOUBLEVECTOR_DATA_INPUT_DIRECTIVE(CMSLorentzVector, std::vector<std::vector<CMSLorentzVector>>)


#define SIMPLE_DATA_OUTPUT_DIRECTIVES \
SIMPLE_DATA_OUTPUT_DIRECTIVE(bool, bool) \
/*SIMPLE_DATA_OUTPUT_DIRECTIVE(uchar, unsigned char)*/ \
/*SIMPLE_DATA_OUTPUT_DIRECTIVE(char, char)*/ \
/*SIMPLE_DATA_OUTPUT_DIRECTIVE(ushort, unsigned short)*/ \
SIMPLE_DATA_OUTPUT_DIRECTIVE(short, short) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(uint, unsigned int) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(int, int) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(ulong, unsigned long) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(long, long) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(ulonglong, unsigned long long) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(longlong, long long) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(float, float) \
SIMPLE_DATA_OUTPUT_DIRECTIVE(double, double) \
/*SIMPLE_DATA_OUTPUT_DIRECTIVE(CMSLorentzVector, CMSLorentzVector)*/

#define VECTOR_DATA_OUTPUT_DIRECTIVES \
VECTOR_DATA_OUTPUT_DIRECTIVE(bool, std::vector<bool>) \
/*VECTOR_DATA_OUTPUT_DIRECTIVE(uchar, std::vector<unsigned char>)*/ \
/*VECTOR_DATA_OUTPUT_DIRECTIVE(char, std::vector<char>)*/ \
/*VECTOR_DATA_OUTPUT_DIRECTIVE(ushort, std::vector<unsigned short>)*/ \
VECTOR_DATA_OUTPUT_DIRECTIVE(short, std::vector<short>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(uint, std::vector<unsigned int>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(int, std::vector<int>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(ulong, std::vector<unsigned long>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(long, std::vector<long>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(ulonglong, std::vector<unsigned long long>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(longlong, std::vector<long long>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(float, std::vector<float>) \
VECTOR_DATA_OUTPUT_DIRECTIVE(double, std::vector<double>) \
/*VECTOR_DATA_OUTPUT_DIRECTIVE(CMSLorentzVector, std::vector<CMSLorentzVector>)*/

#define DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVES \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(bool, std::vector<std::vector<bool>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(uchar, std::vector<std::vector<unsigned char>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(char, std::vector<std::vector<char>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(ushort, std::vector<std::vector<unsigned short>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(short, std::vector<std::vector<short>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(uint, std::vector<std::vector<unsigned int>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(int, std::vector<std::vector<int>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(ulong, std::vector<std::vector<unsigned long>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(long, std::vector<std::vector<long>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(ulonglong, std::vector<std::vector<unsigned long long>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(longlong, std::vector<std::vector<long long>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(float, std::vector<std::vector<float>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(double, std::vector<std::vector<double>>)*/ \
/*DOUBLEVECTOR_DATA_OUTPUT_DIRECTIVE(CMSLorentzVector, std::vector<std::vector<CMSLorentzVector>>)*/


#endif
