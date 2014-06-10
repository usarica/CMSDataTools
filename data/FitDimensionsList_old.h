/*

 */

const int kFitTotalD = 6;
char* strFitDim[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"ZZMass",
	"D_bkg_dec",
	"D_Gamma_c10",
	"D_Gamma_c15",
	"D_Gamma_c20",
	"D_Gamma_c25"
};
char* str1DFit_title[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"T_1D_ZZMass",
	"T_1D_D_bkg_dec",
	"T_1D_D_Gamma_c10",
	"T_1D_D_Gamma_c15",
	"T_1D_D_Gamma_c20",
	"T_1D_D_Gamma_c25"
};
char* strFitDim_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"m_{4l}",
	"D^{dec}_{Bkg}",
	"D_{#Gamma10}",
	"D_{#Gamma15}",
	"D_{#Gamma20}",
	"D_{#Gamma25}"
};
char* str1DFit_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"1D Fit to m_{4l}",
	"1D Fit to D^{dec}_{Bkg}",
	"1D Fit to D_{#Gamma10}",
	"1D Fit to D_{#Gamma15}",
	"1D Fit to D_{#Gamma20}",
	"1D Fit to D_{#Gamma25}"
};

char* str2DFit_title[kFitTotalD-1] = {
	"T_2D_ZZMass_vs_D_bkg_dec",
	"T_2D_ZZMass_vs_D_Gamma_c10",
	"T_2D_ZZMass_vs_D_Gamma_c15",
	"T_2D_ZZMass_vs_D_Gamma_c20",
	"T_2D_ZZMass_vs_D_Gamma_c25"
};
char* str2DFit_label[kFitTotalD-1] = {
	"2D Fit to m_{4l} and D^{dec}_{Bkg}",
	"2D Fit to m_{4l} and D_{#Gamma10}",
	"2D Fit to m_{4l} and D_{#Gamma15}",
	"2D Fit to m_{4l} and D_{#Gamma20}",
	"2D Fit to m_{4l} and D_{#Gamma25}"
};
