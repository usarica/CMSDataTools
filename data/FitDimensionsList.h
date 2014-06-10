/*

 */

const int kFitTotalD = 17;
char* strFitDim[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"ZZMass",
	"D_bkg_kin",

	"D_Gamma",
	"D_Gamma_int",

	"D_Gamma_r10",
	"D_Gamma_r25",
	"D_Gamma_gg_r10",
	"D_Gamma_gg_r25",

	"D_ggZZ",

	"D_Gamma_r5",
	"D_Gamma_r15",
	"D_Gamma_r20",
	"D_Gamma_gg_r5",
	"D_Gamma_gg_r15",
	"D_Gamma_gg_r20",

	"D_Gamma_r1",
	"D_Gamma_gg_r1"
};
char* str1DFit_title[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"ZZMass",
	"T_1D_D_bkg_kin",

	"T_1D_D_Gamma",
	"T_1D_D_Gamma_int",

	"T_1D_D_Gamma_r10",
	"T_1D_D_Gamma_r25",
	"T_1D_D_Gamma_gg_r10",
	"T_1D_D_Gamma_gg_r25",

	"T_1D_D_ggZZ",

	"T_1D_D_Gamma_r5",
	"T_1D_D_Gamma_r15",
	"T_1D_D_Gamma_r20",
	"T_1D_D_Gamma_gg_r5",
	"T_1D_D_Gamma_gg_r15",
	"T_1D_D_Gamma_gg_r20",

	"T_1D_D_Gamma_r1",
	"T_1D_D_Gamma_gg_r1"
};
char* strFitDim_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"m_{ZZ}",
	"D^{kin}_{Bkg}",
	"D_{#Gamma}",
	"D_{#Gamma, int}",

	"D_{#Gamma, 10}",
	"D_{#Gamma, 25}",
	"D_{gg}",
	"D_{gg25}",

	"D_{ggZZ}",

	"D_{#Gamma, 5}",
	"D_{#Gamma, 15}",
	"D_{#Gamma, 20}",
	"D_{gg5}",
	"D_{gg15}",
	"D_{gg20}",

	"D_{#Gamma, 1}",
	"D_{gg1}"
};
char* str1DFit_label[kFitTotalD] = { // Dimension should not exceed kFitTotalD!!!
	"1D Fit to m_{ZZ}",
	"1D Fit to D^{kin}_{Bkg}",
	"1D Fit to D_{#Gamma}",
	"1D Fit to D_{#Gamma, int}",

	"1D Fit to D_{#Gamma, 10}",
	"1D Fit to D_{#Gamma, 25}",
	"1D Fit to D_{gg10}",
	"1D Fit to D_{gg25}",

	"1D Fit to D_{ggZZ}",

	"1D Fit to D_{#Gamma, 5}",
	"1D Fit to D_{#Gamma, 15}",
	"1D Fit to D_{#Gamma, 20}",
	"1D Fit to D_{gg5}",
	"1D Fit to D_{gg15}",
	"1D Fit to D_{gg20}",

	"1D Fit to D_{#Gamma, 1}",
	"1D Fit to D_{gg1}"
};

char* str2DFit_title[kFitTotalD-1] = {
	"T_2D_ZZMass_vs_D_bkg_kin",

	"T_2D_ZZMass_vs_D_Gamma",
	"T_2D_ZZMass_vs_D_Gamma_int",

	"T_2D_ZZMass_vs_D_Gamma_r10",
	"T_2D_ZZMass_vs_D_Gamma_r25",
	"T_2D_ZZMass_vs_D_Gamma_gg_r10",
	"T_2D_ZZMass_vs_D_Gamma_gg_r25",

	"T_2D_ZZMass_vs_D_ggZZ",

	"T_2D_ZZMass_vs_D_Gamma_r5",
	"T_2D_ZZMass_vs_D_Gamma_r15",
	"T_2D_ZZMass_vs_D_Gamma_r20",
	"T_2D_ZZMass_vs_D_Gamma_gg_r5",
	"T_2D_ZZMass_vs_D_Gamma_gg_r15",
	"T_2D_ZZMass_vs_D_Gamma_gg_r20",

	"T_2D_ZZMass_vs_D_Gamma_r1",
	"T_2D_ZZMass_vs_D_Gamma_gg_r1"
};
char* str2DFit_label[kFitTotalD-1] = {
	"2D Fit to m_{ZZ} and D^{kin}_{Bkg}",
	"2D Fit to m_{ZZ} and D_{#Gamma}",
	"2D Fit to m_{ZZ} and D_{#Gamma, int}",

	"2D Fit to m_{ZZ} and D_{#Gamma, 10}",
	"2D Fit to m_{ZZ} and D_{#Gamma, 25}",
	"2D Fit to m_{ZZ} and D_{gg10}",
	"2D Fit to m_{ZZ} and D_{gg25}",

	"2D Fit to m_{ZZ} and D_{ggZZ}",

	"2D Fit to m_{ZZ} and D_{#Gamma, 5}",
	"2D Fit to m_{ZZ} and D_{#Gamma, 15}",
	"2D Fit to m_{ZZ} and D_{#Gamma, 20}",
	"2D Fit to m_{ZZ} and D_{gg5}",
	"2D Fit to m_{ZZ} and D_{gg15}",
	"2D Fit to m_{ZZ} and D_{gg20}",

	"2D Fit to m_{ZZ} and D_{#Gamma, 1}",
	"2D Fit to m_{ZZ} and D_{gg1}"
};
