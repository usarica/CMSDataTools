// Configuration options
{
	"inputDirectory":"<DIR>LHC_<ENERGY>TeV/Templates/<TODAYSDATE>/",
	"outputFile":"<DIR>LHC_<ENERGY>TeV/Templates/<TODAYSDATE>/HtoZZ<CHANNEL>_ConditionalSmoothTemplatesForCombine_<SYST><DJET>.root",
	// template definitions
	"templates":[
		// T_2D_Sig
		{
			"name":"T_2D_Sig",
			"files":[
				"./HtoZZ<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_Sig_Tree",
			"variables":["ZZMass","KD"],
			"weight":"weight*reweight",
			"conserveSumOfWeights":true,
			"selection":"(weight*reweight)>0. && (weight*reweight)<200.",
			"binning":{
				"type":"fixed",
				"bins":[<XBINNING>,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":200, "rescalewidth":1.0}
			]
		},
		// T_2D_ggBkg
		{
			"name":"T_2D_ggBkg",
			"files":[
				"./HtoZZ<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_ggBkg_Tree",
			"variables":["ZZMass","KD"],
			"weight":"weight*reweight",
			"conserveSumOfWeights":true,
			"selection":"(weight*reweight)>0. && (weight*reweight)<200.",
			"binning":{
				"type":"fixed",
				"bins":[<XBINNING>,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":200, "rescalewidth":1.0}
			]
		},
		// T_2D_qqBkg
		{
			"name":"T_2D_qqBkg",
			"files":[
				"./HtoZZ<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_qqBkg_Tree",
			"variables":["ZZMass","KD"],
			"weight":"weight*reweight",
			"conserveSumOfWeights":true,
			"selection":"(weight*reweight)>0. && (weight*reweight)<200.",
			"binning":{
				"type":"fixed",
				"bins":[<XBINNING>,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":200, "rescalewidth":1.0}
			]
		},
		// T_2D_ZX
		{
			"name":"T_2D_ZX",
			"files":[
				"./HtoZZ4mu_ConditionalTemplatesForCombine_ZX_Nominal<DJET>.root",
				"./HtoZZ4e_ConditionalTemplatesForCombine_ZX_Nominal<DJET>.root",
				"./HtoZZ2e2mu_ConditionalTemplatesForCombine_ZX_Nominal<DJET>.root"
				],
			"tree":"T_2D_ZX_Tree",
			"variables":["ZZMass","KD"],
			"weight":"weight",
			"conserveSumOfWeights":true,
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[<XBINNING>,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":200, "rescalewidth":1.0},
				{"type":"reweight", "axes":[0]}
			]
		}

	]
}
