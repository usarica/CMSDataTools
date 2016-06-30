// Configuration options
{
	"inputDirectory":"<DIR>LHC_<ENERGY>TeV/Templates/<TODAYSDATE>/gg/",
	"outputFile":"<DIR>LHC_<ENERGY>TeV/Templates/<TODAYSDATE>/gg/HtoZZ4l_ggTo<CHANNEL>_ConditionalSmoothTemplatesForCombine_<SYST><DJET>.root",
	// template definitions
	"templates":[
		// T_1 
		{
			"name":"T_2D_1",
			"files":[
				"./HtoZZ4l_ggTo<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_1_Tree",
			"variables":["GenHMass","KD"],
			"weight":"weight",
			"conserveSumOfWeights":true,
			"assertion":"1",
			"binning":{
				"type":"adaptive",
				"bins":[<XBINNING>,30,0.,1.],
				"entriesperbin":30
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50, "rescalewidth":0.6},
				{"type":"reweight", "axes":[0,1]}
			]
		},
		// T_2
		{
			"name":"T_2D_2",
			"files":[
				"./HtoZZ4l_ggTo<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_2_Tree",
			"variables":["GenHMass","KD"],
			"weight":"weight",
			"conserveSumOfWeights":true,
			"assertion":"1",
			"binning":{
				"type":"adaptive",
				"bins":[<XBINNING>,30,0.,1.],
				"entriesperbin":30
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50, "rescalewidth":0.6},
				{"type":"reweight", "axes":[0,1]}
			]
		},
		// T_124
		{
			"name":"T_2D_124",
			"files":[
				"./HtoZZ4l_ggTo<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_124_Tree",
			"variables":["GenHMass","KD"],
			"weight":"weight",
			"conserveSumOfWeights":true,
			"assertion":"1",
			"binning":{
				"type":"adaptive",
				"bins":[<XBINNING>,30,0.,1.],
				"entriesperbin":30
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50, "rescalewidth":0.6},
				{"type":"reweight", "axes":[0,1]}
			]
		},
		// T_124_perp
		{
			"name":"T_2D_124_perp",
			"files":[
				"./HtoZZ4l_ggTo<CHANNEL>_ConditionalTemplatesForCombine_<SYST><DJET>.root"
				],
			"tree":"T_2D_124_perp_Tree",
			"variables":["GenHMass","KD"],
			"weight":"weight",
			"conserveSumOfWeights":true,
			"assertion":"1",
			"binning":{
				"type":"adaptive",
				"bins":[<XBINNING>,30,0.,1.],
				"entriesperbin":30
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50, "rescalewidth":0.6},
				{"type":"reweight", "axes":[0,1]}
			]
		}

	]
}
