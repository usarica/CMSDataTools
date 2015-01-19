// Configuration options
{
	"inputDirectory":"<DIR>/LHC_8TeV/",
	"outputFile":"<DIR>/LHC_8TeV/4mu/HtoZZ4l_MCFM_125p6_SmoothTemplates_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root",
	// template definitions
	"templates":[
		// T_1 
		{
			"name":"T_2D_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50}
			]
		},
		// T_1_mZZ2_1 
		{
			"name":"T_2D_1_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_1_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50}
			]
		},
		// T_1_mZZ2_2 
		{
			"name":"T_2D_1_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_1_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":50}
			]
		},
		// T_2
		{
			"name":"T_2D_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}
			]
		},
		// T_2_mZZ_1
		{
			"name":"T_2D_2_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_2_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}
			]
		},
		// T_2_mZZ2_2
		{
			"name":"T_2D_2_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_2_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}
			]
		},
		// T_4
		{
			"name":"T_2D_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_4_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}

			]
		},
		// T_4_mZZ2_1
		{
			"name":"T_2D_4_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_4_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}

			]
		},
		// T_4_mZZ2_2
		{
			"name":"T_2D_4_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_4_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":100}

			]
		},
		// T_qqZZ
		{
			"name":"T_2D_qqZZ",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_qqZZ_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":75},
				{"type":"floor"}
			]
		},
		// T_ZX
		{
			"name":"T_2D_ZX",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_ZX_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":3},
				{"type":"reweight", "axes":[0],
				       	"rebinning":[
						[220,280,360,400,460,520,1140,1400,1600]
						]
				},
				{"type":"floor"}
			]
		},
		// T_ZX_merged
		{
			"name":"T_2D_ZX_merged",
			"files":[
				"../LHC_8TeV/4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root",
				"../LHC_8TeV/4e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root",
				"../LHC_8TeV/2mu2e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root",
				"../LHC_7TeV/4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root",
				"../LHC_7TeV/4e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root",
				"../LHC_7TeV/2mu2e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_Nominal.root"
				],
			"tree":"T_2D_ZX_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":20},
				{"type":"floor"}
			]
		},
		// T_VBF_1 
		{
			"name":"T_2D_VBF_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_1_mZZ2_1 
		{
			"name":"T_2D_VBF_1_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_1_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_1_mZZ2_2 
		{
			"name":"T_2D_VBF_1_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_1_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_1_mZZ2_3 
		{
			"name":"T_2D_VBF_1_mZZ2_3",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_1_mZZ2_3_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_1_mZZ2_4 
		{
			"name":"T_2D_VBF_1_mZZ2_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_1_mZZ2_4_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_2
		{
			"name":"T_2D_VBF_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_2_mZZ2_1
		{
			"name":"T_2D_VBF_2_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_2_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_2_mZZ2_2
		{
			"name":"T_2D_VBF_2_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_2_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_2_mZZ2_3
		{
			"name":"T_2D_VBF_2_mZZ2_3",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_2_mZZ2_3_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_2_mZZ2_4
		{
			"name":"T_2D_VBF_2_mZZ2_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_2_mZZ2_4_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_4
		{
			"name":"T_2D_VBF_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_4_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_4_mZZ2_1
		{
			"name":"T_2D_VBF_4_mZZ2_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_4_mZZ2_1_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_4_mZZ2_2
		{
			"name":"T_2D_VBF_4_mZZ2_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_4_mZZ2_2_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_4_mZZ2_3
		{
			"name":"T_2D_VBF_4_mZZ2_3",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_4_mZZ2_3_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		},
		// T_VBF_4_mZZ2_4
		{
			"name":"T_2D_VBF_4_mZZ2_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_fLQAdded__GenLevelVBF_D_Gamma_gg_r10_<SYST><JET>.root"
				],
			"tree":"T_2D_VBF_4_mZZ2_4_Tree",
			"variables":["ZZMass","D_Gamma_gg_r10"],
			"weight":"templateWeight",
			"conserveSumOfWeights":true,
			"selection":"ZZMass>=220 && ZZMass<=1600",
			"assertion":"1",
			"binning":{
				"type":"fixed",
				"bins":[69,220.,1600.,30,0.,1.]
			},
			"postprocessing":[
				{"type":"smooth", "kernel":"adaptive", "entriesperbin":30}
			]
		}
	]
}
