// Configuration options
{
	"inputDirectory":"<DIR>/LHC_7TeV/",
	"outputFile":"<DIR>/LHC_7TeV/4mu/HtoZZ4l_MCFM_125p6_SmoothTemplates__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root",
	// template definitions
	"templates":[
		// T_1 
		{
			"name":"T_2D_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		// T_2
		{
			"name":"T_2D_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		// T_4
		{
			"name":"T_2D_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		// T_qqZZ
		{
			"name":"T_2D_qqZZ",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
						[220,240,280,320,380,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1340,1360,1380,1400,1420,1440,1460,1480,1500,1520,1540,1560,1580,1600]
						
						]
				},
				{"type":"floor"}

			]
		},
		// T_ZX_merged
		{
			"name":"T_2D_ZX_merged",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root",
				"./4e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root",
				"./2mu2e/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
						[220,240,280,320,380,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1040,1060,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1340,1360,1380,1400,1420,1440,1460,1480,1500,1520,1540,1560,1580,1600]
						
						]
				},
				{"type":"floor"}
				

			]
		},
		// T_VBF_1 
		{
			"name":"T_2D_VBF_1",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		// T_VBF_2
		{
			"name":"T_2D_VBF_2",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		// T_VBF_4
		{
			"name":"T_2D_VBF_4",
			"files":[
				"./4mu/HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_D_Gamma_gg_r10_SysDown_ggPDF.root"
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
		}
	]
}
