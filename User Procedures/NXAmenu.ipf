#pragma rtGlobals=1		// Use modern global access method.
#pragma IndependentModule=NXPmenu

Menu "Macros"
	"NEXAFS Analysis", /Q, Execute/P "INSERTINCLUDE \"NXA_Init\"";Execute/P "COMPILEPROCEDURES ";Execute/P/Q "NXA_init()"
	help={"Analyzes line scan NEXAFS spectra taken at an ALS STXM beamline"}
End
