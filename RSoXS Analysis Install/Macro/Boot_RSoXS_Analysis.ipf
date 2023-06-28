#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//Edited to run JVC_Init_Ver2 5/20/2017

Menu "Macros"
	"RSoXS Analysis Tools Ver 1.0", /Q, Execute/P "INSERTINCLUDE \"NI1_Loader\"";Execute/P "COMPILEPROCEDURES ";Execute/P "INSERTINCLUDE \"NRS_Functions\"";Execute/P "COMPILEPROCEDURES ";Execute/P/Q "NRS_1101Panel_Init()"
	help={"Opens interface for Beamline 11.0.1.2 "}
End
