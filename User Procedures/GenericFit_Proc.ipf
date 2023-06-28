#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "CollinsProcs_Igor7"



//Command function to manage the fit 
FUNCTION Fit_Generic(yw,uw,xw [noFit, Init])
	WAVE yw, uw, xw
	Variable noFit //optionally just make the model based on the current parameters
	Variable init // reinitialize the background fitting waves (do this after altering the model & # parameters, for example)

	STRUCT Generic_FitStruct s
	Variable nParams=2 //number of parameters possible in the fit (change as needed)
	Variable timer=startMStimer //time the entire operation
	String CurrentFolder=GetDataFolder(1), wnote=""
	NewDataFolder/O/S :Fit //make Fit subfolder & put all the stuff related to it in there
	
	//Prepare Data


	////////////////////
	//Prep some things//
	////////////////////
	WAVE s.yw=yw, s.uw=uw, s.xw=xw //if doing simultaneous fit, this should be the concatenated data
	String xwn=NameOfWave(s.xw) //stupid Igor bug
	Duplicate/d/o s.yw, s.fw, s.rw
	s.DataN=GetWavesDataFolder(s.yw,0) //assumes the data is in a folder with a descriptive name for the data	
//	wave s.fw = fw
//	wave s.rw = rw

	//Prepare Model: Save the model parameters to the structure and initialize the background waves
	WAVE/Z pw
	If( !WaveExists(pw) || init ) //initialize fit waves
		Make/o/n=(nParams) pw, hw, W_Sigma, pwub,pwlb
		Make/o/t/n=(nParams) pNames
		WAVE/T s.pNames
		Duplicate/d/o s.yw, s.mw
		WAVE s.pw, s.hw, s.mw, s.pwub, s.pwlb
		W_Sigma=0
		s.pw={0,1}  //add initial parameters values here.
		s.pNames={"b","m"}
		s.mw=1 //initially fitt all points
		hw=0 //initialize hold wave with open parameters
		s.pwub = s.pw + 0.2
		s.pwlb = s.pw - 0.2
	endif
	
	WAVE s.pw=pw, s.hw=hw, s.sw=W_Sigma, s.mw=mw, s.pwub=pwub, s.pwlb=pwlb
	WAVE/T s.pNames = pNames
	
	Generic_ParamTable(s) //display the parameters in a table
	Generic_PlotFitResults(s) //displays the data & fit if not already done
	
	//Prepare Fitting Info
	String allConst="" //add all possible constraints for all parameters here
	s.hold=""
	Variable i, j
	For( i=0; i<nParams; i+=1 )
		s.hold+=num2str(s.hw[i])  //read hold wave into the hold string for the fit
	endfor
	Make/t/o/n=0 s.cw //assemble the constraints
	For( i=0; i<nParams; i+=1 )
		If( str2num(s.hold[i])==0 ) //If the parameter is not held add the constraints
			For( j=0; j<ItemsInList(allConst); j+=1 )  //cycle through list of constraints to find those matching the parameter
				If( StringMatch(StringFromList(j,allConst),"K"+num2str(i)+"*") )
					InsertPoints 0, 1, s.cw
					s.cw[0]=StringFromList(j,allConst)
				endif
			endfor
		endif
	endfor
	KillWaves/Z M_iterates
	
	//Do the Fit
	Variable/G V_FitError=0, V_FitOptions=8 //save the iterates for debugging
	sprintf wnote "DataName:%s;DateTime:%1.8e;", s.DataN, dateTime //optionally add model information here
	Note/K/NOCR s.fw, wnote
	If( !noFit )
		//Actual fit operation:
		wave mask
		FuncFit/H=s.hold/M=2/Q Generic_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /M=s.mw /R=s.rw /W=s.uw /I=1 /C=s.cw /STRC=s ///D=s.fw
//		FuncFit Generic_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /STRC=s ///D=s.fw
	
		WAVE s.fw, s.pw, s.xw=$("::"+xwn) //stupid IGOR thing that loses the pointer after executing FuncFit
		If( V_FitError ) //for debugging
			Wave M_iterates
			DoWindow/F $(s.DataN+"Iterates")
			If( V_Flag==0 && WaveExists(M_iterates) )
				edit/N=$(s.DataN+"Iterates") M_iterates as s.DataN+" Iterates"
				modifytable size=8, width=50
			endif
			printf "FitError %d, IterNo.%d, ParamPerturbed %d\r",V_FitError,s.fi.IterNumber,s.fi.ParamPerturbed
			sprintf wnote "V_FitError:%d;IterNumber:%d;ParamPerturbed:%d;",V_FitError,s.fi.IterNumber,s.fi.ParamPerturbed
			NOTE/NOCR s.fw, wnote
		else
			WAVE M_Covar, s.sw=W_sigma
			Duplicate/d/o M_Covar, s.CorrMat
			s.CorrMat=M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q]) //generate correlation matrix
			s.CorrMat= p==q? nan : s.CorrMat
			Generic_DisplayCorr(s) //display the correlations
			s.nchisq = sqrt(V_chisq/(V_npnts-(V_nterms-V_nheld))) //reduced chi squared value for the fit
			s.rw/=s.uw //residuals in sigmas
			s.rw= s.mw==0 ? nan : s.rw //don't show residuals for masked points
			sprintf wnote, "Nchisq:%g;HoldStr:%s;InterNumber:%d;",s.nchisq,s.hold, s.fi.iternumber
			NOTE/NOCR s.fw, wnote
			Generic_FitFunc(s)
		endif
	else
		//fill the waves with the model (minimized operations to be a fit func)
		Generic_FitFunc(s)
		s.rw=(s.yw-s.fw)/s.uw
		s.sw=0
	endif
	
	SetDataFolder $CurrentFolder
	Variable elapsedTime=stopMStimer(timer)/1e6 //in seconds
	sprintf wnote "ElapsedTime:%g;", elapsedTime
	Note/NOCR s.fw, wnote
	wnote=SelectString(NoFit,"Fit Completed","Model Calculated")
	If( elapsedTime<120 )
		printf "%s: Nchisq:%1.2f, Took %1.2f sec. with %d iterations\r", wnote, s.nChisq, elapsedTime, s.fi.iternumber
	else
		printf "%s: Nchisq:%1.2f, Took %1.2f min. with %d iterations\r", wnote, s.nChisq, elapsedTime/60, s.fi.iternumber
	endif
END


//structure to contain all the information for the fit, potentially including data libraries
STRUCTURE Generic_FitStruct
	WAVE pw //parameters for the fit
	WAVE fw //model data to fit to the experimental data
	WAVE xw //independent parameter
	STRUCT WMFitInfoStruct fi //fit intormation
	WAVE sw //parameter uncertainties
	WAVE/T pNames //parameter names
	String dataN //name of the current data to be fit
	WAVE yw //experimental data
	WAVE uw //uncertainties in the experimental data
	WAVE rw //residuals divided by the uncertainty
	WAVE mw //mask wave to mask specific datapoints
	WAVE corrMat //Correlation matrix
	WAVE/T cw //constraint equations on the parameters
	String hold //hold string for the fit
	WAVE hw //wave displaying the hold string values & allowing user to change in a table
	Variable nXpts //for simultaneous fitting of data, tells how many unique x values there are
	Variable nchisq //normed chi squared value for the fit
	
	//Add data libraries here
	Wave pwlb //lower bound for fitting wave
	Wave pwub //Uppwer bound for fitting wave
	variable x0
	
ENDSTRUCTURE


//Actual Fit Function
FUNCTION Generic_FitFunc(s) : FitFunc
	STRUCT Generic_FitStruct &s
	Variable b = s.pw[0]
	Variable m = s.pw[1]
	
	s.fw = s.pw[0] + s.pw[1]*s.xw
	
END



//////////////////////////////////////////////////
//Nothing Fit related below...Only visualization//
//////////////////////////////////////////////////

//Displays the parameters in a table
FUNCTION Generic_ParamTable(s)
	STRUCT Generic_FitStruct &s
	DoWindow/F $(s.DataN+"_Parameters")
	IF( V_flag )
		return 0
	endif
	Edit/N=$(s.DataN+"_Parameters") /W=(375.75,445.25,592.5,611.75) s.pNames,s.pw, s.pwlb,s.pwub, s.sw, s.hw, as s.DataN+" Parameters"
	ModifyTable size=8,format(Point)=1,width(Point)=12,width(s.pNames)=50,width(s.pw)=50, width(s.sw)=50, width(s.hw)=20
END


//Displays the fit in a graph
FUNCTION Generic_PlotFitResults(s)
	STRUCT Generic_FitStruct &s
	String ywN=NameOfWave(s.yw), fwN=NameOfWave(s.fw), lgndTxt
	DoWindow/F $(s.DataN+"_Fitting")
	If( V_Flag )
		return 0
	endif
	Display/N=$(s.DataN+"_Fitting") /W=(475.5,106.25,733.5,415.25) s.yw vs s.xw as s.DataN+" Fitting"
	AppendToGraph s.fw vs s.xw
	ModifyGraph rgb($fwN)=(0,0,0), rgb($ywN)=(65280,0,0)
	//Add Residuals to a graph in top 20% of window
	String rwN=NameOfWave(s.rw)
	AppendToGraph/L=Res s.rw vs s.xw
	ModifyGraph zero(Res)=1, lblPosMode(Res)=1, freePos(Res)=0, rgb($rwN)=(65280,0,0)
	Label Res "Res [sig]"
	ModifyGraph axisEnab(left)={0,0.8}
	ModifyGraph axisEnab(Res)={0.8,1}

	ModifyGraph margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph grid=2, tick=2, mirror=1, minor=1, standoff=0, lblPos(left)=42, msize=1.5
	Label left "Scattering Intensity [au]"
	Label bottom "Photon Energy [eV]"
	ErrorBars $ywN Y,wave=(s.uw,s.uw)
	Cursor/P/H=2 A $(fwN) 0 //place a cursor with a vertical line on the graph
	ShowInfo
	sprintf lgndTxt, "%s\r\\s(%s) Data\r\\s(%s) Fit",s.DataN,ywN,fwN
	Legend/C/N=text0/J/A=MC/X=23.02/Y=37.08 lgndTxt

END


//Displays the correlation matrix
FUNCTION Generic_DisplayCorr(s)
	STRUCT Generic_FitStruct &s
	String DataN=GetWavesDataFolder(s.yw,0)
	DoWindow/F $(DataN+"_Correlations")
	IF( V_flag )
		return 0
	endif
	Display/N=$(DataN+"_Correlations") /W=(566.25,453.5,781.5,610.25) as DataN+" Correlations"
	AutoPositionWindow/R=$(s.DataN+"_Fitting")
	WAVE/Z corrColorTab
	If( !WAveExists(corrColorTab) )
		MakeCorrelationColorTable()
		WAVE/Z corrColorTab
	endif
	AppendImage/T s.corrMat
	ModifyImage corrMat cindex= corrColorTab
	ModifyGraph margin(left)=22,margin(bottom)=14,margin(top)=22,margin(right)=72,height={Plan,1,left,top}
	ModifyGraph grid=1, mirror=2, nticks=4, minor=1, sep=10, fSize=8, standoff=0, tkLblRot(left)=90
	ModifyGraph btLen=3, tlOffset=-2, manTick(left)={1,1,0,0},manMinor(left)={0,0}
	ModifyGraph manTick(top)={1,1,0,0},manMinor(top)={0,0}
	Label left "Parameter No."
	Label top "Parameter No."
	SetAxis/A/R left
	Cursor/P/I A corrMat 0,0
	ColorScale/C/N=text0/F=0/A=MC/X=80.25/Y=26.54 image=corrMat, heightPct=50
	ColorScale/C/N=text0 tickLen=3, lblMargin=0, axisRange={0.75,NaN,0}
	AppendText "Pos. Corr"
	ColorScale/C/N=text1/F=0/A=MC/X=81.48/Y=-22.22 image=corrMat, heightPct=50
	ColorScale/C/N=text1 tickLen=3, lblMargin=0, axisRange={NaN,-0.75,0}
	AppendText "Neg. Corr."
END
