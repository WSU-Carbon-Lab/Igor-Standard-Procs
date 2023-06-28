#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//#include "CollinsProcs"
//#include "Ferron_Proc"

//Ferron code to get a structure-based fit function started. Containes code for modeling, fitting, debugging, displaying, and assesing
//Replace the word "Ferron" with a unique name, and add the model details


//Command function to manage the fit 
FUNCTION Fit_Ferron(yw,uw,xw, m1N, m2N, BKG, SubBkg, SubBkgE [ noFit, Init])
	WAVE yw, uw, xw
	String m1N, m2N, BKG // Molecule 1 & 2 Detailed NEXAFS file name, BKG is Beamline 11 sample NEXAFS name
	Variable noFit //optionally just make the model based on the current parameters
	Variable init // reinitialize the background fitting waves (do this after altering the model & # parameters, for example)
	Wave SubBkg , SubBkgE //Background wave you want to subtract from your data before fitting (or add during the fit...unsure so far)
	STRUCT Ferron_FitStruct s
	Variable nParams=9 //number of parameters possible in the fit (change as needed)
	Variable timer=startMStimer //time the entire operation
	String CurrentFolder=GetDataFolder(1), wnote=""
	NewDataFolder/O/S :Fit //make Fit subfolder & put all the stuff related to it in there
	
	//Prepare Data
	Duplicate/d/o uw combinedU
	WAVE s.yw=yw, s.uw=combinedU, s.xw=xw
	String xwn=NameOfWave(s.xw) //stupid Igor bug
	Duplicate/d/o s.yw, s.fw, s.rw, s.SubBkg//, s.comb //combined contrast wave outputted from IsoContrast
	Make/o/c/n=(numpnts(s.xw)) s.OC1wkb, s.OC2wkb //WKB-corrected optical constants
	s.DataN=GetWavesDataFolder(s.yw,0) //assumes the data is in a folder with a descriptive name for the data	
	
	//Collect NEXAFS SPectra **NEW** 11/12/2015 ---Updated 11/13/2015 Because I realized I am a dumbass
	WAVE s.RawAbs=$("root:NEXAFS:"+BKG+":"+BKG+"_absSub"), s.BkgE=$("root:NEXAFS:"+BKG+":Energy")
	
	//Setup the background wave
	s.SubBkg = interp(s.xw,SubBkgE,SubBkg)
		
	
	
	//Prepare Model: Save the model parameters to the structure and initialize the background waves
	WAVE/Z pw
	If( !WaveExists(pw) || init ) //initialize fit waves
		Make/o/n=(nParams) pw, hw, W_Sigma
		Make/o/t/n=(nParams) pNames
		Make/o/c/n=(101) WKBth //WKB correction factor as a function of film thickness
		Make/o/n=(101) WKBthReal, WKBthImag //real and imaginary components of WKBth
		WAVE/T s.pNames
		Duplicate/d/o s.yw, s.mw
		WAVE s.pw, s.hw, s.mw
		W_Sigma=0
		s.pw={1,0.1,0.1,100,50,0, 0.5,1E-3,1E-3}//,1E-6}  //add initial parameters values here.
//		s.pNames={"Amplitude","Eoffset","Eshift","Thickness","m2 Conc.","Fract. Rough Scatt.","Dense Diff.","P1","P2"}
		s.pNames={"Amplitude",m1N+" Eoff",m2N+" Eoff","Thickness","m2 Conc.","Svac","Smat","Re[Svm]","Im[Svm]"}//,"Fluor"}
		s.mw=1 //initially fitt all points
		hw={0,0,0,0,0,0,0,0,0}//,0}//,0} //initialize hold wave with open parameters
	endif
	s.m1N=m1N, s.m2N=m2N
	WAVE s.pw, s.hw, s.sw=W_Sigma, s.mw//, s.comb
	WAVE/T s.pNames
	WAVE/C s.OC1=$("root:OpticalConstants:"+s.m1N), s.OC2=$("root:OpticalConstants:"+s.m2N), s.WKBth
	WAVE/C s.OC1u=$("root:OpticalConstants:"+s.m1N+"_u"), s.OC2u=$("root:OpticalConstants:"+s.m2N+"_u")
	WAVE s.en1=$("root:OpticalConstants:"+s.m1N+"_en"), s.en2=$("root:OpticalConstants:"+s.m2N+"_en"), s.WKBthReal, s.WKBthImag
	If( !WaveExists(s.en1) )
		WAVE s.en1=root:OpticalConstants:Energy
	endif
	If( !WaveExists(s.en2) )
		WAVE s.en2=root:OpticalConstants:Energy
	endif
	Duplicate/d/o s.xw OCru
	//OCru= (1-s.pw[4]/100)*imag(cinterp(s.xw,s.en1,s.OC1u))/imag(cinterp(s.xw,s.en1,s.OC1)) + s.pw[4]/100*imag(cinterp(s.xw,s.en1,s.OC2u))/imag(cinterp(s.xw,s.en1,s.OC2))
	OCru= (1-50/100)*imag(cinterp(s.xw,s.en1,s.OC1u))/imag(cinterp(s.xw,s.en1,s.OC1)) + 50/100*imag(cinterp(s.xw,s.en1,s.OC2u))/imag(cinterp(s.xw,s.en1,s.OC2))

	s.uw= s.yw*sqrt( (s.uw/s.yw)^2 + OCru^2 ) //combine the scattering uncertainties with that of the OCs via quadrature of relative uncertainties
	Ferron_ParamTable(s) //display the parameters in a table
	Ferron_PlotFitResults(s) //displays the data & fit if not already done
	
	//Prepare Fitting Info
	String allConst="K4>0;K4<100;" //add all possible constraints for all parameters here
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
		FuncFit/NTHR=0/H=s.hold/M=2/Q Ferron_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /M=s.mw /R=s.rw /W=s.uw /I=1 /C=s.cw /STRC=s //Remove WKB for Born approximation in fitfunc
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
			Ferron_DisplayCorr(s) //display the correlations
			s.nchisq = sqrt(V_chisq/(V_npnts-(V_nterms-V_nheld))) //reduced chi squared value for the fit
			s.rw/=s.uw //residuals in sigmas
			s.rw= s.mw==0 ? nan : s.rw //don't show residuals for masked points
			sprintf wnote, "Nchisq:%g;HoldStr:%s;InterNumber:%d;",s.nchisq,s.hold, s.fi.iternumber
			NOTE/NOCR s.fw, wnote
			
			//
	//		s.xw = s.en1
	//		redimension/N=(numpnts(s.en1)) s.fw
			
			//
			Ferron_FitFunc(s)
		endif
	else
		//fill the waves with the model (minimized operations to be a fit func)
		Ferron_FitFunc(s)
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
STRUCTURE Ferron_FitStruct
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
	String m1N, m2N //name of two materials in sample
	WAVE/C OC1, OC2 //optical constants of the two materials
	WAVE/C OC1u, OC2u //uncertainties in the optical constants (these are measured after all)
	WAVE/C OC1WKB, OC2WKB //corrected optical constants based on the WKB approximation
	WAVE/C WKBth //WKB correction factor as a function of film thickness
	WAVE WKBthReal, WKBthImag //real and imaginary components of WKBth
	WAVE en1, en2 //energy waves of the two materials
//	WAVE Comb //combined contrast wave (non-scaled)

	//New *11/12/2015* Call NEXAFS file for sample
	WAVE yNXFS
	WAVE xNXFS
	WAVE BkgE
	WAVE RawABS
	WAVE AbsEff
	Wave SubBkg
	Wave SubBkgE
	
	//New things for WKB
	Wave RealWKBTemp
	wave ImagWKBTemp
	wave OCReal1
	wave OCImag1
	wave OCReal2
	wave OCImag2
	Wave Slabmodel
	
	
ENDSTRUCTURE


//Actual Fit Function
FUNCTION Ferron_FitFunc(s) : FitFunc
	STRUCT Ferron_FitStruct &s
//	IsoContrastMulti(s.OC1,s.en1,s.OC2,s.en2,s.xw,s.pw)

//	Ferron_WKBOC(s)
	Variable Amp=s.pw[0]
	Variable Eoff1=s.pw[1], Eoff2=s.pw[2]
	Variable th=s.pw[3]//467//s.pw[3]
	Variable Conc2=s.pw[4]
	Variable s1=s.pw[5]
	Variable s2=s.pw[6], s3=s.pw[7], s4=s.pw[8] //change 5 to 6!!!
	Variable Fluor = 0.0025//
	Conc2= Conc2>1 ? Conc2/100 : Conc2 //if User inputs # greater than 1, assume it is in % & convert it to fraction

	Duplicate/o s.xw, hvOff1, hvOff2, Vac1, Vac2, Comp, AbsBlend, Energy, preFact
	Duplicate/o s.xw XRF, MatScattering
	Duplicate/o/c s.OC1WKB CrossPhase
	hvOff1=s.xw+Eoff1
	hvOff2=s.xw+Eoff2
	Variable hc=1239.84 //eV*nm
	Variable re= 2.81794e-6 //classical electron radius [nm]
	Variable alpha= 2*pi/(hc^2)
	
	//Call NEXAFS Info
	

//Setting up for Absorption
	Duplicate/d/o hvOff1, AbsEff,Transmission
	AbsEff = 1-exp(-interp(s.xw,s.BkgE, s.RawAbs))
	Transmission =  exp(-interp(s.xw,s.BkgE, s.RawAbs))
//	variable i
//	For (i=0;hvOff1[i] < 284.1;i+=1)
//		AbsEff[i] = 0
//	endfor
//	
	//dowindow/K AbsEff
	//Display/N=AbsEff AbsEff vs hvOff1
	
//	




//************************************************************************************************
//Born Approximation
	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*th*s.xw/hc) //NOT USED
	//AbsBlend = 1
	Vac1 = magsqr(cinterp(hvOff1,s.en1,s.OC1))
	Vac2 = magsqr(cinterp(hvOff2,s.en2,s.OC2))
	Comp = magsqr( cinterp(hvOff1,s.en1,s.OC1) - cinterp(hvOff2,s.en2,s.OC2) )
	
	CrossPhase = cinterp(hvOff1,s.en1,s.OC1) * Conj(cinterp(hvOff2,s.en2,s.OC2) - cinterp(hvOff1,s.en1,s.OC1))
	//Forcing s1 as positive
	s.fw=abs(s1)*Vac1 + s2*Comp + 2*(s3* Real(CrossPhase) + s4*Imag(crossPhase))

	s.fw *= alpha^2 * s.xw^4 * Transmission // AbsBlend //Took out Amp from everything and only put it in Fluorescence to see what happens
	
	s.fw += AbsEff*Fluor*(1/(8*PI^3)) *0.160729/th///*(1/(4*PI))AbsEff //The 0.107 is the integral of q^2 over the range due to the splitting fluor from scattering / New Q2 value June 2019 0.160729
//	s.fw += (1-AbsBLend)*Fluor*(1/(8*PI^3)) *0.160729/th///*(1/(4*PI))AbsEff //The 0.107 is the integral of q^2 over the range due to the splitting fluor from scattering / New Q2 value June 2019 0.160729

	//1/(8*PI^2) is the conversion to isoctropic scattering 4*PI and the integration factor of 2*PI^2

	s.fw *= Amp
	
	s.fw += s.pw[4]*s.SubBkg/th
	//  0.00275442
	// 0.72849
	
	XRF = Amp * AbsEff*Fluor*(1/(8*PI^3)) *0.160729/th
	MatScattering = Amp * alpha^2 * s.xw^4 * Transmission * s2*Comp

//************************************************************************************************
//NEW WKB Approximation for Fit Function, 06/07/2019
//	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*th*s.xw/hc)
//	Vac1= magsqr(s.OC1WKB)
//	Vac2= magsqr(s.OC2WKB)
//	Comp=magsqr( s.OC1WKB - s.OC2WKB )
//
//	CrossPhase= s.OC1WKB * conj( s.OC2WKB - s.OC1WKB )
// 
//	s.fw= abs(s1)* Vac1 + s2*Comp + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )
//	s.fw *= alpha^2 * s.xw^4 * (AbsBlend) //Transmission ////AbsBlend 
//	s.fw += AbsEff*Fluor*(1/(8*PI^3))*0.160729/th
//
// 	s.fw *= Amp
////	 s.fw += s.pw[4]*s.SubBkg/th
//
//	 XRF = Amp*AbsEff*Fluor*(1/(8*PI^3))*0.160729/th //instead of thickness use 467...the normalized number
//	 MatScattering = Amp * alpha^2 * s.xw^4 * Transmission * s2*Comp




//************************************************************************************************
//OLD WKB Approximation for Fit Function, Not sure if necessary
//	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*th/2*s.xw/hc)
//	Vac1= magsqr(s.OC1WKB)
//	Vac2= magsqr(s.OC2WKB)
//	Comp=magsqr( s.OC1WKB - s.OC2WKB )
//	CrossPhase= s.OC1WKB * conj( s.OC2WKB - s.OC1WKB ) 
//	s.fw= s1* Vac1 + s2*Comp + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )
//	s.fw *= Amp *alpha^2 * s.xw^4 * AbsBlend 
//	s.fw += AbsEff*Fluor
////
////	If( anchor==0 )
////////*******************************************
//	CrossPhase= s.OC1WKB * conj( s.OC2WKB - s.OC1WKB ) 
//	s.fw= s1* Vac1 + s2*Comp + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )
//////*******************************************
//	elseif( anchor==1 )
//		CrossPhase= conj( (cInterp(hvOff1,hv1,OC1)-cInterp(hvOff2,hv2,OC2)) ) * cinterp(hvOff2,hv2,OC2)
//		s.fw= s1* Comp + s2*Vac2 + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )
//	elseif( anchor==2 )
//		CrossPhase= cinterp(hvOff1,hv1,OC1)*conj(cInterp(hvOff2,hv2,OC2))
//		s.fw= s1* Vac1 + s2*Vac2 + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )
//	endif
//
//////*******************************************
//	s.fw *= Amp *alpha^2 * s.xw^4 * AbsBlend  + AbsEff*Fluor
//*******************************************
//	s.fw=Amp * comb

//************************************************************************************************
//Born Approximation
//	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*th*s.xw/hc)
//	//AbsBlend = 1
//	Vac1= magsqr(cinterp(hvOff1,s.en1,s.OC1))
//	Vac2= magsqr(cinterp(hvOff2,s.en2,s.OC2))
//	Comp=magsqr( cinterp(hvOff1,s.en1,s.OC1) - cinterp(hvOff2,s.en2,s.OC2) )
//	
//	CrossPhase = cinterp(hvOff1,s.en1,s.OC1) * Conj(cinterp(hvOff2,s.en2,s.OC2) - cinterp(hvOff1,s.en1,s.OC1))
//	//Forcing s1 as positive
//	s.fw= abs(s1)* Vac1 + s2*Comp + 2*(s3* Real(CrossPhase) + s4*Imag(crossPhase))	
//
//	s.fw *= Amp *alpha^2 * s.xw^4*AbsBlend 
//
//	s.fw += AbsEff*Fluor*127*(1/2*PI^2)*0.72849*10^(-9)
//************************************************************************************************
END

//Actual Fit Function
FUNCTION Ferron_FitFunc_WKB(s) : FitFunc
	STRUCT Ferron_FitStruct &s
//	IsoContrastMulti(s.OC1,s.en1,s.OC2,s.en2,s.xw,s.pw)

//	Ferron_WKBOC(s)
	Variable Amp=s.pw[0]
	Variable Eoff1=s.pw[1], Eoff2=s.pw[2]
	Variable th=467//s.pw[3]
	Variable Conc2=s.pw[4]
	Variable s1=s.pw[5]
	Variable s2=s.pw[6], s3=s.pw[7], s4=s.pw[8] //change 5 to 6!!!
	Variable Fluor = 0.0025//
	Conc2= Conc2>1 ? Conc2/100 : Conc2 //if User inputs # greater than 1, assume it is in % & convert it to fraction
	
	
	Duplicate/o s.xw, hvOff1, hvOff2, Vac1, Vac2, Comp, AbsBlend, Energy, preFact
	Duplicate/o s.xw XRF, MatScattering
	Duplicate/o/c s.OC1WKB CrossPhase
	//WKB Waves
	Duplicate/o s.xw koff1,koff2
	
	//Setup initial information
	hvOff1=s.xw+Eoff1
	hvOff2=s.xw+Eoff2
	Variable hc=1239.84 //eV*nm
	Variable re= 2.81794e-6 //classical electron radius [nm]
	Variable alpha= 2*pi/(hc^2)
	koff1 = 2*PI*(hvoff1)/hc
	koff2 = 2*PI*(hvoff2)/hc
	
	//WKB Parameters
	
	Variable Thick= th//467 //s.pw[3] //Be careful if pw is reorganized!  This will not be correct anymore!
	Variable i,j//, en //wave number 2pi/lambda
	//Make the slab wave that we integrate over
	Make/o/N=(Thick+1) s.SlabModel = p
	//Make temporary arrays before integration
	Make/o/N=((numpnts(hvoff1)),Thick+1) WKB_ImagArray_1, WKB_RealArray_1
	Make/o/N=((numpnts(hvoff2)),Thick+1) WKB_ImagArray_2, WKB_RealArray_2
	Make/o/N=((numpnts(s.xw)),Thick+1) WKB_CompArray, WKB_VacArray, WKB_FullArray
	Make/C/o/N=((numpnts(s.xw)),Thick+1) WKB_CrossArray
	
	duplicate/o hvoff1, RealWKBOC1, ImagWKBOC1
	duplicate/o hvoff2, RealWKBOC2, ImagWKBOC2
	duplicate/o s.xw, WKB_Comp,WKB_RCross,WKB_Vac,WKB_ICross

	duplicate/o s.en1, s.OCReal1, s.OCImag1
	duplicate/o s.en2, s.OCReal2, s.OCImag2

	
	//Call NEXAFS Info
//	CrossPhase= s.OC1WKB * conj( s.OC2WKB - s.OC1WKB )

	
//	Setscale/I x, 0, th, s.WKBth, s.WKBthReal, s.WKBthImag
//Setup the first optical constant set
	s.OCReal1 = real(s.OC1)
	s.OCImag1 = Imag(s.OC1)
	
	WKB_RealArray_1[][] = interp(hvoff1[p],s.en1,s.OCReal1)*COS(koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCReal1)) - interp(hvoff1[p],s.en1,s.OCimag1)*SIN(koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCReal1))
	WKB_RealArray_1[][] *= exp(-koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCimag1))
	
	WKB_ImagArray_1[][] = interp(hvoff1[p],s.en1,s.OCimag1)*COS(koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCReal1)) + interp(hvoff1[p],s.en1,s.OCReal1)*SIN(koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCReal1))
	WKB_ImagArray_1[][] *= exp(-koff1[p]*s.slabmodel[q]*interp(hvoff1[p],s.en1,s.OCimag1))
	
	//Setup the second optical constant set
	s.OCReal2 = real(s.OC2)
	s.OCImag2 = Imag(s.OC2)
	
	WKB_RealArray_2[][] = interp(hvoff2[p],s.en2,s.OCReal2)*COS(koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCReal2)) - interp(hvoff2[p],s.en2,s.OCImag2)*SIN(koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCReal2))
	WKB_RealArray_2[][] *= exp(-koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCImag2))
	
	WKB_ImagArray_2[][] = interp(hvoff2[p],s.en2,s.OCImag2)*COS(koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCReal2)) + interp(hvoff2[p],s.en2,s.OCReal2)*SIN(koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCReal2))
	WKB_ImagArray_2[][] *= exp(-koff2[p]*s.slabmodel[q]*interp(hvoff2[p],s.en2,s.OCimag2))
	
	
//	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*(thick-slabmodel[q])*s.xw/hc)

	WKB_CompArray[][] = (WKB_RealArray_1[p][q] - WKB_RealArray_2[p][q])^2 + (WKB_ImagArray_1[p][q] - WKB_ImagArray_2[p][q])^2// * exp(-4*pi*((1-Conc2)*(Interp(hvOff1[p],s.en1,s.OCImag1))+Conc2*(Interp(hvOff2[p],s.en2,s.OCImag2)))*(thick-s.slabmodel[q])*s.xw[p]/hc)
	WKB_CrossArray[][] = cinterp(hvoff1[p],s.en1,s.OCReal1) * conj(cinterp(hvoff2[p],s.en2,s.OCReal2) - cinterp(hvoff1[p],s.en1,s.OCReal1))
	WKB_VacArray[][] = (WKB_RealArray_1[p][q])^2 + (WKB_ImagArray_1[p][q])^2
	
	WKB_FullArray[][] = (s1*WKB_VacArray[p][q] + s2*WKB_CompArray[p][q] + 2*(s3*Real(WKB_CrossArray[p][q]) + s4*Imag(WKB_CrossArray[p][q]))) * exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1[p],s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2[p],s.en2,s.OC2)))*(thick-s.slabmodel[q])*s.xw[p]/hc)
	//Now Integrate the Thickness side of things using trapezoidal integration dt comes from spacing in wave 'slabmodel' right now its 1
	//So just summing the rows gives the same trapezoidal integration
//	MatrixOP/O RealWKBOC1 = sumrows(WKB_RealArray_1)/Thick
//	MatrixOP/O RealWKBOC2 = sumrows(WKB_RealArray_2)/Thick
//	MatrixOP/O ImagWKBOC1 = sumrows(WKB_ImagArray_1)/Thick
//	MatrixOP/O ImagWKBOC2 = sumrows(WKB_ImagArray_2)/Thick
	
	MatrixOP/O WKB_Comp = s2*sumrows(WKB_CompArray)/Thick
	MAtrixOP/O WKB_RCross = s3*Real(Sumrows(WKB_CrossArray))/Thick
	MatrixOP/O WKB_ICross = s4*Imag(Sumrows(WKB_CrossArray))/Thick
	MatrixOP/O WKB_Vac = s1*sumrows(WKB_VacArray)/Thick
	
	MatrixOP/O s.fw = sumrows(WKB_FullArray)/thick
	
//	s.OC1WKB = cmplx(interp(s.xw+EOff1,s.en1,RealWKBOC1),interp(s.xw+EOff1,s.en1,ImagWKBOC1))
//	s.OC2WKB = cmplx(interp(s.xw+EOff2,s.en2,RealWKBOC2),interp(s.xw+EOff2,s.en2,ImagWKBOC2))

	
	


//Setting up for Absorption
	Duplicate/d/o hvOff1, AbsEff,Transmission
	AbsEff = 1-exp(-interp(s.xw,s.BkgE, s.RawAbs))
	Transmission =  exp(-interp(s.xw,s.BkgE, s.RawAbs))
//	variable i
//	For (i=0;hvOff1[i] < 284.1;i+=1)
//		AbsEff[i] = 0
//	endfor
//	
	//dowindow/K AbsEff
	//Display/N=AbsEff AbsEff vs hvOff1
	
//	

////************************************************************************************************
////NEW WKB Approximation for Fit Function, 06/07/2019
////	AbsBlend= exp(-4*pi*((1-Conc2)*imag(cInterp(hvOff1,s.en1,s.OC1))+Conc2*imag(cInterp(hvOff2,s.en2,s.OC2)))*th*s.xw/hc)
//	Vac1= magsqr(s.OC1WKB)
//	Vac2= magsqr(s.OC2WKB)
//	Comp=magsqr( s.OC1WKB - s.OC2WKB )
//
//	CrossPhase= s.OC1WKB * conj( s.OC2WKB - s.OC1WKB )
// 
//	s.fw= abs(s1)* Vac1 + s2*WKB_Comp + 2*( s3* Real(CrossPhase) + s4*Imag(crossPhase) )

	s.fw *= alpha^2 * s.xw^4 *1// (AbsBlend) //Transmission ////AbsBlend 
	s.fw += AbsEff*Fluor*(1/(8*PI^3))*0.160729/th

 	s.fw *= Amp
	s.fw += abs(s.pw[3])*s.SubBkg/th

	 XRF = Amp*AbsEff*Fluor*(1/(8*PI^3))*0.160729/th //instead of thickness use 467...the normalized number
//	 MatScattering = Amp * alpha^2 * s.xw^4 * s2*WKB_Comp// * Transmission

END




//Correct the optical constants with the WKB approximation
Function Ferron_WKBOC(s)
	STRUCT Ferron_FitStruct &s
	Variable Eoff1=s.pw[1], Eoff2=s.pw[2]
	Variable Thick= 467 //s.pw[3] //Be careful if pw is reorganized!  This will not be correct anymore!
	Variable i,j//, en //wave number 2pi/lambda

	//Make the slab wave that we integrate over
	Make/o/N=(Thick+1) s.SlabModel = p
	//Make temporary arrays before integration
	Make/o/N=((numpnts(s.en1)),Thick+1) WKB_ImagArray_1, WKB_RealArray_1
	Make/o/N=((numpnts(s.en2)),Thick+1) WKB_ImagArray_2, WKB_RealArray_2

	duplicate/o s.en1, koff1, RealWKBOC1, ImagWKBOC1
	duplicate/o s.en2, koff2, RealWKBOC2, ImagWKBOC2

	duplicate/o s.en1, s.OCReal1, s.OCImag1
	duplicate/o s.en2, s.OCReal2, s.OCImag2
	
	Variable hc=1239.84 //eV*nm
	koff1 = 2*PI*(s.en1)/hc
	koff2 = 2*PI*(s.en2)/hc
	
//	Setscale/I x, 0, th, s.WKBth, s.WKBthReal, s.WKBthImag
//Setup the first optical constant set
	s.OCReal1 = real(s.OC1)
	s.OCImag1 = Imag(s.OC1)
	
	WKB_RealArray_1[][] = s.OCReal1[p]*COS(koff1[p]*s.slabmodel[q]*s.OCReal1[p]) - s.OCImag1[p]*SIN(koff1[p]*s.slabmodel[q]*s.OCReal1[p])
	WKB_RealArray_1[][] *= exp(-koff1[p]*s.slabmodel[q]*s.OCImag1[p])
	
	WKB_ImagArray_1[][] = s.OCImag1[p]*COS(koff1[p]*s.slabmodel[q]*s.OCReal1[p]) + s.OCReal1[p]*SIN(koff1[p]*s.slabmodel[q]*s.OCReal1[p])
	WKB_ImagArray_1[][] *= exp(-koff1[p]*s.slabmodel[q]*s.OCImag1[p])
	
	//Setup the second optical constant set
	s.OCReal2 = real(s.OC2)
	s.OCImag2 = Imag(s.OC2)
	WKB_RealArray_2[][] = s.OCReal2[p]*COS(koff2[p]*s.slabmodel[q]*s.OCReal2[p]) - s.OCImag2[p]*SIN(koff2[p]*s.slabmodel[q]*s.OCReal2[p])
	WKB_RealArray_2[][] *= exp(-koff2[p]*s.slabmodel[q]*s.OCImag2[p])
	
	WKB_ImagArray_2[][] = s.OCImag2[p]*COS(koff2[p]*s.slabmodel[q]*s.OCReal2[p]) + s.OCReal2[p]*SIN(koff2[p]*s.slabmodel[q]*s.OCReal2[p])
	WKB_ImagArray_2[][] *= exp(-koff2[p]*s.slabmodel[q]*s.OCImag2[p])
	
	//Now Integrate the Thickness side of things using trapezoidal integration dt comes from spacing in wave 'slabmodel' right now its 1
	//So just summing the rows gives the same trapezoidal integration
	MatrixOP/O RealWKBOC1 = sumrows(WKB_RealArray_1)/Thick
	MatrixOP/O RealWKBOC2 = sumrows(WKB_RealArray_2)/Thick
	MatrixOP/O ImagWKBOC1 = sumrows(WKB_ImagArray_1)/Thick
	MatrixOP/O ImagWKBOC2 = sumrows(WKB_ImagArray_2)/Thick
	
	s.OC1WKB = cmplx(interp(s.xw+EOff1,s.en1,RealWKBOC1),interp(s.xw+EOff1,s.en1,ImagWKBOC1))
	s.OC2WKB = cmplx(interp(s.xw+EOff2,s.en2,RealWKBOC2),interp(s.xw+EOff2,s.en2,ImagWKBOC2))


//	
//	Make/Free/N=(th+1) TempBeta
//	Make/Free/N=(th+1) TempDelta
//	
//	//worst for loop ever made
//	for(j=1;j<=2;j+=1)
//		wave loop = $"s.en"+num2str(j)
//		for(i=0;i<(numpnts(loop));i+=1)
//		
//			TempBeta[] = WKBArray_Beta[i][p]
//			TempDelta[] = WKBArray_Delta[i][p]
//		
//			WKB_Beta[i] = AreaXY(slabmodel,TempBeta) 
//			WKB_Delta[i] = AreaXY(slabmodel,TempDelta)
//		
//		
//		endfor	
//	endfor
//	
//	WKB_Beta /= (thickness+1)
//	WKB_Delta /= (Thickness+1)
//	
//	
		
//	WKBArray_Beta[][] = (OC_Beta[p]*cos(2*PI*(Energy[p]/1240)*Slabmodel[q]*OC_Delta[p]) + OC_Delta[p]*Sin(2*PI*(Energy[p]/1240)*SlabModel[q]*OC_Delta[p]))*exp(-2*PI*(Energy[p]/1240)*SlabModel[q]*OC_Beta[p])
//	WKBArray_Delta[][] = (OC_Delta[p]*cos(2*PI*(Energy[p]/1240)*slabmodel[q]*OC_Delta[p]) - OC_Beta[p]*sin(2*PI*(Energy[p]/1240)*slabmodel[q]*OC_Delta[p]))*exp(-2*PI*(Energy[p]/1240)*slabmodel[q]*OC_Beta[p])


////////////////////SAVE THIS STUFF!!!!!!!!IT WORKS////////////////////////////////////////////////////
	//Setup the first optical constant set
//	s.OCReal1 = Real(s.OC1)
//	s.OCImag1 = Imag(s.OC1)
//	
//	RealWKBOC1 = s.OCReal1*COS(koff1*th*s.OCReal1) - s.OCImag1*SIN(koff1*th*s.OCReal1)
//	RealWKBOC1 *= exp(-koff1*th*s.OCImag1)
//	
//	ImagWKBOC1 = s.OCImag1*COS(koff1*th*s.OCReal1) + s.OCReal1*SIN(koff1*th*s.OCReal1)
//	ImagWKBOC1 *= exp(-koff1*th*s.OCImag1)
//	
//	s.OC1WKB = cmplx(interp(s.xw+EOff1,s.en1,RealWKBOC1),interp(s.xw+EOff1,s.en1,ImagWKBOC1))
//		//Setup the second optical constant set
//	s.OCReal2 = Real(s.OC2)
//	s.OCImag2 = Imag(s.OC2)
//	
//	RealWKBOC2 = s.OCReal2*COS(koff2*th*s.OCReal2) - s.OCImag2*SIN(koff2*th*s.OCReal2)
//	RealWKBOC2 *= exp(-koff2*th*s.OCImag2)
//	
//	ImagWKBOC2 = s.OCImag2*COS(koff2*th*s.OCReal2) + s.OCReal2*SIN(koff2*th*s.OCReal2)
//	ImagWKBOC2 *= exp(-koff2*th*s.OCImag2)
//	
//	s.OC2WKB = cmplx(interp(s.xw+EOff2,s.en2,RealWKBOC2),interp(s.xw+EOff2,s.en2,ImagWKBOC2))
	
////////////////////////////////////////////////////////////////////////////////////////////////////////
	
//	s.ImagWKBTemp = interp(s.xw+Eoff1,s.en1,s.OCImag)*COS(koff1*th*interp(s.xw+Eoff1,s.en1,s.RealWKBTemp))
	
//	s.OC1WKB = interp(s.xw+Eoff1,s.en1,s.imagWKBTemp)*COS(koff1*th*interp(s.xw+Eoff1,s.en1,s.RealWKBTemp)) )
//	s.OC1WKB += Cmplx(-interp(s.xw+Eoff1,s.en1,s.ImagWKBTemp)*COS(koff1*th*interp(s.xw+Eoff1,s.en1,s.realWKBTemp)),interp(s.xw+Eoff1,s.en1,s.RealWKBTemp)*COS(koff1*th*interp(s.xw+Eoff1,s.en1,s.RealWKBTemp)) )
//	s.OC1WKB *= cmplx(exp(-koff1*th*interp(s.xw+EOff1,s.en1,s.ImagWKBTemp)),0)
//	
	//Setup the second optical constant set
//	s.RealWKBTemp = Real(s.OC2)
//	s.ImagWKBTemp = Imag(s.OC2)
//	
//	
//	s.OC2WKB = cmplx(interp(s.xw+Eoff2,s.en2,s.RealWKBTemp)*COS(koff2*th*interp(s.xw+Eoff2,s.en2,s.realWKBTemp)),interp(s.xw+Eoff2,s.en2,s.imagWKBTemp)*COS(koff2*th*interp(s.xw+Eoff2,s.en2,s.RealWKBTemp)) )
//	s.OC2WKB += Cmplx(-interp(s.xw+Eoff2,s.en2,s.ImagWKBTemp)*COS(koff2*th*interp(s.xw+Eoff2,s.en2,s.realWKBTemp)),interp(s.xw+Eoff2,s.en2,s.RealWKBTemp)*COS(koff2*th*interp(s.xw+Eoff2,s.en2,s.RealWKBTemp)) )
//	s.OC2WKB *= cmplx(exp(-koff2*th*interp(s.xw+EOff2,s.en2,s.ImagWKBTemp)),0)
		
//	s.OC1WKB=exp( cmplx(0,1)*2*pi*s.xw/1240*th*cinterp(s.xw+Eoff1,s.en1,s.OC1) )*cinterp(s.xw+Eoff1,s.en1,s.OC1)
//	s.OC2WKB=exp( cmplx(0,1)*2*pi*s.xw/1240*th*cinterp(s.xw+Eoff2,s.en2,s.OC2) )*cinterp(s.xw+Eoff2,s.en2,s.OC2)

//	For( i=0; i<numpnts(s.xw); i+=1 ) //cycle through energies
//		en=s.xw[i]
//		k=2*pi*en/1240 // [1/nm]
//		//Material 1
//		s.WKBth= exp( cmplx(0,1)*k*x*cinterp(en+Eoff1,s.en1,s.OC1) )*cinterp(en+Eoff1,s.en1,s.OC1)
//		s.WKBthReal=real(s.WKBth)
//		s.WKBthImag=imag(s.WKBth)
//		s.OC1WKB[i]=cmplx(faverage(s.WKBthReal),faverage(s.WKBthImag))
//		//Material 2
//		s.WKBth= exp( cmplx(0,1)*k*x*cinterp(en+Eoff2,s.en2,s.OC2) )*cinterp(en+Eoff2,s.en2,s.OC2)
//		s.WKBthReal=real(s.WKBth)
//		s.WKBthImag=imag(s.WKBth)
//		s.OC2WKB[i]=cmplx(faverage(s.WKBthReal),faverage(s.WKBthImag))
//	endfor
end


//Displays the parameters in a table
FUNCTION Ferron_ParamTable(s)
	STRUCT Ferron_FitStruct &s
	DoWindow/F $(s.DataN+"_Parameters")
	IF( V_flag )
		return 0
	endif
	Edit/N=$(s.DataN+"_Parameters") /W=(375.75,445.25,592.5,611.75) s.pNames,s.pw,s.sw, s.hw, as s.DataN+" Parameters"
	ModifyTable size=8,format(Point)=1,width(Point)=12,width(s.pNames)=50,width(s.pw)=50, width(s.sw)=50, width(s.hw)=20
END


//Displays the fit in a graph
FUNCTION Ferron_PlotFitResults(s)
	STRUCT Ferron_FitStruct &s
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
FUNCTION Ferron_DisplayCorr(s)
	STRUCT Ferron_FitStruct &s
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

//combine the beta and delta waves into a complex OC wave
Function MakeCmplxOC(del,bet,name)
	WAVE del, bet
	String name
	String CurrentFolder=GetDataFolder(1)
	WAVE delU=$( nameOfWave(del) + "U"), betU=$( nameOfWave(bet)+"U")
	WAVE en=$(removeEnding( nameOfWave(del), "_delta")+"_e")
	Variable npts=numpnts(bet)
	NewDataFolder/O/S root:OpticalConstants
	Make/o/d/n=(npts)/C $name, $(name+"_u")
	WAVE/C OC=$name, OCu=$(name+"_u")
	string wnote=Note(del)
	OC=cmplx(del,bet)
	OCu=cmplx(delU,betU)
	Duplicate/d/o en $(name+"_en")
	Note/K/NOCR OC, wnote
	Note/K/NOCR OCu, wnote
	SetDataFolder $CurrentFolder
end

Function FitCompare(name, Escans [All, Svac, Smat, Thick, Conc, Fluor, inter, OOP,d,reset])
	string name // Name of whatever you want to call your folder where the waves will be stored
	string Escans // Should be a list of the different folders in Escans that you want to compare fit parameters
	variable All, Svac, Smat, Thick, Conc,OOP,Fluor,inter //Set equal to 1 if you want to graph them otherwise ignore.
	variable d , reset //Set if you want to display things and or reset everything if changes are made
	string CurrentFolder = GetDataFOlder(1), nfldr = "root:Escans:", namefldr = "root:FitCompare:" + name
	
	//Locate Folders that hold fit data
	SetDataFolder $nfldr
	String FldrList=ReplaceString(",",StringByKey("FOLDERS",DataFolderdir(1)),";") //Full List of files in Escans Folder
	//print FldrList
	variable i, j
	if(All ==1)
		Smat = 1
		Thick = 1
		Conc = 1
		Svac = 1
		OOP = 1
		Fluor = 1
		Inter = 1
	endif
	For(i=0;i<ItemsinList(Escans);i+=1)
		string Scan = ListMatch(FldrList,StringFromList(i,Escans))
		string Fit = ReplaceString(";", Scan, "")
		string scanfldr = nfldr + Fit +":Fit"
		if (!DataFolderExists(scanfldr))
			print "Fit TSI on " + fit + " Before Continuing"
			return 0
		endif
		SetDataFolder $scanfldr
		wave pwU = W_sigma
		wave pwfit = pw
		NewDataFolder/o/s $("root:FitCompare")
		NewDataFolder/o/s $namefldr
		//print ItemsinList(Escans)
		//Make Text Wave Listing Names for later plots
		Make/o/t/n = (itemsinList(Escans)) Samples
		Samples[i] = Fit
		//Make Svac wave
		if (Svac ==1)
			Make/o/n=(Itemsinlist(Escans)) Svac_C
			Make/o/n=(Itemsinlist(Escans)) Svac_U
			Svac_C[i] = pwfit[5]
			Svac_U[i] = pwU[5]
		endif
		//Make Smat wave
		if (Smat ==1)
			make/o/n=(ItemsinList(Escans)) Smat_C
			make/o/n=(ItemsinList(Escans)) Smat_U
			Smat_C[i] = pwfit[6]
			Smat_U[i] = pwU[6]
		endif
		//Make Thickness Wave
		if (Thick ==1)
			make/o/n=(ItemsinList(Escans)) Thick_C
			make/o/n=(ItemsinList(Escans)) Thick_U
			Thick_C[i] = pwfit[3]
			Thick_U[i] = pwU[3]
		endif		
		//Make Concentration Wave
		if (Conc ==1)
			make/o/n=(ItemsinList(Escans)) Conc_C
			make/o/n=(ItemsinList(Escans)) Conc_U
			Conc_C[i] = pwfit[4]
			Conc_U[i] = pwU[4]
		endif
		if (Fluor ==1)
			make/o/n=(ItemsinList(Escans)) Fluor_C
			make/o/n=(ItemsinList(Escans)) Fluor_U
			Fluor_C[i] = pwfit[9]
			Fluor_U[i] = pwU[9]
		endif
		//Always looks at Energy Offsets
		make/o/n=(ItemsinList(Escans)) PSOff_C
		make/o/n=(ItemsinList(Escans)) PSOff_U
		PSOff_C[i] = pwfit[1]
		PSOff_U[i] = pwU[1]
		make/o/n=(ItemsinList(Escans)) PMOff_C
		make/o/n=(ItemsinList(Escans)) PMOff_U
		PMOff_C[i] = pwfit[2]
		PMOff_U[i] = pwU[2]
	endfor
	
	wave knownconc
	if (!waveexists(knownconc)==1)
		Make/o/n=(itemsinlist(Escans)) knownconc
		edit/N=Concentrations Samples, knownconc
		NewPanel /W=(1000,50,250,140) /N=temppanel
		DrawText 40,30,"Please Input known concentrations"
		Button button1,pos={60,45},size={100,50},title="Continue",proc = KillTempPanel
		Pauseforuser temppanel, Concentrations	
		Killwindow Concentrations

	endif
	
	////////////////////////////////////////////////////
	//Gather OOP correction Factors
	If (OOP==1)
		string OOPFldr = namefldr+":" //"root:Escans:"+name+":"
		If (!DataFolderExists("root:q2D")==1)
			print "Run the Out of Plane Analysis on at least 1 sample to continue"
			return 0
		endif
		NewDataFolder/o/s root:q2D
		String QFldrList=ReplaceString(",",StringByKey("FOLDERS",DataFolderdir(1)),";") //Full List of files in Q2D Folder
		variable l
		Make/o/n=(itemsinlist(QFldrList)) $(OOPFldr+"OOPs")
		Make/o/t/n=(itemsinlist(QFldrList)) $(OOPFldr+"OOPSamples")

		Wave OOPCorr = $(OOPFldr+"OOPs")
		Wave/t OOPSamp = $(OOPFldr+"OOPSamples")

		for (l=0;l<itemsinlist(QFldrList);l+=1)
			string Qscan = StringFromList(l, QFldrList)
			SetDataFolder $("root:Q2D:"+Qscan +":"+Qscan+"_Integrate")
			variable/g ratio
			OOPCorr[l] = ratio
			OOPSamp[l] = Qscan		
		endfor
		
	endif
	//Gather Thicknesses from NEXAFS load files
	If (Thick ==1)
		string ThickFldr = namefldr +":"//"root:Escans:"+name+":"
		If (!DataFolderExists("Root:NEXAFS")==1)
			print "Import BL11 NEXAFS files"
			return 0
		endif
		NewDataFolder/o/s root:NEXAFS
		String TFldrList=ReplaceString(",",StringByKey("FOLDERS",DataFolderdir(1)),";") //Full List of files in Q2D Folder
		Make/o/n=(itemsinlist(TFldrList)) $(ThickFldr+"Thick_NXFS")
		Make/o/n=(itemsinlist(TFldrList)) $(ThickFldr+"Thick_NXFSU")
		
		Wave TNEXAFS = $(ThickFldr+"Thick_NXFS")
		Wave TNEXAFSU = $(ThickFldr+"Thick_NXFSU")
		for (i=0;i<itemsinlist(TFldrList);i+=1)
			string Tscan = StringFromList(i, TFldrList)
			SetDataFolder $("root:NEXAFS:"+Tscan)
			variable/g Thickness
			variable/g ThicknessU
				for (j=0; j<numpnts(Samples) ; j+=1)
					if( StringMatch((Samples[j]),Tscan+"*") ==1)
						TNEXAFS[j] = Thickness
						TNEXAFSU[j] = ThicknessU	
					endif
				endfor
		endfor
		
	endif
			
	///////////////////////////////////////////////////////////////////////////////////////////////	
	//Make a wave of the Smat values scaled based on the OOP correction AND INTERPHASE 5/31/2016
	NewDataFolder/o/s $namefldr
	if (OOP==1)
		duplicate/o Smat_C Smat_OOP
	
		variable m,n,o
		string SampleName
		for(m=0 ; m<numpnts(Samples) ; m+=1)//Cycle through Sample text wave
			for(n=0 ; n<(strlen(Samples[m]));n+=1) //Finds the Specific name of each sample (First 3 digits before underscore)
				if (cmpstr ((Samples[m])[n], "_") == 0)
					SampleName = (Samples[m])[0,n-1]
					//print SampleName
					break
				endif
			endfor
			//Compare Sample code (Ex. B43) with all OOP folders
			//print "here"
			for(o=0; o<numpnts(OOPSamp) ; o+=1)
				variable done = 0
				if( StringMatch(OOPSamp[o],(SampleName +"*")) ==1)
					//print SampleName
					Smat_OOP[m] *= OOPCorr[o]
					//print m, o 
					done =1
				endif
				if (done==1) // This will stop the loop if it finds a match already
					break
				endif
			endfor		
		endfor
	endif
	
	////////////////////////////////////////////////////////////
	//Make the ideal conditions in which to compare the product of volume fractions between
	SetDataFolder $namefldr	
	wave IdealSmat
	if(!waveexists(IdealSmat)==1)	
		Make/N=100 IdealConc
		Make/N=100 IdealSmatC
		Make/N=(numpnts(Samples)) IdealSmat
		if (knownconc[0] > 1)
			knownconc /= 100
		endif
		IdealConc = p/100
		IdealSmatC =  idealConc*(1 - Idealconc)
		IdealSmat = knownConc *(1-KnownConc)
		IdealConc *= 100
	endif
	
	
//Calculate the interphase correction for each sample given its morphology
	If (Inter==1)
		SetDataFolder $namefldr
		NewDataFolder/o/s InterPhase
		for(i=0;i<itemsinList(Escans);i+=1)
			string IScan = ReplaceString(";", (ListMatch(FldrList,StringFromList(i,Escans))), "")
			string IScanFldr = nfldr + Iscan
			Make/o/N=(itemsinList(Escans)) Interphase
			wave Morphology
			if (!waveexists(Morphology)==1) //Create Morphology Wave for future Use....Gives the user the chance to guess morphology
				Make/o/N=(itemsinList(escans)) Morphology
				edit/N=Morph  $("root:FitCompare:"+name+":samples"), Morphology
				NewPanel /W=(1283,386,1575,536) /N= TempPanel
				SetDrawLayer UserBack
				SetDrawEnv fsize= 14
				DrawText 23.5543105614553,88.1326905587186,"Please Record Expected Morphology \r 1 - Lamelle \r 2 - Cylindrical \r 3 - Spherical"
				Button button1,pos={148,59},size={84,62},proc=KillTempPanel,title="Continue"
				Pauseforuser temppanel, Morph	
				Killwindow Morph
			endif
			
			SetDataFOlder IScanFldr
			wave qvals, intProf
			wavestats/Q IntProf
			variable qmax = qvals[V_maxrowloc]
			
			SetDataFolder $(nameFldr+":Interphase")
			//Gathers Interphase Correction
			variable TempConc
			if(knownconc[i] >= 0.5)
				TempConc = 1-knownconc[i]
			else
				TempConc = knownconc[i]
			endif
			
			Interphase[i] = Calcinterphase(Morphology[i],qmax,TempConc)
			
			
			
		endfor
	//Scale the Ideal Smat Values by the calculated Interphase
		SetDataFolder $namefldr
		duplicate/o IdealSmat IdealSmat_Interphase
		IdealSmat_Interphase -= Interphase

	endif
		

	
	
	
	//////////////////////////////////////////////////////
	//Convert the Fluorescence value into an efficiency that should all be equal (WORK IN PROGRESS!!!)
	SetDataFolder $namefldr
	if (Fluor==1)
		NewDataFolder/o/s FluorCalc
		//print FldrList
		For(i=0;i<ItemsinList(Escans);i+=1)
			string FScan =ReplaceString(";", (ListMatch(FldrList,StringFromList(i,Escans))), "")
			string Fscanfldr = nfldr +FScan
			Make/o/N=(itemsinlist(Escans)) q2Int

			wave qtest = $(FscanFldr+":qvals"), TempQ, TempQ2
			duplicate/o qtest TempQ
			duplicate/o TempQ TempQ2
			TempQ2 = TempQ^2
			q2int[i] = AreaXY(TempQ, TempQ2)
		endfor
		SetDataFolder $namefldr
		duplicate/o Fluor_C Efficiency
		Efficiency /= 4*PI*q2int*Thick_C*10^(-9)
		//Tabulated PS-b-PMMA efficiecny 0.002575


	endif	
	
	//Now Display All the data Only a layout will show up
	SetDataFolder $namefldr
	if (d==1)
		if (Smat ==1)
			string SmatTitle = "Smat_" + name
			if (reset==1)
				DoWindow/K $SmatTitle
			endif
			Dowindow $SmatTitle
			if (V_flag == 0)
				display/hide=1/n=$smattitle Smat_C vs knownconc as "Concentration Comparison"
				appendtograph IdealSmatC vs IdealConc
				appendtograph IdealSmat vs Knownconc
				if (OOP==1)
					appendtograph Smat_OOP vs knownconc
					ModifyGraph marker(Smat_OOP)=17,rgb(Smat_OOP)=(65280,43520,0), mode(Smat_OOP) = 3
					Legend/C/N=text0/J/A=MC "\\s(Smat_C) Smat\r\\s(IdealSmat) Ideal Smat\r\\s(Smat_OOP) Smat OOP Corrected"
				endif
				if(Inter==1)
					appendtograph IdealSmat_Interphase vs knownconc
					ModifyGraph marker(IdealSmat_Interphase)=16, rgb(IdealSmat_Interphase)=(13056,13056,13056), mode(IdealSmat_Interphase) = 3
				endif
				ModifyGraph mode(Smat_C)=3,marker(Smat_C)=19
				ModifyGraph tick=2,mirror=1
				ModifyGraph lowTrip(left)=0.01
				ModifyGraph mode(IdealSmatC)=0,rgb(IdealSmatC)=(0,0,0)
				ModifyGraph marker(IdealSmat)=16, rgb(IdealSmat)=(26112,26112,26112), mode(IdealSmat)=3
				SetAxis left 0,0.3;SetAxis bottom 0,100
				Label left "Smat"
				Label bottom "M2 Conc"
				ErrorBars/W=$Smattitle/T=0.5/L=0.5/Y=3/x=3 Smat_C Y wave=(smat_u,smat_u)
			endif
		endif
		if (Thick ==1)
			string ThickTitle = "Thick_" + name
			if (reset==1)
				DoWindow/K $ThickTitle
			endif
			Dowindow $ThickTitle
			if (V_flag == 0)
				wave Thick_NXFS, Thick_NXFSU
				display/hide=1/n=$Thicktitle Thick_C vs Samples as "Thickness Comparison"
				appendtograph Thick_NXFS vs Samples
				ModifyGraph mode=3,marker=16
				ModifyGraph tick=2,mirror=1
				ModifyGraph rgb=(0,15872,65280)
				ModifyGraph rgb(Thick_NXFS)=(65280,0,0)
				Legend/C/N=text0/J/A=MC "\\s(Thick_C) Fit Thickness\r\\s(Thick_NXFS) NEXAFS Thickness"
				Label left "Thickness [nm]"
				Label bottom "Sample"
				ErrorBars/W=$thicktitle/T=0.5/L=0.5/Y=3/x=3 thick_C Y,wave=(Thick_U,Thick_U)
				ErrorBars/W=$thicktitle/T=0.5/L=0.5/Y=3/x=3 thick_C Y,wave=(Thick_NXFSU,Thick_NXFSU)

			endif
		endif
		
		if (Conc ==1)
			string ConcTitle = "Concentration_" + name
			if (reset==1)
				DoWindow/K $ConcTitle
			endif
			Dowindow $ConcTitle
			wave Knownperc = knownconc
			if (knownconc[2]<1)
				knownperc = knownconc*100
			endif
			if (V_flag == 0)
				display/hide=1/n=$Conctitle Conc_C vs Samples as "Concentration Comparison"
				appendtograph knownperc vs Samples
				ModifyGraph mode=3,marker=16
				ModifyGraph tick=2,mirror=1
				ModifyGraph rgb=(0,15872,65280)
				ModifyGraph rgb(Conc_C)=(65280,0,0)
				Label left "Concentration [%PMMA]"
				Label bottom "Sample"
				Legend/C/N=text0/J "\\s(Conc_C) Parameter Concentration\r\\s(knownconc) Real Concentration"		
				ErrorBars/W=$Conctitle/T=0.5/L=0.5/Y=3/x=3 Conc_C Y,wave=(Conc_U,Conc_U)
			endif
		endif
		
		if (Svac ==1)
			string SvacTitle = "Svac_" + name
			if (reset==1)
				DoWindow/K $SvacTitle
			endif
			Dowindow $SvacTitle
			if (V_flag == 0)
				display/hide=1/n=$Svactitle Svac_C vs Samples as "Roughness Comparison"
				ModifyGraph mode=3,marker=16
				ModifyGraph tick=2,mirror=1
				ModifyGraph rgb=(65280,43520,0)
				ModifyGraph prescaleExp(left)=6
				Label left "Roughness [\\u]"
				Label bottom "Sample"
				ErrorBars/W=$Svactitle/T=0.5/L=0.5/Y=3/x=3 Svac_C Y,wave=(Svac_U,Svac_U)
			endif
		endif
		if (OOP ==1)
			string OOPTitle = "OOP_" + name
			if (reset==1)
				DoWindow/K $OOPTitle
			endif
			Dowindow $OOPTitle
			if (V_flag == 0)
				display/Hide=1/n=$OOPtitle OOPCorr vs OOPSamp as "OOP Correction"
				ModifyGraph mode=3,marker=16
				ModifyGraph tick=2,mirror=1
				ModifyGraph rgb=(26368,0,52224)	
				ModifyGraph prescaleExp(left)=6
				Label left "OOP Correction Factor [\\u]"
				Label bottom "Sample"
				SetAxis left 0.7,1.3
				//ErrorBars/W=$Svactitle/T=0.5/L=0.5/Y=3/x=3 Svac_C Y,wave=(Svac_U,Svac_U)
			endif
		endif
		
		if (Fluor ==1)
			string FluorTitle = "Fluor_" + name
			if (reset==1)
				DoWindow/K $FluorTitle
			endif
			Dowindow $FluorTitle
			if (V_flag == 0)
				display/hide=1/n=$Fluortitle Fluor_C vs Samples as "Fluorescence Comparison"
				ModifyGraph mode=3,marker=16
				ModifyGraph tick=2,mirror=1
				ModifyGraph rgb=(0,52224,0)				
				Label left "Fluorescence [au]"
				Label bottom "Sample"
				ErrorBars/W=$Fluortitle/T=0.5/L=0.5/Y=3/x=3 Fluor_C Y,wave=(Fluor_U,Fluor_U)
			endif
		endif
		
		string ECTitle = "EC_" + name
			if (reset==1)
				DoWindow/K $ECTitle
			endif
		Dowindow $ECTitle
		if (V_flag == 0)
			display/hide=1/n=$ECtitle PSOff_C vs Samples as "Energy Offsets"
			appendtograph PMOff_C vs Samples
				
			ModifyGraph mode=3,marker=16
			ModifyGraph tick=2,mirror=1
			ModifyGraph rgb=(0,15872,65280)
			ModifyGraph rgb(PSOff_C)=(65280,0,0)
			Legend/C/N=text0/J/A=MC "\\s(PSOff_C) PS Offset\r\\s(PMOff_C) PMMA Offset"
			ModifyGraph lowTrip(left)=0.001,prescaleExp(left)=-1
			Label left "Energy Offset [eV]"
			Label bottom "Sample"
			ErrorBars/W=$ECtitle/T=0.5/L=0.5/Y=3/x=3 PSOff_C Y,wave=(PSOff_U,PSOff_U)
			ErrorBars/W=$ECtitle/T=0.5/L=0.5/Y=3/x=3 PMOff_C Y,wave=(PMOff_U,PMOff_U)

		endif
		Dowindow/F $Name
			if (reset==1)
				DoWindow/K $Name
			endif
		if (All ==1 & V_Flag==0)
			String LayoutName = "Parameter Comparison Layout for "+Name
			NewLayout/n=$Name as LayoutName
			appendlayoutobject/w=$name graph $Svactitle
			appendlayoutobject/w=$name graph $thicktitle
			appendlayoutobject/w=$name graph $Smattitle
			appendlayoutobject/w=$name graph $Conctitle
			appendlayoutobject/w=$name graph $OOPtitle
			appendlayoutobject/w=$name graph $FluorTitle
			appendlayoutobject/w=$name graph $ECTitle
			Execute "Tile"
		endif
	endif
		

end

Function Calcinterphase(m,q,phi)
	variable m,q,phi //m relates the type of morphology, 1 for lam, 2 for cyl, and 3 for sphere
				// q is the long period
				//phi is the concentration
	variable SV,Iphase, L

	If (phi > 1)
		phi /=100
	endif
	L = 2*PI/q
	//print L
	if (m==1)
		SV = 2*(1/L)
	endif
	if (m==2)
		SV = Sqrt((8/Sqrt(3)) * PI * phi)/L
	endif
	if(m==3)
		SV = 2*(9*PI*Phi^2)^(1/3)*(1/L)
	endif

	//return SV	
	//print SV
	//iphase = (1/6)*(SV)*5
	//print "Interphase =", Iphase
	//return Iphase
end


/////////////////////////////////////////////////////////////BUTTON CONTROLS FOR THE ABOVE

//Function KillTempPanel(name) : ButtonControl
//	string name
//	KillWIndow temppanel
//end