#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "CollinsProcs_Igor7"



//Command function to manage the fit 
FUNCTION Fit_BCP_SF(yw,uw,xw [noFit, Init])
	WAVE yw, uw, xw
	Variable noFit //optionally just make the model based on the current parameters
	Variable init // reinitialize the background fitting waves (do this after altering the model & # parameters, for example)

	STRUCT BCP_SF_FitStruct s
	Variable nParams=6 //number of parameters possible in the fit (change as needed)
	Variable timer=startMStimer //time the entire operation
	String CurrentFolder=GetDataFolder(1), wnote=""
	NewDataFolder/O/S :Fit //make Fit subfolder & put all the stuff related to it in there
	
	//Prepare Data
	
	//Information on the Polymer chosen
	//Inputs by the User
	
	Variable rho1 = 1.05 // [g/cm^3] Density of monomer for block 1
	Variable mw1 = 104.15 //[g/mol] Molecular weight of repeat unit 
	Variable Mw1_mol = 13E3 // 14.9E3 //[g/mol] MOlecular weight of the full segment (for N calculation)
	
	Variable rho2 = 1.17 // [g/cm^3] Density of monomer for block 2
	Variable Mw2 = 100.12 //g/mol] Molecular weight of repeat unit 
	Variable Mw2_mol = 13E3 // 13.1E3 //[g/mol] MOlecular weight of the full segment (for N calculation)
	
//Inputs into the model
	s.PDI = 1.1 //Polydispersity from the manufacturer
	s.volu = 0.118 //nm^3 From Bates ACS Macro Letters Viewpoint from 2015 (Taken from Table 1 as standard reference volume)

	//Contants 
	Variable Na = (6.02214076E23) //Avogadros number
	
	//Calculations from Inputs
	s.N1_Avg = Floor(Mw1_mol/Mw1) //Average Degree of Polymerization of block 1
	s.N2_Avg = Floor(Mw2_mol/Mw2) //Average Degree of Polymerization of block 2

	s.N_Avg = s.N1_Avg+s.N2_Avg //Full molecule N used in the Nbar calculations from the equation

	//Calculate the volume fraction distributions for reference
	Variable/G vol1_mol = (Mw1/(Na*rho1)) * (1E21)
	Variable/G vol2_mol = (Mw2/(Na*rho2)) * (1E21)
	Variable/G vol1 = s.N1_Avg * Vol1_Mol //Volume of block 1 segment
	Variable/G vol2 = s.N2_Avg * vol2_mol //Volume of block 2 segment
	
	Variable/G vol = vol1 + vol2 //Volume of the total block copolymer
	
	Variable/G f1 = vol1/vol //calculate volume fraction of block 1
	Variable/G f2 = vol2/vol //calculate volume fraction of block 2
	
	
	s.f1 = f1
	s.f2 = f2
	
	s.vol1 = vol1
	s.vol2 = vol2

	s.vol = vol

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
		s.pw={0.22,0.66,0.03,-208,175,0.17}  //add initial parameters values here.
		s.pNames={"A","Rg","chi","A","B","p"}
		s.mw=1 //initially fitt all points
		hw={1,0,0,0,0,0} //initialize hold wave with open parameters
		s.pwub = s.pw + 0.2
		s.pwlb = s.pw - 0.2
	endif
	
	WAVE s.pw=pw, s.hw=hw, s.sw=W_Sigma, s.mw=mw, s.pwub=pwub, s.pwlb=pwlb
	WAVE/T s.pNames = pNames
	
	BCP_SF_ParamTable(s) //display the parameters in a table
	BCP_SF_PlotFitResults(s) //displays the data & fit if not already done
	
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
	
	//Set the bounds for masking based on cursers
	//NOT TESTING WITH MORE THAN ONE ENERGY
	//Quick check for cursers...mask out everything outside of cursers
	DoWindow/F $(s.DataN+"_Fitting")
	Variable MaskedPoints = 0
	if( numtype(hcsr(A))==0 && numtype(hcsr(B))==0)
		Variable MinVal, MaxVal, NumData
		Make/O/N=0 xw_fit,yw_fit,uw_Fit
		
		MinVal = min(pcsr(B),pcsr(A))
		MaxVal = max(pcsr(B),pcsr(A))
		NumData = MaxVal - MinVal
		
			Duplicate/O/R=[minVal,MaxVal] s.xw, xwTemp
			Duplicate/O/R=[minVal,MaxVal] s.yw, ywTemp
			Duplicate/O/R=[minVal,MaxVal] s.uw, uwTemp
			Concatenate/NP/O {xw_fit, xwTemp}, TempWv
			Duplicate/O TempWv,xw_fit
			Concatenate/NP/O {yw_fit, ywTemp}, TempWv
			Duplicate/O TempWv,yw_fit
			Concatenate/NP/O {uw_Fit, uwtemp}, TempWv
			Duplicate/O TempWv,uw_Fit

	//	wave s.yw_Initial = s.yw , s.xw_initial = s.xw , s.uw_initial = s.uw
		wave s.yw = yw_fit , s.xw = xw_fit , s.uw = uw_fit
		Redimension/d/n=(numpnts(yw_fit)) s.fw, s.rw, s.mw //remake waves
		s.mw=1
		MaskedPoints = 1
		//Plot a new graph of the masked data ONLY
		BCP_SF_PlotReducedFit(s)
		
	else
		MaskedPoints = 0
		
	endif
	
	//Do the Fit
	Variable/G V_FitError=0, V_FitOptions=8 //save the iterates for debugging
	sprintf wnote "DataName:%s;DateTime:%1.8e;", s.DataN, dateTime //optionally add model information here
	Note/K/NOCR s.fw, wnote
	If( !noFit )
		//Actual fit operation:
		wave mask
		FuncFit/H=s.hold/M=2/Q BCP_SF_FitFunc1, s.pw, s.yw /X=s.xw /D=s.fw /M=s.mw /R=s.rw /W=s.uw /I=1 /C=s.cw /STRC=s ///D=s.fw
//		FuncFit BCP_SF_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /STRC=s ///D=s.fw
	
		WAVE s.fw, s.pw, s.xw=$("::"+xwn) //stupid IGOR thing that loses the pointer after executing FuncFit
		wave s.xw = xw
		if(MaskedPoints)
			wave s.xw= xw_fit
		endif
		
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
			BCP_SF_DisplayCorr(s) //display the correlations
			s.nchisq = sqrt(V_chisq/(V_npnts-(V_nterms-V_nheld))) //reduced chi squared value for the fit
			s.rw/=s.uw //residuals in sigmas
			s.rw= s.mw==0 ? nan : s.rw //don't show residuals for masked points
			sprintf wnote, "Nchisq:%g;HoldStr:%s;InterNumber:%d;",s.nchisq,s.hold, s.fi.iternumber
			NOTE/NOCR s.fw, wnote
			BCP_SF_FitFunc1(s)
		endif
	else
		//fill the waves with the model (minimized operations to be a fit func)
		BCP_SF_FitFunc1(s)
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
STRUCTURE BCP_SF_FitStruct
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
	
	//Polymer Variables
	Variable rho1
	Variable rho2
	Variable PDI
	
	Variable f1
	Variable f2
	Variable vol1
	Variable vol2
	Variable vol
	Variable volu
	
	//Parameter Distributions
	Wave N1_Prob
	Wave N2_Prob
	Wave N_Prob
	Wave N1_Vals
	Wave N2_Vals
	Wave N_Vals
	
	//Median Values
	Variable N1_Avg
	Variable N2_Avg
	Variable N_Avg
	
	Wave Rg2DMap
	Wave SZ2DMap
	
ENDSTRUCTURE

FUNCTION BCP_SF_FitFunc(s) : FitFunc
	STRUCT BCP_SF_FitStruct &s
	Variable Scale = s.pw[0]
	Variable KuhnL = s.pw[1]
	Variable chi = s.pw[2]
	
	//Bkg, might not be needed afterall
	Variable A = s.pw[3]
	Variable B = s.pw[4]
	Variable pow = s.pw[5]
	
		
	
	///////////////////////////////Calculate the Debye Functions for the Form-Factor, D(Rg1),D(Rg2),D(Rg)
	//Make Debye Functions
	Make/D/O/N=(numpnts(s.xw)) Debye_Rg1, Debye_Rg2, Debye_Rg
	if(s.f1 > 0.5)
		Debye_Rg1 = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f2,2)
		Debye_Rg2 = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f2,1)
		Debye_Rg = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f2,0)

	else
		Debye_Rg1 = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f1,1)
		Debye_Rg2 = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f1,2)
		Debye_Rg = Calc_Debye_PD(s.xw[p],s.PDI,KuhnL,s.N1_Avg,s.N2_Avg,s.f1,0)
	
	endif	
			

	////Calculate the FF
	
	s.fw = Debye_Rg
	s.fw /= ( (Debye_Rg1*Debye_Rg2) - (1/4)*(Debye_Rg - Debye_Rg1 - Debye_Rg2)^2 )

	//Calculate the SF
	s.fw -= 2*chi*s.N_Avg//s.vol
	s.fw = s.N_Avg/s.fw //Not s.N_Avg

	//
	s.fw *= 4*PI^2
	s.fw /= (1239.8)^4
	s.fw *= 74960 //Delta n^2
	s.fw *= 1E7 //Convert from nm^-1 to cm^-1
	
	s.fw *= 0.69 //Reduction in Transmission
	
//	Calc_BCP_Scattering(s.fw,s.xw,KuhnL,s.f1,s.vol,chi)
	s.fw *= Scale
	
	s.fw += A*s.xw^pow + B //Power Law

//	s.fw += A*s.xw + B //Line
//	s.fw += A*s.xw^p + B //Power Law
	
	
END

Function Calc_Debye_PD(qValue,PDI,b,N1,N2,Vf,mode)

	Variable qValue //Qvalue
	Variable PDI //Polydispersity
	Variable b //Kuhn Length
	
	Variable N1, N2 //Degree of polymerization
	Variable Vf //Volume fraction of the low VF block
	
	Variable mode //thing you want to calculate //phi, 1-phi, 1
	
	///////
	Variable Result //Outvalue
	///////
	
	///////
	Variable N = N1+N2
	///////
	Wave Debuff = root:Data:Fit:zpp
	////
	
	//Calculate the Rg variable
	Variable Rg = Sqrt(b^2*N/6) //Outside the brackets
	Variable x = qValue^2*Rg^2
	
	//Constant variables for the different hypergeometric variables
	Variable k = (2-PDI)/(PDI-1) //The same for all of them
	Variable z = (1- (2*Vf) )/(1-Vf) //The same for everything
	Variable zp, zpp //Two parameters that change based on x1 / x2
	
	//Check what mode we are calculating
	If(mode == 0) //The calculation for Rg
		Result = (2*N/x^2) * (x + 1/(2-PDI)*( (1+x*(PDI-1))^(-k) -1) )
		Result /= N
		
	elseif(Mode == 1) //The calculation for Rg1 where phi <= 0.5	
		zp = 1 - (Vf/(1-Vf))*( 1 + Vf*x*(PDI-1) )^(-1)
	
		Result = ( (1 + Vf*x*(PDI-1))^(-(k+1))*HyperG2F1(k+1,1,(2*k)+2,zp) ) //Square brackets first
		Result -= HyperG2F1(1,k+1,(2*k)+2,z)
		Result *= (1/((1-Vf)*(3-PDI)))
		Result += (x/2)*(Vf/(1-Vf))*HyperG2F1(1,k+2,(2*k)+3,z)
		Result *= (2*N/x^2)
		Result /= N

	elseif(Mode == 2)
		zpp = 1 - (Vf/(1-Vf))*( 1 + (1-Vf)*x*(PDI-1) )
		////////////////
				
		Result = ( (1 + (1-Vf)*x*(PDI-1) )^(-k) )*HyperG2F1(1,k+1,(2*k)+2,zpp)
		Result -= HyperG2F1(k+1,1,(2*k)+2,z)
		Result *= (1/((1-Vf)*(3-PDI)))
		Result += (x/2)*HyperG2F1(1,k+1,(2*k)+3,z)
		Result *= (2*N/x^2)
		Result /= N
	
	endif
	Return Result

end



FUNCTION BCP_SF_FitFunc1(s) : FitFunc
	STRUCT BCP_SF_FitStruct &s
	Variable Scale = s.pw[0]
	Variable KuhnL = s.pw[1]
	Variable chi = s.pw[2]
	
	//Bkg, might not be needed afterall
	Variable A = s.pw[3]
	Variable B = s.pw[4]
	Variable pow = s.pw[5]
	
		
	
	///////////////////////////////Calculate the Debye Functions for the Form-Factor, D(Rg1),D(Rg2),D(Rg)
	//Make Debye Functions
	Make/D/O/N=(numpnts(s.xw)) Debye_RgA1, Debye_RgA2, Debye_RgB1, Debye_RgB2

	Variable x1 = s.N1_Avg*KuhnL^2/6 //Not including q^2
	Variable x2 = s.N2_Avg*KuhnL^2/6 //not including q^2
	Variable v0 = (s.vol1*s.vol2)^(1/2)
	Variable rcn = (s.vol1/v0)*s.N1_Avg + (s.vol2/v0)*s.N2_Avg

	Debye_RgA1 = Calc_Debye_PD_1(s.xw[p],s.PDI,x1) 
	Debye_RgA2 = Calc_Debye_PD_2(s.xw[p],s.PDI,x1) 

	Debye_RgB1 = Calc_Debye_PD_1(s.xw[p],s.PDI,x2) 
	Debye_RgB2 = Calc_Debye_PD_2(s.xw[p],s.PDI,x2) 
	
	s.fw = ( rcn*s.f1^2*Debye_RgA2 + 2*rcn*s.f1*s.f2*Debye_RgA1*Debye_RgB1 + rcn*s.f2^2*Debye_RgB2 )
	s.fw /= (  (rcn*s.f1^2*Debye_RgA2*rcn*s.f2^2*Debye_RgB2) - (rcn*s.f1*s.f2*Debye_RgA1*Debye_RgB1)^2 )

	s.fw -= 2*chi
	
	s.fw = 1 / s.fw
	
		//
	s.fw *= 4*PI^2
	s.fw /= (1239.8)^4
	s.fw *= 74960 //Delta n^2
	s.fw *= 1E7 //Convert from nm^-1 to cm^-1
	
	s.fw *= 0.69 //Reduction in Transmission
	
//	Calc_BCP_Scattering(s.fw,s.xw,KuhnL,s.f1,s.vol,chi)
	s.fw *= Scale
	
	s.fw += A*s.xw^pow + B //Power Law
//	s.fw += A*s.xw + B //Line




end



Function Calc_Debye_PD_1(qValue,PDI,xVal)

	Variable qValue //Qvalue
	Variable PDI //Polydispersity
	Variable xVal //either x1 or x2
	
	Variable xx = xVal*qValue^2
	Variable result
	
	result = ( xx*(PDI-1)+1 )^( -(PDI-1)^(-1) )
	result = (1/xx)*(1 - result)
	
	return result


end
Function Calc_Debye_PD_2(qValue,PDI,xVal)

	Variable qValue //Qvalue
	Variable PDI //Polydispersity
	Variable xVal //either x1 or x2
	
	Variable xx = xVal*qValue^2
	Variable result
	
	result = ( xx*(PDI-1)+1 )^( -(PDI-1)^(-1) )
	result = (2/xx^2)*(-1 + xx + result)
	
	return result


end
//////////////////////////////////////////////////////
//DISTRIBUTION FUNCTIONS FOR INCLUDED POLYDISPERSITY//
//////////////////////////////////////////////////////
Function Initialize_Distribution(OutputName,NumberofPoints,Precision,AveMw,PDI)

	String OutputName
	
	Variable NUmberofPoints //NUmber of points in the distribution
	Variable Precision //The lowest % tolerable in distribution
	Variable AveMw //The average Molecular weight (number average)
	Variable PDI //Polydispersity
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFolder/o/s TempDistributions
	
	//MakeDistribution Wave
	
	Make/O/N=(NumberofPoints) $(OutputName + "_N"), $(OutputName+"_Dist")
	Wave Prob = $(OutputName+"_Dist")
	Wave PolyDeg = $(OutputName + "_N")
	
	Calc_N_Distribution(PolyDeg,NUmberofPOints,Precision,AveMw,PDI)
	
	Prob = CalcSZ_Distribution(PolyDeg[p],AveMw,PDI)
	
	Variable ScaleDist = 1/AreaXY(PolyDeg,Prob,-inf,inf)
	
	Prob *= ScaleDist
	
	SetDataFOlder $CurrentFOlder


end

//
Function Calculate_CombinedDistributions(N1Dist,N2Dist,N1Vals,N2Vals,[Mw1,Mw2,rho1,rho2])

	Wave N1Dist, N2Dist //Probabilty distributions
	Wave N1Vals, N2Vals //The probability distri ution x-waves basically
		
	Variable Mw1,Mw2 //Optional Molecular weight parametrs for calculate volume fraction
	Variable rho1,rho2
	
	Mw1 = ParamisDefault(Mw1) ? 1 : Mw1 
	Mw2 = ParamisDefault(Mw2) ? 1 : Mw2 
	rho1 = ParamisDefault(rho1) ? 1 : rho1 
	rho2 = ParamisDefault(rho2) ? 1 : rho2 
	Variable Na = 6.02214076E23
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFolder/o/s TempDistributions
	
	//Make the 2D map waves
	
	Make/N=(numpnts(N1Dist),numpnts(N2Dist))/O/D VolumeFraction1_2D,VolumeFraction2_2D, Rg2D, Probability2D


	//First make the Probability2D Wave
	Probability2D[][] = N1Dist[p]*N2Dist[q]
	//Calculate the Rg Value at each location given an R1 and R2
	Rg2D[][] = Sqrt(N1Vals[p]^2 + N2Vals[q]^2)/6 //DOes not include the Kuhn Length Squared***
	//Calculate the Volume fraction
	VolumeFraction1_2D[][] = ( N1Vals[p] * (Mw1/(Na*rho1)) ) /( (N1Vals[p]) * (Mw1/(Na*rho1)) + (N2Vals[q]) * (Mw2/(Na*rho2)))
	VolumeFraction2_2D = 1 - VolumeFraction1_2D


end

	Wave N1Dist, N2Dist //The two distributions
///

Function Calc_N_Distribution(OutputWave, NumberofPoints,Precision,AveMw,PDI)
	Wave OutputWave //The final wave that you want the distribution to be set in
	
	Variable NUmberofPoints //NUmber of points in the distribution
	Variable Precision //The lowest % tolerable in distribution
	Variable AveMw //The average Molecular weight (number average)
	Variable PDI //Polydispersity
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFolder/o/s TempDistributions
	
	Variable startx, endx, guess, Step, mode, TempVal, tempresult
	Step = 5
	Mode = AveMw
	
	//look for Min
	Variable MinimumPossible = 1 //1 monomer per chain is the limit
	
	TempVal = mode
	
	do
		TempVal = TempVal - Step
		
		if (tempval < MinimumPossible)
			TempVal = MinimumPossible
		endif
		
		TempResult = SZ_Cumulative(tempval,AveMw,PDI)
		
	while ((tempresult > Precision) && (tempval > MinimumPossible))
	startx = tempval
	
	//Look for Max
	Variable MaximumPossible = 10000 //Not going to have this many monomers
	do
		TempVal = TempVal + Step
		
		if (tempval > MaximumPossible)
			TempVal = MaximumPossible
		endif
		
		TempResult = SZ_Cumulative(tempval,AveMw,PDI)
		
	while ((tempresult < (1-Precision)) && (tempval < MaximumPossible))
	endx = tempval
	
	//Now make the data

	Make/O/D/N=(numberofpoints) Temp_CumulTargets
	Make/O/D/N=(3*Numberofpoints) Temp_CumulativeWave, Temp_Polymerization
	
	Temp_Polymerization = startx + p*(endx-startx)/(3*numberofpoints - 1)
	Temp_CumulTargets=Precision+p*(1-Precision-Precision)/(numberOfPoints-1) //this puts equally spaced values between myprecision and (1-myprecision) in this wave

	//Build the temp wave
	
	Temp_CumulativeWave = SZ_Cumulative(Temp_Polymerization,AveMw,PDI)
	//Make the last wave
	
	OutputWave = interp(Temp_CumulTargets, Temp_CumulativeWave, Temp_Polymerization)
	
	SetDataFolder $CurrentFOlder
	
end

Function SZ_Cumulative(xpos,AveMw,PDI)

	Variable xpos, AveMw, PDI
	
	Variable result
	Make/O/N=500/Free TempSzDist
	SetScale/I x 0,xpos,"",TempSzDist
	
	TempSzDist = CalcSZ_Distribution(x,AveMw,PDI)
	
	Result = Area(TempSZDist)
	return result

end

Function CalcSZ_Distribution(x,AveMw,PDI)

	Variable x,AveMw, PDI
	
	Variable Result, a, b //b is equal to k for the PDI version
	
	b = (PDI - 1)^(-1)
	a = b / AveMw

	Result = (  (a^(b+1))/gamma(b+1) * x^b / exp(a*x)  )
	
	if (numtype(result)!=0)
		result = 0
	endif
	
	return result
end
///


/////////////////////////////////////////
//Fit stuff needed to calculate Gmatrix//
/////////////////////////////////////////

Function StartBinN(NWave,i)
	Wave NWave
	Variable i
	
	variable start
	variable Imax=numpnts(NWave)
	
	if (i==0)
		start=NWave[0] - (NWave[1]-NWave[0])/2
		if (start<0)
			start=1		//we will enforce minimum size of the scatterer as 1 molecule
		endif
	elseif (i==Imax-1)
		start=NWave[i]-(NWave[i]-NWave[i-1])/2
	else
		start=NWave[i]-((NWave[i]-NWave[i-1])/2)
	endif
	return start


end

Function EndBinN(NWave,i)

	Wave NWave
	Variable i
	
	variable endL
	variable Imax=numpnts(NWave)
	
	if (i==0)
		endL=NWave[0]+(NWave[1]-NWave[0])/2
	elseif (i==Imax-1)
		endL=NWave[i]+((NWave[i]-NWave[i])/2)
	else
		endL=NWave[i]+((NWave[i+1]-NWave[i])/2)
	endif
	return endL
	
	
end

Function BinWidth_N(N_Distribution,i)			//calculates the width in diameters by taking half distance to point before and after
	variable i								//returns number in A
	Wave N_distribution
	
	variable width
	variable Imax=numpnts(N_distribution)
	
	if (i==0)
		width=N_distribution[1]-N_distribution[0]
		if ((N_distribution[0]-(N_distribution[1]-N_distribution[0])/2)<0)
			width=N_distribution[0]+(N_distribution[1]-N_distribution[0])/2
		endif
	elseif (i==Imax-1)
		width=N_distribution[i]-N_distribution[i-1]
	else
		width=((N_distribution[i]-N_distribution[i-1])/2)+((N_distribution[i+1]-N_distribution[i])/2)
	endif
	return abs(width)		//9/17/2010, fix for user models when bins are sorted from large to small
end





//////////////////////////////////////////////////
//Nothing Fit related below...Only visualization//
//////////////////////////////////////////////////

//Displays the parameters in a table
FUNCTION BCP_SF_ParamTable(s)
	STRUCT BCP_SF_FitStruct &s
	DoWindow/F $(s.DataN+"_Parameters")
	IF( V_flag )
		return 0
	endif
	Edit/N=$(s.DataN+"_Parameters") /W=(375.75,445.25,592.5,611.75) s.pNames,s.pw, s.pwlb,s.pwub, s.sw, s.hw, as s.DataN+" Parameters"
	ModifyTable size=8,format(Point)=1,width(Point)=12,width(s.pNames)=50,width(s.pw)=50, width(s.sw)=50, width(s.hw)=20
END


//Displays the fit in a graph
FUNCTION BCP_SF_PlotFitResults(s)
	STRUCT BCP_SF_FitStruct &s
	String ywN=NameOfWave(s.yw), fwN=NameOfWave(s.fw), lgndTxt
	DoWindow/F $(s.DataN+"_Fitting")
	If( V_Flag )
		return 0
	endif
	Display/N=$(s.DataN+"_Fitting") /W=(475.5,106.25,733.5,415.25) s.yw vs s.xw as s.DataN+" Fitting"
	AppendToGraph s.fw vs s.xw
	ModifyGraph rgb($fwN)=(0,0,0), rgb($ywN)=(65280,0,0)
	ModifyGraph mode($ywn)=3, marker($ywn)=19
	//Add Residuals to a graph in top 20% of window
	String rwN=NameOfWave(s.rw)
	AppendToGraph/L=Res s.rw vs s.xw
	ModifyGraph zero(Res)=1, lblPosMode(Res)=1, freePos(Res)=0, rgb($rwN)=(0,0,65280)
	Label Res "Res [sig]"
	ModifyGraph axisEnab(left)={0,0.8}
	ModifyGraph axisEnab(Res)={0.8,1}

	ModifyGraph margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph grid=2, tick=2, mirror=1, minor=1, standoff=0, lblPos(left)=42, msize=1.5
	Modifygraph log(bottom)=1
	Label left "Scattering Intensity [au]"
	Label bottom "q [nm\S-1\M]"
	ErrorBars/T=0 $ywN Y,wave=(s.uw,s.uw)
	Cursor/P/H=2 A $(fwN) 0 //place a cursor with a vertical line on the graph
	ShowInfo
	sprintf lgndTxt, "%s\r\\s(%s) Data\r\\s(%s) Fit",s.DataN,ywN,fwN
	Legend/C/N=text0/J/A=MC/X=23.02/Y=37.08 lgndTxt

END

//Display masked fit results
Function BCP_SF_PlotReducedFit(s)
	STRUCT BCP_SF_FitStruct &s
	Variable update
	String ywN=NameOfWave(s.yw), fwN=NameOfWave(s.fw), lgndTxt
	DoWindow/F $(s.DataN+"_Fitting_MaskedData")
	If( V_Flag && update == 0)
		return 0
	endif
	//Add the initial non-cut data to pick from
	Display/N=$(s.DataN+"_Fitting_MaskedData") /W=(475.5,106.25,733.5,415.25) s.yw vs s.xw as s.DataN+" Fitting Masked Data"
	AppendToGraph s.fw vs s.xw
	ModifyGraph rgb($fwN)=(0,0,0), rgb($ywN)=(65280,0,0)
	ModifyGraph mode($fwN)=3,marker($fwN)=19,msize($fwN)=3
	ModifyGraph mode($ywN)=3,msize($ywN)=1.5
	
	//Optionally append other waves here if this is a multi wave fit
	
	ModifyGraph margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph grid=2, tick=2, mirror=1, minor=1, standoff=0, lblPos(left)=42, msize=1.5
	ModifyGraph log=1
	Label left "Scattering Intensity [au]"
	Label bottom "Q [Å\S-1\M]"
	ErrorBars $ywN Y,wave=(s.uw,s.uw)
	Cursor/P/H=2 A $(fwN) 0 //place a cursor with a vertical line on the graph
	ShowInfo
	sprintf lgndTxt, "%s\r\\s(%s) Data\r\\s(%s) Fit",s.DataN,ywN,fwN
	Legend/C/N=text0/J/A=MC/X=23.02/Y=37.08 lgndTxt

	//Add Residuals to a graph in top 20% of window
	String rwN=NameOfWave(s.rw)
	AppendToGraph/L=Res s.rw vs s.xw
	ModifyGraph zero(Res)=1, lblPosMode(Res)=1, freePos(Res)=0, rgb($rwN)=(1,16019,65535), mode($rwn)=3, msize($rwn)=1.5
	ModifyGraph tick=2,mirror=1,standoff=0
	Label Res "Res [sig]"
	ModifyGraph axisEnab(left)={0,0.8}
	ModifyGraph axisEnab(Res)={0.8,1}
END


//Displays the correlation matrix
FUNCTION BCP_SF_DisplayCorr(s)
	STRUCT BCP_SF_FitStruct &s
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
