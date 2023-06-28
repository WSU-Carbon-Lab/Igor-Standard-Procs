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
	
	Variable rho1 = 1.04 // [g/cm^3] Density of monomer for block 1
	Variable mw1 = 104.15 //[g/mol] Molecular weight of repeat unit 
	Variable Mw1_mol = 14.9E3 //[g/mol] MOlecular weight of the full segment (for N calculation)
	
	Variable rho2 = 1.18 // [g/cm^3] Density of monomer for block 2
	Variable Mw2 = 100.12 //g/mol] Molecular weight of repeat unit 
	Variable Mw2_mol = 13.1E3 //[g/mol] MOlecular weight of the full segment (for N calculation)
	
	//Distribution Info
	Variable PDI = 1.05
	Variable NumberofPoints = 50
	

	Variable Precision = 0.01
	
	//Contants 
	Variable Na = (6.02214076E23) //Avogadros number
	s.volu = 0.118 //nm^3 From Bates ACS Macro Letters Viewpoint from 2015 (Taken from Table 1 as standard reference volume)
	s.PDI = PDI
	//Calculations from Inputs
	Variable N1 = Floor(Mw1_mol/Mw1) //Average Degree of Polymerization of block 1
	Variable N2 = Floor(Mw2_mol/Mw2) //Average Degree of Polymerization of block 2

	s.N1_Avg = N1
	s.N2_Avg = N2
	s.N_Avg = N1+N2
	//Calculate Distributions based on Polydispersity
	Initialize_Distribution("Block1", NumberofPoints , Precision , N1, PDI) //For block 1
	Initialize_Distribution("Block2", NumberofPoints , Precision , N2, PDI) //For block 2
	Initialize_Distribution("BCP" , NumberofPOints, Precision, s.N_Avg, PDI)
	
	Wave s.N1_Prob = $(CurrentFolder+"fit:TempDistributions:Block1_Dist")
	Wave s.N2_Prob = $(CurrentFolder+"Fit:TempDistributions:Block2_Dist")
	Wave s.N_Prob = $(CurrentFolder+"Fit:TempDistributions:BCP_Dist")
	
	Wave s.N1_Vals = $(CurrentFolder+"Fit:TempDistributions:Block1_N")
	Wave s.N2_Vals = $(CurrentFolder+"Fit:TempDistributions:Block2_N")
	Wave s.N_Vals = $(CurrentFolder+"Fit:TempDistributions:BCP_N")
	
//	Calculate_CombinedDistributions(s.N1_Prob,s.N2_Prob,s.N1_Vals,s.N2_Vals,Mw1=Mw1,Mw2=Mw2,rho1=rho1,rho2=rho2)
	
	Wave s.Rg2DMap = $(CurrentFolder+"Fit:TempDistributions:Rg2D")
	Wave s.SZ2DMap = $(CurrentFolder+"Fit:TempDistributions:Probability2D")

	//Calculate the volume fraction distributions for reference
	Variable/G vol1 = N1 * (Mw1/(Na*rho1)) * (1E21) //Volume of block 1 segment
	Variable/G vol2 = N2 * (Mw2/(Na*rho2)) * (1E21) //Volume of block 2 segment
	
	Variable/G vol = vol1 + vol2 //Volume of the total block copolymer
	
	Variable/G f1 = vol1/vol //calculate volume fraction of block 1
	Variable/G f2 = vol2/vol //calculate volume fraction of block 2

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
		s.pw={1,6.8,0.3,1,0,-4}  //add initial parameters values here.
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
		FuncFit/H=s.hold/M=2/Q BCP_SF_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /M=s.mw /R=s.rw /W=s.uw /I=1 /C=s.cw /STRC=s ///D=s.fw
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
			BCP_SF_FitFunc(s)
		endif
	else
		//fill the waves with the model (minimized operations to be a fit func)
		BCP_SF_FitFunc(s)
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
	
	//Make the RgMap
	
	//s.Rg2DMap *= s.pw[1]^2 //Finalize with the fit parameter b such that Rg^2 = N*b^2/6
	
	///////////////////////////////Calculate the Debye Functions for the Form-Factor, D(Rg1),D(Rg2),D(Rg)
	//Make Debye Functions
	Make/D/O/N=(numpnts(s.xw)) Debye_Rg1, Debye_Rg2, Debye_Rg
	//////////////////Make GMatrix for each thing
	
	Variable m = numpnts(s.xw)
	Variable N = numpnts(s.N1_Vals) //Same as N2_Vals
	Variable NwAvg = (s.N1_Avg + s.N2_Avg)*s.PDI
	Variable RgAvg = Sqrt(NwAvg * KuhnL^2 / 6) 
	
	Make/D/O/N=(M,N) GMatrix_Rg1, GMatrix_Rg2, GMatrix_Rg //Each GMatrix is a 2D matrix with the calculated scattering given N1 and N2 for each location int he distribution
	
	Calc_GMatrix_Rgx_2D(Gmatrix_Rg1,s.xw,s.N1_Vals,s.N2_Avg,KuhnL) //(N1+N2)D(Rg1)
	Calc_GMatrix_Rgx_2D(Gmatrix_Rg2,s.xw,s.N2_Vals,s.N1_Avg,KuhnL) //(N1+N2)D(Rg2)
	Calc_GMatrix_Rg_2D(Gmatrix_Rg,s.xw,s.N_Vals,KuhnL)   //(N1+N2)D(Rg)
	
	//Finalize the integral 
	//Duplicate/o s.Sz2DMap Temp2DIntegrand
	Duplicate/o S.N_Prob TempIntegrand
	
	//Temp2DIntegrand[][] = s.Sz2DMap[p][q]*BinWidth_N(s.N1_Vals,p)*BinWidth_N(s.N2_Vals,q) //O(N1)*O(N2)*DeltaN1*DeltaN2
	
	TempIntegrand = s.N1_Prob * Binwidth_N(s.N1_Vals,p)
//	MatrixOP/o TempGMatrix_Rg1 = GMatrix_Rg1 x TempIntegrand //s.N1Prob * BinWidth_N(s.N1_Vals,//Temp2DIntegrand
//	MatrixOP/o TempDebyeScattering_Rg1 = Sum(TempGMAtrix_Rg1)
//	Debye_Rg1[] = TempDebyeScattering_Rg1[0][p]
	MatrixOP/o Debye_Rg1 = GMatrix_Rg1 x TempIntegrand //s.N1Prob * BinWidth_N(s.N1_Vals,//Temp2DIntegrand
	Debye_Rg1 /= NwAvg
	
	TempIntegrand = s.N2_Prob * Binwidth_N(s.N2_Vals,p)
//	MatrixOP/o TempGMatrix_Rg2 = GMatrix_Rg2 x TempIntegrand // * Temp2DIntegrand
//	MatrixOP/o TempDebyeScattering_Rg2 = Sum(TempGMAtrix_Rg2)
//	Debye_Rg2[] = TempDebyeScattering_Rg2[0][p]
	MatrixOP/o Debye_Rg2 = GMatrix_Rg2 x TempIntegrand // * Temp2DIntegrand
	Debye_Rg2 /= NwAvg
	
	TempIntegrand = s.N_Prob * Binwidth_N(s.N_Vals,p)
//	MatrixOP/o TempGMatrix_Rg = GMatrix_Rg x TempIntegrand // * Temp2DIntegrand
//	MatrixOP/o TempDebyeScattering_Rg = Sum(TempGMAtrix_Rg)
//	Debye_Rg[] = TempDebyeScattering_Rg[0][p]	
	MatrixOP/o Debye_Rg = GMatrix_Rg x TempIntegrand // * Temp2DIntegrand
	Debye_Rg /= NwAvg
	
	////Calculate the FF
	
	s.fw = Debye_Rg
	s.fw /= ( (Debye_Rg1*Debye_Rg2) - (1/4)*(Debye_Rg - Debye_Rg1 - Debye_Rg2)^2 )

	//Calculate the SF
	
	s.fw -= 2*chi*s.vol
	s.fw = 1/s.fw
	

	//
	s.fw *= 4*PI^2
	s.fw /= (1239.8)^4
	s.fw *= 74960 //Delta n^2
	
//	Calc_BCP_Scattering(s.fw,s.xw,KuhnL,s.f1,s.vol,chi)
//	s.fw *= Scale
	
//	s.fw += A*s.xw + B //Line
//	s.fw += A*s.xw^p + B //Power Law
	
	
END

/////////////////////
Function Calc_GMatrix_Rgx_2D(OutWave,qw,N1_Vals,N2_Avg,bVal)

	Wave Outwave
	Wave qw
	Wave N1_Vals
	
	Variable N2_Avg
	Variable bVal //How to calculate Rg from N
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFOlder/o/s ScatteringCalc
	
	Make/D/O/N=(numpnts(qw)) TempWave //The temp scattering wave to populate the GMatrix during the loop
	
	Variable i, j //Loop variables for the two distributions
	
	Variable N1_Current //The current values for the Degree of Polymerization
	Variable N1_TempLow, N1_TempHigh //The upper and lower bins to average over
	
	//Start loop
	For(i=0;i < numpnts(N1_Vals) ; i += 1)
	
		N1_Current = N1_Vals[i]
		N1_TempLow = StartBinN(N1_Vals,i)
		N1_TempHigh = EndBinN(N1_Vals,i)
		

		Multithread Tempwave =  CalcDebyePoints_Rgx_2D(qw[p],N1_Current,N1_Templow,N1_TempHigh,N2_AVg,bval)
		
		Outwave[][i] = Tempwave[p]
			
		
	endfor
	
	SetDataFOlder $CurrentFOlder
end

ThreadSafe Function CalcDebyePoints_Rgx_2D(QValue, N1, N1min, N1Max, N2Avg, b)

	Variable QValue
	Variable N1, N1min, N1max
	Variable N2Avg
	Variable b

	//Calculate Rg from N values
	Variable Rg1 = Sqrt(N1*b^2/6)
	Variable Rg1min = Sqrt(N1min*b^2/6)
	Variable Rg1max = Sqrt(N1max*b^2/6)
	
	//Variable Rg2 = Sqrt(N2*b^2/6) (Used to calculate Rg and nothing else)
	Variable Rg2 = Sqrt(N2Avg*b^2/6)
	
	//Calculate the Average value of Rg for the denominator based on the location
	Variable Rgmin = Sqrt(Rg1min^2 + Rg2^2)
	Variable Rgmax = Sqrt(Rg1max^2 + Rg2^2)
	Variable Rg = (Rgmax + Rgmin) / 2

	//Calc things needed for the loops
	Variable NumberofSteps = Floor(3+abs(N1max - N1min)) //In incremends of N=1 with a min of 3
	if(NumberofSteps > 60) //Just not that many
		NumberofSteps = 60
	endif
	
	Variable step = (N1max - N1min)/(NumberofSteps-1)
	Step = Sqrt(step*b^2/6) //Convert Delta N into Delta Rg
	Variable i
	
	//Stuff to store the value in
	Variable result = 0 //Final answer
	Variable TempResult=0 //Each step that will sum to result
	
	For(i=0 ; i < numberofsteps ; i+=1)
		Rg1 = Rg1min + i*step
		
		Tempresult = 2*(exp(-(QValue^2 * Rg1^2)) + QValue^2 * Rg1^2 - 1)
		TempResult /= QValue^4
		TempResult /= Rg^4
		
		result += TempResult
	
	endfor
	Result /= NUmberofSteps //D(Rg1)
	
	Result *= (N1+N2Avg) //Inside the integrand
	
	Return result

end

Function Calc_GMatrix_Rg_2D(OutWave,qw,N_Vals,bVal)

	Wave Outwave
	Wave qw
	Wave N_Vals
	
	Variable bVal //How to calculate Rg from N
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFOlder/o/s ScatteringCalc
	
	Make/D/O/N=(numpnts(qw)) TempWave //The temp scattering wave to populate the GMatrix during the loop
	
	Variable i, j //Loop variables for the two distributions
	
	Variable N_Current //The current values for the polydispersities
	Variable N_TempLow, N_TempHigh //The upper and lower bins to average over
	
	//Start with N1 (I guess)
	For(i=0;i < numpnts(N_Vals) ; i += 1)
	
		N_Current = N_Vals[i]
		N_TempLow = StartBinN(N_Vals,i)
		N_TempHigh = EndBinN(N_Vals,i)
		
		
		//Now that we have the range of N that we are running, we calculate the scattering gi
	
		Multithread Tempwave = CalcDebyePoints_Rg_2D(qw[p],N_Current,N_Templow,N_TempHigh,bval)
		
		Outwave[][i] = Tempwave[p]
			
	endfor
	
	SetDataFOlder $CurrentFOlder
end

ThreadSafe Function CalcDebyePoints_Rg_2D(QValue, N, Nmin, NMax, b)

	Variable QValue
	Variable N, Nmin, Nmax

	Variable b

	//Calculate Rg from N values
	Variable Rgmin = Sqrt(Nmin*b^2/6)
	Variable Rgmax = Sqrt(Nmax*b^2/6)
	Variable Rg = 0
//	Variable Rg = (Rgmax + Rgmin) / 2

	//Calc things needed for the loops
	Variable NumberofSteps = Floor(3+abs((Nmax) - (Nmin))) //In incremends of N=1 with a min of 3
	if(NumberofSteps > 60) //Just not that many
		NumberofSteps = 60
	endif
	
	Variable step = ((Nmax) - (Nmin))/(NumberofSteps-1)
	Step = Sqrt(step*b^2/6) //Convert Delta N into Delta Rg
	Variable i
	
	//Stuff to store the value in
	Variable result = 0 //Final answer
	Variable TempResult //Each step that will sum to result
	
	For(i=0 ; i<numberofsteps ; i+=1)
		Rg = Rgmin + i*step
		
		Tempresult = 2*(exp(-(QValue^2 * Rg^2)) + QValue^2 * Rg^2 - 1)
		TempResult /= QValue^4
		TempResult /= Rg^4
		
		result += TempResult
	
	endfor
	Result /= NUmberofSteps
	Result *= (N)
	
	Return result

end









///////////////////////////
///////////////////////////
////////////////////////
/////////////////////////
//////////////////////////////
////////////////////////////////


//Actual Fit Function
FUNCTION BCP_SF_FitFunc1(s) : FitFunc
	STRUCT BCP_SF_FitStruct &s
	Variable Scale = s.pw[0]
	Variable KuhnL = s.pw[1]
	Variable chi = s.pw[2]
	
	//Bkg, might not be needed afterall
	Variable A = s.pw[3]
	Variable B = s.pw[4]
	Variable pow = s.pw[5]
	
	//Make the RgMap
	
	s.Rg2DMap *= s.pw[1]^2 //Finalize with the fit parameter b such that Rg^2 = N*b^2/6
	
	///////////////////////////////Calculate the Debye Functions for the Form-Factor, D(Rg1),D(Rg2),D(Rg)
	//Make Debye Functions
	Make/D/O/N=(numpnts(s.xw)) Debye_Rg1, Debye_Rg2, Debye_Rg
	//////////////////Make GMatrix for each thing
	
	Variable m = numpnts(s.xw)
	Variable N = numpnts(s.N1_Vals) //Same as N2_Vals
	Variable NwAvg = (s.N1_Avg + s.N2_Avg)*s.PDI
	
	Make/D/O/N=(N,N,M) GMatrix_Rg1, GMatrix_Rg2, GMatrix_Rg //Each GMatrix is a 3D matrix with the calculated scattering given N1 and N2 for each location int he distribution
	
	Calc_GMatrix_Rgx(Gmatrix_Rg1,s.xw,s.N1_Vals,s.N2_Vals,KuhnL) //(N1+N2)D(Rg1)
	Calc_GMatrix_Rgx(Gmatrix_Rg2,s.xw,s.N2_Vals,s.N1_Vals,KuhnL) //(N1+N2)D(Rg2)
	Calc_GMatrix_Rg(Gmatrix_Rg,s.xw,s.N1_Vals,s.N2_Vals,KuhnL)   //(N1+N2)D(Rg)
	
	//Finalize the integral 
	Duplicate/o s.Sz2DMap Temp2DIntegrand
	
	Temp2DIntegrand[][] = s.Sz2DMap[p][q]*BinWidth_N(s.N1_Vals,p)*BinWidth_N(s.N2_Vals,q) //O(N1)*O(N2)*DeltaN1*DeltaN2
	
	MatrixOP/o TempGMatrix_Rg1 = GMatrix_Rg1 * Temp2DIntegrand
	MatrixOP/o TempDebyeScattering_Rg1 = Sum(TempGMAtrix_Rg1)
	Debye_Rg1[] = TempDebyeScattering_Rg1[0][0][p]
	Debye_Rg1 /= NwAvg
	
	MatrixOP/o TempGMatrix_Rg2 = GMatrix_Rg2 * Temp2DIntegrand
	MatrixOP/o TempDebyeScattering_Rg2 = Sum(TempGMAtrix_Rg2)
	Debye_Rg2[] = TempDebyeScattering_Rg2[0][0][p]
	Debye_Rg2 /= NwAvg
	
	MatrixOP/o TempGMatrix_Rg = GMatrix_Rg * Temp2DIntegrand
	MatrixOP/o TempDebyeScattering_Rg = Sum(TempGMAtrix_Rg)
	Debye_Rg[] = TempDebyeScattering_Rg[0][0][p]
	Debye_Rg /= NwAvg
	
	////Calculate the FF
	
	s.fw = Debye_Rg
	s.fw /= ( (Debye_Rg1*Debye_Rg2) - (1/4)*(Debye_Rg - Debye_Rg1 - Debye_Rg2)^2 )

	//Calculate the SF
	
	s.fw -= 2*chi*s.vol
	s.fw = 1/s.fw
	

	//
	s.fw *= 4*PI^2
	s.fw /= (1239.8)^4
	s.fw *= 74960 //Delta n^2
	
//	Calc_BCP_Scattering(s.fw,s.xw,KuhnL,s.f1,s.vol,chi)
//	s.fw *= Scale
	
//	s.fw += A*s.xw + B //Line
//	s.fw += A*s.xw^p + B //Power Law
	
	
END

////////////////////////
//GMAtrix calculations//
////////////////////////

Function Calc_GMatrix_Rgx(OutWave,qw,N1_Vals,N2_Vals,bVal)

	Wave Outwave
	Wave qw
	Wave N1_Vals, N2_Vals
	
	Variable bVal //How to calculate Rg from N
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFOlder/o/s ScatteringCalc
	
	Make/D/O/N=(numpnts(qw)) TempWave //The temp scattering wave to populate the GMatrix during the loop
	
	Variable i, j //Loop variables for the two distributions
	
	Variable N1_Current, N2_Current //The current values for the polydispersities
	Variable N1_TempLow, N1_TempHigh, N2_TempLow,N2_TempHigh //The upper and lower bins to average over
	
	//Start with N1 (I guess)
	For(i=0;i < numpnts(N1_Vals) ; i += 1)
	
		N1_Current = N1_Vals[i]
		N1_TempLow = StartBinN(N1_Vals,i)
		N1_TempHigh = EndBinN(N1_Vals,i)
		
		For(j=0 ; j < numpnts(N2_Vals) ; j+=1)
			N2_Current = N2_Vals[j]
			N2_TempLow = StartBinN(N2_Vals,j)
			N2_TempHigh = EndBinN(N2_Vals,j)
			
			//Now that we have the range of N that we are running, we calculate the scattering gi
		
			Multithread Tempwave =  CalcDebyePoints_Rgx(qw[p],N1_Current,N1_Templow,N1_TempHigh,N2_Current,N2_Templow,N2_Temphigh,bval)
		
			Outwave[i][j][] = Tempwave[r]
			
			endfor
	endfor
	
	SetDataFOlder $CurrentFOlder
end

ThreadSafe Function CalcDebyePoints_Rgx(QValue, N1, N1min, N1Max, N2, N2min, N2Max, b)

	Variable QValue
	Variable N1, N1min, N1max
	Variable N2, N2min, N2max
	Variable b

	//Calculate Rg from N values
	Variable Rg1 = Sqrt(N1*b^2/6)
	Variable Rg1min = Sqrt(N1min*b^2/6)
	Variable Rg1max = Sqrt(N1max*b^2/6)
	
	//Variable Rg2 = Sqrt(N2*b^2/6) (Used to calculate Rg and nothing else)
	Variable Rg2min = Sqrt(N2min*b^2/6)
	Variable Rg2max = Sqrt(N2max*b^2/6)
	
	//Calculate the Average value of Rg for the denominator based on the location
	Variable Rgmin = Sqrt(Rg1min^2 + Rg1min^2)
	Variable Rgmax = Sqrt(Rg1max^2 + Rg2max^2)
	Variable Rg = (Rgmax + Rgmin) / 2

	//Calc things needed for the loops
	Variable NumberofSteps = Floor(3+abs(N1max - N1min)/9) //In incremends of N=1 with a min of 3
	if(NumberofSteps > 60) //Just not that many
		NumberofSteps = 60
	endif
	
	Variable step = (N1max - N1min)/(NumberofSteps-1)
	Step = Sqrt(step*b^2/6) //Convert Delta N into Delta Rg
	Variable i
	
	//Stuff to store the value in
	Variable result = 0 //Final answer
	Variable TempResult=0 //Each step that will sum to result
	
	For(i=0 ; i < numberofsteps ; i+=1)
		Rg1 = Rg1min + i*step
		
		Tempresult = 2*(exp(-(QValue^2 * Rg1^2)) + QValue^2 * Rg1^2 - 1)
		TempResult /= QValue^4
		TempResult /= Rg^4
		
		result += TempResult
	
	endfor
	Result /= NUmberofSteps //D(Rg1)
	
	Result *= (N1+N2) //Inside the integrand
	
	Return result

end
Function Calc_GMatrix_Rg(OutWave,qw,N1_Vals,N2_Vals,bVal)

	Wave Outwave
	Wave qw
	Wave N1_Vals, N2_Vals
	
	Variable bVal //How to calculate Rg from N
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFOlder/o/s ScatteringCalc
	
	Make/D/O/N=(numpnts(qw)) TempWave //The temp scattering wave to populate the GMatrix during the loop
	
	Variable i, j //Loop variables for the two distributions
	
	Variable N1_Current, N2_Current //The current values for the polydispersities
	Variable N1_TempLow, N1_TempHigh, N2_TempLow,N2_TempHigh //The upper and lower bins to average over
	
	//Start with N1 (I guess)
	For(i=0;i < numpnts(N1_Vals) ; i += 1)
	
		N1_Current = N1_Vals[i]
		N1_TempLow = StartBinN(N1_Vals,i)
		N1_TempHigh = EndBinN(N1_Vals,i)
		
		For(j=0 ; j < numpnts(N2_Vals) ; j+=1)
			N2_Current = N2_Vals[j]
			N2_TempLow = StartBinN(N2_Vals,j)
			N2_TempHigh = EndBinN(N2_Vals,j)
			
			//Now that we have the range of N that we are running, we calculate the scattering gi
		
			Multithread Tempwave = CalcDebyePoints_Rg(qw[p],N1_Current,N1_Templow,N1_TempHigh,N2_Current,N2_Templow,N2_Temphigh,bval)
		
			Outwave[i][j][] = Tempwave[r]
			
			endfor
	endfor
	
	SetDataFOlder $CurrentFOlder
end

ThreadSafe Function CalcDebyePoints_Rg(QValue, N1, N1min, N1Max, N2, N2min, N2Max, b)

	Variable QValue
	Variable N1, N1min, N1max
	Variable N2, N2min, N2max
	Variable b

	//Calculate Rg from N values
	Variable Rg1min = Sqrt(N1min*b^2/6)
	Variable Rg1max = Sqrt(N1max*b^2/6)
	
	//Variable Rg2 = Sqrt(N2*b^2/6) (Used to calculate Rg and nothing else)
	Variable Rg2min = Sqrt(N2min*b^2/6)
	Variable Rg2max = Sqrt(N2max*b^2/6)
	
	//Calculate the Average value of Rg for the denominator based on the location
	Variable Rgmin = Sqrt(Rg1min^2 + Rg1min^2)
	Variable Rgmax = Sqrt(Rg1max^2 + Rg2max^2)
	Variable Rg = (Rgmax + Rgmin) / 2

	//Calc things needed for the loops
	Variable NumberofSteps = Floor(3+abs((N1max+N2max) - (N1min+N1min))/9) //In incremends of N=1 with a min of 3
	if(NumberofSteps > 60) //Just not that many
		NumberofSteps = 60
	endif
	
	Variable step = ((N1max+N2max) - (N1min+N2min))/(NumberofSteps-1)
	Step = Sqrt(step*b^2/6) //Convert Delta N into Delta Rg
	Variable i
	
	//Stuff to store the value in
	Variable result = 0 //Final answer
	Variable TempResult //Each step that will sum to result
	
	For(i=0 ; i<numberofsteps ; i+=1)
		Rg = Rgmin + i*step
		
		Tempresult = 2*(exp(-(QValue^2 * Rg^2)) + QValue^2 * Rg^2 - 1)
		TempResult /= QValue^4
		TempResult /= Rg^4
		
		result += TempResult
	
	endfor
	Result /= NUmberofSteps
	Result *= (N1+N2)
	
	Return result

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
