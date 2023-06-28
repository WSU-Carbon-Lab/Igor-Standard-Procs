#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "CollinsProcs_IGOR7"
#Include "MEF_Distributions"
#Include "MEF_SCatteringFunc"
#Include "MEF_FitFunctions"

//Generic code to get a structure-based fit function started. Containes code for modeling, fitting, debugging, displaying, and assesing
//Replace the word "Generic" with a unique name, and add the model details

//Command function to manage the fit
FUNCTION Fit_MultiEnScattering(yw,uw,xw, [noFit, init,test]) //optionally add parameters indicating libraries to use or other datasets to simultaneously fit
	WAVE yw, uw, xw
	Variable noFit //optionally just make the model based on the current parameters
	Variable init //reinitialize the background fitting waves (do this after altering the model in some way)
	Variable test
	STRUCT MultiEnScatt_FitStruct s
	Variable nParams=0//16 //number of parameters possible in the fit (change as needed)
	Variable timer=startMStimer //time the entire operation
	String CurrentFolder=GetDataFolder(1), wnote=""
	NewDataFolder/O/S :Fit //make Fit subfolder & put all the stuff related to it in there
	String FitDataFolder = GetDataFolder(1) //allows me to return to the fit folder whenever I want

	Variable i,j //loop parameters
	//Create Folders Required
	NewDataFolder/o root:Packages:MEF_Vars

	//Prepare Data
	//If doing a simultaneous fit, concatenate the data here into a long 1D wave
/////////////////////
//Fit Functionality//
/////////////////////
	Variable NumIterations = 80 //Set the number of iterations you want (Default 40)
//	Wave FitTolerance
	Variable FitTolerance = 1E-4
		
//*****************************************************************
//Let the code know what kind of fit you want to have done
//INitialize the number of pops and give them appropriate names
//*****************************************************************
	s.numpops = 2//HOW MANY POPULATIONS DO YOU WANT //I(q) = I(pop1) + I(pop2) + .....
	S.En = 5//HOW MANY ENERGIES DID YOU STITCH TOGETHER IN YOUR IMPORT
//I am current running under the assumption that the number of points per energy is the same

//********************************************
//Fill in types of models you want to fit
//Form factor is the basic model (may include things that are not technically form factors)
//Structure Factors are additional scales...See IRENA manual for possible Structure Factors

	//Available Form Factor Models
	//CoreSHell      - This is the same Core_Shell model as IRENA but with complex optical constants
	//CoreShell_Mod1 - Distribution of the core and shell width...basi function copied from IRENA
	//CoreShell_Mod2 - Include the SF into the distribution popuation, better connect SF to size of object
	//Power_Law      - Decay function based on power with included linear term f(x) = Ax^p + Bx + C
	//Sphere         - Its a sphere, because of course it is
	
//Assign populations to available models
//Will only uses models up to 's.numpops' all else ignored
	s.PopModel_FF_1 = "CoreShell"
	s.PopModel_SF_1 = "HardSpheres"
	s.PopModel_FF_2 = "Power_Law"
	s.PopModel_SF_2 = ""
	s.PopModel_FF_3 = ""
	s.PopModel_SF_3 = ""
	s.PopModel_FF_4 = ""
	s.PopModel_SF_4 = ""
	s.PopModel_FF_5 = ""
	s.PopModel_SF_5 = ""
//Sets all of these selections to become globals	
	MultiEnScatt_SetModelFits(s) 
	
	//Initialize some fit requirements here
	//Collect the end point for each energy dataset
	Make/N=(s.En+1)/O EndEnLoc
	wave  s.EndenLoc = EndEnLoc
	s.EndEnLoc[0] = 0
	//Assumes equal number of points in each scan
	s.EndEnloc = p*numpnts(xw)/s.En
	
	s.numdatasets = s.En //I don't think this does anything, but I don't want to uncomment it....
	//Datasets should equal the number of energies...
	
	//Check for correct initialization parameters
	if(s.numpops < 1 || s.numpops > 5)
		SetDataFolder $CUrrentFolder
		Print "FIT STOPPED"
		Abort "ERROR: Set number of populations between 1 and 5 -- See line 27"
	elseif(s.En < 1)
		SetDataFolder $CurrentFOlder
		Print "FIT STOPPED"
		Abort "ERROR: Number of energies is set below 0...please fix -- See line 28"
	endif
		
	////////////////////
	//Prep some things//
	////////////////////
	WAVE s.yw=yw, s.uw=uw, s.xw=xw //if doing simultaneous fit, this should be the concatenated data
	String xwn=NameOfWave(s.xw) //stupid Igor bug
	Duplicate/d/o s.yw, s.fw, s.rw
	s.DataN=GetWavesDataFolder(s.yw,0) //assumes the data is in a folder with a descriptive name for the data	
	
	//Prepare Model: Save the model parameters to the structure and initialize the background waves
	WAVE/Z pw
	If( !WaveExists(pw) || init ) //initialize fit waves
		Make/o/n=(nParams) pw, pwlb, pwub
		Redimension/D pw, pwlb, pwub
		Make/o/t/n=(nParams) pNames, pGlobal	
		WAVE/T s.pNames, s.pGlobal
		Duplicate/d/o s.yw, s.mw
		WAVE s.pw, s.mw, s.pwlb, s.pwub
		s.nParams = MultiEnScatt_InitializePW(s)
		duplicate/d/o s.pw s.hw, W_Sigma
		W_Sigma=0
		s.mw=1 //initially fitt all points
		s.hw=1 //initialize hold wave with open parameters		
	endif
	WAVE s.pw=pw, s.hw=hw, s.sw=W_Sigma, s.mw=mw, s.pwub=pwub, s.pwlb=pwlb, s.EndPwPoploc
	WAVE/T s.pNames, s.pGlobal
	//DOUBLE CHECK THE PARAMETER NUMBER HERE
	s.nParams = numpnts(s.pw)
	nParams = numpnts(s.pw)

//**********************************************************************
//IF ANY DISTRIBUTION INITIALIZATION IS REQUIRED FOR A MODEL PUT IT IN THIS FUNCTION
//**********************************************************************
	MultiEnScatt_InitDistribution(s)


//////////////////////////////
//Initialize Parameter Table//
//////////////////////////////
	MultiEnScatt_ParamTable(s) //display the parameters in a table
	MultiEnScatt_PlotFitResults(s) //displays the data & fit if not already done
	//MultiEnScatt_PlotDistributions(pop) //Displays the distributions used in the  fit....should update with fit
	
	
////////////////////////
//Prepare Fitting Info//
////////////////////////

//*************************
//Create Constraint for fit
//*************************

	//String allConst="K6>0;K6<1" //add all possible constraints for all parameters here
	String CustomConst = ""//"K14<K16;K12<K16" //Any custom constraints that you want to add that do not include > or <
	String allconst =  MultiEnScatt_Parameterbounds(s) //function that will create bounds based on the waves that give that information
//	print Allconst
	s.hold="" //initialize the hold string

	For( i=0; i<nParams; i+=1 )
		s.hold+=num2str(s.hw[i]) //read hold wave into the hold string for the fit
	endfor
	Make/t/o/n=0 s.cw //assemble the constraints
	For( i=0; i<nParams; i+=1 )
		If( str2num(s.hold[i])==0 ) //If the parameter is not held add the constraints
			For( j=0; j<ItemsInList(allConst); j+=1 )  //cycle through list of constraints to find those matching the parameter
				If( StringMatch(StringFromList(j,allConst),"K"+num2str(i)+" <*")||StringMatch(StringFromList(j,allConst),"K"+num2str(i)+" >*") )
					InsertPoints inf, 1, s.cw
					s.cw[inf]=StringFromList(j,allConst)
				endif
			endfor
		endif
	endfor
	//Additional parameters
//	For(i=0;i<itemsinlist(CustomConst,";");i+=1)
//		insertpoints inf, 1, s.cw
//		s.cw[inf] = stringfromlist(i,CustomConst)
//	endfor
	//
	KillWaves/Z M_iterates
	
	//Set the bounds for masking based on cursers
	//NOT TESTING WITH MORE THAN ONE ENERGY
	//Quick check for cursers...mask out everything outside of cursers
	DoWindow/F $(s.DataN+"_Fitting")
	Variable MaskedPoints = 0
	if( numtype(hcsr(A))==0 && numtype(hcsr(B))==0)
		Variable MinVal, MaxVal,numpnts_En
		Make/O/N=0 xw_fit,yw_fit,uw_Fit
		
		MinVal = min(pcsr(B),pcsr(A))
		MaxVal = max(pcsr(B),pcsr(A))
		Numpnts_En = MaxVal - MinVal
		
		for(i=1;i<numpnts(s.EndEnLoc);i+=1) 
			Duplicate/O/R=[minVal+s.EndEnLoc[i-1],MaxVal+s.EndEnLoc[i-1]] s.xw, xwTemp
			Duplicate/O/R=[minVal+s.EndEnLoc[i-1],MaxVal+s.EndEnLoc[i-1]] s.yw, ywTemp
			Duplicate/O/R=[minVal+s.EndEnLoc[i-1],MaxVal+s.EndEnLoc[i-1]] s.uw, uwTemp
			Concatenate/NP/O {xw_fit, xwTemp}, TempWv
			Duplicate/O TempWv,xw_fit
			Concatenate/NP/O {yw_fit, ywTemp}, TempWv
			Duplicate/O TempWv,yw_fit
			Concatenate/NP/O {uw_Fit, uwtemp}, TempWv
			Duplicate/O TempWv,uw_Fit

		endfor
		s.EndEnLoc = p*numpnts(xw_fit)/s.En
		//Redefine structure variables
	//	wave s.yw_Initial = s.yw , s.xw_initial = s.xw , s.uw_initial = s.uw
		wave s.yw = yw_fit , s.xw = xw_fit , s.uw = uw_fit
		Redimension/d/n=(numpnts(yw_fit)) s.fw, s.rw, s.mw //remake waves
		s.mw=1
		MaskedPoints = 1
		//Plot a new graph of the masked data ONLY
		MultiEnScatt_PlotReducedFit(s)
		
	else
		MaskedPoints = 0
		
	endif

//	
//Change #points based on masking....takes the area of mw in the range....either a 1 or a 0 so it easy
//	duplicate/o/Free s.EndEnLoc TempEnLoc
//	for(i=1;i<numpnts(s.EndEnLoc);i+=1)
//		TempEnLoc[i] = round(area(s.mw,s.Endenloc[i-1],s.Endenloc[i])) + TempEnloc[i-1]
//	endfor
//	s.endEnloc = tempEnloc

	//Do the Fit
	Variable/G V_FitError=0, V_FitOptions=8//save the iterates for debugging
	Variable/G V_FitMaxIters=NumIterations //Sets the number of Iterations
	Variable/G Var1=0,Var2=0,Var3=0,Var4=0,Var5=0
	sprintf wnote "DataName:%s;DateTime:%1.8e;", s.DataN, dateTime //optionally add model information here
	Note/K/NOCR s.fw, wnote
	Wave/Z Epwave
	Make/N=(numpnts(s.pw))/O EpWave = FitTolerance*s.pw
	If( !noFit )
		//Actual fit operation:
		SetDataFolder $FitDataFolder //Make sure everything is created in the correct folder
		FuncFit/NTHR=0/H=s.hold/M=2/Q MultiEnScatt_FitFunc, s.pw, s.yw /X=s.xw /D=s.fw /M=s.mw /R=s.rw /W=s.uw /I=1 /C=s.cw /STRC=s /E=EpWave
		WAVE s.fw, s.pw, s.xw=$("::"+xwn) //stupid IGOR thing that loses the pointer after executing FuncFit
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
			WAVE M_Covar, s.sw=W_sigma, s.fw, s.pw
			Duplicate/d/o M_Covar, s.CorrMat
			s.CorrMat=M_Covar[p][q]/sqrt(M_Covar[p][p]*M_Covar[q][q]) //generate correlation matrix
			s.CorrMat= p==q? nan : s.CorrMat
			MultiEnScatt_DisplayCorr(s) //display the correlations
			s.nchisq = sqrt(V_chisq/(V_npnts-(V_nterms-V_nheld))) //reduced chi squared value for the fit
			s.rw/=s.uw //residuals in sigmas
			s.rw*=s.mw //don't show residuals for masked points
			sprintf wnote, "Nchisq:%g;HoldStr:%s;InterNumber:%d;",s.nchisq,s.hold, s.fi.iternumber
			NOTE/NOCR s.fw, wnote

			MultiEnScatt_FitFunc(s) //recalculate the model to fillin model if user masked the data
		endif
	else
		//fill the waves with the model (minimized operations to be a fit func)
		MultiEnScatt_FitFunc(s)
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
STRUCTURE MultiEnScatt_FitStruct
	WAVE pw //parameters for the fit
	WAVE fw //model data to fit to the experimental data
	WAVE xw //independent parameter
	STRUCT WMFitInfoStruct fi //fit intormation
	WAVE sw //parameter uncertainties
	wave pwlb //lower bound of the fit parameters for the fit
	wave pwub //upper bound of the fit parameter for the fit
	wave EPWave //Epsilon values for each parameter
	WAVE/T pNames //parameter names
	WAVE/T pDim //name the parameter units/dimensions here
	Wave/T pGlobal //GLobal names of the parameters for updating purposes //SCATTERING THING
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
	Variable x
	//Add data libraries here
	Variable En //number of energies to be fit
	Variable numpops //number of pops in the desired fit
	Variable Pop_Params //Number of Params for each population...depends on the model
	Variable En_Params //Number of Params for each energy //Depends on the model as well
	
	//I THINK THESE ARE OLD AND I DON"T NEED THEM...NOT GOING TO DELETE JUST IN CASE
	Variable numDataSets //Old? Don't think I need this anymore
	Variable nParams //Same as above
	Variable TotalEnParams
	Variable TotalPopParams

	//For the Energy combining
	Wave EndEnLoc
	wave EndPwPopLoc
	
	//INFO FOR EACH POPULATION
	String PopModel_FF_1
	String PopModel_SF_1
	String PopModel_FF_2
	String PopModel_SF_2	
	String PopModel_FF_3
	String PopModel_SF_3	
	String PopModel_FF_4
	String PopModel_SF_4
	String PopModel_FF_5
	String PopModel_SF_5
	
ENDSTRUCTURE




//Given the set parameters in the Fit function this sets the globals to be called in the scattering function
Function MultiEnScatt_SetModelFits(s)
	STRUCT MultiEnScatt_FitStruct &s
	
	String/G $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(1)) = s.PopModel_FF_1
	String/G $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(1)) = s.PopModel_SF_1
	
	String/G $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(2)) = s.PopModel_FF_2
	String/G $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(2)) = s.PopModel_SF_2
	
	String/G $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(3)) = s.PopModel_FF_3
	String/G $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(3)) = s.PopModel_SF_3
	
	String/G $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(4)) = s.PopModel_FF_4
	String/G $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(4)) = s.PopModel_SF_4
		
	String/G $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(5)) = s.PopModel_FF_5
	String/G $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(5)) = s.PopModel_SF_5
	
end

Function MultiEnScatt_InitializePW(s)
	STRUCT MultiEnScatt_FitStruct &s

	//Setup the parameter wave since fuck everything
	//1/22/2019 - Fuck everything again as I add in functionality
	//12/11/2019 - Less upset now that things have been working, but adding more functionality
	//Non-Energy Parameters
		//Param1 - Bkg
	//Pop specific
		Variable Pop_params=0 //Running tally of parameters added to the wave for population control
	//En Specific
		Variable En_Params
		//Core Delta
		//Core Beta
		//Shell Delta
		//Shell Beta
		//Solvent Delta
		//Solvent Beta
		
		//Make a wave that tracks how many parameters exist per popultion
	Make/O/N=1 s.EndPwPopLoc = s.En //Start @ number of background shit
	Variable i,j //loop var

		
	//Add bkg param (number based on Energies!)
	
	insertpoints 0, s.En, s.pw //Bkg parameter
	insertpoints 0, s.En, s.pNames //Bkg
	insertpoints 0, s.En, s.pGlobal
	insertpoints 0, s.En, s.pwlb //lower bound
	insertpoints 0, s.En, s.pwub //upper bound
	for(i=0;i<s.En;i+=1)
		Variable/G $("root:Packages:MEF_Vars:Bkg_En"+num2str(i+1)) /N=var = 1E-14
		s.pGLobal = "Bkg_En"+num2str(i+1)
	endfor
	
	s.pw = 1E-14
	s.pNames = "Bkg_En"+num2str(p)

	//Setup Populations
	
	for(j=0;j<s.numpops;j+=1)	
		SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(j+1))
		SVAR Pop_SF = $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(j+1))
		//*************************
		//Setup the FormFactor
		//*************************
		if(!cmpstr(Pop_FF,"CoreShell_Mod1"))
			s.Pop_Params = 6 //Might not want ot change these for now
			InsertPoints inf, s.Pop_Params, s.pw //Pops as described above
			insertpoints inf, s.Pop_Params, s.pNames
			insertpoints inf, s.Pop_Params, s.pGlobal
			InsertPoints inf, s.Pop_Params, s.pwlb //
			InsertPoints inf, s.Pop_Params, s.pwub //
			
			//Append points and make them as global variables for the fit functions to use....All will be identical..need to manually change after fit begins
			//Core Scale
			Variable/G $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(j+1)) = 1E-16
			NVAR var =  $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+1] = var
			s.pNames[pop_params+(s.En-1)+1] = "CScale_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+1] = "Core_DistScale_pop"+num2str(j+1)

			//Core Mean Size
			Variable/G $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(j+1)) = 70
			NVAR var = $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+2] = var
			s.pNames[pop_params+(s.En-1)+2] = "CMS_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+2] = "Core_DistMeanSize_pop"+num2str(j+1)

			//Core WIdth
			Variable/G $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(j+1)) = 25
			NVAR var = $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+3] = var
			s.pNames[pop_params+(s.En-1)+3] = "CW_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+3] = "Core_DistWidth_pop"+num2str(j+1)

			//Shell Scale
			Variable/G $("root:Packages:MEF_Vars:Shell_DistScale_pop"+num2str(j+1)) = 1
			NVAR var = $("root:Packages:MEF_Vars:Shell_DistScale_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+4] = var
			s.pNames[pop_params+(s.En-1)+4] = "SScale_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+4] = "Shell_DistScale_pop"+num2str(j+1)
			//Shell Mean Size
			Variable/G $("root:Packages:MEF_Vars:Shell_DistMeanSize_pop"+num2str(j+1)) = 15
			NVAR var = $("root:Packages:MEF_Vars:Shell_DistMeanSize_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+5] = var
			s.pNames[pop_params+(s.En-1)+5] = "SMS_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+5] = "Shell_DistMeanSize_pop"+num2str(j+1)
			//Shell Width
			Variable/G $("root:Packages:MEF_Vars:Shell_DistWidth_pop"+num2str(j+1)) = 5
			NVAR var = $("root:Packages:MEF_Vars:Shell_DistWidth_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+6] = var
			s.pNames[pop_params+(s.En-1)+6] = "SW_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+6] = "Shell_DistWidth_pop"+num2str(j+1)
			
		elseif(!cmpstr(Pop_FF,"CoreSHell") || !cmpstr(Pop_FF,"CoreShell_Mod2"))	
			s.Pop_Params = 7 //Might not want ot change these for now
			InsertPoints inf, s.Pop_Params, s.pw //Pops as described above
			insertpoints inf, s.Pop_Params, s.pNames
			insertpoints inf, s.Pop_Params, s.pGlobal
			InsertPoints inf, s.Pop_Params, s.pwlb //
			InsertPoints inf, s.Pop_Params, s.pwub //
			
			//Append points and make them as global variables for the fit functions to use....All will be identical..need to manually change after fit begins
			//Core Scale
			Variable/G $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(j+1)) = 1E-16
			NVAR var =  $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+1] = var
			s.pNames[pop_params+(s.En-1)+1] = "CScale_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+1] = "Core_DistScale_pop"+num2str(j+1)

			//Core Mean Size
			Variable/G $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(j+1)) = 70
			NVAR var = $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+2] = var
			s.pNames[pop_params+(s.En-1)+2] = "CMS_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+2] = "Core_DistMeanSize_pop"+num2str(j+1)

			//Core WIdth
			Variable/G $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(j+1)) = 25
			NVAR var = $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+3] = var
			s.pNames[pop_params+(s.En-1)+3] = "CW_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+3] = "Core_DistWidth_pop"+num2str(j+1)

			//Shell Scale
			Variable/G $("root:Packages:MEF_Vars:Shell_Size_pop"+num2str(j+1)) = 1
			NVAR var = $("root:Packages:MEF_Vars:Shell_Size_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+4] = var
			s.pNames[pop_params+(s.En-1)+4] = "SSize_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+4] = "Shell_Size_pop"+num2str(j+1)
			
			
			//Shell Conc Scale
			Variable/G $("root:Packages:MEF_Vars:Phi_pop"+num2str(j+1)) = 1
			NVAR var = $("root:Packages:MEF_Vars:Phi_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+5] = var
			s.pNames[pop_params+(s.En-1)+5] = "Phi_pop"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+5] = "Phi_pop"+num2str(j+1)
			
			//Attenuation Attenuation thickness
			Variable/G $("root:Packages:MEF_Vars:Attenu_T_pop"+num2str(j+1)) = 0
			NVAR var = $("root:Packages:MEF_Vars:Attenu_T_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+5] = var
			s.pNames[pop_params+(s.En-1)+6] = "Attenu_T_pop"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+6] = "Attenu_T_pop"+num2str(j+1)
			
			//Shell Conc Scale
			Variable/G $("root:Packages:MEF_Vars:Attenu_phi_pop"+num2str(j+1)) = 1
			NVAR var = $("root:Packages:MEF_Vars:Attenu_phi_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+5] = var
			s.pNames[pop_params+(s.En-1)+7] = "Attenu_phi_pop"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+7] = "Attenu_phi_pop"+num2str(j+1)
		
		
		elseif(!cmpstr(pop_FF,"Power_Law")) // f(x) = A/x^p +Bx
			s.pop_Params = 3
			InsertPoints inf, s.Pop_Params, s.pw //Pops as described above
			insertpoints inf, s.Pop_Params, s.pNames
			insertpoints inf, s.Pop_Params, s.pGlobal
			InsertPoints inf, s.Pop_Params, s.pwlb //
			InsertPoints inf, s.Pop_Params, s.pwub //
			//Scale for power law decay
			Variable/G $("root:Packages:MEF_Vars:Porod_Scale_pop"+num2str(j+1)) /N=var = 1E-17
			s.pw[pop_params+(s.En-1)+1] = var
			s.pNames[pop_params+(s.En-1)+1] = "PLScale_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+1] = "Porod_Scale_pop"+num2str(j+1)
			//Linear Term
			Variable/G $("root:Packages:MEF_Vars:Porod_line_pop"+num2str(j+1)) /N=var = 0
			s.pw[pop_params+(s.En-1)+2] = var
			s.pNames[pop_params+(s.En-1)+2] = "PLLine_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+2] = "Porod_line_pop"+num2str(j+1)
			//Power term
			Variable/G $("root:Packages:MEF_Vars:Porod_power_pop"+num2str(j+1)) /N=var = -4
			s.pw[pop_Params+(s.En-1)+3] = var
			s.pNames[pop_params+(s.En-1)+3] = "PLPwr_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+3] = "Porod_power_pop"+num2str(j+1)
			
		elseif(!cmpstr(pop_FF,"Sphere")) // sphere with distribution
		
			s.pop_Params = 3
			InsertPoints inf, s.Pop_Params, s.pw //Pops as described above
			insertpoints inf, s.Pop_Params, s.pNames
			insertpoints inf, s.Pop_Params, s.pGlobal
			InsertPoints inf, s.Pop_Params, s.pwlb //
			InsertPoints inf, s.Pop_Params, s.pwub //
			//Scale Parameter for the sphere
			Variable/G $("root:Packages:MEF_Vars:Sphere_Scale_pop"+num2str(j+1)) /N=var = 1
			s.pw[pop_params+(s.En-1)+1] = var
			s.pNames[pop_params+(s.En-1)+1] = "Sphere_Scale"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+1] = "Sphere_Scale_pop"+num2str(j+1)
			//Average size of aggregate
			Variable/G $("root:Packages:MEF_Vars:Sphere_Rad_pop"+num2str(j+1)) /N=var = 75
			s.pw[pop_params+(s.En-1)+2] = var
			s.pNames[pop_params+(s.En-1)+2] = "Sphere_rad"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+2] = "Sphere_Rad_pop"+num2str(j+1)
			//Width of aggregate distribution
			Variable/G $("root:Packages:MEF_Vars:Sphere_Width_pop"+num2str(j+1)) /N=var = 25
			s.pw[pop_params+(s.En-1)+3] = var
			s.pNames[pop_params+(s.En-1)+3] = "Sphere_width"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+3] = "Sphere_width_pop"+num2str(j+1)
		

		endif
			Pop_Params += s.Pop_Params //End of Form Factor
		
		//*************************
		//Setup the Structure Factor
		//*************************

		if(!cmpstr(Pop_SF,""))
			s.Pop_Params = 0
		
		elseif(!cmpstr(Pop_SF,"HardSpheres"))
			s.Pop_Params = 2
			InsertPoints inf, s.Pop_Params, s.pw //Pops as described above
			insertpoints inf, s.Pop_Params, s.pNames
			insertpoints inf, s.Pop_Params, s.pGlobal
			InsertPoints inf, s.Pop_Params, s.pwlb //
			InsertPoints inf, s.Pop_Params, s.pwub //
		
			//Structure Factor Phi
			Variable/G $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(j+1)) = 0.2
			NVAR var = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+1] = var
			s.pNames[pop_params+(s.En-1)+1] = "SF_Phi_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+1] = "SF_Phi_Param_pop"+num2str(j+1)
			//Structure Factor Radius
			Variable/G $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(j+1)) = 100
			NVAR var = $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(j+1))
			s.pw[pop_params+(s.En-1)+2] = var	
			s.pNames[pop_params+(s.En-1)+2] = "SF_Rad_p"+num2str(j+1)
			s.pGlobal[pop_params+(s.En-1)+2] = "SF_Eta_Param_pop"+num2str(j+1)
		endif
	
		Pop_Params += s.Pop_Params
//	endfor
	
		Variable numpw = numpnts(s.pw) //total number of parameters so far in the thing
		s.TotalPopParams = numpw-s.En

	
	//Add in energy parameters based on Pops //NEW ADDITION 12/11/2019
	//*************************
	//Setup the Energy Dependence for each potential Form factor pop
	//*************************
	
	//Add in energy dependence for a core shell - 2 parameters per Core/SHell making a total of 4 parameters
		if(!cmpstr(Pop_FF,"CoreSHell"))
			s.En_Params = 6 //DOn't change this ever, how many parameters per energy
			En_Params = s.En_Params
			for(i=0;i<s.En;i+=1)
				Insertpoints inf, En_Params, s.pw
				insertpoints inf, En_Params, s.pNames
				insertpoints inf, En_Params, s.pGlobal
				InsertPoints inf, En_Params, s.pwlb //
				InsertPoints inf, En_Params, s.pwub //
				
				//Param8 - Core Delta
				Variable/G $("root:Packages:MEF_Vars:Delta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 10
				NVAR var = $("root:Packages:MEF_Vars:Delta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[numpw+i*En_Params] = var
				s.pNames[numpw+i*En_Params] = "Delta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGLobal[numpw+i*En_Params] = "Delta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)
			
				//Param9 - Core Beta
				Variable/G $("root:Packages:MEF_Vars:Beta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 0
				NVAR var = $("root:Packages:MEF_Vars:Beta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1))	
				s.pw[1+numpw+i*En_Params] = var
				s.pNames[1+numpw+i*En_Params] = "Beta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGLobal[1+numpw+i*En_Params] = "Beta_Core_En"+num2str(i+1)+"_pop"+num2str(j+1)
				
				//Param10 - Shell Delta
				Variable/G $("root:Packages:MEF_Vars:Delta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 15
				NVAR var = $("root:Packages:MEF_Vars:Delta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[2+numpw+i*En_Params] = var
				s.pNames[2+numpw+i*En_Params] = "Delta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGLobal[2+numpw+i*En_Params] = "Delta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)			
		
				//Param11 - Shell Beta
				Variable/G $("root:Packages:MEF_Vars:Beta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 0
				NVAR var = $("root:Packages:MEF_Vars:Beta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[3+numpw+i*En_Params] = var
				s.pNames[3+numpw+i*En_Params] = "Beta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGLobal[3+numpw+i*En_Params] = "Beta_Shell_En"+num2str(i+1)+"_pop"+num2str(j+1)	
						
		
				//Param12 - Solvent Delta
				Variable/G $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 1
				NVAR var = $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[4+numpw+i*En_Params] = var
				s.pNames[4+numpw+i*En_Params] = "Delta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGLobal[4+numpw+i*En_Params] = "Delta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)							

				//Param13 - Solvent Beta
				Variable/G $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 0
				NVAR var = $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[5+numpw+i*En_Params] = var
				s.pNames[5+numpw+i*En_Params] = "Beta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)
				s.pGlobal[5+numpw+i*En_Params] = "Beta_Solvent_En"+num2str(i+1)+"_pop"+num2str(j+1)
			
			endfor
			
		elseif(!cmpstr(Pop_FF,"Sphere"))
			s.En_Params = 1
			En_Params = s.En_Params
			
			for(i=0;i<s.En;i+=1)
				Insertpoints inf, En_Params, s.pw
				insertpoints inf, En_Params, s.pNames
				insertpoints inf, En_Params, s.pGlobal
				InsertPoints inf, En_Params, s.pwlb //
				InsertPoints inf, En_Params, s.pwub //
				
				//Param - Core Delta_effective
				Variable/G $("root:Packages:MEF_Vars:Sphere_Delta_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 1
				NVAR var = $("root:Packages:MEF_Vars:Sphere_Delta_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[numpw+i*En_Params] = var
				s.pNames[numpw+i*En_Params] = "Sphere_Delta_En"+num2str(i+1)
				s.pGLobal[numpw+i*En_Params] = "Sphere_Delta_En"+num2str(i+1)+"_pop"+num2str(j+1)
				
			
			endfor	
			
		elseif(!cmpstr(Pop_FF,"Power_Law"))
			s.En_Params = 2
			En_Params = s.En_Params
			
			for(i=0;i<s.En;i+=1)
				Insertpoints inf, En_Params, s.pw
				insertpoints inf, En_Params, s.pNames
				insertpoints inf, En_Params, s.pGlobal
				InsertPoints inf, En_Params, s.pwlb //
				InsertPoints inf, En_Params, s.pwub //
				
				//Param - Core Delta_effective
				Variable/G $("root:Packages:MEF_Vars:PorodScale_En"+num2str(i+1)+"_pop"+num2str(j+1)) = 1
				NVAR var = $("root:Packages:MEF_Vars:PorodScale_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[numpw+i*En_Params] = var
				s.pNames[numpw+i*En_Params] = "PorodScale_En"+num2str(i+1)
				s.pGLobal[numpw+i*En_Params] = "PorodScale_En"+num2str(i+1)+"_pop"+num2str(j+1)
				
				//Param - Core Delta_effective
				Variable/G $("root:Packages:MEF_Vars:PorodPower_En"+num2str(i+1)+"_pop"+num2str(j+1)) = -3
				NVAR var = $("root:Packages:MEF_Vars:PorodPower_En"+num2str(i+1)+"_pop"+num2str(j+1))
				s.pw[1+numpw+i*En_Params] = var
				s.pNames[1+numpw+i*En_Params] = "PorodPower_En"+num2str(i+1)
				s.pGLobal[1+numpw+i*En_Params] = "PorodPower_En"+num2str(i+1)+"_pop"+num2str(j+1)
			
			endfor	
		
		endif
		
		////////////////////////////////////
		//End Adding Population Parameters//
		////////////////////////////////////
		//Recalculate total number per population
		Pop_Params += En_Params*s.En
	
		//Update the end of the population
		insertpoints inf, 1, s.EndPwPopLoc
		s.EndPwPopLoc[j+1] = Pop_Params + s.En //because otherwise I have to recode a lot I don't want too
		
	endfor
	
	//////////////////////////////////////////////////////////////
	//Add Solvent parameters for ALL POPULATIONS TO COMPARE WITH//
	//////////////////////////////////////////////////////////////
//	Variable Total_ModelParams = numpnts(s.pw)
//	for(i=0;i<s.En;i+=1)
//		Insertpoints inf, 2, s.pw
//		insertpoints inf, 2, s.pNames
//		insertpoints inf, 2, s.pGlobal
//		InsertPoints inf, 2, s.pwlb //
//		InsertPoints inf, 2, s.pwub //
//		
//		//Param12 - Solvent Delta
//		Variable/G $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(i+1)) = 1
//		NVAR var = $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(i+1))
//		s.pw[Total_ModelParams + i*2] = var
//		s.pNames[Total_ModelParams+ i*2] = "Delta_Solvent_En"+num2str(i+1)
//		s.pGLobal[Total_ModelParams+ i*2] = "Delta_Solvent_En"+num2str(i+1)							
//		
//		//Param13 - Solvent Beta
//		Variable/G $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(i+1)) = 0
//		NVAR var = $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(i+1))
//		s.pw[1+Total_ModelParams+ i*2] = var
//		s.pNames[1+Total_ModelParams+ i*2] = "Beta_Solvent_En"+num2str(i+1)
//		s.pGlobal[1+Total_ModelParams+ i*2] = "Beta_Solvent_En"+num2str(i+1)
//
//	endfor


	////Make boundary waves
	s.pwlb = s.pw * 0.8
	s.pwub = s.pw * 1.2

	return numpnts(s.pw)
end

//Initialize any distributions you may require based on pops

Function MultiEnScatt_InitDistribution(s)
	STRUCT MultiEnScatt_FItStruct &s
	
	Variable i //Loop variable

	String DistType = "Schulz-Zimm" //I think this is the only distribution we want
	for(i=0;i<s.numpops;i+=1)
		SVAR Pop_Dist = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(i+1))
		
		if(!cmpstr(Pop_Dist,"CoreShell_Mod1"))
		
			Variable/G $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1)) = 50
			Variable/G $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1)) = 10
			Variable/G $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1)) = 0.01
			Variable/G $("root:Packages:MEF_Vars:NumDistributions_pop"+num2str(i+1)) = 2 //NUmber of distributions required
	

			//Core DIstribution #points
			NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1))
			//Shell Distribution #points
			NVAR ShellPnts = $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1))
			//Distribution precision
			NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1))
			//Create the initial distributions

//	s.EndPwPopLoc
			MEF_InitializeDistribution(i+1,"Core",Corepnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+1],s.pw[s.EndPwPopLoc[i]+2],s.pw[s.EndPwPopLoc[i]])
			MEF_InitializeDistribution(i+1,"Shell",Shellpnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+4],s.pw[s.EndPwPopLoc[i]+5],s.pw[s.EndPwPopLoc[i]+3])
			MEF_CoreShellCompositeDistribution(i+1,"Core","Shell")
			MultiEnScatt_PlotDistributions(i+1) //Displays the distributions used in the  fit....should update with fit

		
		elseif(!cmpstr(Pop_Dist,"CoreShell") || !cmpstr(Pop_Dist,"CoreShell_Mod2"))
			Variable/G $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1)) = 50
			Variable/G $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1)) = 10
			Variable/G $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1)) = 0.01	
			Variable/G $("root:Packages:MEF_Vars:NumDistributions_pop"+num2str(i+1)) = 1 //NUmber of distributions required


			//Core DIstribution #points
			NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1))
			//Shell Distribution #points
			//Distribution precision
			NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1))
			//Create the initial distributions
			MEF_InitializeDistribution(i+1,"Core",Corepnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+1],s.pw[s.EndPwPopLoc[i]+2],s.pw[s.EndPwPopLoc[i]])
			MultiEnScatt_PlotDistributions(i+1) //Displays the distributions used in the  fit....should update with fit
		
		elseif(!cmpstr(Pop_Dist, "Sphere"))
			Variable/G $("root:Packages:MEF_Vars:Spherepnts_pop"+num2str(i+1)) = 50
			Variable/G $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1)) = 0.01	
			Variable/G $("root:Packages:MEF_Vars:NumDistributions_pop"+num2str(i+1)) = 1 //NUmber of distributions required
		
			
			//Sphere DIstribution #points
			NVAR SpherePnts = $("root:Packages:MEF_Vars:Spherepnts_pop"+num2str(i+1))
			//Shell Distribution #points
			//Distribution precision
			NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1))
			//Create the initial distributions
			MEF_InitializeDistribution(i+1,"Sphere",SpherePnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+1],s.pw[s.EndPwPopLoc[i]+2],s.pw[s.EndPwPopLoc[i]])
			MultiEnScatt_PlotDistributions(i+1)
		endif
	endfor


end


//Actual Fit Function
FUNCTION MultiEnScatt_FitFunc(s) : FitFunc
	STRUCT MultiEnScatt_FitStruct &s
	//Section on naming the parameters
	String CUrrentFOlder=GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
//	SaveGlobal(s) //Quick save 
	//Loop Variables
	Variable i,j
//****************************************

//	s.pw[8] = s.pw[7]/4


//*****************************************
	//Update global variables based on s.pw....Written by Jan! Great idea, Line 1521 IR2L_NLSQFFunction.ipf
	Variable NumParams = numpnts(s.pw)
	String ParamName
	
	For(i=0;i<NumParams;i+=1)
		ParamName = s.pGlobal[i]
		NVAR TempVar = $(ParamName)
		TempVar=s.pw[i]
	endfor
//****************************************

	//Reinitialize distributions assuming parameters have been adjusted during fit
	MultiEnScatt_InitDistribution(s)

	//Create Pop models for each set //Gonna be a clusterfuck
	//duplicate the x-wave to create Qmodel_Final for later solving
	String ModelFinalQ = "Qmodel_Final"
	Duplicate/o s.xw $ModelFinalQ
	//Now we cycle through the fitting
	for(j=1;j<s.En+1;j+=1)
	//Create appropriate xw
		string ModelQs = "Qmodel_Set"+num2str(j) //Model Q is literally the Qvalues...x wave for displaying
		Duplicate/O/R=[s.EndEnLoc[j-1],s.EndEnLoc[j]-1] s.xw $ModelQs //Set up the Qwave to have the appropriate number of points as split from the datawave
		wave debuggpointer = $modelQs
		for(i=1;i<s.numpops+1;i+=1)
			//Check what kind of fit we want to run
			SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(i))
			SVAR Pop_SF = $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(i))
			
			//Check what the population is, either run the scattering function or another function with simpler functions
			if(!cmpstr(Pop_FF,"CoreShell_Mod1") || !cmpstr(Pop_FF,"CoreShell")|| !cmpstr(Pop_FF,"CoreShell_mod2"))
				MEF_CalcScatteringInt(i,j) //J is the dataset thing
			elseif(!cmpstr(Pop_FF,"Power_Law"))
				MEF_CalcInt(i,j)
			endif
			//Wave result = $("root:Packages:MEF_VArs:IntensityModel_Set"+num2str(j)+"_pop"+num2str(1))
			DoUpdate /W=$("DistributionGraph_pop"+num2str(i))
		endfor
	endfor
	
	MultiEnScatt_SumModel(s)
	Wave result = $("root:Packages:MEF_VArs:IntensitySumModel")

	s.fw = result	
	DoUpdate /W=$(s.DataN+"_Fitting")
	DoUpdate /W=$(s.DataN+"_Fitting_MaskedData")
	DoUpdate /W=$(s.DataN+"_Parameters")
	SetDataFolder $CurrentFolder
END

Function MultiEnScatt_SumModel(s)
	
	STRUCT MultiEnScatt_FitStruct &s
	String CurrentFolder = GetDataFOlder(1)
	SetDataFolder root:Packages:MEF_Vars
	//Setup the final wave BEFORE we loop so it won't get written over and left with a lot of zeros
	Wave/Z FinalModelQ = $("Qmodel_Final")
	duplicate/o FinalMOdelQ $("IntensitySumModel")
	Wave ModelCombinedSum = $("IntensitySumModel")

//
////Grab some waves because fuck me
////
	Wave/Z F127Mu = root:SimulateAbs_122019:F127_Mu
	Wave/Z H2OMu = root:SimulateAbs_122019:H2O_Mu
	Wave/Z EnMu = root:SimulateAbs_122019:OC_En
	Wave/Z AttenuationWave = root:SimulateAbs_122019:Sample_Abs_10
	Wave/Z SampleEn = root:Data_FullFit_Res:EN_Wave
////


	Make/N=0/O Tempwave
	Variable i,j,k //Datasets and populations //Although I only have 1 dataset
	for(j=1;j<s.En+1;j+=1)
//		Wave/Z FinalModelQ = $("Qmodel_Final")
		wave/Z CurrentPopQ = $("Qmodel_set"+num2str(j))
		if(!WaveExists(CurrentPopQ))
			print "Some issue happened in Summing the model, Can't find the Qwave...Line 527 of Fit Procedure"
			return 1
		endif
		duplicate/o CurrentPopQ $("IntensitySumModel_set"+num2str(j))
//		duplicate/o FinalMOdelQ $("IntensitySumModel")
		Wave ModelIntSum = $("IntensitySumModel_set"+num2str(j))
//		Wave ModelCombinedSum = $("IntensitySumModel")
		ModelIntSum = 0		
		for(i=1;i<s.numpops+1;i+=1)
			wave CurrentPopModel = $("IntensityModel_set"+num2str(j)+"_pop"+num2str(i))
			ModelIntSum += CurrentPopModel
		endfor
		
		//Add in Resonant scaling based on energy
		
		Variable hc = 1239.9 //nm*ev
		
		
		ModelIntSum *= 4*PI^2*(SampleEn[j-1])^4
		ModelIntSum /= (hc)^4
		
		
//		//Add in attenuation based on energy
////		
		NVAR Attenu_Phi = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(1)) //$("root:Packages:MEF_Vars:Attenu_phi_pop"+num2str(1))
		NVAR Attenu_t = $("root:Packages:MEF_Vars:Attenu_t_pop"+num2str(1))
		
		AttenuationWave = e^(-((1-Attenu_Phi)*H2Omu + Attenu_Phi*F127mu)*Attenu_t)
		
		Variable Transmission = interp(SampleEn[j-1],EnMu,AttenuationWave)
		
		ModelIntSum *= Transmission
	//	ModelIntSum += s.pw[j-1]


//		concatenate/np/O {ModelCombinedSum, ModelIntSum}, ModelCombinedSum
		
	//	ModelCombinedSum[s.EndEnLoc[j-1],s.EndEnLoc[j]-1] = ModelIntSum[p]
		for(k=0;k<s.EndEnloc[1];k+=1)
			ModelCombinedSum[s.EndEnLoc[j-1]+k] = ModelIntSum[k]+s.pw[j-1] //Adds the backgrounds based on energy now.
		endfor	
	endfor
	
	
	
	SetDataFolder $currentfolder
end



//Displays the parameters in a table
FUNCTION MultiEnScatt_ParamTable(s)
	STRUCT MultiEnScatt_FitStruct &s
	DoWindow/F $(s.DataN+"_Parameters")
	IF( V_flag )
		return 0
	endif
	Edit/N=$(s.DataN+"_Parameters") /W=(375.75,445.25,592.5,611.75) s.pNames,s.pw,s.pwlb,s.pwub,s.sw, s.hw, as s.DataN+" Parameters"
	ModifyTable size=8,format(Point)=1,width(Point)=18,width(s.pNames)=50,width(s.pw)=50,width(s.pwlb)=50,width(s.pwub)=50,width(s.sw)=50, width(s.hw)=20
END

Function/S MultiEnScatt_Parameterbounds(s)
	STRUCT MultiEnScatt_FitStruct &s
	
	Variable numparams = numpnts(s.pw) //number of total constraints
	String allconst = "" //the final constraint wave to pass to the fit function

	if(!waveexists(s.pwub) || !waveexists(s.pwlb)) //double check to make sure the bound waves exist
		return allconst
	endif
	
	Variable i //loop variable
	for(i=0;i<numparams;i+=1) //loop through parameters and assign to the constraint wave
		//lower bound
		if(numtype(s.pwlb[i])==0)
			allconst += "K"+num2str(i)+" > "+num2str(s.pwlb[i]) + ";"
		endif
		//upper bound
		if(numtype(s.pwub[i])==0)
			allconst += "K"+num2str(i)+" < "+num2str(s.pwub[i]) + ";"
		endif
	endfor
	return allconst
end





//Displays the fit in a graph
FUNCTION MultiEnScatt_PlotFitResults(s)
	STRUCT MultiEnScatt_FitStruct &s
	Variable update
	String ywN=NameOfWave(s.yw), fwN=NameOfWave(s.fw), lgndTxt
	DoWindow/F $(s.DataN+"_Fitting")
	If( V_Flag && update == 0)
		return 0
	endif
	//Add the initial non-cut data to pick from
	Display/N=$(s.DataN+"_Fitting") /W=(475.5,106.25,733.5,415.25) s.yw vs s.xw as s.DataN+" Fitting"
	AppendToGraph s.fw vs s.xw
	ModifyGraph rgb($fwN)=(0,0,0), rgb($ywN)=(65280,0,0)
	ModifyGraph mode($fwN)=3,marker($fwN)=8,msize($fwN)=3
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



//Display masked fit results
Function MultiEnScatt_PlotReducedFit(s)
	STRUCT MultiEnScatt_FitStruct &s
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

Function MultiEnScatt_PlotDistributions(pop)

	Variable pop
	//FLDR
	String CurrentFOlder = GetdataFOlder(1)
	SetDataFolder root:packages:MEF_Vars
	//Initialization
	String lgndTxt
	
	//Gather distribution waves to fit
	//Identify the type of thing we want to plot
	SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(pop))
	//Local names to plot
	String Dist1 = ""
	String Dist2 = ""
	String Rad1 = ""
	String Rad2 = ""
	
	//Grab the name of the Dist and Radius
	if(!cmpstr(Pop_FF,"Sphere"))
		Dist1 = "Sphere_VolumeDist_Pop"+num2str(pop)
		Dist2 = ""
		Rad1 = "Sphere_Radius_Pop"+num2str(pop)
		Rad2 = ""
		
	elseif(!cmpstr(Pop_FF,"CoreShell"))
		Dist1 = "Core_VolumeDist_pop"+num2str(pop)
		Dist2 = "Shell_VolumeDist_pop"+num2str(pop)
		Rad1 = "Core_Radius_pop"+num2str(pop)
		Rad2 = "Shell_Radius_pop"+num2str(pop)
	//Others
	endif
	
	
	
//	String CoreDistS = "Core_VolumeDist_pop"+num2str(pop)
//	String ShellDistS = "Shell_VolumeDist_pop"+num2str(pop)
//	String CoreRadS = "Core_Radius_pop"+num2str(pop)
//	String ShellRads = "Shell_Radius_pop"+num2str(pop)
	
	//Gater waves
	Wave/Z Dist1_w = $Dist1
	Wave/Z Dist2_w = $Dist2
	Wave/Z Rad1_w = $Rad1
	Wave/Z Rad2_w = $Rad2
	//Wave CompositeDist = $("CompositeDist_pop"+num2str(pop))
	//Check for windows
	String DistName = "DistributionGraph_pop"+num2str(pop)
	
	Dowindow/F $DistName
	If(!V_Flag)
		Display/N=$(DistName) as "Distribution Graphs"
		
		Appendtograph/L/B Dist1_w vs Rad1_w
		//Dist1 Stuff
		ModifyGraph rgb($NameofWave(Dist1_w))=(65535,0,0)
		Label Left "Probability Density [% \U]";Label bottom "Dist 1 Radius [Å]"
		ModifyGraph axRGB(left)=(65535,0,0),tlblRGB(left)=(65535,0,0),alblRGB(left)=(65535,0,0)
		ModifyGraph axRGB(bottom)=(65535,0,0),tlblRGB(bottom)=(65535,0,0),alblRGB(bottom)=(65535,0,0)


		//Maybe Dist2 Stuff
		if(waveexists(Dist2_w))
			Appendtograph/R/T Dist2_w vs Rad2_w
			ModifyGraph rgb($NameofWave(Dist2_w))=(1,16019,65535)
			Label Right "Probability Density [% \U]";Label Top "Dist 2 Radius [Å]"
			ModifyGraph axRGB(right)=(1,16019,65535),tlblRGB(right)=(1,16019,65535),alblRGB(right)=(1,16019,65535)
		ModifyGraph axRGB(top)=(1,16019,65535),tlblRGB(top)=(1,16019,65535),alblRGB(top)=(1,16019,65535)
		endif



	endif
	
	ModifyGraph tick=2
	ModifyGraph Mode=4,marker=41,msize=3
	Legend/C/N=lgnd/J/A=MC/X=-4.00 /Y=-14.00
	SetDataFolder $CurrentFolder

end



//Displays the correlation matrix
FUNCTION MultiEnScatt_DisplayCorr(s)
	STRUCT MultiEnScatt_FitStruct &s
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


//////////////////////Functions from IRENA adapted for multiEn fitting

Function SaveGlobal(s)
	STRUCT MultiEnScatt_FitStruct &s

	///////////All of these need to be created elsewhere to start.
		
	
	//Setup the parameter wave since fuck everything
	//Non-Energy Parameters
		//Param1 - Bkg
	//Pop specific
		Variable Pop_params = s.Pop_Params
		//CoreScale pop X
		//CoreMeanSize pop X
		//CoreWidth pop X
		//ShellScale pop X
		//SHellMeanSize pop X
		//ShellWidth pop X
		//SF_Radius pop X
		//SF_phi pop X
	//En Specific
		Variable En_Params = 6// s.En_Params
		//Core Delta
		//Core Beta
		//Shell Delta
		//Shell Beta
		//Solvent Delta
		//Solvent Beta


	Variable i,j //loop var
	Variable NumParam_FF
	//Setup non-energy things
	for(j=0;j<s.numpops;j+=1)
		SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(j+1))
		SVAR Pop_SF = $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(j+1))

		
		//Append points and make them as global variables for the fit functions to use....All will be identical..need to manually change after fit begins
		//FORMFACTOR FIRST
		if(!cmpstr(Pop_FF, "CoreShell_Mod1"))
			//Core Scale
			NVAR var =  $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(j+1))
			var = s.pw[s.EndPwPopLoc[j]]
			//Core Mean Size
			NVAR var = $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(j+1))
			var = s.pw[1+s.EndPwPopLoc[j]]
			//Core WIdth
			NVAR var = $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(j+1))
			var = s.pw[2+s.EndPwPopLoc[j]]
			//Shell Scale
				NVAR var = $("root:Packages:MEF_Vars:Shell_DistScale_pop"+num2str(j+1))
			var = s.pw[3+s.EndPwPopLoc[j]]
			//Shell Mean Size
			NVAR var = $("root:Packages:MEF_Vars:Shell_DistMeanSize_pop"+num2str(j+1))
			var = s.pw[4+s.EndPwPopLoc[j]]
			//Shell Width
			NVAR var = $("root:Packages:MEF_Vars:Shell_DistWidth_pop"+num2str(j+1))
			var = s.pw[5+s.EndPwPopLoc[j]]
			NumParam_FF = 6
	
		elseif(!cmpstr(Pop_FF,"Power_Law"))
			NVAR var = $("root:Packages:MEF_Vars:Porod_Scale_pop"+num2str(j+1))
			var = s.pw[s.EndPwPopLoc[j]]
			
			NVAR var = $("root:Packages:MEF_Vars:Porod_line_pop"+num2str(j+1))
			var = s.pw[s.EndPwPopLoc[j]+1]
			
			Variable/G $("root:Packages:MEF_Vars:Porod_power_pop"+num2str(j+1)) /N=var = -4
			var = s.pw[s.EndPwPopLoc[j]+2]
			NumParam_FF = 3
			
		endif
		
		
		if(!cmpstr(pop_SF,"HardSpheres"))
			//Structure Factor Phi
			NVAR var = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(j+1))
			var = s.pw[NumParam_FF+s.EndPwPopLoc[j]]
			//Structure Factor Radius
			NVAR var = $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(j+1))
			var = s.pw[1+NumParam_FF+s.EndPwPopLoc[j]]
		endif
	endfor
	
	Variable TotalPopParams = (s.EndPwPopLoc[inf]-2)+1 //Legacy, but it still works somehow
	
	//Add in energy parameters
	for(i=0;i<s.En;i+=1)
		
		//Param8 - Core Delta
		NVAR var = $("root:Packages:MEF_Vars:Delta_Core_En"+num2str(i+1))
		var = s.pw[TotalPopParams+i*En_Params+(s.En-1)]
		//Param9 - Core Beta
		NVAR var = $("root:Packages:MEF_Vars:Beta_Core_En"+num2str(i+1))
		var = s.pw[1+TotalPopParams+i*En_Params+(s.En-1)]

		//Param10 - Shell Delta
		NVAR var = $("root:Packages:MEF_Vars:Delta_Shell_En"+num2str(i+1))
		var = s.pw[2+TotalPopParams+i*En_Params+(s.En-1)]

		//Param11 - Shell Beta
		NVAR var = $("root:Packages:MEF_Vars:Beta_Shell_En"+num2str(i+1))
		var = s.pw[3+TotalPopParams+i*En_Params+(s.En-1)]

		//Param12 - Solvent Delta
		NVAR var = $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(i+1))
		var = s.pw[4+TotalPopParams+i*En_Params+(s.En-1)]

		//Param13 - Solvent Beta
		NVAR var = $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(i+1))
		var = s.pw[5+TotalPopParams+i*En_Params+(s.En-1)]
		
		
	endfor

	
end


//OLD WAY TO ASSIGN PARAMETER WAVES
//	If( !WaveExists(pw) || init ) //initialize fit waves
//		Make/o/n=(nParams) pw, hw, W_Sigma, pwlb,pwub
//		Make/o/t/n=(nParams) pNames		
//		WAVE/T s.pNames
//		Duplicate/d/o s.yw, s.mw
//		WAVE s.pw, s.hw, s.mw, s.pwlb, s.pwub
//		MultiEnScatt_InitializePW(s)
//
//		W_Sigma=0
//		s.pw={0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1E-14}  //add initial parameters values here.
//		s.pwlb=NaN
//		s.pwub=NaN
//		s.pNames={"CScale","CMS","CW","SScale","SMS","SW","VF","Rad","CD","CB","SD","SB","SolD","SolB","Scale","Bkg"} //name the parameters here.
//	//	s.pDim={""} //name the parameter units here.
//		s.mw=1 //initially fitt all points
//		hw=0 //initialize hold wave with open parameters
//		
//			
//		////////////////////////Setup the parameter wave intitial conditions here
//		//Parameter wave information for Core-Shell
//		//Form Factor Info*******************************
//		//Param0 - Core scale
//		NVAR var = $("root:Packages:MEF_Vars:Core_DistScale_pop"+num2str(pop))
//		s.pw[0] = var
//		//param1 - Core Mean Size
//		NVAR var = $("root:Packages:MEF_Vars:Core_DistMeanSize_pop"+num2str(pop))
//		s.pw[1] = var
//		//Param2 - Core width
//		NVAR var = $("root:Packages:MEF_Vars:Core_DistWidth_pop"+num2str(pop))
//		s.pw[2] = var
//		//Param3 - Shell Scale
//		NVAR var = $("root:Packages:MEF_Vars:Shell_DistScale_pop"+num2str(pop))
//		s.pw[3] = var
//		//Param4 - Shell Mean Size
//		NVAR var = $("root:Packages:MEF_Vars:Shell_DistMeanSize_pop"+num2str(pop))
//		s.pw[4] = var
//		//Param5 - Shell Width
//		NVAR var = $("root:Packages:MEF_Vars:Shell_DistWidth_pop"+num2str(pop))
//		s.pw[5] = var
//		//Structure factor info**************************
//		//Param6 - Volume Fraction
//		NVAR var = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(pop))
//		s.pw[6] = var
//		//Param7 - Radius thing
//		NVAR var = $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(pop))
//		s.pw[7] = var
//		//Index of refractions***************************
//		//Param8 - Core Delta
//		NVAR var = $("root:Packages:MEF_Vars:Delta_Core_pop"+num2str(pop))
//		s.pw[8] = var
//		//Param9 - Core Beta
//		NVAR var = $("root:Packages:MEF_Vars:Beta_Core_pop"+num2str(pop))
//		s.pw[9] = var
//		//Param10 - Shell Delta
//		NVAR var = $("root:Packages:MEF_Vars:Delta_Shell_pop"+num2str(pop))
//		s.pw[10] = var
//		//Param11 - Shell Beta
//		NVAR var = $("root:Packages:MEF_Vars:Beta_Shell_pop"+num2str(pop))
//		s.pw[11] = var
//		//Param12 - Solvent Delta
//		NVAR var = $("root:Packages:MEF_Vars:Delta_Solvent_pop"+num2str(pop))
//		s.pw[12] = var
//		//Param13 - Solvent Beta
//		NVAR var = $("root:Packages:MEF_Vars:Beta_Solvent_pop"+num2str(pop))
//		s.pw[13] = var
//		///////////Al of these need to be created elsewhere to start.
//		//set initial boudns
//		s.pwub = s.pw*10
//		s.pwlb = s.pw/10
//		//Non Fit parameter things we will need later
//	endif
//	WAVE s.pw, s.hw, s.sw=W_Sigma, s.mw, s.pwub, s.pwlb
//	WAVE/T s.pNames
//	//Could prepare model in a subroutine here
//
//	//Distribution type to fit
//	String DistType = "Schulz-Zimm"
//	//Core DIstribution #points
//	NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(pop))
//	//Shell Distribution #points
//	NVAR ShellPnts = $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(pop))
//	//Distribution precision
//	NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(pop))
////Create the initial distributions


//Function RunScatteringFit(yw,uw,xw [noFit,init,pops,Energies])
//
//	Wave yw,uw,xw
//	Variable noFit
//	Variable Init
//	Variable pops //number of populatons we want to fit
//	Variable Energies //number of energies that we are fitting
//	//We will now basically setup all initial conditions that we would want in order to fit
//	
//	NewDataFolder/o/s root:Packages:MEF_Vars
//	
//	Variable i,j //loop variables for pops and energies
//	For(i=0 ;i<pops ;i+=1)
//		For(j=0;j<Energies;j+=1)
//		
//		
//		//Set the index values
//		Variable CoreDelta = 10
//		Variable CoreBeta = 0
//	
//		Variable SHellDelta = 20
//		Variable ShellBeta = 0
//			
//		Variable SolventDelta = 1
//		Variable SolventBeta = 0 
//	
//	//Set distribution values
//		String DistType = "Schulz-Zimm"
//		Variable/G DistPrecision_pop1 = 0.01 //Same for both
//		//Core
//		Variable/G Corepnts_pop1 = 50
/////////////////////////////////
//		Variable CoreMean = 150
//		Variable CoreWidth = 50
//		Variable CoreScale = 1
/////////////////////////////////
//
//		//Shell
//		Variable/G Shellpnts_pop1 = 10
/////////////////////////////////
//		Variable ShellMean = 25
//		Variable ShellWidth = 1
//		Variable ShellScale = 1
/////////////////////////////////
//	//Structure Factor
//		Variable/G SF_Phi_Param_pop1 = 0
//		Variable/G SF_Eta_Param_pop1 = 1
//	
//		String/G StructureFactor_pop1 = "HardSpheres"
//	
//		wave nothing = $("root:Packages:MEF_VArs:Qmodel_set1")
//		endfor
//	endfor
//	//Do everything
//	SetIndex("Core",CoreDelta,CoreBeta,1)
//	SetIndex("Shell",ShellDelta,ShellBeta,1)
//	SetIndex("SOlvent",SOlventDelta,SOlventBeta,1)
//	
//end

//Code for distributions.....not used anymore
	//Create the distributions for the fitting
	//Distribution type to fit
//	String DistType = "Schulz-Zimm"
//	for(i=0;i<s.numpops;i+=1)
//		SVAR Pop_Dist = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(i+1))
//		
//		if(!cmpstr(Pop_Dist,"CoreShell_Mod1"))
//		
//			Variable/G $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1)) = 50
//			Variable/G $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1)) = 10
//			Variable/G $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1)) = 0.01	
//
//			//Core DIstribution #points
//			NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1))
//			//Shell Distribution #points
//			NVAR ShellPnts = $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1))
//			//Distribution precision
//			NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1))
//			//Create the initial distributions
//
////	s.EndPwPopLoc
//			MEF_InitializeDistribution(i+1,"Core",Corepnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+1],s.pw[s.EndPwPopLoc[i]+2],s.pw[s.EndPwPopLoc[i]])
//			MEF_InitializeDistribution(i+1,"Shell",Shellpnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+4],s.pw[s.EndPwPopLoc[i]+5],s.pw[s.EndPwPopLoc[i]+3])
//			MEF_CoreShellCompositeDistribution(i+1,"Core","Shell")
//			MultiEnScatt_PlotDistributions(i+1) //Displays the distributions used in the  fit....should update with fit
//		endif
//	endfor
//	//Core DIstribution #points
//	NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(pop))
//	//Shell Distribution #points
//	NVAR ShellPnts = $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(pop))
//	//Distribution precision
//	NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(pop))
////Create the initial distributions
//
//	MEF_InitializeDistribution(1,"Core",Corepnts,DistPrecision,DistType,s.pw[1],s.pw[2],s.pw[0])
//	MEF_InitializeDistribution(1,"Shell",Shellpnts,DistPrecision,DistType,s.pw[4],s.pw[5],s.pw[3])
//	MEF_CoreShellCompositeDistribution(1,"Core","Shell")

//FROM FIT FUNCTION
//	String DistType = "Schulz-Zimm"
//	Variable i, j
//	For(i=0;i<s.numpops;i+=1)
//		//Core DIstribution #points
//		NVAR CorePnts = $("root:Packages:MEF_Vars:Corepnts_pop"+num2str(i+1))
//		//Shell Distribution #points
//		NVAR ShellPnts = $("root:Packages:MEF_Vars:Shellpnts_pop"+num2str(i+1))
//		//Distribution precision
//		NVAR DistPrecision = $("root:Packages:MEF_Vars:DistPrecision"+num2str(i+1))
//		//Create the initial distributions
//
//		MEF_InitializeDistribution(i+1,"Core",Corepnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+1],s.pw[s.EndPwPopLoc[i]+2],s.pw[s.EndPwPopLoc[i]])
//		MEF_InitializeDistribution(i+1,"Shell",Shellpnts,DistPrecision,DistType,s.pw[s.EndPwPopLoc[i]+4],s.pw[s.EndPwPopLoc[i]+5],s.pw[s.EndPwPopLoc[i]+3])
//		MEF_CoreShellCompositeDistribution(i+1,"Core","Shell")	
//	endfor
//	
//	//Sort energy stuff
////	for(i=0;i<s.En;i+=1)
////		SetIndex("Core",s.pw[(s.TotalPopParams+1)+i*s.En_Params],s.pw[1+(s.TotalPopParams+1)+i*s.En_Params],i)
////		SetIndex("Shell",s.pw[2+(s.TotalPopParams+1)+i*s.En_Params],s.pw[3+(s.TotalPopParams+1)+i*s.En_Params],i)
//		SetIndex("SOlvent",s.pw[4+(s.TotalPopParams+1)+i*s.En_Params],s.pw[5+(s.TotalPopParams+1)+i*s.En_Params],i)
//	endfor
	
//	//print num2str(s.pw[1])
//	//May need to remake optical constants if we are fitting these eventually.
////	print "Setting Index" + num2str(var1)
//	SetIndex("Core",s.pw[8],s.pw[9],1)
//	SetIndex("Shell",s.pw[10],s.pw[11],1)
//	SetIndex("SOlvent",s.pw[12],s.pw[13],1)
////	print "Initialize Distributions" + num2str(var2)
//	MEF_InitializeDistribution(1,"Core",Corepnts,DistPrecision,DistType,s.pw[1],s.pw[2],s.pw[0])
//	MEF_InitializeDistribution(1,"Shell",Shellpnts,DistPrecision,DistType,s.pw[4],s.pw[5],s.pw[3])
//	MEF_CoreShellCompositeDistribution(1,"Core","Shell")
////	print "Calc Scattering" + num2str(var3)