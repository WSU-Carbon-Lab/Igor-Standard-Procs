#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//Calculate the scattering from the coreshell thing

Function MEF_CalcScatteringInt(pop,dataset)
	Variable pop,dataset //DATASET ACTUALLY MEANS ENERGY, IT SELECTS WHAT OPTICAL CONSTANTS TO USE!!!!!!!
	
	//Folder stuff
	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
	
	//I am currently assuming we only have a single dataset and multiple pops should be able to run...
	
	//Start going through Gmatrix stuff initialization.
	string ModelQs = "Qmodel_Set"+num2str(Dataset) //Model Q is literally the Qvalues...x wave for displaying
	//Duplicate/o xw $ModelQs
	Wave ModelQ = $ModelQS
	if (!waveexists(modelQ))
		Abort "Select original data first"
	endif
	//Initialize the modelint? unsure of where ModelQ comes from
	Duplicate/o ModelQ, $("IntensityModel_Set"+num2str(Dataset)+"_pop"+Num2str(pop)) //Makes intensity wave of the same size
	Wave ModelInt = $("IntensityModel_Set"+num2str(Dataset)+"_pop"+Num2str(pop))
	ModelInt = 0//initializes the wave so it equals zero and not the Q value....duh
	
	//Grab our fit parameters that include the optical constants of the components
	NVAR Delta_Core1 = $("root:Packages:MEF_Vars:Delta_Core_En"+num2str(dataset)+"_pop"+num2str(pop))
	NVAR Delta_Shell1 = $("root:Packages:MEF_Vars:Delta_Shell_En"+num2str(dataset)+"_pop"+num2str(pop))
	NVAR Delta_Solvent1 = $("root:Packages:MEF_Vars:Delta_Solvent_En"+num2str(dataset)+"_pop"+num2str(pop))
	
	NVAR Beta_Core1 = $("root:Packages:MEF_Vars:Beta_Core_En"+num2str(dataset)+"_pop"+num2str(pop))
	NVAR Beta_Shell1 = $("root:Packages:MEF_Vars:Beta_Shell_En"+num2str(dataset)+"_pop"+num2str(pop))
	NVAR Beta_Solvent1 = $("root:Packages:MEF_Vars:Beta_Solvent_En"+num2str(dataset)+"_pop"+num2str(pop))
	
	//NEW THING
	NVAR Phi = $("root:Packages:MEF_Vars:Phi_pop"+num2str(pop))
	
//	//ADDITIONAL FIXES THAT FUCK UP THE MODEL BUT MIGHT WORK FOR THE FUTURE
	Variable Delta_Core = Delta_Core1
	Variable Delta_Shell = Delta_Shell1*Phi + (1-Phi)*Delta_Solvent1
	Variable Delta_Solvent = Delta_Solvent1
//	
	Variable Beta_Core = Beta_Core1
	Variable Beta_Shell = Beta_Shell1*Phi + (1-Phi)*Beta_Solvent1
	Variable Beta_Solvent = Beta_Solvent1
	
//	Wave EnergyList = root:Data_FullEzFit_Result:En_Wave
	
//	Variable Delta_Shell = CalcRatio(phi,EnergyList[dataset-1]) * Delta_Core
	
	//Update GLobal Here,
//	
//	NVAR UpdateGlobal = $("root:Packages:MEF_Vars:Delta_Shell_En"+num2str(dataset)+"_pop"+num2str(pop))
//	UpdateGlobal = Delta_Shell
	
	//For non-CoreShell things
	Variable/C ComplexIndex_Object = cmplx(1-delta_core,beta_core)
	Variable/C ComplexIndex_Matrix = cmplx(1-delta_Solvent,Beta_solvent)
	Variable LocalContrast = magsqr(ComplexIndex_Object - ComplexINdex_Matrix)
	/////////////////////**Now we break up into different possible gmatrix calculations based on number of distributions
	SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(pop)) //Grab the FF
	SVAR Pop_SF = $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(pop)) //Grab the SF
	NVAR Pop_Dist = $("root:Packages:MEF_Vars:NumDistributions_pop"+num2str(pop)) //NUmber of distributions required
	
	//Variahbles for defining number of points
	Variable M,N,O //numpnts in Q, numnpts in Core, nuimnpts in Shell (for Gmatrix dimension calculation)

	if(Pop_Dist == 1)
		if(!cmpstr(Pop_FF,"CoreShell") || !cmpstr(Pop_FF,"CoreShell_mod2"))
			Wave Core_Radius = $("root:Packages:MEF_Vars:Core_Radius_pop"+num2str(pop))
			Wave NumDist_Sum = $("root:Packages:MEF_Vars:NumberDistSum_pop"+num2str(pop))
			Wave VolumeDist = $("root:Packages:MEF_Vars:core_VolumeDist_pop"+num2str(pop))
			NVAR ShellThickness = $("root:Packages:MEF_Vars:Shell_Size_pop"+num2str(pop))
				
			Wave/Z Gmatrix=$("root:Packages:MEF_VARS:Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
			M=numpnts(ModelQ)
			N=numpnts(Core_Radius)
		//	if(!WaveExists(Gmatrix))
				Make/D/O/N=(M,N) $("Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
				Wave Gmatrix=$("root:Packages:MEF_Vars:Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
		//	endif	
			//calculate G matrix...
			MEF_CalculateGMatrix_2D(Gmatrix,ModelQ,Core_Radius,ShellThickness,pop,0,0,1,Delta_Core,Delta_Shell,Delta_Solvent,Beta_Core,Beta_Shell,Beta_Solvent,Pop_FF)
		//elseifs.,....
		
		endif			

		
	//Round out the single distribution waves
		if(cmpstr(Pop_FF,"CoreShell")==0 || cmpstr(Pop_FF,"CoreShell_mod2")==0 || cmpstr(Pop_FF,"CoreShellCylinder")==0 || cmpstr(Pop_FF,"CoreShellShell")==0|| stringmatch(Pop_FF,"Janus CoreShell Micelle*"))
			MatrixOP/O GmatrixTemp = Gmatrix * 1e20			//this shape contains contrast already in...
		else
			MatrixOP/O GmatrixTemp = Gmatrix * LocalContrast*1e20		//this multiplies by scattering contrast
		endif
		
		duplicate/O VolumeDist, TepVolumeDist
		TepVolumeDist=VolumeDist[p]* IR1_BinWidthInDiameters(Core_Radius,p)
		MatrixOp/O ModelInt =GmatrixTemp x TepVolumeDist 
		Killwaves/Z GmatrixTemp
		
	elseif(Pop_Dist == 2)		
		//FUNCTIONS THAT CONTAIN TWO DISTRIBUTIONS!
		if(!cmpstr(Pop_FF,"CoreShell_Mod1"))
			//Lets grab all of our distribution waves that we created earlier.
			Wave Core_Radius = $("root:Packages:MEF_Vars:Core_Radius_pop"+num2str(pop))
			Wave Shell_radius = $("root:Packages:MEF_Vars:Shell_Radius_pop"+num2str(pop))
			Wave NumDist_Sum = $("root:Packages:MEF_Vars:NumberDistSum_pop"+num2str(pop))
			Wave VolumeDist_Sum = $("root:Packages:MEF_Vars:VolumeDistSum_pop"+num2str(pop))
			Wave VolumeDist = $("root:Packages:MEF_VArs:CompositeDist_pop"+num2str(pop))
		
			//debugg
			wave core_volumedist = $("root:Packages:MEF_Vars:core_VolumeDist_pop"+num2str(pop))
			
			Wave/Z Gmatrix=$("root:Packages:MEF_VARS:Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
			M=numpnts(ModelQ)
			N=numpnts(Core_Radius)
			O=numpnts(Shell_Radius)
		//	if(!WaveExists(Gmatrix))
				Make/D/O/N=(N,O,M) $("Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
				Wave Gmatrix=$("root:Packages:MEF_Vars:Gmatrix_set"+num2str(DataSet)+"_pop"+num2str(pop))
		//	endif	
			//calculate G matrix...
		
			MEF_CalculateGMAtrix_3D(Gmatrix,modelQ,Core_Radius,Shell_Radius,1,Delta_Core,Delta_Shell,Delta_Solvent,Beta_Core,Beta_Shell,Beta_Solvent)
		
			//Right now is whn we would covert to cm^3...not sure I care (its a factor of 1E20)
			//Since we include contast within the Gmatrix calculation we do not need to nclude it here....
			//Multiple Gmatrix by contrast here
			
			//Gap for when we add more form factors
			
			//
			
			//This is where we would add contrast...because the Coreshell is used...contrast has always occured, this is now just a unit scale
			GMatrix = GMatrix * 1E20
			//Incorporating volume dist	
			Duplicate/o VolumeDist TempVolumeDist
			//duplicate/o core_volumedist tempdebugg
			TempVolumeDist[][] = VolumeDist[p][q]*BinWidthInDiameters(Core_Radius,p)*BinWidthInDiameters(Shell_Radius,q) //The full Propwave times the two Deltas r_core and r_Shell
			//tempdebugg = core_volumedist[p]*BinWidthInDiameters(Core_Radius,p)
			SetDataFolder root:Packages:MEF_Vars	
			
			//MatrixOP
			MatrixOP/o TempGMatrix = GMatrix * TempVolumeDIst //layer by layer multiplication of teh 2D map onto the gmatrix...each point in a scattering pattern is scaled by its probability
			MatrixOP/o TempModelInt = Sum(TempGMatrix) // Sum all values in each layer...creats a 1D matrix in wrong dimension....
			ModelInt[] = TempModelInt[0][0][p] //Transpose the matrix into a row
		
		endif
	endif
		//Get and use Structure factor stuff
	
	//Quickly Save some Wave and variables
	if(cmpstr(Pop_FF,"CoreShell_Mod2")) //List of functions with variable structure factors
	
		String FormFactorSave = ("FormFactorSave_Pop"+num2str(pop)+"_En"+num2str(Dataset))
		String StructureFactorSave = ("StructureFactorSave_pop"+num2str(pop)+"_En"+num2str(Dataset))
		duplicate/o ModelInt $FOrmFactorSave
		duplicate/o ModelInt $StructureFactorSave
		wave StructureFactor = $StructureFactorSave
		
		NVAR Phi = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(pop))
		NVAR Eta = $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(pop))
		
		SVAR StrFac = $("root:Packages:MEF_Vars:StructureFactor_pop"+num2str(pop)) 
		
		StructureFactor = IR2S_CalcStructureFactor(StrFac,ModelQ,Eta,Phi,1,1,1,1)
		ModelInt *= IR2S_CalcStructureFactor(StrFac,ModelQ,Eta,Phi,1,1,1,1)
	endif

//		NVAR Index
//	duplicate/o ModelInt $("ModelInt"+num2str(INdex))
//	Wave ModelIntAppend = $("ModelInt"+num2str(INdex))
//	Index += 1
//	Appendtograph/W=Graph4 ModelIntAppend
	SetDataFolder $CurrentFolder
	
	
end

Function MEF_CalculateGMatrix_2D(Gmatrix,Model_Q,Core_Radius,FFParam1,FFPAram2,FFParam3,FFParam4,VolumePower,Delta_Core,Delta_Shell,Delta_Solvent,Beta_Core,Beta_Shell,Beta_Solvent,Pop_FF)
	Wave GMatrix //3D gmatrix to fill up with scattering data
	Wave Model_Q //The q values we want to use
	Wave Core_Radius //The radius wave for the Core given by the fit parameters and the distribution wave calculation
	Variable FFParam1,FFPAram2,FFParam3,FFPAram4 //up to 4 possible form factor params based on what we want to do (in the future)
	
	Variable VolumePower //Volume to scale at the end of the function....equals 1!!!
	
	Variable Delta_Core, Beta_Core //Core OC
	Variable Delta_Shell, Beta_Shell //Shell OC
	Variable Delta_Solvent, Beta_Solvent //SOlvent OC
	
	String Pop_FF //What form factor are we running
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFolder/o/s root:Packages:MEF_FFCalc
		
	//check the volume multiplier, shoudl be either 1 or 2 or dissasters can happen
	if(!(VolumePower==1) &&!(VolumePower==0) && !(VolumePower==2))
		Abort "Wrong input for volume muliplier in  IR1T_GenerateGMatrix, can be only 0, 1 or 2"
	endif
	
	variable M=numpnts(Model_Q)
	variable N=numpnts(Core_Radius)
	Make/D/O/N=(M) TempWave


	variable Recalculate=0
	variable i, currentR, j //loop variables and the current core/shell radius
	variable tempVal1, tempVal2,tempval3//Start and end bin values for the core and shell respectively...

	//Any other initialization can be made here
	
	
	if(!cmpstr(Pop_FF,"CoreShell"))
		Variable ShellThickness = FFParam1 //shell thickness
				
		For(i=0;i<N;i+=1)
			CurrentR = Core_Radius[i]
			TempVal1 = StartBinRadius(Core_Radius,i)
			TempVal2 = EndBinRadius(Core_Radius,i)
			
			Multithread Tempwave = MEF_CalcCoreSHellComplexIRENA(Model_Q[p],CurrentR,tempval1,tempval2,Volumepower,Shellthickness,delta_Core,beta_core,delta_shell,beta_shell,delta_solvent,beta_solvent)
			//note, the above calculated form factor contains volume^1 in it... So we need to multiply by volume^(power-1) here. Also we use volume of the core for particle volume!!!

			Multithread Tempwave *= ((4/3)*PI*(CurrentR+ShellThickness)*(CurrentR+ShellThickness)*(CurrentR+ShellThickness))^(volumepower-1)
	//		tempwave *= IR2S_CalcStructureFactor(StrFac,Model_Q,Eta,Phi,1,1,1,1)
			Gmatrix[][i]=TempWave[p]		//and here put it into G wave
		endfor

	elseif(!cmpstr(Pop_FF,"CoreShell_Mod2"))
		Variable ShellThickness_SF = FFParam1 //shell thickness
		Variable pop = FFParam2 //Population for SF extraction stuff
		
		NVAR Phi = $("root:Packages:MEF_Vars:SF_Phi_Param_pop"+num2str(2))
		NVAR Eta = $("root:Packages:MEF_Vars:SF_Eta_Param_pop"+num2str(2))
		Make/O/N=2 SFWave = {eta, phi}

		For(i=0;i<N;i+=1)
			CurrentR = Core_Radius[i]
			TempVal1 = StartBinRadius(Core_Radius,i)
			TempVal2 = EndBinRadius(Core_Radius,i)
			SFWave[0] = eta+CurrentR
			Multithread Tempwave = MEF_CalcCoreSHellComplexIRENA_SF(Model_Q[p],CurrentR,tempval1,tempval2,Volumepower,Shellthickness_SF,delta_Core,beta_core,delta_shell,beta_shell,delta_solvent,beta_solvent,SFWave)
			//note, the above calculated form factor contains volume^1 in it... So we need to multiply by volume^(power-1) here. Also we use volume of the core for particle volume!!!

			Multithread Tempwave *= ((4/3)*PI*(CurrentR+ShellThickness)*(CurrentR+ShellThickness)*(CurrentR+ShellThickness))^(volumepower-1)
	//		tempwave *= IR2S_CalcStructureFactor(StrFac,Model_Q,Eta,Phi,1,1,1,1)
			Gmatrix[][i]=TempWave[p]		//and here put it into G wave
		endfor
		
	endif
	//Convert units to match what Jan does...this puts it in terms of cm^3...later scaled by 1E20 with the contrast value
	GMatrix = Gmatrix*(1E-24)^VolumePower //Volumepower will always be 1.
//	GMatrix_DoneTime = StopMSTimer(GMatrixCalcTimer)
	
//	print "G-matrix calculated: "+num2str(GMatrix_DoneTime/(1E6)) +" seconds"

	SetDataFolder $CUrrentFOlder
end
Function MEF_CalculateGMatrix_3D(Gmatrix,Model_Q,Core_Radius,Shell_Radius,VolumePower,Delta_Core,Delta_Shell, Delta_Solvent,Beta_Core,Beta_shell,Beta_Solvent)

	Wave GMatrix //3D gmatrix to fill up with scattering data
	Wave Model_Q //The q values we want to use
	Wave Core_Radius //The radius wave for the Core given by the fit parameters and the distribution wave calculation
	Wave Shell_Radius //Same but for the shell radius
	
	Variable VolumePower //unsure of what this does so far
	
	Variable Delta_Core, Beta_Core //Core OC
	Variable Delta_Shell, Beta_Shell //Shell OC
	Variable Delta_Solvent, Beta_Solvent //SOlvent OC
	
	String CurrentFolder = GetDataFolder(1)
	NewDataFolder/o/s root:Packages:MEF_FFCalc

	
	//check the volume multiplier, shoudl be either 1 or 2 or dissasters can happen
	if(!(VolumePower==1) &&!(VolumePower==0) && !(VolumePower==2))
		Abort "Wrong input for volume muliplier in  IR1T_GenerateGMatrix, can be only 0, 1 or 2"
	endif
	
	variable M=numpnts(Model_Q)
	variable N=numpnts(Core_Radius)
	Variable O=Numpnts(Shell_Radius)
	
	Make/D/O/N=(M) TempWave
	
	variable Recalculate=0
	variable i, currentR_C,currentR_S, j //loop variables and the current core/shell radius
	variable tempVal1_C, tempVal2_C, tempval1_S,tempval2_S //Start and end bin values for the core and shell respectively...
//	string OldNote=note(Gmatrix)
//	string NewNote = "", VolDefL=""
//	string reason=""
	
	//I don't think I am going to keep any of the recalculate questions...just remake it every time...tough shit
	
	//Running code Line 581 or IR1_FormFactors.ipf to setup stuff...this is gonna be loooong

	//Starting Coreshell info
	//Initialization things are already done earlier with Deltas and Betas
	For(i=0;i<N;i+=1) //Calculate in chunks of the gmatrix...this is core loop
		CurrentR_C = Core_Radius[i]
		TempVal1_C = StartBinRadius(Core_Radius,i)
		TempVal2_C = EndBinRadius(Core_Radius,i)
		//print TempVal1_C, TempVal2_C
	//	print "RUNNING A NEW CORE RADII " + num2str(i)
		//New timer for loop tie
//		GMatrix_RadTimer = StartMSTimer
		
		for(j=0;j<O;j+=1) //Calculate now the shell loop
			CurrentR_S = Shell_Radius[j]
			TempVal1_S = StartBinRadius(Shell_Radius,j)
			TempVal2_S = EndBinRadius(Shell_Radius,j)
			
			//This should be everything we want before we begin the looping
			Multithread tempwave = MEF_CalcCoreShellFFPoints(Model_Q[p],CurrentR_C,CurrentR_S,Volumepower,tempval1_C,tempval2_C,tempval1_S,tempval2_S,delta_Core,beta_core,delta_shell,beta_shell,delta_solvent,beta_solvent)
			//tempwave = MEF_CalcCoreShellFFPoints(Model_Q[p],Core_Radius,Shell_Radius,Volumepower,tempval1_C,tempval2_C,tempval1_S,tempval2_S,delta_Core,beta_core,delta_shell,beta_shell,delta_solvent,beta_solvent)

			//want to take out the factor of volume that is still sitting around
			Multithread Tempwave *= ((4/3)*PI*(CurrentR_C+CurrentR_S)*(CurrentR_C+CurrentR_S)*(CurrentR_C+CurrentR_S))^(volumepower-1)
			//print "NOT YET FROZEN"
			Gmatrix[i][j][] = Tempwave[r]
		
		endfor
		
		//Print timer
//		Gmatrix_RadTimeCheck = StopMSTimer(GMatrix_RadTimer)
	//	print "Radii "+num2str(i) + "Done: "+num2str(GMatrix_RadTimeCHeck/(1E6)) +" seconds"
		
	endfor
	//Convert units to match what Jan does...this puts it in terms of cm^3...later scaled by 1E20 with the contrast value
	GMatrix = Gmatrix*(1E-24)^VolumePower //Volumepower will always be 1.
//	GMatrix_DoneTime = StopMSTimer(GMatrixCalcTimer)
	
//	print "G-matrix calculated: "+num2str(GMatrix_DoneTime/(1E6)) +" seconds"

	SetDataFolder $CurrentFOlder

end


//DIAGNOSTICS
Function ReadGMatrix(GMatrix,corept,shellpt)
	Wave GMatrix
	Variable CorePt // What location in Core do you want
	Variable ShellPt //What location in Shell do you want
	
	Variable NumQPts = dimsize(GMatrix,2)
	
	Make/N=(numQpts)/O GMat_Wave =0
	
	GMat_Wave[] = GMatrix[Corept][ShellPt][p] 



end

//Standard CoreShell Model with complex values....adopted from IRENA by Jan Ilavsky
ThreadSafe Function MEF_CalcCoreSHellComplexIRENA(QValue,radius,radiusmin,radiusmax,VOlumePower,ShellRadius,Delta_C,beta_C,Delta_S,Beta_S,Delta_Sol,Beta_Sol)

	Variable Qvalue //The qvalue we are calculating all this stuff from
	Variable radius
	Variable RadiusMin //Core radius min for bin
	Variable Radiusmax //Core Radius max for bin,
	Variable VolumePower //Equal to 1
	Variable ShellRadius //The shell radius (fit parameter)
	
	Variable Delta_C, Beta_C //Optical constants for core
	Variable Delta_S, Beta_S //Optical constants for Shell
	Variable Delta_Sol, Beta_Sol //Optical constants for Solvent
	
	////////Set up full complex optical constants for contrast
	
	Variable/C n_Core = cmplx(1-Delta_C,Beta_C)
	Variable/C n_Shell = cmplx(1-Delta_S,Beta_S)
	Variable/C n_Solvent = cmplx(1-Delta_sol,Beta_sol)
	//INITIAL CALCULATION OF THE CORE SCATTERING
	Variable QR = Qvalue*Radius
	Variable QRMin = Qvalue*Radiusmin
	Variable QRMax = QValue*Radiusmax

	//other parameters
	variable tempResult
	variable result=0							//here we will stash the results in each point and then divide them by number of points
	Variable/C ComplexResult =0
	variable tempRad
	
	variable numbOfSteps=floor(3+abs((10*(QRMax-QRMin)/pi)))		//depending on relationship between QR and pi, we will take at least 3
	if (numbOfSteps>60)											//steps in QR - and maximum 60 steps. Therefore we will get reaasonable average 
		numbOfSteps=60											//over the QR space. 
	endif
	variable step=(QRMax-QRMin)/(numbOfSteps-1)					//step in QR
	variable stepR=(radiusMax-radiusMin)/(numbOfSteps-1)			//step in R associated with above, we need this to calculate volume
	variable i
	
	For (i=0;i<numbOfSteps;i+=1)									//here we go through number of points in QR (and R)
		QR=QRMin+i*step
		tempRad=radiusMin+i*stepR

		tempResult=(3/(QR*QR*QR))*(sin(QR)-(QR*cos(QR)))				//calculate sphere scattering factor 
	
		result+=tempResult//* (IR1T_SphereVolume(tempRad))				//scale by volume add the values together...
	endFor
	result=result/numbOfSteps											//this averages the values obtained over the interval....
	result=result*((4/3)*PI*(Radius*Radius*Radius))						//multiply by volume of sphere)
	ComplexResult = result*(n_core - n_shell)
//Now the shell

	QRMin=Qvalue*(radiusMin+ShellRadius)
	QRMax=Qvalue*(radiusMax+ShellRadius)
	step=(QRMax-QRMin)/(numbOfSteps-1)	
	stepR=((radiusMax+ShellRadius)-(radiusMin+ShellRadius))/(numbOfSteps-1)
	Variable result1 = 0
	Variable/C ComplexResult1 = 0
	
	For (i=0;i<numbOfSteps;i+=1)									//here we go through number of points in QR (and R)
		QR=QRMin+i*step
		tempRad=radiusMin+ShellRadius+i*stepR

		tempResult=(3/(QR*QR*QR))*(sin(QR)-(QR*cos(QR)))				//calculate sphere scattering factor 
	
		result1+=tempResult//*(IR1T_SphereVolume(tempRad)) 			//and add the values together...
	endFor
	result1 = result1/numbofsteps
	result1 = result1 * ((4/3)*PI*(Radius+ShellRadius)*(Radius+ShellRadius)*(Radius+ShellRadius))
	ComplexResult1 = result1*(n_shell-n_solvent)
//	print (n_shell - n_solvent)
	Variable FinalResult = magsqr(ComplexResult + ComplexResult1) //Uncommented on 12/12/2019...trying to better incorporate the fit parameters

//	Variable FinalResult = (result^2 * abs(Delta_C)) + (result1^2 * abs(Delta_s)) + 2*(result*result1 * (Delta_sol)) //Manually do the magnitude square stuff

	FinalResult = Finalresult / ((4/3)*PI*(Radius+ShellRadius)*(Radius+ShellRadius)*(Radius+ShellRadius)) //normalizes to volume of the core only? May want to play with this later.


//
//	if(abs(n_shell - n_solvent) == 0)
//		FinalResult = Finalresult / ((4/3)*PI*(Radius)*(Radius)*(Radius)) //normalizes to volume of the core only? May want to play with this later.
//	else
//		FinalResult = Finalresult / ((4/3)*PI*(Radius+ShellRadius)*(Radius+ShellRadius)*(Radius+ShellRadius)) //normalizes to volume of the core only? May want to play with this later.
//	endif
	return Finalresult



end





//Coreshell model with variable structure factor based on core size....adopted from IRENA by Jan Ilavsky

ThreadSafe Function MEF_CalcCoreSHellComplexIRENA_SF(QValue,radius,radiusmin,radiusmax,VOlumePower,ShellRadius,Delta_C,beta_C,Delta_S,Beta_S,Delta_Sol,Beta_Sol,SFWave)

	Variable Qvalue //The qvalue we are calculating all this stuff from
	Variable radius
	Variable RadiusMin //Core radius min for bin
	Variable Radiusmax //Core Radius max for bin,
	Variable VolumePower //Equal to 1
	Variable ShellRadius //The shell radius (fit parameter)
	
	Variable Delta_C, Beta_C //Optical constants for core
	Variable Delta_S, Beta_S //Optical constants for Shell
	Variable Delta_Sol, Beta_Sol //Optical constants for Solvent
	
	wave SFWave //to grab SF info
	
	////////Set up full complex optical constants for contrast
	
	Variable/C n_Core = cmplx(1-Delta_C,Beta_C)
	Variable/C n_Shell = cmplx(1-Delta_S,Beta_S)
	Variable/C n_Solvent = cmplx(1-Delta_sol,Beta_sol)

	//INITIAL CALCULATION OF THE CORE SCATTERING
	Variable QR = Qvalue*Radius
	Variable QRMin = Qvalue*Radiusmin
	Variable QRMax = QValue*Radiusmax

	//other parameters
	variable tempResult
	variable result=0							//here we will stash the results in each point and then divide them by number of points
	Variable/C ComplexResult =0
	variable tempRad
	
	variable numbOfSteps=floor(3+abs((10*(QRMax-QRMin)/pi)))		//depending on relationship between QR and pi, we will take at least 3
	if (numbOfSteps>60)											//steps in QR - and maximum 60 steps. Therefore we will get reaasonable average 
		numbOfSteps=60											//over the QR space. 
	endif
	variable step=(QRMax-QRMin)/(numbOfSteps-1)					//step in QR
	variable stepR=(radiusMax-radiusMin)/(numbOfSteps-1)			//step in R associated with above, we need this to calculate volume
	variable i
	
	For (i=0;i<numbOfSteps;i+=1)									//here we go through number of points in QR (and R)
		QR=QRMin+i*step
		tempRad=radiusMin+i*stepR

		tempResult=(3/(QR*QR*QR))*(sin(QR)-(QR*cos(QR)))				//calculate sphere scattering factor 
	
		result+=tempResult//* (IR1T_SphereVolume(tempRad))				//scale by volume add the values together...
	endFor
	result=result/numbOfSteps											//this averages the values obtained over the interval....
	result=result*((4/3)*PI*(Radius*Radius*Radius))						//multiply by volume of sphere)
	ComplexResult = result*(n_core - n_shell)
//Now the shell

	QRMin=Qvalue*(radiusMin+ShellRadius)
	QRMax=Qvalue*(radiusMax+ShellRadius)
	step=(QRMax-QRMin)/(numbOfSteps-1)	
	stepR=((radiusMax+ShellRadius)-(radiusMin+ShellRadius))/(numbOfSteps-1)
	Variable result1 = 0
	Variable/C ComplexResult1 = 0
	
	For (i=0;i<numbOfSteps;i+=1)									//here we go through number of points in QR (and R)
		QR=QRMin+i*step
		tempRad=radiusMin+ShellRadius+i*stepR

		tempResult=(3/(QR*QR*QR))*(sin(QR)-(QR*cos(QR)))				//calculate sphere scattering factor 
	
		result1+=tempResult//*(IR1T_SphereVolume(tempRad)) 			//and add the values together...
	endFor
	result1 = result1/numbofsteps
	result1 = result1 * ((4/3)*PI*(Radius+ShellRadius)*(Radius+ShellRadius)*(Radius+ShellRadius))
	ComplexResult1 = result1*(n_shell-n_solvent)
//	print (n_shell - n_solvent)
	Variable FinalResult = magsqr(ComplexResult + ComplexResult1)

	if(abs(n_shell - n_solvent) == 0)
		FinalResult = Finalresult / ((4/3)*PI*(Radius)*(Radius)*(Radius)) //normalizes to volume of the core only? May want to play with this later.
	else
		FinalResult = Finalresult / ((4/3)*PI*(Radius+ShellRadius)*(Radius+ShellRadius)*(Radius+ShellRadius)) //normalizes to volume of the core only? May want to play with this later.
	endif
	
	FinalResult *= IR2S_HardSphereStruct(SFWave,Qvalue)
	return Finalresult



end

///Core shell calculation with SHell Radius as a distribution
Threadsafe Function MEF_CalcCoreShellFFPoints(Qvalue,radius_C,radius_S,VolumePower,radiusmin_C,radiusmax_C,radiusmin_S,radiusmax_S,Delta_C,beta_C,Delta_S,Beta_S,Delta_Sol,Beta_Sol)

	Variable Qvalue //The qvalue we are calculating all this stuff from
	Variable Radius_C,Radiusmax_C,radiusmin_C //Core radius information
	Variable Radius_S,RadiusMax_S,RadiusMin_S //Shell radius information
	
	Variable VolumePower //no idea yet
	
	Variable Delta_C, Beta_C //Optical constants for core
	Variable Delta_S, Beta_S //Optical constants for Shell
	Variable Delta_Sol, Beta_Sol //Optical constants for Solvent
	
	////////Set up full complex optical constants for contrast
	
	Variable/C n_Core = cmplx(1-Delta_C,Beta_C)
	Variable/C n_Shell = cmplx(1-Delta_S,Beta_S)
	Variable/C n_Solvent = cmplx(1-Delta_sol,Beta_sol)

	//INITIAL CALCULATION OF THE CORE SCATTERING
	Variable QR = Qvalue*Radius_C
	Variable QRMin = Qvalue*Radiusmin_C
	Variable QRMax = QValue*Radiusmax_C
	//print QRmin, QRMax
	
//	//Shell limiting values
//	
//	Variable QR_S = Qvalue*Radius_S
//	Variable QRMin_S = Qvalue*Radiusmin_S
//	Variable QRMax_S = QValue*Radiusmax_S
//	
	//other parameters
	variable tempResult
	variable result=0							//here we will stash the results in each point and then divide them by number of points
	Variable/C ComplexResult =0
	variable tempRad
	
	Variable numbOfSteps = floor(3+abs((10*(QRMax - QRMin)/pi)))
//	Variable numbOfSteps_S = floor(3+abs((10*(QRMax_S - QRMin_S)/pi))) //no idfea what the point of this is..but whatevs
	if(numbOfSteps > 60)
		numbOfSteps = 60
//	elseif(numbOfSteps_S > 60)
//		numbOfSteps_S = 60
	endif //don't go to crazy with step sizes.
	//print numbofsteps
	//print ""
	Variable Step = (QRmax - QRMin)/(numbofsteps - 1)
	Variable StepR = (RadiusMax_C - RadiusMin_C)/(numbofsteps-1)
//	Variable Step_S = (QRmax_S - QRMin_S)/(numbofsteps_S - 1)
//	Variable Step_RS = (RadiusMax_S - RadiusMin_S)/(numbofsteps_S-1)

	Variable i,j //loops I am sure I will ned rifght now
	
	For(i=0;i<numbofsteps;i+=1)
		
		QR = QRMin + i*step
		TempRad = Radiusmin_C + i*StepR
		
		TempResult = (3/(QR*QR*QR))*(sin(QR) - (QR*cos(QR)))
		
		result += Tempresult	
	endfor
	result = result/numbofsteps
	result = result * ((4/3)*PI*(Radius_C*RAdius_C*Radius_C))
	ComplexResult = result*(n_core - n_shell)
	
	//Now complete the shell component
	QRmin = QValue*(Radiusmin_C + Radiusmin_S)
	QRMax = Qvalue*(Radiusmax_C + Radiusmax_S)
	step = (QRMax - QRmin)/(numbofsteps-1)
	StepR = ((Radiusmax_C+RadiusMax_S) - (Radiusmin_C +Radiusmin_S))/(numbofsteps-1)
	Variable result1 = 0
	Variable/C ComplexResult1 = 0
	
	for(i=0;i<numbofsteps;i+=1)
		QR = QRmin + i*step
		TempRad = Radiusmin_C + radiusmin_S + i*StepR
		
		Tempresult = (3/(QR*QR*QR))*(sin(QR) - (QR*cos(QR)))
		
		result1 += tempresult	
	endfor
	
	result1 = result1/numbofsteps
	result1 = result1 * ((4/3)*PI*(Radius_C+Radius_S)*(Radius_C+Radius_S)*(Radius_C+Radius_S))
	ComplexResult1 = result1*(n_shell-n_solvent)
	
	Variable FinalResult = magsqr(ComplexResult + ComplexResult1)
	FinalResult = Finalresult / ((4/3)*PI*(Radius_C+Radius_S)*(Radius_C+Radius_S)*(Radius_C+Radius_S)) //normalizes to volume of the core only? May want to play with this later.
	
	return Finalresult
end

//Function MEF_Structurefactor(SFName,QValue,Param1,Param2,Param3,Param4,Param5.Param6)
//	String SFName //pick the name of the structure factor
//	Variable Qvalue,Param1,Param2,Param3,Param4,Param5,Param6
//
//	Variable result
//	String CurrentFOlder = GetDataFOlder(1)
//	
//	NewDataFolder/o/s root:Packages:StructureFactorCalc
//	Make/O/N=6 ParWv //Parameter wave
//	
//	if(cmpstr(SFName,"HardSpheres")==0
//		result = IR2S_HardSphereStruct(ParWv,QValue)
//	else
//		result = 1
//	endif
//	
//	SetDataFolder $CurrentFolder
//	return result
//endif
//	
//
//end