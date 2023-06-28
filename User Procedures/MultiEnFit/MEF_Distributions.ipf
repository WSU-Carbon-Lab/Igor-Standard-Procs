#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

///Initialize all the waves required for later in calculating the distribution

Function MEF_InitializeDistribution(pop,NameofDist,DistPnts,DistPrecision,DistType,Meansize,width,scale)
	Variable pop //what population are we working with
	String NameofDist
	String DistType
	Variable DistPnts //Number of points in the distribution
	Variable DistPrecision //Eactly what it sounds like

	
	//FIT PARAMETER INITIALIZATION
	Variable Meansize // Mean size of distribution
	Variable Width		//Width of whatever distribution this is
	Variable Scale		//Optional Scale factor that ill most likely be held @ 1
	
	
	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
	
		
	//Debugging
	
	//Parameter initializations (might not want here)
	
//	//Var Names
//	String MeanSizeName = NameofDist + "_DistMeanSize_pop"+num2str(pop)
//	String WidthName = NameofDist+ "_DistWidth_pop"+num2str(pop)
//	String ScaleName = NameofDist+ "_DistScale_pop"+Num2str(pop)
//	//Create globals for later use
//	Variable/G $MeanSizeName = MeanSize
//	Variable/G $WidthName = Width
//	Variable/G $ScaleName = Scale
		
	//Make the radius and other distribution waves
	Make/O/N=(Distpnts) $(NameofDist+"_Radius_Pop"+num2str(pop)),$(NameofDist+"_Diameter_Pop"+num2str(pop)),$(NameofDist+"_VolumeDist_Pop"+num2str(pop)),$(NameofDist+"_NumberDist_Pop"+num2str(pop))
	Wave Radius=$(NameofDist+"_Radius_Pop"+num2str(pop))
	Wave Diameter=$(NameofDist+"_Diameter_Pop"+num2str(pop))
	Wave VolumeDist = $(NameofDist+"_VolumeDist_Pop"+num2str(pop))
	Wave NumberDist = $(NameofDist+"_NumberDist_Pop"+num2str(pop))
	//Use Jans code to build the radius distribution because it works really well
//	print "MeanSize"+num2str(meansize)
//	print "Width_"+num2str(width)
	
//	if(width > 300)
//		width = 300
//	elseif(meansize > 500)
//		meansize = 500
//	endif
	IR2L_GenerateRadiiDist(DistType,Radius,DistPnts,DistPrecision,meansize,width,width) //Last two params for SZ are the same....fine for now
	//print "end"
	Diameter = 2*Radius
	//Creates the distributions based o the recent creation of the Radius Wave
	Duplicate/Free VolumeDist TempDist, TempDistL
	Redimension/D TempDist,TempDistL
	//This is outlined in the function MEF_CalcDistributions...I want to incorporate it here so all distribution creation is done simultaneiously
	if(Stringmatch(DistType,"Schulz-Zimm")) //Calc radius distributions
		TempDist = MES_Calc_ShultzZimm(Radius[p],MeanSize,Width,0)
	endif
	
	//Quick normalization if the area under the Shultz zimm is not quite 1...sometime it is around 0.99...this fixes that problem
	//Also includes the scale factor
	Variable ScaleDist = scale / AreaXY(Radius,TempDist,-inf,inf)
	TempDist *= ScaleDist
	
	VolumeDist = TempDist //These will be changed later for the composite volume calculation...for now I am just saving the data
	NumberDist = TempDist //Again, not complete information till later.

	SetDataFolder $CurrentFOlder

end

Function MEF_CoreShellCompositeDistribution(pop,CoreDistName,ShellDistName)
	
	Variable pop
	String CoreDistName //Name of the Core Distribution stuff
	String ShellDistName //Name of teh Shell Distribution Stuff
	
	//Distinction is made here because the core and shell volume calculations are different based on what is the core/shell
	
	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
	
	//Grab the parametrs that I need 
	NVAR CoreScale = $("root:Packages:MEF_Vars:"+CoreDistName+"_DistScale_Pop"+num2str(pop))
	NVAR ShellScale = $("root:Packages:MEF_Vars:"+ShellDistName+"_DistScale_Pop"+num2str(pop))


	//Grab waves for calculations
	Wave Core_RadiusWave = $(CoreDistName + "_Radius_pop"+num2str(pop))
	Wave Core_VolumeDist = $(CoreDistName + "_VolumeDist_pop"+num2str(pop)) //Currently the full distribution for the core
	Wave Core_NumberDist = $(CoreDistName + "_NumberDIst_pop"+num2str(pop))
	
	Wave Shell_RadiusWave = $(ShellDistName + "_Radius_pop"+num2str(pop))
	Wave Shell_VolumeDist = $(ShellDistName + "_VolumeDist_pop"+num2str(pop)) //Curently the full distribution for the shell
	Wave Shell_NumberDist = $(ShellDistName + "_NumberDist_pop"+num2str(pop))

//Make temp waves that we can update and later change starting waves
	Duplicate/Free Core_VolumeDist, Core_TempDist, Core_TempVolDistL
	Redimension/D Core_TempDist, Core_TempVolDistL //Makes them higher precision waves for calculations as to not lose data
	//And again for shell distributions
	Duplicate/Free Shell_VolumeDist, Shell_TempDist, Shell_TempVolDistL
	Redimension/D Shell_TempDist, Shell_TempVolDistL //Makes them higher precision waves for calculations as to not lose data
	//Check for the type of distribution we want....can add more later
	
	//Create composite distributions
	String TotalDistName = ("CompositeDist_pop" + num2str(pop))
	String VOlumeDistSum = "VolumeDistSum_pop"+num2str(pop)
	String NumDistSum = "NumDistSum_pop"+num2str(pop)
	Make/N=(numpnts(Core_RadiusWave),numpnts(Shell_RadiusWave))/O $VolumeDistSum, $NumDistSum
	Make/N=(numpnts(Core_RadiusWave),numpnts(Shell_RadiusWave))/O $TotalDistName

	Wave VolumeDistSum_w = $VolumeDistSum
	Wave NumDistSum_w = $NumDistSum
	Wave TotalDist = $TotalDistName
	//Create the composite Distribution
	TotalDist[][] = Core_VolumeDist[p] * Shell_VolumeDist[q]

	
	MEF_CoreShell_AveVolumeDist(VolumeDistSum_w,Core_TempVolDistL,SHell_TempVolDistL,Core_RadiusWave,Shell_radiusWave,pop)
	
	//Calculate the Num/VolumeDistributions
	
	Core_VolumeDist = Core_TempDist
	Shell_VolumeDist = Shell_TempDist
	Core_NumberDist = Core_VolumeDist / Core_TempVolDistL
	Shell_NumberDist = Shell_VolumeDist / Shell_TempVolDistL

	//NumDistSum_w = TotalDist / VolumeDistSum_w
	
	//Normalize to the scale factor and the area under the distribution which is almost (not always 1...sometimes 0.99....)

	Variable ScaleCore,ScaleShell
	ScaleCore = CoreScale / AreaXY(Core_RadiusWave,Core_VolumeDist,-inf,inf)
	ScaleShell = ShellScale / AreaXY(Shell_RadiusWave,Shell_VolumeDist,-inf,inf)
	
	Core_VolumeDist *= ScaleCore
	Core_NumberDist *= ScaleCore
	Shell_VolumeDist *= ScaleShell
	SHell_NumberDist *= ScaleShell
	
	//Finally add together the possible Distributions
	NumDistSum_w = TotalDist / VolumeDistSum_w
	
	//All distributions should be done now, Now we want to calculate the scattering with newfound volumes.
	SetDataFolder $CUrrentFOlder
end





Function MES_Calc_ShultzZimm(x,meanpos,width,shape)
	Variable x //This is the x-wave of the distribution that is made earlier...NOT A FIT PARAM
	Variable meanpos //This is the meansize of the core radius we want to fit //FIT PARAMETER!!!!
	Variable width //This is the width of the distribution that we want to fit //FIT PARAMETER
	Variable shape //This seems to just be 0....No idea why I have it here //MIGHT JUST REMOVE
	//Taken from Line 130 of IR1_Functions.ipf in IRENA...all credit to Jan, adapted here so I can change things without impacting Irena
	//no Typos as of 9:21am 12/19/2018 (works)
	Variable result, a, b
	
	b = 1/(width/(2*Meanpos))^2
	a = b/meanpos
	
	if(b<70)
		result = ( (a^(b+1))/gamma(b+1) * x^b / exp(a*x) )
	else //run it with log numbers to avoid large numbers //???unsure of what this means
		result = exp( (b+1)*ln(a) - gammln(b+1) + B*ln(x) - (a*x) )
	endif
	
	if (numtype(result)!=0)
		result = 0
	endif
	
	return result
end




Function MEF_CoreShell_AveVolumeDist(ResultsWave,Core_VolumeDist,Shell_VolumeDist,Core_RadiusWave,Shell_RadiusWave,pop)
	
	Wave ResultsWave //2D map for the full averaged volume distributions together as the full particle
	Wave Core_VolumeDist //The volume distributio that you want to update with your newfound cool radius distribution
	Wave Core_RadiusWave //The same radius wave that the distribution will be created
	Wave Shell_VolumeDist //Same for shell
	Wave Shell_RadiusWave //Same for Shell
	Variable pop //What population is this being created fotr
	
	String CurrentFolder = GetDataFolder(1)
	
	//Parameters that will look into each statistical bin of the radius distribution and average the volume around each point
	variable i,j ,m,n//loop variables
	Variable C_StartValue //Start bin value
	Variable C_endValue //End bin value (These values surround each data point in the distribution to be averaged)
	Variable C_TempVolume //Exactly what it soudns like
	Variable C_TempRadius //See above
	
	Variable S_StartValue //Start bin value
	Variable S_endValue //End bin value (These values surround each data point in the distribution to be averaged)
	Variable S_TempVolume //Exactly what it soudns like
	Variable S_TempRadius //See above
	
	
	String cmd2, infostr //might not need these
	
	SetDataFolder root:Packages
	NewDataFolder/o/S root:Packages:FormFactorCalc
	Variable/G TempVolCalc // Temp volume that will be updated for each bin
	String VolDefL //Definition of the core shell volume will always be the same for now
	
	VolDefL = "Whole Particle" //allows for later calculations of the volume distribution of the shell or radius solo...
	//NOT IMPLEMENTED YET!!!
	
	//Begin cylcing through the Radius of the Core..Will calculate both core and shell independent..then add....
	for(i=0 ; i<numpnts(Core_RadiusWave); i+=1)
	
		C_StartValue = StartBinRadius(Core_RadiusWave,i)
		C_EndValue = EndBinRadius(Core_RadiusWave,i)
		
		C_TempVolume = 0
	//	TempVolCalc =0
		
		For(j=0;j<50;j+=1)
			C_TempRadius = (C_StartValue + j*(C_Endvalue-C_StartValue)/50)
			C_TempVolume += (4/3)*(pi*C_TempRadius*C_TempRadius*C_TempRadius)		
		endfor
		C_TempVolume /= 50
		//C_TempVolume *= 10^(-24) //convert A to cm
		Core_VolumeDist[i] = C_TempVolume
		C_TempRadius = ((3/4)*(1/PI)*C_TempVolume)^(1/3) //THIS IS NOW TEH RADIUS CORRESPONDING TO AVE VOLUME
	
	//Start calculating possible Shell Radii around this new core radius average
		For(m=0; m <numpnts(SHell_RadiusWave); m+=1)
		
			S_StartValue = StartBinRadius(Shell_RadiusWave,m)
			S_EndValue = EndBinRadius(Shell_RadiusWave,m)
		
			S_TempVolume = 0
			//TempVolCalc =0
			
			For(n=0;n<50;n+=1)
				S_TempRadius = (S_StartValue + n*(S_Endvalue-S_StartValue)/50)
				S_TempVolume += (4/3)*(pi*(C_TempRadius+S_TempRadius)*(C_TempRadius+S_TempRadius)*(C_TempRadius+S_TempRadius)) - C_TempVolume	
			endfor
		
			S_TempVolume /= 50
			//S_TempVolume *= 10^(-24) //convert A to cm
			Shell_VolumeDist[m] = S_TempVolume
		
		endfor
		
	
	
	endfor

	
	//Bring the distributions together by looking at the Entire Particle Radii
	
	ResultsWave[][] = Core_VolumeDist[p] + Shell_VolumeDist[q]
	
	
	
	
	
	
end


Function StartBinRadius(RadiusWave,i)
	Wave radiusWave
	Variable i
	
	variable start
	variable Imax=numpnts(RadiusWave)
	
	if (i==0)
		start=RadiusWave[0]-(RadiusWave[1]-RadiusWave[0])/2
		if (start<0)
			start=1		//we will enforce minimum size of the scatterer as 1 A
		endif
	elseif (i==Imax-1)
		start=RadiusWave[i]-(RadiusWave[i]-RadiusWave[i-1])/2
	else
		start=RadiusWave[i]-((RadiusWave[i]-RadiusWave[i-1])/2)
	endif
	return start


end

Function EndBinRadius(RadiusWave,i)

	Wave RadiusWave
	Variable i
	
	variable endL
	variable Imax=numpnts(RadiusWave)
	
	if (i==0)
		endL=RadiusWave[0]+(RadiusWave[1]-RadiusWave[0])/2
	elseif (i==Imax-1)
		endL=RadiusWave[i]+((RadiusWave[i]-RadiusWave[i])/2)
	else
		endL=RadiusWave[i]+((RadiusWave[i+1]-RadiusWave[i])/2)
	endif
	return endL
	
	
end

Function BinWidthInDiameters(D_distribution,i)			//calculates the width in diameters by taking half distance to point before and after
	variable i								//returns number in A
	Wave D_distribution
	
	variable width
	variable Imax=numpnts(D_distribution)
	
	if (i==0)
		width=D_distribution[1]-D_distribution[0]
		if ((D_distribution[0]-(D_distribution[1]-D_distribution[0])/2)<0)
			width=D_distribution[0]+(D_distribution[1]-D_distribution[0])/2
		endif
	elseif (i==Imax-1)
		width=D_distribution[i]-D_distribution[i-1]
	else
		width=((D_distribution[i]-D_distribution[i-1])/2)+((D_distribution[i+1]-D_distribution[i])/2)
	endif
	return abs(width)		//9/17/2010, fix for user models when bins are sorted from large to small
end





/////****************************OLD CODE

//Trimmed down version of the Irena function IR2L_CalculateDistrubionts (line 1612 IR2L_NLSQFCalc.ipf)
//Want the option to change things if I need too
////OLD CODE AS OF 1/3/2019
Function MEF_CalcDistributions(pop,Core_RadiusWave,Shell_RadiusWave,Core_NumDist,Shell_NumDist,Core_VolumeDist,Shell_VolumeDist)
	Variable pop //Current pop value..will be cycled outside this func
	Wave Core_RadiusWave, Shell_RadiusWave // xw for Radii in your distribution
	Wave Core_NumDist,Shell_NumDist //Final Wave for number distribution
	Wave Core_VolumeDist,Shell_VolumeDist //Final wave for Volume Distribution
	
	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
	
	//Grab global variables that we set at the start of the function
	SVAR DistShape = $("root:Packages:MEF_vars:PopDistShape_pop"+num2str(pop))
	NVAR Core_SZMeanSize = $("root:Packages:MEF_Vars:Core_SZDistMeanSize_pop"+num2str(pop)) //specific pop meansize for ShultzZimm
	NVAR Core_SZWidth = $("root:packages:MEF_Vars:Core_SZDistWidth_pop"+num2str(pop)) //Specific pop width for Shultzzimm
	NVAR Core_SZScale = $("root:Packages:MEF_Vars:Core_Scale_pop"+num2str(pop))
	NVAR Shell_SZMeanSize = $("root:Packages:MEF_Vars:Shell_SZDistMeanSize_pop"+num2str(pop)) //specific pop meansize for ShultzZimm
	NVAR Shell_SZWidth = $("root:packages:MEF_Vars:Shell_SZDistWidth_pop"+num2str(pop)) //Specific pop width for Shultzzimm
	NVAR Shell_SZScale = $("root:Packages:MEF_Vars:Shell_Scale_Pop"+num2str(pop))
	
	//Make temp waves that we can update and later change starting waves
	Duplicate/Free Core_VolumeDist, Core_TempDist, Core_TempVolDistL
	String TotalDistName = ("TotalDist_pop" + num2str(pop))
	Make/N=(numpnts(Core_RadiusWave),numpnts(Shell_RadiusWave))/O $TotalDistName
	Wave TotalDist = $TotalDistName
	Redimension/D Core_TempDist, Core_TempVolDistL,TotalDist //Makes them higher precision waves for calculations as to not lose data
	//And again for shell distributions
	Duplicate/Free Shell_VolumeDist, Shell_TempDist, Shell_TempVolDistL
	Redimension/D Shell_TempDist, Shell_TempVolDistL //Makes them higher precision waves for calculations as to not lose data
	//Check for the type of distribution we want....can add more later
	
	//CALCULATES ONLY THE PROBABILITY FUNCTION...WILL BE SCALED BY VOLUMES LATER>...AREA OF THIS CURVE is ~1
	if(StringMatch(DistShape,"SchultzZimm")) //Calc Radius Distributions
		Core_TempDist = MES_Calc_ShultzZimm(Core_RadiusWave[p],Core_SZMeanSize,Core_SZWidth,0) // Use functyion below...no diea what the last value is
		Shell_TempDist = MES_Calc_ShultzZimm(Shell_RadiusWave[p],Shell_SZMeanSize,Shell_SZWidth,0) // Use functyion below...no diea what the last value is
		
		TotalDist[][] = Core_TempDist[p] * Shell_TempDist[q]
	
	//elseif()//Add other distributions here
	
	endif
	String VOlumeDistSum = "VolumeDistSum_pop"+num2str(pop)
	String NumDistSum = "NumDistSum_pop"+num2str(pop)
	Make/N=(numpnts(Core_RadiusWave),numpnts(Shell_RadiusWave))/O $VolumeDistSum, $NumDistSum
	Wave VolumeDistSum_w = $VolumeDistSum
	Wave NumDistSum_w = $NumDistSum
	
	
	MEF_CoreShell_AveVolumeDist(VolumeDistSum_w,Core_TempVolDistL,SHell_TempVolDistL,Core_RadiusWave,Shell_radiusWave,pop)
	
	//Calculate the Num/VolumeDistributions
	
	Core_VolumeDist = Core_TempDist
	Shell_VolumeDist = Shell_TempDist
	Core_NumDist = Core_VolumeDist / Core_TempVolDistL
	Shell_NumDist = Shell_VolumeDist / Shell_TempVolDistL

	//NumDistSum_w = TotalDist / VolumeDistSum_w
	
	//Normalize to the scale factor and the area under the distribution which is almost (note always 1)

	Variable ScaleCore,ScaleShell
	ScaleCore = Core_SZScale / AreaXY(Core_RadiusWave,Core_VolumeDist,-inf,inf)
	ScaleShell = Shell_SZScale / AreaXY(Shell_RadiusWave,Shell_VolumeDist,-inf,inf)
	
	Core_VolumeDist *= ScaleCore
	Core_NumDist *= ScaleCore
	Shell_VolumeDist *= ScaleShell
	SHell_NumDist *= ScaleShell
	
	//Finally add together the possible Distributions
	NumDistSum_w = TotalDist / VolumeDistSum_w
	
	//All distributions should be done now, Now we want to calculate the scattering with newfound volumes.


end

