#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

////////////////////////////////////////////////////////////////////////////////////////////
//
// Written by Devin Grabner, Washington State University, March 2024
//
//	If you have any further questions please contact Devin Grabner (devin.grabner@wsu.edu) or 
//		Brian Collins (brian.collins@wsu.edu) for more information.
//
// This panel is made to facilitate the generation of scattering contrast functions between phases
// within a system. It uses the complex indices of refraction for all the materials within the system
// and appropriately weights them by their mass fraction within a phase. Once the complex indices of
// refraction have been calulated for each phase, contrast functions for each phase pair are calulated
// and plotted.
//
// The initial parameters and each step are saved to the specificed material folder for later use if needed. 
//
////////////////////////////////////////////////////////////////////////////////////////////

Menu "RSoXS Tools"
	"Complex Contrast Calculator Panel", Contrast_Calc_BuildPanel()
End

Function Contrast_Calc_BuildPanel()
	NewDataFolder/O/S root:MaterialContrast
	
	Variable/G NumMat = 2 //Number of different materials
	Variable/G NumPhase = 3 //Number of phases chosen by user
	string/G SystemName = ""
	Make/T/O OCwaveType = {"Energy_Wave","Delta_Wave","Beta_Wave"}
	
	//Build Panel
	DoWindow/K Contrast_Calc_MainPanel
	NewPanel/N=$("Contrast_Calc_MainPanel")/K=1 /W=(50,50,1000,770) as "Complex Refractive Index Contrast Calculator"
	
	SetVariable SystemName, pos={5,5},size={300,0},fsize=12,title="Name of your System:",help={"This will be the name of the folder/graph in which your contrast data will be stored/displayed."}
	SetVariable SystemName, value = SystemName
	TitleBox SyntaxStatment, pos={310,5},size={250,0},fsize=12,frame=0,title="(Follow IGOR syntax for naming waves.)"
	
	SetVariable NumberOfMaterials,pos={5,35},size={210,0},fsize=12,title="Number of Different Materials:"
	SetVariable NumberOfMaterials,limits={2,6,1},proc=NumberUpdate,value = NumMat
	MatWaveUpdate()
	TitleBox MaxMaterials, pos={220,37},size={100,0},fsize=10,frame=0,title="(Max 6 Materials)"
	
	PopupMenu OldData, pos={600,5},size={100,0},fsize=12,mode=0,title="Load Prior System"
	PopupMenu OldData, proc=LoadOldData, value=DataFolderList("*",";")
	
	TitleBox SetParentFolder,pos={5,85},size={750,0},fsize=12,fstyle=1,frame=2,title="Put in the full wave path for the associated material optical constants!!"
	
	GroupBox OptConst,pos={5,120},fsize=16,size={930,110},title="Material Optical Constants"
	TitleBox MaterialDesc,pos={10,145},size={750,0},frame=0,fstyle=2,fsize=14,title="Material Name\nEnergy Wave\nDelta Wave\nBeta Wave"
	
	ListBox OCwaveMatrix,pos={114,145},size={810,80},fsize=12,frame=2
	ListBox OCwaveMatrix,listWave=root:MaterialContrast:OCwaveMatrix,selWave=root:MaterialContrast:MatSelWave
	
	SetVariable NumberOfPhases,pos={5,245},size={335,0},fsize=12,title="How many different phases are in your system?"
	SetVariable NumberOfPhases,limits={NumMat,126,1},proc=NumberUpdate,value = NumPhase
	CompUpdate()
	PhasePairUpdate()
	
	TitleBox FillTable,pos={5,275},fsize=12,frame=0,title="\t\t***Attention: Fill out the following tabel with the mass fraction of each material within the given phase. Enter it is as a decimal. (e.g. 45% -> 0.45)***"
	
	GroupBox PhaseComp,pos={5,300},fsize=16,size={930,205},title="Phase Composition"
	
	ListBox MaterialNames,pos={10,325},size={910,40},fsize=14,frame=2,listWave=root:MaterialContrast:matNames
	
	ListBox PercMassFrac,pos={10,375},size={910,120},fsize=14,frame=2
	ListBox PercMassFrac,listWave=root:MaterialContrast:CompositionMatrix,selWave=root:MaterialContrast:CompSelWave
	
	TitleBox PhasePairNote,pos={5,515},fsize=14,frame=0,title="What Phase_i_j pairs do you want to have plotted? (Shift + Click) to select multiple"
	
	GroupBox PhasePair,pos={5,540},fsize=16,size={930,88},title="Phase Contrast Pairs"
	
	ListBox PhasePairs,pos={10,565},size={910,55},fsize=14,frame=2,mode=8
	ListBox PhasePairs,listWave=root:MaterialContrast:PhasePairs,selWave=root:MaterialContrast:PhaseSelWave
	
	Button GenerateContrastFunc,pos={5,655},size={930,45},frame=3,fcolor=(8738,35723,20000),clickEventModifiers=4
	Button GenerateContrastFunc, proc=GenContrastFunc, title="Calculate Contrast Functions"
END

Function NumberUpdate(sv) : SetVariableControl
	STRUCT WMSetVariableAction &sv
	MatWaveUpdate()
	CompUpdate()
	PhasePairUpdate()
	return 0
End

Function LoadOldData(PU_Struct) : PopupMenuControl
	STRUCT WMPopupAction &PU_Struct
	
	If(PU_Struct.eventcode!=2)
		return 0
	elseif(stringmatch(PU_Struct.popStr,"_none_") == 1)
		abort "You do not have any prior data to load."
	endif
	
	SVAR System_Name_main = root:MaterialContrast:SystemName
	NVAR Num_Mat_main = root:MaterialContrast:NumMat
	NVAR Num_Phase_main = root:MaterialContrast:NumPhase
	WAVE/T OC_wavePathMatrix_main = root:MaterialContrast:OCwaveMatrix
	WAVE/T Composition_Matrix_main = root:MaterialContrast:CompositionMatrix
	
	SVAR System_Name_sub = root:MaterialContrast:$(PU_Struct.popStr):SystemName
	NVAR Num_Mat_sub = root:MaterialContrast:$(PU_Struct.popStr):NumMat
	NVAR Num_Phase_sub = root:MaterialContrast:$(PU_Struct.popStr):NumPhase
	WAVE/T OC_wavePathMatrix_sub = root:MaterialContrast:$(PU_Struct.popStr):OCwaveMatrix
	WAVE/T Composition_Matrix_sub = root:MaterialContrast:$(PU_Struct.popStr):CompositionMatrix
	
	setdataFolder root:MaterialContrast
	System_Name_main = System_Name_sub
	Num_Mat_main = Num_Mat_sub
	Num_Phase_main = Num_Phase_sub
	
	MatWaveUpdate()
	CompUpdate()
	PhasePairUpdate()
	
	OC_wavePathMatrix_main = OC_wavePathMatrix_sub
	Composition_Matrix_main = Composition_Matrix_sub
	
	abort //This stops all other subroutines from running that would get triggered by the data being changed
End

Function MatWaveUpdate()
	setdataFolder root:MaterialContrast
	NVAR NumMat = root:MaterialContrast:NumMat
	
	Make/T/O/N=(4,NumMat) OCwaveMatrix
	Make/O/N=(4,NumMat) MatSelWave = 2
	Make/T/O/N=(1,NumMat+1) matNames
	matNames[0][0] = "Material Name"
	Variable i
	for(i=0;i<NumMat;i++)
		matNames[0][i+1]=OCwaveMatrix[0][i]
	endfor
End

Function CompUpdate()
	setdataFolder root:MaterialContrast
	NVAR NumPhase = root:MaterialContrast:NumPhase
	NVAR NumMat = root:MaterialContrast:NumMat
	WAVE/T matNames = root:MaterialContrast:matNames
	
	if(waveExists(root:MaterialContrast:CompositionMatrix) ==1)
		DeletePoints 0,1000, root:MaterialContrast:CompositionMatrix
	endif
	Make/O/T/N=(NumPhase,NumMat+1) CompositionMatrix
	Variable i
	for(i=0;i<NumPhase;i++)
		CompositionMatrix[i][0] = "Phase_" + num2str(i+1)
	endfor
	
	Make/O/N=(NumPhase,NumMat+1) CompSelWave = 2
	CompSelWave[*][0] = 8
End

Function PhasePairUpdate()
	WAVE/T Composition_Matrix = root:MaterialContrast:CompositionMatrix
	NVAR NumPhase = root:MaterialContrast:NumPhase
	variable PhaseComb //Number of Phase combinations
	variable i, j, tempVal
	
	PhaseComb = binomial(NumPhase, 2)
	
	Make/O/T/N=(PhaseComb) PhasePairs
	
	variable k=1
	variable s=0
	for(i=0; i<NumPhase; i++)		
		for(j=k; j<NumPhase; j++)		
			PhasePairs[s] = Composition_Matrix[i][0] + "_" + Composition_Matrix[j][0] //Contrast between Phase i and Phase j
			s=s+1
		endfor
		k=k+1
	endfor
	
	tempVal = ceil(PhaseComb/6)
	
	Redimension/N=(tempVal,6) PhasePairs
	Make/O/N=(tempVal,6) PhaseSelWave = 4
End

Function GenContrastFunc(ctrlName) : ButtonControl
	String ctrlName
	SVAR System_Name = root:MaterialContrast:SystemName
	NVAR Num_Mat = root:MaterialContrast:NumMat
	NVAR Num_Phase = root:MaterialContrast:NumPhase
	
	if(stringmatch(System_Name,"") == 1)
		abort "Please name your system before continuing."
	endif
	
	WAVE/T OC_wavePathMatrix = root:MaterialContrast:OCwaveMatrix
	WAVE/T Composition_Matrix = root:MaterialContrast:CompositionMatrix
	
	Variable i,j,tempVar1,tempVar2
	string tempName1,tempName2,tempName3,tempName4
	string windowName
	
	string PlotWindowName	= "Contrast_Functions_" + System_Name
	
	KillWindow/Z $PlotWindowName
	
	NewDataFolder/O root:MaterialContrast:$(System_Name)
	NewDataFolder/O root:MaterialContrast:$(System_Name):ComplexOCs
	
	SetDataFolder root:MaterialContrast:$(System_Name)
	
	//Build Energy Wave for the Contrast Function
	Make/Free/N=(Num_Mat) EnRangeMin //List of the first energy for each component 
	Make/Free/N=(Num_Mat) EnRangeMax //List of the highest energy for each component
	for(i=0; i<Num_Mat; i++)
		tempName1 = OC_wavePathMatrix[1][i]
		WAVE tempEnergy = $(tempName1)
		EnRangeMin[i] = tempEnergy[0]
		EnRangeMax[i] = waveMax(tempEnergy)		
	endfor
	
	//These will provide the range to which the contrast function can be make. All component NEXAFS will fall within this range.
	variable EnMin = ceil(waveMax(EnRangeMin))
	variable EnMax = floor(waveMin(EnRangeMax))
	variable points = (EnMax-EnMin)*10
	Make/O/N=(points) Contrast_En =  EnMin + x/10
	
	SetDataFolder root:MaterialContrast:$(System_Name):ComplexOCs
	for(i=0; i<Num_Mat; i++)
		tempName1 = "C_" + OC_wavePathMatrix[0][i]
		make/C/O/N=(points) $tempName1 = cmplx((1-(interp(Contrast_En[x], $(OC_wavePathMatrix[1][i]), $(OC_wavePathMatrix[2][i])))),(interp(Contrast_En[x], $(OC_wavePathMatrix[1][i]), $(OC_wavePathMatrix[3][i]))))
	endfor
	
	for(i=0; i<Num_Phase; i++)
		tempName1 = Composition_Matrix[i][0]
		
		Make/C/O/N=(points) $tempName1 = 0 // Set the whole wave to zero
		Wave/C nP = $tempName1
		  
		for(j=0; j<Num_Mat; j++)
			tempName2 = "C_" + OC_wavePathMatrix[0][j]
			tempVar1 = str2num(Composition_Matrix[i][j+1])
			WAVE/C COC = root:MaterialContrast:$(System_Name):ComplexOCs:$tempName2 //Get the complex index for a material (Complex Optical Constants)
			nP = nP + tempVar1*COC
		endfor
	endfor
	SetDataFolder root:MaterialContrast:$(System_Name)
	
	WAVE/T PhasePairs = root:MaterialContrast:PhasePairs
	WAVE PhaseSelWave = root:MaterialContrast:PhaseSelWave
	
	variable s, k=1
	for(i=0; i<Num_Phase; i++)
		tempName1 = Composition_Matrix[i][0]
		WAVE/C np_i = root:MaterialContrast:$(System_Name):ComplexOCs:$tempName1
		
		for(j=k; j<Num_Phase; j++)
			tempName2 = Composition_Matrix[j][0]
			WAVE/C np_j = root:MaterialContrast:$(System_Name):ComplexOCs:$tempName2
			
			tempName3 = "CmplxContrast_" + Composition_Matrix[i][0] + "_" + Composition_Matrix[j][0] //Complex contrast between Phase i and Phase j
			make/C/O/N=(points) $tempName3 = np_i - np_j
			WAVE/C CCP = $tempName3
			
			tempName4 = Composition_Matrix[i][0] + "_" + Composition_Matrix[j][0] //Contrast between Phase i and Phase j
			make/O/N=(points) $tempName4 = magsqr(CCP)*Contrast_En^4
			
			SetDataFolder root:MaterialContrast
			tempName1 = Composition_Matrix[i][0] + "_" + Composition_Matrix[j][0]
			FindValue/TEXT=tempName1 PhasePairs
			tempVar1 = PhaseSelWave[V_row][V_col]
			SetDataFolder root:MaterialContrast:$(System_Name)
			
			if(tempVar1 == 5)
				PlotFunc($tempName4,Contrast_En,PlotWindowName)
			endif
		endfor
		k=k+1
	endfor
	
	Duplicate/O OC_wavePathMatrix,OCwaveMatrix
	Duplicate/O Composition_Matrix,CompositionMatrix
	Variable/G	NumMat = Num_Mat
	Variable/G NumPhase = Num_Phase
	String/G SystemName = System_Name
	
	SetDataFolder root:MaterialContrast
	CommonColorsButtonProc2()
End

//Displays the Contrast Functions in a graph
FUNCTION PlotFunc(yWave, xWave, windowSuffix)

	wave yWave
	wave xWave
	string windowSuffix
	string windowName	= windowSuffix
	string leftName		= "Contrast * Energy^4"
	string bottomName	= "Energy (eV)"
	
	DoWindow/F $windowName
	if(V_flag == 0)
		Display/W=(500,50,1000,500) yWave vs xWave
		ModifyGraph grid=2, tick=2, mirror=1, minor=1, standoff=0, msize=1.5e
		Label left leftName
		Label bottom bottomName
		DoWindow/C/T $windowName, windowName
		Legend/C/N=ContrastPlots/A=RT
		SetAxis bottom 270,350
	else
		CheckDisplayed yWave
		if(V_flag == 0)
			AppendToGraph yWave vs xWave
		endif
	endif
	SetAxis/A
	ModifyGraph lsize=2
	ModifyGraph log(left)=1
END

///////////////////Below Functions are Copied from KBColorizeTraces.ipf
Function CommonColorsButtonProc2()
	
	String graphName = KBActiveGraph()
	if (strlen(graphName) == 0)
		return -1
	endif
	
	Variable numTraces = KBTracesInGraph(graphName)

	if (numTraces <= 0)
		return -1
	endif

	Variable red, green, blue
	Variable i, index
	for(i=0; i<numTraces; i+=1)
		index = mod(i, 10)				// Wrap after 10 traces.
		switch(index)
			case 0:
				red = 0; green = 0; blue = 0;
				break

			case 1:
				red = 65535; green = 16385; blue = 16385;
				break
				
			case 2:
				red = 2; green = 39321; blue = 1;
				break
				
			case 3:
				red = 0; green = 0; blue = 65535;
				break
				
			case 4:
				red = 39321; green = 1; blue = 31457;
				break
				
			case 5:
				red = 48059; green = 48059; blue = 48059;
				break
				
			case 6:
				red = 65535; green = 32768; blue = 32768;
				break
				
			case 7:
				red = 0; green = 65535; blue = 0;
				break
				
			case 8:
				red = 16385; green = 65535; blue = 65535;
				break
				
			case 9:
				red = 65535; green = 32768; blue = 58981;
				break
		endswitch
		ModifyGraph/W=$graphName rgb[i]=(red, green, blue)
	endfor
End

static Function/S KBActiveGraph()
	// Active graphs can be embedded in panels and layouts
	String topWin= WinName(0,1+4+64)	// Graphs+Layouts+Panels, usually KBColorizePanel
	String nextWin= WinName(1,1+4+64)// Graphs+Layouts+Panels
	String graphName=WinName(0,1)
	if( CmpStr(topWin,"KBColorizePanel") == 0 )
		topWin= nextWin
	endif
	if( strlen(topWin) )
		GetWindow $topWin activeSW
		if( WinType(S_Value) == 1 )
			graphName= S_value
		endif
	endif
	return graphName
End

// Find the number of traces on the top graph
static Function KBTracesInGraph(win)	// "" for top graph
	String win
	
	if( strlen(win) == 0 )
		win= KBActiveGraph()
		if( strlen(win) == 0 )
			return 0
		endif
	endif
	return ItemsInList(TraceNameList(win,";",3))
End
