#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions and panel support for improved RSoXS workflow.
//All panel support originally designed by Elliot Gann (now at NSLS-II) utilizing the NIKA Data reduction package by Jan Ilavsky.

//The following panel was designed following the December 2017 NIKA Beta that includes RSoXS support.
//All primary data reduction functions are written and maintained by Jan Ilavsky http://usaxs.xray.aps.anl.gov/staff/ilavsky/nika.html
//Jan Ilavsky, "Nika - software for 2D data reduction", J. Appl. Cryst. (2012), vol. 45, pp. 324-328. DOI:10.1107/S0021889812004037
//

//Adapted to IGOR 7.05 by Thomas Ferron (Washington State University) Contact: thomas.ferron@wsu.edu


//Version Updates
//V 1.0 -
//Fixed a bug that would not update the popupmenu upon changing fits file
//Fixed a bug that would invert the image in the display window without inverting the sectors as well.
	//No change to performance, simply visually confusing.
//Install file created for more user friendly experience.
//Install file subsequentially fixed for an even more friendly experience.
//Fixed a bug that would override the loaded I0 upon importing a single image
//Allows loading even if no matching dark exists
//Fixed a bug that allows loaded darks to contain decimal points

//V 0.2 - 
//Fixed an issue that would have the main panel pop up every time the macro initializes...Now it will only pop up upon compiling.
//	//If no editing occurs then it will only pop up if you rewrite some code and compile.
//Fixed a bug that would open the mask panel/Beam centering panel with nothing selected. (otherwise it would Open with no image ready)
//Fixed a bug that would not let you change folders with less data than previously
//Fixed a bug that would error if you click on any whitespace in the listboxes (It would error otherwise)
//Fixed 'Open to Beam Centering' to work even though 'Open for Mask' had not been clicked first.
//Fixed Browse/Update button to properly change paths and not just add more options to the same list.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Initialization function for the panel. 
Function NRS_1101Panel_Init()

	String CurrentFolder = GetDataFolder(1)
	
	
	//INITIALIZE NIKA 2D SAS Macros. Menu will be built but no main panel construction (should not be needed ever)
	LoadNika2DSASMacros() //This function is currently built for Igor ver. 7.05, THIS WILL NOT LOAD WITH A LESSER VERSION.
	//Some later functions require the panel to 'exist' so we will make it.....then hide it.

	
	//Variables and Strings to be used later on
	
	//MOVED THIS TO ANY BUTTON ON THE MAIN PANEL
	//Create the main NIKA panel to allow access to buttons and panel support. WIll immediately hide it so it will not get in the way.	
//	NI1A_Convert2Dto1DMainPanel()
//	DoWindow/Hide=1 NI1A_Convert2Dto1DPanel
	
	//Initialize NIKA RSoXS Parameters  
//	NI1_RSoXSConfigureNika() //Developed by Jan Ilavsky to set basic parameters that should be used by all users
	//Variables and Strings to be used later on
	NVAR/Z UseRSoXSCodeModifications=root:Packages:Nika_RSoXS:UseRSoXSCodeModifications
	UseRSoXSCodeModifications = 1
//	CheckBox UseRSoXSCodeModifications value=1, win=NI1_RSoXSMainPanel

//	Initialize additional Basic Parameters

//	These can be edited for your preferences in data reduction...Just add or remove preferences in the following function (at the end of the procedure file)
//	Or, rewrite your own function and put it here.
//	NRS_LabPreferences()
	
	//Now we construct the necessary strings and variables for the custom user interface for faster and more efficient workflow
	//Make folders
	//Folder names allowing easier references later
	String NRS_Fldr = "root:Packages:NIKA_RSoXS"
	String Panel_Fldr = NRS_Fldr + ":NRS_Vars"
	NewDataFolder/o/s root:Packages
	NewDataFolder/o/s $NRS_Fldr
	NewDataFolder/o/s $Panel_Fldr
	
	//Create Global Strings...All of these strings and panels will be used to store variables and workflow information.
	//A basic description is given next to each one (If I ended up figuring out what they were suppose to do)
	
	String/G pathtodata //Path to the header .txt filed used to organize the .fits files
	String/G pathtoFits //Path to the .fits files you want to load...This may be different than the headers...(TBD)
	String/G pathname = "Path_1101panel" //Name of the Igor symbolic path that governs 'Pathtodata' this specific name will not change
	String/G pathname_Fits = "Path_1101Fits" //optional parameter if .fits files are different than header files (Not yet developed)
	//String/G loadeddatadir
	String/G header //This is a temp string for the wave header
	string/G headerinfo //This is another temp string? for the wave header? who knows.
	String/G imagekeys //List of keywords searchable from the given .fits header
	string/G imagetime //The time the .fits file was created
	String/G imagevalue //Seems to be an old variable that is no longer used
//	String/G filebasename 
	String/G imagekeypick = "Exposure" //Specific keyword to be displayed in the popup menu governing all header key words
//	String/G test123123
	//Create GLobal Variables
	Variable/G ckautoshow=1
	Variable/G logimage = 1 //Variable to determine if the preview window should be log (default)
//	Variable/G NormI0
	Variable/G FitsTopSelected =0 //The top selected fits file for loading
	Variable/G SeriesSelected = 0//The Highlighted series in listbox 1
	Variable/G n2save
	Variable/G showmask
	Variable/G endplotrun
	Variable/G plotrunsp
	Variable/G plotrunst
	Variable/G plotrunend
	Variable/G NormalizeData=1
	
	Variable/G SetPreferences =0 //Quick monitor to see if preferences were ever ran. Mostly if the main panel is opened before any buttons pressed.

	
	Variable/G maxslider
	Variable/G minslider
	Variable/G HotPixel = 0
	
	Variable/G writedata = 0
	
	//POSSIBLY CHANGING IN THE FUTURE VARIABLES
	Variable/G xpix_dim = 1024
	variable/G ypix_dim = 1024
	//Make Text Waves
	make/T/O/N=0 Basenames
	make/T/O/N=0 Motors
	make/T/O/N=0 Scans //File display for the first listbox window...Updated upon browsing a file or updating.
	Make/T/O/N=0 FitsDirectory
	Make/T/O/N=0 Series_Display
	Make/T/O/N=0 ScanID_Display
	make/T/O/N=0 Scantypes
	make/T/O/N=0 Times
	make/T/O/N=0 Files
	make/T/O/N=0 Motor1Pos 
	make/T/O/N=0 Motor2Pos
	make/T/O/N=0 I0s
	make/T/O/N=0 M1Wave
	Make/T/O/N=0 M2Wave
	Make/T/O/N=0 angwave
	Make/T/O/N=0 Files
	Make/T/O/N=0 Filedisc //File display for the second listbox window...Updated upon clicking an option in box 1
	Make/O/N=0 Filesel //Selected files to import in the fits listbox

	//Make Variable waves / Display waves
	Make/O/N=0 Frames
	Make/O/N=0 nmotors
	Make/O/N=0 ScanNum
	Make/O/N=0 FSindex
	Make/D/O/N=(xpix_dim,ypix_dim) Data = 0 //Data that is imported.....The displayed may be edited...this is no edits,that will display the currently selected folder
	Make/D/O/N=(xpix_dim,ypix_dim) Data_disp = 0 //Data that will display the currently selected folder
	Make/D/O/N=(xpix_dim,ypix_dim) Saturate_Mask = NaN //Mask that will point out saturated data points...
	Make/B/U/O/N=(xpix_dim,ypix_dim) Data_mask = NaN

	Make/D/O/N=(xpix_dim,ypix_dim) Datax
	Make/D/O/N=(xpix_dim,ypix_dim) Datay
	Make/B/U/O/N=(xpix_dim,ypix_dim) Avg_Mask = 1

	//Assign preliminary wave results
	scans=""
	motor1pos=""
	motor2pos=""
	I0s=""
	Basenames = ""
	ImageValue = ""
	
	
	//initialize panel
	//MakePanel Here
	
	NRS_BuildMainPanel()
	SetDataFolder $CurrentFolder
	
	
	
	print "//////////////////////////////////////////////////////////////////////"
	Print "///                                                                ///"
	print "///               11.0.1.2 Image Analysis Panel                    ///"
	Print "///                       Version 1.0                              ///"
	Print "///  Download @ https://labs.wsu.edu/carbon/xray-analysis-tools/   ///"
	Print "///                                                                ///"
	print "//////////////////////////////////////////////////////////////////////"

end

///////////////////////////////////////////////////////////////////////////////////////////////////////
///These are RSoXS Data reduction preferences for the Collins group at WSU/////////////////////////////
///If your group uses other preferences they can be added here/////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

//Setup our labs preferences for the main panel of NIKA. These simply check boxes and set settings. Nothing more
Function NRS_LabPreferences()

	//Global Variables from NIKA that we like to use....these can be adjusted but for now we are just setting them
	NVAR DezingerData = root:Packages:Convert2Dto1D:DezingerCCDData
	NVAR UseSectors = root:Packages:Convert2Dto1D:UseSectors
	NVAR DoCircularAverage = root:Packages:Convert2Dto1D:DoCircularAverage
	NVAR DoSectorAverages = root:Packages:Convert2Dto1D:DoSectorAverages
	NVAR NumberOfSectors = root:Packages:Convert2Dto1D:NumberOfSectors
	NVAR AngleBetweenSectors = root:Packages:Convert2Dto1D:SectorsStepInAngle
	NVAR Display1DGraph = root:Packages:Convert2Dto1D:DisplayDataAfterProcessing
	NVAR StoreDatainIgor = root:Packages:Convert2Dto1D:StoreDataInIgor
	NVAR OverwriteData = root:Packages:Convert2Dto1D:OverwriteDataIfExists
	NVAR UseCustomNameFnct = root:Packages:Convert2Dto1D:UseSampleNameFnct
	NVAR UseInputDataName = root:Packages:Convert2Dto1D:Use2DDataName
	NVAR UseMask = root:Packages:Convert2Dto1D:UseMask
	NVAR Process_Individually = root:packages:Convert2Dto1D:Process_Individually
	NVAR NumPoints = root:Packages:Convert2Dto1D:QvectorNumberPoints
	
	NVAR InvertImage = root:Packages:Convert2Dto1D:InvertImages
		
//	NVAR UsePolarizationCor = root:Packages:Convert2Dto1D:DoPolarizationCorrection
//	NVAR PolarizationCor = root:Packages:Convert2Dto1D:Use2DPolarizationCor

	
	SVAR SampleNameFnct = root:Packages:Convert2Dto1D:SampleNameFnct

	//Just some basic data processing that is mostly up to the user
	DezingerData = 1
	UseSectors = 1
	DoCircularAverage = 1
	NumPoints = 150
//	UsePolarizationCor = 1
//	PolarizationCor = 1
	UseMask = 1
	Process_Individually=1
	InvertImage = 1

	//Setup Sectors
	DoSectorAverages = 1
	NumberOfSectors = 8
	AngleBetweenSectors = 45
	
	//Display Data with no Q scaling to verify import looks good before more in depth processing
	Display1DGraph = 1
	StoreDatainIgor = 1
	OverwriteData = 1
	//Checkbox is also required to prepare to use the custom filename
	UseCustomNameFnct = 1
	UseInputDataName = 0
	
	SampleNameFnct = "NRS_CustomFileName"
end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Im testing where I want to put this...to make sure it happens...but also so the main panel does not pop up...ever//
//This function initializes all the NIKA panels that you would want in order to run the any NIKA function////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function NRS_Initialize_NIKAMainPanel()

	NI1A_Convert2Dto1DMainPanel()
	NI1_RSoXSInitialize()
	
	NVAR/Z UseRSoXSCodeModifications=root:Packages:Nika_RSoXS:UseRSoXSCodeModifications
	UseRSoXSCodeModifications = 1
		
	NI1_RSoXSConfigureNika()
	
	//Need to change this if you want to change your lab preferences
	NRS_LabPreferences()

	DoWindow/Hide=1 NI1A_Convert2Dto1DPanel
	
	
end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function used to create the save name in the manor that we prefer. Format is based on scan ID, energy and polarization//
//Support for more advanced naming support will come at a later version.//////////////////////////////////////////////////
//ScanID_Energy_Pol_Sector////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function/S NRS_CustomFileName(My2DImage, OriginalDataFileName)
	wave My2DImage
	string OriginalDataFileName
	//do something to create a new name, here I will trunkate the name to 15 characters
	
	String wvnote = note(My2DImage)
	String Energy = Stringbykey("BEAMLINE ENERGY", wvnote, "=",";")
	String Pol = Stringbykey("EPU POLARIZATION", wvnote, "=",";")
	String TempScanName, TempScanID, TempScanNum
	Splitstring /E="(.*)_?(\\d{4})-(\\d{5})$" OriginalDataFilename,TempScanName,TempScanID,tempScanNum //grabs the scan name and the last 4 digits of the ID
	
	//New from IGOR 7...Removes specific endings from polarization and energy terms...to be automated later for any keyword 
	Pol = Removeending(Pol,". [Counts]")
	Energy = RemoveEnding(Energy," [eV]")
	//Round Energy to nearest decimal...a bit sloppy and can be done better some other time
	Variable En = str2num(Energy)
	En *= 10; En = Round(En); En /= 10

	//Create the final name string
	string tempName= TempScanID + "_" + num2str(En) + "_" + TempScanNum + "Dark1"
	//string tempName= TempScanID + "_" + num2str(En) + "_" + TempScanNum
	//string tempName= TempScanID + "_" + num2str(En) + "_" + Pol
	//or here look up "SampleTitle" in wave note
	//string Wvnote = note(My2DImage)
	//tempName = stringByKey("SampleTitle", Wvnote)

	return tempName
end

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//The following functions are all panel controls. All button commands run from a single function for simplicity////////////////
//More in-depth functions may be called and will be outlined later/////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//Button functions
Function NRS_ButtonProc(ba) : ButtonControl
	Struct WMButtonAction & ba
	
	If(ba.eventcode!=2)
		return 0
	endif
	
	String CurrentFolder = GetDataFolder(1)
	//Gather parameters and variables that will be used in buttons
	SVAR pathtodata = root:packages:NIKA_RSoXS:NRS_Vars:pathtodata
	SVAR pathname = root:packages:NIKA_RSoXS:NRS_Vars:pathname
	Wave/T Filesel = root:Packages:Nika_RSoXS:NRS_Vars:Filesel
	
	SVAR CCDFileExtension=root:Packages:Convert2Dto1D:CCDFileExtension
	
	NVAR SetPreferences = root:Packages:NIKA_RSoXS:NRS_Vars:SetPreferences

	
	//Gather NIKA variables that we will append values too as we go ///NOTE: Some are called later after they are created through NIKAs creation macros
	Wave/T	ListOf2DSampleData=root:Packages:Convert2Dto1D:ListOf2DSampleData
	Wave	ListOf2DSampleDataNumbers=root:Packages:Convert2Dto1D:ListOf2DSampleDataNumbers
	Wave/T	ListOf2DEmptyData=root:Packages:Convert2Dto1D:ListOf2DEmptyData
	Wave/T	ListOffilenames=root:Packages:Convert2Dto1D:ListOfCCDDataInCCDPath
	Wave/T	ListofBCFilenames=root:Packages:Convert2Dto1D:ListOfCCDDataInBmCntrPath
	Wave/Z M_ROIMask=root:Packages:Convert2Dto1D:M_ROIMask

	//Strings for later
	String filenamelist
	String filename
	
	//Variables for later
	Variable i

	//Buttons commands
	If(StringMatch(ba.win,"NRS_MainPanel")) //Check to see if you are on the main panel (most likely wont be more panels)
	
	//Quick check to see if the main panel exists...its required to be created in order to import data. Although it does not need to be visible.
	
	Dowindow/HIDE=? NI1A_Convert2Dto1DPanel
	if(V_Flag == 0)
		NRS_Initialize_NIKAMainPanel()
		SetPreferences = 1
	endif
	
	//Second check to see if the main preferences were ever established...quick check to verify the main panel was not opened first
	if(SetPreferences == 0 )
		NI1_RSoXSConfigureNika()
		NRS_LabPreferences()
		SetPreferences = 1
	endif
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Prompt for the user to find the folder used for importing data. Setup the leftmost listbox and assign paths to NIKA//
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_Browse"))
		
			NewPath/O Path_1101panel //Create path to the header files
			if(V_Flag != 0) //Quick check to see if the user cancelled searching for a folder
				print "User cancelled action, no data folder selected"
				return 0
			endif
			GetFileFolderInfo/D/Q/P=$pathname
			pathtodata = s_path
			//Set the NIKA paths for internal functions to be the same as ours
			NRS_convertpathtonika(main=1)
			NRS_LoadSeriesList(pathname,pathname)
			
		endif
		//Refreshes listbox if more data has been taken
		if(!cmpstr(ba.ctrlname,"NRS_UpdateBrowse"))
			PathInfo Path_1101panel
			pathtodata = s_path
			//Set the NIKA paths for internal functions to be the same as ours
			NRS_convertpathtonika(main=1)
			NRS_LoadSeriesList(pathname,pathname)
		endif
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//This button will convert the selected data into NIKA identical to the 'Process Image(s)' button//
		///////////////////////////////////////////////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_ConvertData"))
		
			filenamelist = getfilename(selwave=filesel)
			
			NRS_convertpathtonika(main=1)
			if(itemsinlist(filenamelist,";") == 0)
				print "Please select a scan file to import"
				return 0
			endif	
			ListOf2DSampleDataNumbers = 0
			for(i=0;i<itemsinlist(filenamelist);i+=1)
				filename = stringfromlist(i,filenamelist)
				FindValue /TEXT=filename /TXOP=6 /Z ListOf2DSampleData
				if(v_value>=0)
					ListOf2DSampleDataNumbers[v_value] = 1
				endif
			endfor
			print "zoop"
			if(itemsinlist(filenamelist,";") == 1)
				NRS_CheckforDark()
				NRS_CheckforI0()
			endif
			
			
			
			NI1A_CheckParametersForConv()
			//set selections for using RAW/Converted data...
			NVAR LineProfileUseRAW=root:Packages:Convert2Dto1D:LineProfileUseRAW
			NVAR LineProfileUseCorrData=root:Packages:Convert2Dto1D:LineProfileUseCorrData
			NVAR SectorsUseRAWData=root:Packages:Convert2Dto1D:SectorsUseRAWData
			NVAR SectorsUseCorrData=root:Packages:Convert2Dto1D:SectorsUseCorrData
			LineProfileUseRAW=0
			LineProfileUseCorrData=1
			SectorsUseRAWData=0
			SectorsUseCorrData=1
			//selection done
			NI1A_LoadManyDataSetsForConv()
		//	DoWindow/K CCDImageToConvertFig
		endif
		//////////////////////////////////////////////////////////////////////////////////////////
		//This will load the selected images as darks and saved using NIKAs new import funciton.//
		//This is a reduced function of 'NI1A_LoadEmptyOrDark().//////////////////////////////////
		//Originally written by Jan Ilavsky///////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_LoadasDark"))
			setDataFolder root:Packages:Convert2Dto1D
			filenamelist = getfilename(selwave=filesel)
			
			if(itemsinlist(filenamelist,";") == 0)
				print "Please select a scan file to load as dark"
				return 0
			endif	

			NRS_convertpathtonika(dark=1)
			for(i=0;i<itemsinlist(filenamelist);i+=1)
				filename = stringfromlist(i,filenamelist)		
				NI1A_UniversalLoader("Convert2Dto1DEmptyDarkPath",filename,"fits","DarkFieldData")
				wave NewCCDData = $("DarkFieldData")
				//this is modification needed for RSoXS data processing
				NVAR/Z UseRSoXSCodeModifications = root:Packages:Nika_RSoXS:UseRSoXSCodeModifications
				if(NVAR_Exists(UseRSoXSCodeModifications))
					Execute("NI1_RSoXSCopyDarkOnImport()")		//this guarrantees compile if function not present. 
				endif	
			endfor
			print "\r"
			print "Darks have finished loading"
		endif
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//This will open the mask panel with the selected image (first if multiple are selected) and be prepared to make a mask//
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_OpenMask"))
			filenamelist = getfilename(selwave=filesel)
			
			//Quick check to make sure you have clicked on a .fits file.		
			if(itemsinlist(filenamelist,";") == 0)
				print "Please select a scan file to open for mask"
				return 0
			endif			
			
			NRS_convertpathtonika(mask=1)
			doupdate
			FindValue /TEXT=(Stringfromlist(0,filenamelist,";")) /TXOP=6 /Z ListOffilenames
			if(V_Value >= 0)
				listbox CCDDataSelection win=NI1M_ImageROIPanel, selrow=v_value
				doupdate
				NI1M_MaskCreateImage() 
			endif
			popupmenu CCDFileExtension win=NI1M_ImageROIPanel, popmatch=".hdf"
			CCDFileExtension = ".hdf"
			NI1M_UpdateMaskListBox()
		endif
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//This will open the Beam centering panel with the selected image (first if multiple are selected) and be prepared to center//
		//ALso contains some preferences for our lab to make the centering a bit easier///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		if(!cmpstr(ba.ctrlname,"NRS_BeamCenter"))
			filenamelist = getfilename(selwave=filesel)

			//Quick check to make sure you have clicked on a .fits file.
			if(itemsinlist(filenamelist,";") == 0)
				print "Please select a scan file to open for beam centering"
				return 0
			endif
			
			NRS_ConvertpathtoNIKA(Beamcenter=1)
			doupdate
			Wave/T	ListofBCFilenames=root:Packages:Convert2Dto1D:ListOfCCDDataInBmCntrPath //Call a newly formed wave after the creation of the beam center utilities

			FindValue /TEXT=(Stringfromlist(0,filenamelist,";")) /TXOP=6 /Z ListofBCFilenames
			if(v_Value >= 0)
				listbox CCDDataselection win=NI1_CreateBMCntrFieldPanel, selrow=V_Value
				doupdate
				NI1BC_BMCntrCreateImage()
				NRS_BeamCenterPref()
			endif		
		endif
			
		//////////////////////////////////////////////
		//Opens the new NIKA panel to load RSoXS I0s//
		//Function built by Jan Ilavsky///////////////
		//////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_LoadI0"))
			NI1_RSoXSCreateGUI()
			SVAR DefaultI0 = root:Packages:NIKA_RSoXS:I0DataToLoad
			SVAR DefaultPD = root:Packages:NIKA_RSoXS:PhotoDiodeDataToLoad
			DefaultI0 = "AI_3_IZero"
			DefaultPD = "Photodiode"
			//Popupmenu I0DataToLoad, win=NI1_RSoXSMainPanel, 
			
			

		endif
		
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//This simply opens the advanced options panel for more NIKA options. You are able to change sectors and beam center and such//
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"NRS_Advanced"))
			NRS_AdvancedOptions()
		endif		
	endif
	
	if(stringmatch(ba.win,"NRS_MainPanel#NRS_AdvancedOptions"))
				
		/////////////////////////////////////////////////////////
		//Loads a new mask from the listbox next to this button//
		/////////////////////////////////////////////////////////
		if(!cmpstr(ba.ctrlname,"LoadMask"))
			NRS_LoadMask()
			DoWindow/Hide=1 NI1A_Convert2Dto1DPanel

		endif
		/////////////////////////////////////////
		//Selects an appropriate path to a mask//
		/////////////////////////////////////////
		if(!cmpstr(ba.ctrlname, "MaskSelectPath"))
			NewPath/C/O/M="Select path to your data" Convert2Dto1DMaskPath	
			NI1A_UpdateMainMaskListBox()	
			NI1M_UpdateMaskListBox()
			DoWindow/Hide=1 NI1A_Convert2Dto1DPanel
		endif
	endif
	
		if(!cmpstr(ba.ctrlname, "DownloadLink"))
			BrowseURL/Z "https://labs.wsu.edu/carbon/xray-analysis-tools/"
		endif
	
	
	SetDataFolder $CurrentFolder
end

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Main panel check box control used to alter the image display to include beam center sectors etc./////
///////////////////////////////////////////////////////////////////////////////////////////////////////

Function NRS_CheckBoxProc(cba) : CheckBoxControl

	STRUCT WMCHeckBoxAction &cba 
	//Grab variables
	//NVAR LogPreviewImage = root:packages:NIKA_RSoXS:NRS_Vars:LogPreviewImage
	wave Data_Disp = root:Packages:NIKA_RSoXS:NRS_Vars:Data_Disp
	wave Data = root:Packages:NIKA_RSoXS:NRS_Vars:Data
	
	NVAR UseSectors = root:Packages:Convert2Dto1D:UseSectors
	Wave/Z M_ROIMask=root:Packages:Convert2Dto1D:M_ROIMask

	If(cba.eventcode!=2)
		return 0
	endif
	String CurrentFolder = GetDataFolder(1)
	
	
	if(StringMatch(cba.win,"NRS_MainPanel"))
	
		///Sets the image to a log scale
		if(!cmpstr(cba.ctrlname, "NRS_LogPreview"))
			Variable Checked = cba.checked
			if(checked)
				data_disp = log(data)
			else
				data_disp = data	
			endif
			//Update the sliders now that the color scales have changed
			NRS_ResetSlider()
		endif
		//Sets up the sectors to view
		if(!cmpstr(cba.ctrlname,"NRS_SectorPreview"))
			checked = cba.checked
			if(checked && UseSectors==1)
				NRS_DrawSectorsIn2DGraph()
			else
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
			endif
		endif
		
		//Displays the beam center on the graph		
		if(!cmpstr(cba.ctrlname,"NRS_BeamCenterPreview"))
			checked = cba.checked
			if(checked)
				NRS_DrawCenterIn2DGraph()
			else
	    		setDrawLayer/K/W=NRS_MainPanel#disp Overlay
			endif
		endif
		//overlays the current mask on the image
		if(!cmpstr(cba.ctrlname,"NRS_MaskPreview"))
			checked = cba.checked
			
			//Quick check to see if any masks exist
			if(!waveexists(M_ROIMask))
				print "No Mask created. Please create or select a mask and try again"
				checkbox NRS_MaskPreview value=0
				return 0
			endif
			
			
			if(checked)
				RemoveImage /z /W=NRS_MainPanel#Disp M_ROIMask
				APPENDIMAGE/W=NRS_MainPanel#Disp M_ROIMask
				ModifyImage/W=NRS_MainPanel#Disp M_ROIMask ctab= {0,.9,Grays,0}
				ModifyImage/W=NRS_MainPanel#Disp M_ROIMask minRGB=0,maxRGB=NaN
			else
				RemoveImage /W=NRS_MainPanel#Disp M_ROIMask
//				RemoveImage /W=DataReduction#G0 data_mask
			endif
		endif
		
		
		
	endif
	
	SetDataFolder $CurrentFolder
end


/////////////////////////////////////////////////////////////////////////////////
//List box support to display the header files along with the appropriate .fits//
/////////////////////////////////////////////////////////////////////////////////

Function NRS_ListBoxProc(lba) : ListBoxControl
	STRUCT WMListboxaction & lba
	
	//Strings and waves to be used later on
	
	wave/T ImportCheck = root:Packages:Nika_RSoXS:NRS_Vars:Basenames //Wave to check if the given file has imported yet
	Wave/T FitsCheck = root:Packages:NIKA_RSoXS:NRS_Vars:files //Used to make sure you click on an actual file. Just an error check
	wave/T scans = root:Packages:Nika_RSoXS:NRS_Vars:scans //Wave to check if the given file has imported yet
	Wave/T FitsDirectory = root:Packages:NIKA_RSoXS:NRS_Vars:FitsDirectory
	Wave FitsSelect = root:Packages:NIKA_RSoXS:NRS_Vars:filesel


	SVAR pathtoheader = root:Packages:Nika_RSoXS:NRS_Vars:pathname
	SVAR PathtoFits = root:Packages:Nika_RSoXS:NRS_Vars:PathName_FIts

	SVAR header = root:Packages:Nika_RSoXS:NRS_Vars:header
	NVAR/Z FitsTopSelected = root:Packages:Nika_RSoXS:NRS_Vars:FitsTopSelected
	NVAR/Z SeriesSelected = root:packages:NIKA_RSoXS:NRS_Vars:SeriesSelected
	
	
	Variable row = lba.row
	Variable col = lba.col
	Wave/T/Z listwave = lba.listwave
	Wave/Z selwave = lba.selwave
	
	String CurrentFolder = GetDataFolder(1)
	
	if(StringMatch(lba.win,"NRS_MainPanel"))
	
		if(!cmpstr(lba.ctrlname, "NRS_PrimarySeriesList"))
			switch(lba.eventcode)
				case -1:
					break
				case 1: //Mouse down
					break
				case 3: //Double Click
					break
				case 4: //Cell selection
					//If it has not loaded data it will load relevant information 
					if(row >= numpnts(ImportCheck)) //Quick check to make sure you don't click on any whitespace
						return 0
					endif
					SeriesSelected = row
					//Check the current fits file slot for quick import when changing series
					If (strlen(ImportCheck[row]) < 1)
						NRS_LoadandStoreData(Scans[row],pathtoheader,row)
						print "Importing Data"
					endif
											
					//Update the fits path if it is not the same as the header file
					NewPath/O/Q/Z Path_1101Fits, FitsDirectory[row]
					
					NRS_LoadDataList(row)
					
					//Load the currently selected fits file when changing series. Requested Feature
					if(FitsTopSelected < numpnts(FitsCheck) && FitsTopSelected >= 0)
					//QUick reality check to make sure fits exist with that header
						if(StringMatch(FitsCheck[FitsTopSelected],"No Files Found"))
							print "No files found, please selected a different scan file"
							return 0
						endif
						
						header = previewfits(FitsTopSelected)
						NRS_PreviewVals()
						NRS_UpdatePopUpDisplay()
					else
						FitsTopSelected = 0
					endif
					//FitsSelect = 0 //Unselect all fits files
				case 5: // cell selection plus shift key
					break
				case 6: // begin edit
					break
				case 7: // finish edit
					break	
			endswitch
		endif
		if(!cmpstr(lba.ctrlname, "NRS_FitsList"))
		
			switch( lba.eventCode )
				case -1: // control being killed
					break
				case 1: // mouse down
					break
				case 3: // double click
					break
				case 4: // cell selection		
					//Quick error catches
					if(row >= numpnts(FitsCheck)) //Quick check to make sure you don't click on any whitespace...Otherwise it would error
						return 0
					elseif(StringMatch(FitsCheck[FitsTopSelected],"No Files Found"))
							print "No files found, please selected a different scan file"
							return 0
					endif
					FitsTopSelected = row //Save the currently selected row

					
					//String CurrentFolder=GetDataFolder(1)
					SetDataFolder root:Packages:Nika_RSoXS:NRS_Vars
					NewPath/O/Q/Z Path_1101Fits, FitsDirectory[SeriesSelected] //Grab the updated fits path

					header = PreviewFits(row)
					if(strlen(header) < 1 || !Svar_Exists(Header))
						print "No header loaded"
						return 0
					endif
					
					String Tempheader = NRS_UpdateHeaderKeys(header)
					NRS_PreviewVals()
					NRS_UpdatePopUpDisplay()

				case 5: // cell selection plus shift key
					break
				case 6: // begin edit
					break
				case 7: // finish edit
					break
			endswitch		
		endif
	endif
	SetDataFolder $CurrentFOlder
	return 0
end


/////////////////////////////////////////////////////////////////////////
//Quick function that updates the value display upon changing selection//
/////////////////////////////////////////////////////////////////////////

Function NRS_UpdatePopUpDisplay()

	SVAR/Z imagevalue = root:Packages:NIKA_RSoXS:NRS_Vars:imagevalue
	SVAR/Z ImageKeyPick = root:Packages:NIKA_RSoXS:NRS_Vars:imagekeypick
	SVAR/Z Header = root:Packages:NIKA_RSoXS:NRS_Vars:header
	
	
	if(!SVAR_Exists(ImageValue) || !Svar_Exists(ImageKeyPick) || !Svar_Exists(Header))
		print "Unable to update value display for popup menu"
		return 0
	endif

	imagevalue = stringbykey(imagekeypick,header,"=",";")

end

//////////////////////////////////////////////////////////////////////////////////////
//Controls to update the image with any changing variables in the advanced panel//////
//Allows for a scroll wheel to change locations of beam center or sector information//
//////////////////////////////////////////////////////////////////////////////////////
Function NRS_SetVarProc(sv) : SetVariableControl
	STRUCT WMSetVariableAction & sv
	
	NVAR UseSectors = root:Packages:Convert2Dto1D:UseSectors
	
	String CurrentFolder = GetDataFolder(1)
	

	Variable BCSelect

	if(Stringmatch(sv.win, "NRS_MainPanel#NRS_AdvancedOptions"))
	
		if(!cmpstr(sv.ctrlname,"NumberOfSectors"))
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			endif
		
		endif
		if(!cmpstr(sv.ctrlname,"SectorsStartAngle"))
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			endif
		
		endif
		if(!cmpstr(sv.ctrlname,"SectorsHalfWidth"))
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			endif
		
		endif
		if(!cmpstr(sv.ctrlname,"SectorsStepInAngle"))
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			endif
		
		endif
		
		if(!cmpstr(sv.ctrlname,"BeamCenterY"))
			ControlInfo/W=NRS_MainPanel NRS_BeamCenterPreview
			BCSelect = V_Value
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			elseif(BCSelect)
				setDrawLayer/K/W=NRS_MainPanel#disp Overlay
				NRS_DrawCenterIn2DGraph()		
			endif
		
		endif
		
		if(!cmpstr(sv.ctrlname,"BeamCenterX"))
			ControlInfo/W=NRS_MainPanel NRS_BeamCenterPreview
			BCSelect = V_Value
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			elseif(BCSelect)
				setDrawLayer/K/W=NRS_MainPanel#disp Overlay
				NRS_DrawCenterIn2DGraph()
			endif
		
		endif
		
		if(!cmpstr(sv.ctrlname,"HorizontalTilt"))
			ControlInfo/W=NRS_MainPanel NRS_BeamCenterPreview
			BCSelect = V_Value
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			elseif(BCSelect)
				setDrawLayer/K/W=NRS_MainPanel#disp Overlay
				NRS_DrawCenterIn2DGraph()
			endif
		
		endif
		
		if(!cmpstr(sv.ctrlname,"VerticalTilt"))
			ControlInfo/W=NRS_MainPanel NRS_BeamCenterPreview
			BCSelect = V_Value
			ControlInfo/W=NRS_MainPanel NRS_SectorPreview
			if(V_Value && UseSectors==1)
				setDrawLayer/K/W=NRS_MainPanel#disp ProgFront
				NRS_DrawSectorsIn2DGraph()
			elseif(BCSelect)
				setDrawLayer/K/W=NRS_MainPanel#disp Overlay
				NRS_DrawCenterIn2DGraph()
			endif
		
		endif
	endif	
	
	
	SetDataFolder $CurrentFolder
	
	
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Slider commands to change the color scale of the image based on the min and max values present in the data//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function NRS_SliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa	
	NVAR maxslide = root:Packages:Nika_RSoXS:NRS_Vars:maxslider
	NVAR minslide = root:Packages:Nika_RSoXS:NRS_Vars:minslider
	String CurrentFolder = GetDataFolder(1)
	
	if(Stringmatch(sa.win,"NRS_MainPanel"))
		if(!cmpstr(sa.ctrlname, "UpperBoundImage") || !cmpstr(sa.ctrlname,"LowerBoundImage"))
			switch( sa.eventCode )
				case -1: // kill
					break
				default:
					if( sa.eventCode & 1 ) // value set
						ModifyImage /W=NRS_MainPanel#Disp data_disp ctab= {minslide,maxslide,terrain256,0}
					endif
					break
			endswitch			
		endif	
	endif
	SetDataFolder $CurrentFolder
	return 0
End

//Quick function to reset the sliders back to their respective min and max points...Easy to put this here with other slider stuff
Function NRS_ResetSlider()
	NVAR maxslide = root:Packages:Nika_RSoXS:NRS_Vars:maxslider
	NVAR minslide = root:Packages:Nika_RSoXS:NRS_Vars:minslider
	Wave data_disp = root:Packages:Nika_RSoXS:NRS_Vars:data_disp	
	minslide = wavemin(data_disp)
	maxslide = wavemax(data_disp)		
	ModifyImage/W=NRS_mainPanel#disp Data_Disp ctab={*,*,terrain256,0}
	Slider UpperBoundImage, limits={wavemin(data_disp),wavemax(data_disp),0},win=NRS_MainPanel
	Slider LowerboundImage, limits={wavemin(data_disp),wavemax(data_disp),0},win=NRS_MainPanel
end

//////////////////////////////////////////////////////////////////////////////////////////////
//Support for popup menus on the main panel///////////////////////////////////////////////////
//Additional support to reroute the beam centering panel to use calibrations common in RSoXS//
//////////////////////////////////////////////////////////////////////////////////////////////

Function NRS_popupmenuproc(pm) : PopupMenuControl
	STRUCT WMPopupAction &pm
	
	//Variables and strings called
	SVAR Imagekeypick = root:packages:NIKA_RSoXS:NRS_Vars:imagekeypick
	SVAR header = root:packages:NIKA_RSoXS:NRS_Vars:header
	SVAR Imagevalue = root:packages:NIKA_RSoXS:NRS_Vars:imagevalue
	
	String CurrentFolder = GetDataFolder(1)

	
	
	
	if(Stringmatch(pm.win,"NRS_MainPanel"))
		if(!cmpstr(pm.ctrlname,"NRS_ImageKeyPopup"))
			ImageKeyPick = pm.popstr		 //key word you want to investigate on the file header
			imagevalue = stringbykey(imagekeypick,header,"=",";") //value to display
		endif
	
	
	endif
	//Updates the beam centering panel to look for different common calibration factors used by RSoXS.
	if(Stringmatch(pm.win,"NI1_CreateBmCntrFieldPanel"))
		if(!cmpstr(pm.ctrlname,"BmCalibrantName"))
			//All calibrant information
			NVAR BMCalibrantD1=root:Packages:Convert2Dto1D:BMCalibrantD1
			NVAR BMCalibrantD2=root:Packages:Convert2Dto1D:BMCalibrantD2
			NVAR BMCalibrantD3=root:Packages:Convert2Dto1D:BMCalibrantD3
			NVAR BMCalibrantD4=root:Packages:Convert2Dto1D:BMCalibrantD4
			NVAR BMCalibrantD5=root:Packages:Convert2Dto1D:BMCalibrantD5
			NVAR BMCalibrantD6=root:Packages:Convert2Dto1D:BMCalibrantD6
			NVAR BMCalibrantD7=root:Packages:Convert2Dto1D:BMCalibrantD7
			NVAR BMCalibrantD8=root:Packages:Convert2Dto1D:BMCalibrantD8
			NVAR BMCalibrantD9=root:Packages:Convert2Dto1D:BMCalibrantD9
			NVAR BMCalibrantD10=root:Packages:Convert2Dto1D:BMCalibrantD10
			NVAR BMUseCalibrantD1=root:Packages:Convert2Dto1D:BMUseCalibrantD1
			NVAR BMUseCalibrantD2=root:Packages:Convert2Dto1D:BMUseCalibrantD2
			NVAR BMUseCalibrantD3=root:Packages:Convert2Dto1D:BMUseCalibrantD3
			NVAR BMUseCalibrantD4=root:Packages:Convert2Dto1D:BMUseCalibrantD4
			NVAR BMUseCalibrantD5=root:Packages:Convert2Dto1D:BMUseCalibrantD5
			NVAR BMUseCalibrantD6=root:Packages:Convert2Dto1D:BMUseCalibrantD6
			NVAR BMUseCalibrantD7=root:Packages:Convert2Dto1D:BMUseCalibrantD7
			NVAR BMUseCalibrantD8=root:Packages:Convert2Dto1D:BMUseCalibrantD8
			NVAR BMUseCalibrantD9=root:Packages:Convert2Dto1D:BMUseCalibrantD9
			NVAR BMUseCalibrantD10=root:Packages:Convert2Dto1D:BMUseCalibrantD10
			if(!cmpstr(pm.popstr,"PS 100nm Spheres"))
				//added by Thomas Ferron 6/1/2016
				//Spacing calibrated through PS 300nm Spheres from Beamline
				BMCalibrantD1=0
				BMCalibrantD2=490
				BMCalibrantD3=322
				BMCalibrantD4=242
				BMCalibrantD5=195
				BMUseCalibrantD1=0
				BMUseCalibrantD2=1
				BMUseCalibrantD3=1
				BMUseCalibrantD4=1
				BMUseCalibrantD5=1
				BMUseCalibrantD6=0
				BMUseCalibrantD7=0
				BMUseCalibrantD8=0
				BMUseCalibrantD9=0
				BMUseCalibrantD10=0	
	
			elseif(!cmpstr(pm.popstr,"PS 300nm Spheres"))
				//added by Eliot
				//the D Spacing Measured on 11.0.1.2 is 391.4nm, or 3914 Angstroms.  The second Value
				//should be 3914/sqrt(3) = 2260
				BMCalibrantD1=1675
				BMCalibrantD2=955
				BMCalibrantD3=618
				BMCalibrantD4=0
				BMCalibrantD5=0
				BMUseCalibrantD1=1
				BMUseCalibrantD2=1
				BMUseCalibrantD3=1
				BMUseCalibrantD4=0
				BMUseCalibrantD5=0
				BMUseCalibrantD6=0
				BMUseCalibrantD7=0
				BMUseCalibrantD8=0
				BMUseCalibrantD9=0
				BMUseCalibrantD10=0	
			//From NIKA Calibration (Jan Ilavsky)
			elseif(!cmpstr(pm.popstr,"Ag behenate"))
			//The number I use is q = 0.1076 (1/Angstrom), d = 58.380 Angstroms.  The
			//reference is T.C. Huang et al, J. Appl. Cryst. (1993), 26, 180-184.
				BMCalibrantD1=58.380
				BMCalibrantD2=29.185
				BMCalibrantD3=19.46
				BMCalibrantD4=14.595
				BMCalibrantD5=11.676	//fixed form 11.767 on 2-12-2015, typo
				BMCalibrantD6=9.73
				BMCalibrantD7=8.34
				BMCalibrantD8=7.2975
				BMCalibrantD9=6.48667
				BMCalibrantD10=5.838
				BMUseCalibrantD1=1
				BMUseCalibrantD2=1
				BMUseCalibrantD3=1
				BMUseCalibrantD4=1
				BMUseCalibrantD5=1
				BMUseCalibrantD6=1
				BMUseCalibrantD7=1
				BMUseCalibrantD8=1
				BMUseCalibrantD9=1
				BMUseCalibrantD10=1	
			endif
			SetVariable BMCalibrantD1,win=NI1_CreateBmCntrFieldPanel, disable=(BMUseCalibrantD1==0)
			SetVariable BMCalibrantD1LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=(BMUseCalibrantD1==0)
			SetVariable BMCalibrantD2,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD2==0)
			SetVariable BMCalibrantD2LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD2==0)
			SetVariable BMCalibrantD3,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD3==0)
			SetVariable BMCalibrantD3LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD3==0)
			SetVariable BMCalibrantD4,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD4==0)
			SetVariable BMCalibrantD4LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD4==0)
			SetVariable BMCalibrantD5,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD5LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD6,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD6==0)
			SetVariable BMCalibrantD6LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD7,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD7==0)
			SetVariable BMCalibrantD7LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD8,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD8==0)
			SetVariable BMCalibrantD8LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD9,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD9==0)
			SetVariable BMCalibrantD9LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			SetVariable BMCalibrantD10,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD10==0)
			SetVariable BMCalibrantD10LineWidth,win=NI1_CreateBmCntrFieldPanel, disable=( BMUseCalibrantD5==0)
			NI1BC_DisplayCalibrantCircles()				
		endif
	endif
	
	SetDataFolder $CurrentFolder
	
end




////////////////////////////////////////////////////////////////////
//This region is for longer functions within buttons/listbox procs//
////////////////////////////////////////////////////////////////////


//Given a path for the header files and the fits files this will sort them numerically and update the SeriesListBox
//No information will be loaded here. Only the names of the files will be imported. Updating the text wave 'Scans'
//All other waves will be redimensioned to the number of series run
Function NRS_LoadSeriesList(path_header,path_fits)
	string path_header
	String path_fits
	String CurrentFolder = GetDataFolder(1)
	
	//Go to variable folder to grab the pre-created waves/strings required to populate listbox
	//Most of these waves will simply be redimensioned in this function. Things will be appended in a later function
	//Upon initially loading data
	SetDataFolder root:Packages:NIKA_RSoXS:NRS_Vars
	wave/T BaseNames, motors, scans, scantypes, times, motor1pos, motor2pos, I0s, FitsDirectory //grabs files to be used later upon loading data
	wave Frames, nmotors, scannum, FSIndex //grab non-text files for the same purpose
	
	//Makes waves for sorting data
	wave/T Series_Display //Names to display in the listbox
	wave/T Filedisc //Names to display in the fits listbox
	Wave/T ScanID_Display //Used to sort the scan files according to scanID...IGOR 7 changes sorting so previous method does not work
	Make/O/N=(numpnts(ScanID_DIsplay)) ScanID_Index =0 //Used for sorting later down the line

	//Quick check to make sure a path has been specified
	if(strlen(path_header) < 1)
		print "No data folder specified, please 'Browse' to locate data folder"
		return 0
	endif
	
	//Strings and Variables for later
	String Filename_List //List of filenames
	String Filename //Individual file name for later parsing
	String FitsName_List //List of Fits Files found 
	String FitsName //Individual fits name for later parsing
	
	//The following strings will be used for  organization and creating the name for the listbox
	String TempSeriesName //Full Name of series scan file Ex. "B48_27231"
	String TempScanID // Strickly the number of the scan Ex. "27231"
	String TempScanName //Strictly the name of the series Ex. "B48"
	Variable NumFits //Number of fits files under each scan name
	String TempDisplayName //Name to display in the list file


	//String graphname
	//Variable Result
	
	//Indexes used for loops
	Variable SeriesIndex = 0
	Variable Index = 0
	

	//Grab filenames of both the header files and fits files. The fits file will only be a different path if given
	FileName_List = indexedfile($path_header,-1,".txt") //grab the list of text files that are present in the folder given
	FitsName_List = Indexedfile($path_fits,-1,".fits")
	
	//Quick check to make sure each list even contains fits files or txt files. If not you will be asked to pick a different folder
	if(strlen(FileName_List) < 1)
		print "No header files found in selected folder, please select a different DataPath"
		Series_DIsplay = ""
		Filedisc = ""
		FitsDirectory = ""
		redimension/n=(0) Scans, Series_Display, ScanID_Display, ScanID_Index,FitsDirectory,filedisc
		return 0
	endif
//	if(strlen(FitsName_List) < 1)
//		print "No"
//	//	return 0
//	endif
//		
		//Resets everything as necessary
	ResetDataLists()
//Here we want to redimension the display wave since otherwise it will not allow you to change files into one with less files.		
//	redimension /N=0 Series_DIsplay, ScanID_Display , scans//Both the main wave and that for sorting. This simply sets them to zero.
// Loop through each header file and create display for the listbox Format: 4-digit ID: Num Scans" Full name
	Do		
		fileName = stringfromlist(index,filename_List, ";") //grab the first filename

		if (strlen(fileName) <1)			// No more files?
			break									// Break out of loop
		endif
		
		//This is to check to see if files were added to the folder and redimension the wave appropriately. (starts at 0 points in the wave)
		If(Index > (numpnts(Scans)-1))
			redimension/n=(numpnts(Scans)+1) Scans, Series_Display, ScanID_Display, ScanID_Index,FitsDirectory
		endif

		//This finds the matching header files
		if(stringmatch(filename,"*-AI.txt"))
			//Split filename into components used for labeling 
			Splitstring /E="(.*)_?(\\d{4})-AI.txt$" filename,TempScanName,TempScanID //grabs the scan name and the last 4 digits of the ID
			Splitstring /E="(.*)-AI.txt$" filename,TempSeriesName //grabs everything as a full name
			NumFits = LocateFits(path_header,TempSeriesName,index,0) //Counts the number of fits files corresponding to each .txt file
			
//			//CHeck here to see if it found any fits files...If not it will look in any child directory for matching fits files
//			if(Numfits == 0)
//				//Run Some function here
//				NumFits = LocateFits(path_header,TempSeriesName,index)
//			else
//				pathinfo $path_header
//				FitsDirectory[index] = S_Path
//			endif
//			
			//Create name for the listbox and sorting later on
			TempDisplayname = TempScanID + " - " + num2str(NumFits) + " image(s)" + " '"+TempScanName+"'"
			//Gather all information needed from this scan to allow later loading
			Series_Display[index] = TempDisplayName
			ScanID_Display[index] = TempScanID
			Scans[Index] = Filename
						
			SeriesIndex +=1 //number of files to be displayed
		endif
		index += 1 //Go to the next filename
	while (1)

	//After reading the files it will dimension a bunch of waves for saving data later. A bit of an artifact of older import methods. Does not slow things down so I kept it.
	redimension /n=(SeriesIndex) scanNum, FSindex, basenames, times, scantypes, motors, nmotors, frames,FitsDirectory
	redimension /n=(SeriesIndex,-1) motor1pos, motor2pos, I0s
	
	//Sort everything the same way based on ID. Not sure if this is the cleanest way to sort things. But it ensures everything is sorted equally.
	MakeIndex ScanID_Display, ScanID_Index
	IndexSort ScanID_Index, Series_Display, ScanID_Display, Scans, FitsDirectory
	
	SetDataFolder $CurrentFolder

end

Function LocateFits(pathname,SeriesName,index,level)

	String pathname,SeriesName
	Variable Index,level // Index for saving fits directory, and level for recursive checks for subfolders
	//Global Variables
	
	Wave/T/Z FitsDirectory = root:Packages:NIKA_RSoXS:NRS_Vars:FitsDirectory
	
	
	//Local Variables
	String FolderList,Fitslist,DummyFolderList,path
	variable numfits
	Variable FolderIndex
	Variable i,j //for loops
	
	//In the top level folder check for a match
	Fitslist = indexedfile($pathname,-1,".fits")
	numfits = NRS_Countfits(SeriesName,FitsList)
	
	if(NumFits > 0)
		pathinfo $pathname
		FitsDirectory[index] = S_path
		return NumFits
	endif
	
	
	//Find the list of all child folders to look there.
	
	FolderIndex = 0
	do
		path = IndexedDir($pathname,Folderindex,1)
		if(strlen(path) == 0)
			break //run out of folders to run through
		endif
		
		String SubFolderPathName = "Dummypath_"+num2istr(level+1)
		
		//Moving to the new parent folder
		String subFolderPath = path
		
		NewPath/Q/O $SubFolderPathName, SubFolderPath
		NumFits = LocateFits(SubFOlderPathName,SeriesName,Index,level+1)
		Killpath/Z $SubFolderPathName
	
		FolderIndex += 1
		
		if(NUmFits > 0)
			break //quick check to see if you already found the appropriate number of fit functions.
		endif
		
	while(1)
	
	return NumFits //should be nonzero by now
	
end
//	for(i=0;i<itemsinlist(FolderList,";");i+=1)
//		//Make a temp path that will change
//		NewPath/O/Q/Z Dummypath, Stringfromlist(i,Folderlist,";")
//		DummyFolderList = IndexedDir(Dummypath,-1,1)
//		
//		Fitslist = Indexedfile(Dummypath,-1,".fits")
//		numfits = NRS_CountFits(SeriesName,FitsList)
//		
//		if(NumFits > 0)
//			FitsDirectory[index] = Stringfromlist(i,FolderList,";")
//			return NumFits
//		elseif(NumFits == 0 && !Itemsinlist(DUmmyFolderList,";") == 0)
//		endif
//	endfor
//	
//	//If you got this far you were unable to find any fits
//	pathinfo $path
//	FitsDirectory[index] = S_Path
//	
//	return 0
//end





Function ResetDataLists()

	String CurrentFolder = GetDataFolder(1)
	
	//Go to variable folder to grab the pre-created waves/strings required to populate listbox
	//Most of these waves will simply be redimensioned in this function. Things will be appended in a later function
	//Upon initially loading data
	SetDataFolder root:Packages:NIKA_RSoXS:NRS_Vars
	wave/T BaseNames, motors, scans, scantypes, times, motor1pos, motor2pos, I0s //grabs files to be used later upon loading data
	wave Frames, nmotors, scannum, FSIndex //grab non-text files for the same purpose
	redimension/N=0 BaseNames, Motors, Scans, Scantypes, times, motor1pos, motor2pos, I0s
	BaseNames = ""; motors = ""; Scans = "", Scantypes = ""; times = ""; motor1pos = ""; motor2pos = ""; I0s = ""
	
	redimension/N=0 Frames, nmotors, scannum, FSindex
	Frames = 0; nmotors=0; scannum = 0; FSindex = 0
	
	Wave/T Series_Display, ScanID_Display, ScanID_Index
end



Function NRS_LoadDataList(row)
	variable row
	
	String CurrentFolder=GetDataFolder(1)
	SetDataFolder root:Packages:NIKA_RSoXS:NRS_Vars
	string/g filelist
	svar pathname_Fits
	wave/t basenames,motors,motor1pos,motor2pos
	wave nmotors
	
	
	print pathname_fits

//	print IndexedFile($pathname,-1, ".fits")
//	print "Break"
//	Print basenames
	
	filelist = ListMatch(sortlist(IndexedFile($pathname_fits,-1, ".fits"),";",17),basenames[row]+"*") //add png here as well

	//if only png, then the time is unknown
	if(strlen(filelist)<1)
		make /o /t /n=1 files,filedisc
		files = "No files found"
		filedisc = "No files found"
		filelist = ListMatch(IndexedFile($pathname_fits,-1, ".png"),basenames[row]+"*")
		if(strlen(filelist)>0)
//			displaypng(filelist)
		endif
		return 0
	endif
	make /o /t /n=(itemsinlist(filelist,";")) files,filedisc
	make /o /n=(itemsinlist(filelist,";")) filesel
	wave /t files,filedisc
	files = stringfromlist(p,filelist,";")
	If(nmotors[row]==0)
		filedisc = num2str(p+1)+"_"+basenames[row]
	elseif( nmotors[row]==1 )
		filedisc = num2str(p+1)+"_"+stringfromlist(p,motor1pos[row])+" ("+motors[row]+")"
	elseif( nmotors[row]==2 )
		filedisc = num2str(p+1)+"_"+stringfromlist(p,motor1pos[row])+" vs "+stringfromlist(p,motor2pos[row])+" ("+motors[row]+")"
	endif
//	ControlInfo fitslist
////	print row, filedisc //numpnts(s_datafolder)
//
//
////	if(v_value>0 && v_value<numpnts($ (s_datafolder+s_value) ) )
////		listbox /Z fitslist selrow=v_value
////		displayfits(v_value)
////	elseif(v_value>numpnts($ (s_datafolder+s_value) ))
////		listbox /Z fitslist selrow=numpnts($ (s_datafolder+s_value) )-1
////		displayfits(numpnts($ (s_datafolder+s_value) )-1)
////	else
////		listbox /Z fitslist selrow=0
////		displayfits(0)
////	endif
//
//
//	svar header = root:Packages:NIKA_RSoXS:NRS_Vars:headerinfo
//	variable hour,minute,second
//	sscanf stringbykey("DATETIME",header),"%*d/%*d/%*d %d%*[:]%d%*[:]%d",hour,minute,second
//	string /g imagetime =  num2str(hour)+":" + num2str(minute)
//	string/g imagekeys=""
//	string s1
//	variable i
//	for(i=0;i<itemsinlist(header,";");i+=1)
//		splitstring /e="^([^:]*):.*$" stringfromlist(i,header,";"),s1
//		imagekeys+=s1+";"
//	endfor
//	svar imagekeypick
//	string/g imagevalue
//	imagevalue = stringbykey(imagekeypick,header)
//	//PopupMenu Key,mode=1,popvalue=imagekeypick, win=DataReduction//,value=imagevalue
	SetDataFolder $CurrentFolder
	
	
end

//This seems to be a legacy function that loads data from the header file and stores it..
//looks at specific motors that are moved within a scan file and stores it. Common users will often not use this information

Function NRS_LoadandStoreData(Filename,path,index)
	string Filename, Path
	variable index
	variable fileref,nmotors1,len,i,ech
	string name,indexnum,basename,lineread,scantype,dummy,motor1,motor2

	String CurrentFolder=GetDataFolder(1)
	SetDataFolder root:Packages:NIKA_RSoXS:NRS_Vars
	Splitstring /E="(.*)_?(\\d{4})-AI.txt$" filename,name,indexnum
	Splitstring /E="(.*)-AI.txt$" filename,basename
	wave /t basenames,motors,scans,scantypes,times,motor1pos,motor2pos,I0s
	wave frames,nmotors,scanNum,FSindex
	basenames[index]=basename
	
	GetFileFolderInfo/z /q /p=$path filename
	if(v_flag !=0)
		print "File not found,"
		scans[index]="file could not be opened"
		return 0
	endif
	times[index]=secs2time(V_modificationDate,1)
	
	open /z /R /P=$path /T="Text" fileref as filename
	if(v_flag !=0)
		print "File could not be opened"
		scans[index]="File could not be opened"
		return 0
	endif
	
	FReadline /N=200 /T=(num2char(13)) fileref,lineread
	FReadline /N=200 /T=(num2char(13)) fileref,lineread
	FReadline /N=200 /T=(num2char(13)) fileref,lineread
	if(strlen(lineread)==0)
		print "not enough information in file"
		scans[index]="not enough information in file"
		close fileref
		return 0
	endif
	Splitstring /E=".*:(.*)$" lineread,scantype
	scantypes[index]=scantype
	len=strlen(lineread)
	ech=-1
	for(;(len>0)&&(ech==-1);) //Loop here to find the line before the header line...
		FReadline /N=500 /T=(num2char(13)) fileref,lineread
		ech =  strsearch(lineread, "Frame #",0) //As of 2/2/2018 the preceding line is the Data File Path:
		len = strlen(lineread)
	endfor
	//If the header ever reaches 'Frame #' that means its on the current line with all the header info. Might as well just read up till that line
//	FReadline /N=500 /T=(num2char(13)) fileref,lineread //Going to skip this line for now but use the above loop to hit Frame #
	if(strlen(lineread)==0)
		print "motor scan did not complete"
		scans[index]= "motor scan did not complete"
		close fileref
		return 0
	endif
	motor1 = stringfromlist(1,lineread,"	")
	motor2 = stringfromlist(2,lineread,"	")
	///Collins
	motor1= ShortenMotor(motor1)
	motor2= ShortenMotor(motor2)
	
	if(strlen(motor1)==0 || (!cmpstr(motor1,"CCD Temperature") || !cmpstr(motor1,"M3 Mirror Current")) || (!cmpstr(motor1,"Beam Current")))
		motors[index]="N/A"
		scans[index]= ""//"no Motors"
		nmotors1=0
	else 
		if(strlen(motor2)==0 || (!cmpstr(motor2,"I0"))  || (!cmpstr(motor2,"Beam Current"))  || (!cmpstr(motor2,"CCD Temperature")) || (!cmpstr(motor2,"M3 Mirror Current")))
			motors[index]= motor1
			scans[index]= motor1
			nmotors1=1
		else
			motors[index]=motor1 + " & " + motor2
			scans[index]=motor1 + " & " + motor2
			nmotors1=2
		endif
	endif
	variable I0loc,M1loc,M2loc,BCloc
	BCloc = whichlistitem("Beam Current",lineread,"	")
	m1loc = whichlistitem(Lmotor(motor1),lineread,"	")
	m2loc = whichlistitem(Lmotor(motor2),lineread,"	")
	if(m1loc==-1)
		m1loc=nan
	endif
	if(m2loc==-1)
		m2loc=nan
	endif
	I0loc = whichlistitem("I0",lineread,"	")
	nmotors[index]=nmotors1
	len=strlen(lineread)
	For(i=0;len != 0;i+=1)
		FReadline /N=500 /T=(num2char(13)) fileref,lineread
		len = strlen(lineread)
		if(len>0)
			motor1pos[index] += num2str(str2num(StringFromList(M1loc, lineread , "	")))+";"
			motor2pos[index] += num2str(str2num(StringFromList(M2loc, lineread , "	")))+";"
			I0s[index] += StringFromList(I0loc, lineread , "	")+";"
		endif
	endfor
	frames[index]= i-2
	if(i==1)
		frames[index]=1
	endif
	scanNum[index]=str2num(indexnum)
	FSindex[index]=index
	scans[index]=indexnum + " - " + num2str(frames[index]) + " image(s) " + scans[index] + " '" + Basename+"'"
	if(frames[index]>1)
		print scans[index]
	endif
	close fileref
	SetDataFolder $CurrentFolder
	return 0
end	

//Some simple naming changes that are legacy code 
//Not worth changing this
Function/T ShortenMotor(mName)
	string mName
	If( stringmatch(mName,"Beamline Energy Goal") )
		return "BLenGoal"
	elseif( stringmatch(mName,"EPU Polarization") )
		return "EPUpol"
	elseif( stringmatch(mName,"Sample X") )
		return "samX"
	elseif( stringmatch(mName,"Sample Y") )
		return "samY"
	else
		return mName
	endif
end

Function/T Lmotor(mName)
	string mName
	If( stringmatch(mName,"BLenGoal") )
		return "Beamline Energy Goal"
	elseif( stringmatch(mName,"EPUpol") )
		return "EPU Polarization"
	elseif( stringmatch(mName,"samX") )
		return "Sample X"
	elseif( stringmatch(mName,"samY") )
		return "Sample Y"
	else
		return mName
	endif
end



//After selecting a .fits file to load it will load the file and show a preview in the window
//Additionally it will grab the header and commonly desired information such as exposure and sample position
//This will NOT load data into NIKA.

//returns the header
Function/S PreviewFits(row)

	variable row
	
	//Grab necessary information to load the correct .fits file
	wave/T files = root:Packages:NIKA_RSoXS:NRS_Vars:files
	svar path = root:packages:NIKA_RSoXS:NRS_Vars:pathname_Fits
	//Other variables used for fixing the color scale and other nonsense like that
	NVAR minslide = root:packages:NIKA_RSoXS:NRS_Vars:minslider
	NVAR maxslide = root:packages:NIKA_RSoXS:NRS_Vars:maxslider
	NVAR HotPixel = root:packages:NIKA_RSoXS:NRS_Vars:HotPixel

	
	//This is borrowed from Line 1905 from NI1A_ConvProc - 'NI1A_ImportThisOneFile' All I want is to use the universal loader to get the CCDImagetoConvert file.
	//Once that is created we are good to go to display. This singular function should be all we need upon importing, we will grab the header afterwords as well
	//The header will only be imported assuming the RSoXS mods have been checked...Otherwise it will not save it
	variable loadedOK=NI1A_UniversalLoader(path,files[row],"fits","CCDImageToConvert")
	if(LoadedOK==0) //check to see if it loaded ok.
		return ""
	endif
	
	Wave CCDImageToConvert = root:Packages:Convert2Dto1D:CCDImageToConvert
	wave DisplayWave = root:Packages:NIKA_RSoXS:NRS_Vars:Data_Disp
	wave Data = root:Packages:NIKA_RSoXS:NRS_Vars:Data
	
	//Kills logo to not confuse people
	NRS_KillLogo()
	//
	
	Data = CCDImagetoConvert
	DisplayWave = Data
	NRS_CheckSaturatedPixel(DisplayWave)
	HotPixel = wavemax(Data) //reports the highest count in the image....lets you know how close you are to saturation
	
	wavestats/Q Displaywave
	//Displaywave /= V_Max * 0.24 //Makes the display wave a value from 0 to 4
	
	//check to see if you want it log or not
	ControlInfo/W=NRS_MainPanel NRS_LogPreview
	if(V_Value)
		DisplayWave = Log(data)
	endif
	
	NRS_ResetSlider()
	
	
	return note(CCDImageToConvert)
	
end

//Fills in some values to preview when you initially select a fits file
Function NRS_PreviewVals()

	SVAR header = root:Packages:Nika_RSoXS:NRS_Vars:header
	NVAR Energy = root:Packages:Nika_RSoXS:NRS_Vars:Preview_Energy
	NVAR EPUPol = root:Packages:Nika_RSoXS:NRS_Vars:Preview_EPUPol
	NVAR SAMX = root:Packages:Nika_RSoXS:NRS_Vars:Preview_SAMX
	NVAR SAMY = root:Packages:Nika_RSoXS:NRS_Vars:Preview_SAMY
	NVAR SAMZ = root:Packages:Nika_RSoXS:NRS_Vars:Preview_SAMz
	NVAR Exposure = root:Packages:Nika_RSoXS:NRS_Vars:Preview_Exposure

	
	
	Energy = str2num(Stringbykey("Beamline Energy",header,"=",";"))
	
	//Very rough way to round the energy to the nearest decimal point
	Energy *= 10
	Energy = Round(Energy)
	Energy /= 10
	
	EPUPol = round(str2num(Stringbykey("EPU Polarization",header,"=",";")))
	SAMX = (str2num(Stringbykey("Sample X",header,"=",";")))
	SAMY = (str2num(Stringbykey("Sample Y",header,"=",";")))
	SAMZ = (str2num(Stringbykey("Sample Z",header,"=",";")))
	Exposure = (str2num(Stringbykey("Exposure",header,"=",";")))
end


/////////////////////////////////////////////////////////////////
//This is the construction of the main panel


Function NRS_BuildMainPanel()


	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:packages:NIKA_RSoXS:NRS_Vars
	//Create variables that are unique to the panel itself
	Variable/G Saturated_Pixel = 0
	Variable/G LogPreviewImage
	Variable/G Preview_Energy = 0
	Variable/G Preview_EPUPol = 0
	Variable/G Preview_SAMX = 0
	Variable/G Preview_SAMY = 0
	Variable/G Preview_SAMZ = 0
	Variable/G Preview_Exposure = 0

	//Grab values that you will need following initialization
	
	SVAR pathtodata = root:packages:NIKA_RSoXS:NRS_Vars:pathtodata
	SVAR Imagekeys = root:packages:NIKA_RSoXS:NRS_Vars:imagekeys

	NVAR minslide = root:packages:NIKA_RSoXS:NRS_Vars:minslider
	NVAR minlimit = root:packages:NIKA_RSoXS:NRS_Vars:minslider


	NVAR maxslide = root:packages:NIKA_RSoXS:NRS_Vars:maxslider
	NVAR maxlimit = root:packages:NIKA_RSoXS:NRS_Vars:Maxlimit
	
	wave/T PrimarySeriesList = root:packages:NIKA_RSoXS:NRS_Vars:Series_Display
	wave/T FitsScanList= root:packages:NIKA_RSoXS:NRS_Vars:FileDisc
	Wave/T FitsSelectList = root:packages:NIKA_RSoXS:NRS_Vars:Filesel
	Wave Data_Disp = root:packages:NIKA_RSoXS:NRS_Vars:Data_Disp
	wave Saturate_Mask = root:Packages:Nika_RSoXS:NRS_Vars:Saturate_Mask

	//Build Panel
	DoWindow/K NRS_MainPanel
	NewPanel/N=$("NRS_MainPanel")/K=1 /W=(50,50,890,470) as "11.0.1.2 Data Reduction Panel"
	
	SetDrawEnv fillfgc= (39064,7710,12850),fillbgc= (0,65280,0)
	DrawRect 5,370,485,415
	SetDrawEnv fillfgc= (65535,65535,65535)
	DrawRect 270,310,485,365
//	SetDrawEnv fillfgc= (65280,65280,65280),fillbgc= (65280,65280,65280)
//	DrawRect 5,240,265,530
	
	
	//Buttons for data path selection
	Button NRS_Browse,pos={10,10},size={50,20},proc=NRS_Buttonproc,title="Browse",fsize=10,help={"Select a folder containing X-ray data"}
	TitleBox NRS_Pathtitle,pos={60,15},size={250,10},fSize=10,frame=0,fixedsize=1,variable = pathtodata
	Button NRS_UpdateBrowse,pos={270,10},size={50,20},proc=NRS_Buttonproc,title="Update",fsize=10 ,help={"Reload the chosen data folder"}
	
	//List RSoXS Data Series
	ListBox NRS_PrimarySeriesList, pos={10,35},size={250,330},proc=NRS_ListBoxProc,frame=3,editStyle= 1,fsize=10
	Listbox NRS_PrimarySeriesList, listwave = PrimarySeriesList,mode=1,selRow=1
	
	//List of Fits files for the given sample
	ListBox NRS_fitslist,pos={270,35},size={215,270},proc=NRS_Listboxproc,frame=2,fsize=10,editstyle=1,mode=4
	ListBox NRS_fitslist,listWave=FitsScanList,selWave=FitsSelectList
	
	//Convert Data into NIKA button
	Button NRS_ConvertData, fsize=10,pos={390,375},size={90,35},proc=NRS_Buttonproc,title="Convert\r Selection",help={"Load the selected data"}
	Button NRS_OpenMask, fsize=10, pos={105,375},size={90,35},proc=NRS_Buttonproc,title="Open for\r Mask",help={"Open the NIKA mask drawing panel to create a data mask for the selected sample"}
	Button NRS_BeamCenter, fsize=10,pos={200,375},size={90,35},proc=NRS_Buttonproc,title="Open for\r Beam Centering",help={"Open the NIKA beam centering panel to calibrate the geometry on the selected sample"}
	Button NRS_LoadasDark, fsize=10,pos={295,375},size={90,35},proc=NRS_Buttonproc,title="Load as Dark" , help={"Load the selected images for dark background subtraction.\r Each loaded image requires a matching dark with an equal exposure time"}
	Button NRS_LoadI0, fsize=10,pos={10,375},size={90,35},proc=NRS_Buttonproc,title="Load I0",help={"Open the 'I0 Load' panel"}
	Button NRS_Advanced, fsize=10,pos={750,10},size={80,20},proc=NRS_Buttonproc,title="Advanced" ,help={"More options to manually change NIKA parameters common in RSoXS analysis"}

	//Other stuff for value display
	ValDisplay NRS_SaturatePixel title="Saturated Pixels",size={90,15},pos={500,15},value=#"root:packages:NIKA_RSoXS:NRS_Vars:Saturated_Pixel",mode=2,fsize=10,limits={0,1,0},lowcolor=(0,0,0),barmisc={0,0},help={"Indicator if any pixels are saturated. These will show up red on the displayed image."}
	ValDisplay NRS_HotPixel,pos={600,15},size={100.00,15},title="Max Counts",fsize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:HotPixel", help = {"Highest count of a single pixel. Saturation occurs around 65000"}

	Popupmenu NRS_ImageKeyPopup proc=NRS_popupmenuproc,value=#"root:Packages:NIKA_RSoXS:NRS_Vars:imagekeys",pos={345,10},bodywidth=70
	Titlebox NRS_ImageKeyPreview pos={400,10},size={85,20},fsize=10,fixedsize=1,Variable=root:packages:NIKA_RSoXS:NRS_Vars:imagevalue

	//Options to display other data reduction parameters
	CheckBox NRS_LogPreview title="Log Image",side=0,fSize=10,pos={738,385},proc=NRS_CheckBoxProc,value=1,help={"Option to view image on a log scale"}
	CheckBox NRS_MaskPreview title="Display Mask",side=0,fSize=10,pos={738,400},proc=NRS_CheckBoxProc,value=0,help = {"Option to view the data mask"}
	CheckBox NRS_SectorPreview title="Display Sectors",side=0,fSize=10,pos={655,400},proc=NRS_CheckBoxProc,value=0, help={"Option to view the chosen sectors. No sectors will preview if you have them disabled"}
	CheckBox NRS_BeamCenterPreview title="Display Center",side=0,fSize=10,pos={655,385},proc=NRS_CheckBoxProc,value=0,help = {"Option to view the beam center"}

	//Just reading the header and displaying the most common values
	ValDisplay NRS_PreviewEnergy,pos={275,315},size={90,15},title="Energy",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_Energy"
	ValDisplay NRS_PreviewPol,pos={275,330},size={90,15},title="Polarization",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_EPUPol"
	ValDisplay NRS_PreviewPolarization,pos={275,345},size={90,15},title="Exposure",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_Exposure"
	ValDisplay NRS_PreviewSAMX,pos={380,315},size={100,15},title="Sample X",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_SAMX"
	ValDisplay NRS_PreviewSAMY,pos={380,330},size={100,15},title="Sample Y",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_SAMY"
	ValDisplay NRS_PreviewSAMZ,pos={380,345},size={100,15},title="Sample Z",fSize=10,value= #"root:packages:NIKA_RSoXS:NRS_Vars:Preview_SAMZ"

	//Preview Image commands and associated actions
	Display/W=(500,35,830,365)/HOST=#
	Appendimage Data_Disp	
	Appendimage Saturate_Mask

//	Appendimage /W=NRS_MainPanel#Disp Saturate_Mask
//	ModifyImage /W=NRS_MainPanel#Disp Saturate_Mask explicit=1,eval={1,65535,0,0},eval={255,-1,-1,-1}
	ModifyImage Data_Disp ctab={*,*,terrain256,0} //Display data
	ModifyImage Saturate_Mask explicit=1,eval={1,65535,0,0},eval={255,-1,-1,-1}	//Append mask

	
	ModifyGraph margin(left)=1,margin(bottom)=1,margin(top)=1,margin(right)=1,frameInset=2,mirror=2 //Some preferences
	
	//Add color control to graph...these will change based on log preferences
	Slider UpperBoundImage,pos={500,390},size={150,20},proc=NRS_SliderProc,limits={0,5,0},variable= maxslide,vert= 0,ticks= 0
	Slider LowerBoundImage,pos={500,370},size={150,20},proc=NRS_SliderProc,limits={0,5,0},variable= minslide,vert= 0,ticks= 0

	//Rename the subwindow (This worked and I don't want to change it)
	RenameWindow #,Disp
	SetActiveSubwindow ##
	
	NRS_DisplayLogo()
	SetActiveSubWindow NRS_MainPanel
	SetDataFolder $CurrentFolder


end

//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////ADVANCED PANEL OPTIONS/////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

	//Boxes to edit NIKA base variables such as geometry and sectors and such
	//These are exact duplicate buttons and check boxes to interact with those of the NIKA main panel. This removes the need to play around
	//With an additional panel...All credit goes to Jan Ilavsky at APS
	////GEOMETRY
	
Function NRS_AdvancedOptions()

	
	DoWindow /W=NRS_AdvancedOptions NRS_MainPanel
	Button NRS_Advanced, title="Hide Advanced",win=NRS_MainPanel,fcolor=(52428,1,1)
	
	//Quick check to see if the panel exists....if it does it kill it instead.
	if(V_flag)
		KillWindow/Z NRS_AdvancedOptions
		Button NRS_Advanced, title="Advanced",win=NRS_MainPanel,fcolor=(0,0,0)
		return 0
	endif
	NewPanel/HOST=NRS_MainPanel/EXT=0/N=$("NRS_AdvancedOptions")/K=1 /W=(0,0,310,400) as "Advanced"

	GroupBox NRS_AdvanedPanel,pos={5.00,5.00},size={300,350},title="Beamline Options"
	GroupBox NRS_AdvanedPanel,labelBack=(65535,65535,65535),fSize=14,frame=0

	TitleBox NRS_GeometryControl title="Sample Geometry",fsize=10,labelback=(65535,65535,65535),pos={10,25},frame=2
	
	SetVariable SampleToDetectorDistance,pos={100,27},size={200,15},fsize=10,proc=NI1A_PanelSetVarProc,title="Sample to CCD distance [mm]"
	SetVariable SampleToDetectorDistance,limits={0,inf,1},value= root:Packages:Convert2Dto1D:SampleToCCDDistance
	
	TitleBox GeometryDesc,pos={10,50},size={215.00,12.00},title="Direction          X (horizontal)                      Y (vertical)              "
	TitleBox GeometryDesc,labelBack=(65535,65535,65535),fSize=10,frame=2
	
//	SetVariable PixleSizeX,pos={13,300},size={120,10},bodyWidth=50,proc=NI1A_PanelSetVarProc,title="CCD pixel size [mm]"
//	SetVariable PixleSizeX,limits={0,inf,1},value= root:Packages:Convert2Dto1D:PixelSizeX
//	SetVariable PixleSizeY,pos={140,300},size={120,10},bodyWidth=50,proc=NI1A_PanelSetVarProc,title="CCD pixel size [mm]"
//	SetVariable PixleSizeY,limits={0,inf,1},value= root:Packages:Convert2Dto1D:PixelSizeY
	
	SetVariable BeamCenterX,pos={10,75},size={130,15},proc=NRS_SetVarProc,title="Beam center [pix]",fsize=10
	SetVariable BeamCenterX,value= root:Packages:Convert2Dto1D:BeamCenterX
	SetVariable BeamCenterY,pos={150,75},size={130,15},proc=NRS_SetVarProc,title="Beam center [pix]",fsize=10
	SetVariable BeamCenterY,value= root:Packages:Convert2Dto1D:BeamCenterY
	
	SetVariable HorizontalTilt,pos={10,95},size={130,15},proc=NRS_SetVarProc,title="Horizontal Tilt [deg]",fsize=10
	SetVariable HorizontalTilt,help={"Tilt of the image in horizontal plane (around 0 degrees)"}
	SetVariable HorizontalTilt,limits={-90,90,0},value= root:Packages:Convert2Dto1D:HorizontalTilt
	
	SetVariable VerticalTilt,pos={150,95},size={130,15},proc=NRS_SetVarProc,title="Vertical Tilt [deg]",fsize=10
	SetVariable VerticalTilt,help={"Tilt of the image in vertical plane (around 90 degrees)"}
	SetVariable VerticalTilt,limits={-90,90,0},value= root:Packages:Convert2Dto1D:VerticalTilt
	
	/////SECTORS
	TitleBox NRS_SectorControl,pos={10,120},size={70,20.00},title="Sector Options",frame=2
	TitleBox NRS_SectorControl,labelBack=(65535,65535,65535),fSize=10
	
	CheckBox DoSectorAverages,pos={90,123},size={95.00,13.00},title="Make sector averages?",side=1,fsize=10//,proc=NI1A_checkProc
	CheckBox DoSectorAverages,help={"Create data with sector average?"}
	CheckBox DoSectorAverages,variable= root:Packages:Convert2Dto1D:DoSectorAverages
	
	SetVariable NumberOfSectors,pos={10,145},size={140,20},proc=NRS_SetVarProc,title="Number of sectors",fsize=10
	SetVariable NumberOfSectors,help={"Number of sectors you want to create"}
	SetVariable NumberOfSectors,limits={0,inf,1},value= root:Packages:Convert2Dto1D:NumberOfSectors
	
	SetVariable SectorsStartAngle,pos={160,145},size={140,20},proc=NRS_SetVarProc,title="Start angle of sectors",fsize=10
	SetVariable SectorsStartAngle,help={"Angle around which first sectors is centered"}
	SetVariable SectorsStartAngle,limits={0,inf,1},value= root:Packages:Convert2Dto1D:SectorsStartAngle
	
	SetVariable SectorsHalfWidth,pos={10,165},size={140,20},proc=NRS_SetVarProc,title="Width of sector +/- ",fsize=10
	SetVariable SectorsHalfWidth,help={"Half width of sectors in degrees"}
	SetVariable SectorsHalfWidth,limits={0,inf,1},value= root:Packages:Convert2Dto1D:SectorsHalfWidth
	
	SetVariable SectorsStepInAngle,pos={160,165},size={140,20},proc=NRS_SetVarProc,title="Angle between sectors",fsize=10
	SetVariable SectorsStepInAngle,help={"Angle between center directions of sectors"}
	SetVariable SectorsStepInAngle,limits={0,inf,1},value= root:Packages:Convert2Dto1D:SectorsStepInAngle
	/////Data Importing Options
	
	TitleBox NRS_1DImportData, pos={10,190},size={70,20},title="1D Import Options"
	TitleBox NRS_1DImportData, labelBack=(65535,65535,65535),fSize=10,frame=2
	
	SetVariable QbinPoints,pos={10,215},size={140,14.00},title="Number of points",fsize=10
	SetVariable QbinPoints,limits={0,inf,10},value= root:Packages:Convert2Dto1D:QvectorNumberPoints

	CheckBox DisplayDataAfterProcessing,pos={160,215},size={78.00,13.00},title="Create 1D graph? ",side=1
	CheckBox DisplayDataAfterProcessing,variable= root:Packages:Convert2Dto1D:DisplayDataAfterProcessing,fsize=10


	/////MASK OPTIONS
	TitleBox NRS_MaskControl,pos={10,240},size={70,20.00},title="Mask Options"
	TitleBox NRS_MaskControl,labelBack=(65535,65535,65535),fSize=10,frame=2
	
	ListBox MaskListBoxSelection,pos={110,240},size={190,110},proc=NI1_MaskListBoxProc
	ListBox MaskListBoxSelection,help={"Select 2D data set for mask"}
	ListBox MaskListBoxSelection,listWave=root:Packages:Convert2Dto1D:ListOf2DMaskData
	ListBox MaskListBoxSelection,mode= 1,selRow= 0,fsize=10
		
	Checkbox NRS_UseMask, pos={13,265},size={78,13},title="Use Mask?",side=1
	Checkbox NRS_UseMask, Variable = root:Packages:Convert2Dto1D:UseMask,fsize=10
	
	Button MaskSelectPath,pos={10,285},size={90,30},proc=NRS_ButtonProc,title="Select mask\r data path",fsize=10
	Button MaskSelectPath,help={"Select path to mask file"}
	
	Button LoadMask,pos={10,320},size={90,30},proc=NRS_Buttonproc,title="Load mask",fsize=10
	Button LoadMask,help={"Load the mask file "}
	
	////Download Link
	BUtton DownloadLink, pos={95,360},size={120,30}, proc=NRS_ButtonProc, title="Support Website", fsize=10
	Button DownloadLink, Help={"Link to website to download analysis panel"},fcolor=(52428,1,1)



end
//
//	


/////////////////////////////////////////////////////////////////////////////////
//This section is for smaller functions that are called throughout
//Each one is fairly small with a singular function
/////////////////////////////////////////////////////////////////////////////////


Function NRS_BeamCenterPref()

	NVAR Logimage = root:Packages:Convert2Dto1D:BmCntrDisplayLogImage
	NVAR DezingImage = root:Packages:Convert2Dto1D:BMDezinger
	NVAR UseMask = root:Packages:Convert2Dto1D:BMUseMask
	
	
//	LogImage = 1
	DezingImage = 1
	//UseMask = 1
	
	NI1BC_DisplayHelpCircle()
	NI1BC_DisplayMask()
	showinfo /w=CCDImageForBmCntr
	
	//Sets up PS100 as the calibrant...This may need to be changed later on
end

Function NRS_DisplayLogo()

	String LogoPath = SpecialDirPath("Igor Pro User files",0,0,0) + "User Procedures:NRS_Procedures:"
	NewPath/q/z/o Path2Logo, LogoPath
	if(!waveexists(root:packages:NIKA_RSoXS:NRS_Vars:NRS_Logo))
		ImageLoad/Z/P=Path2Logo/T=tiff/Q/N=NRS_Logo "NRS_Logo.tif"
		if(!V_Flag)
			return 0 //No logo found, no hard feelings
		endif
	endif
	
	Display/W=(501,50,829,300)/N=Logo/HOST=NRS_MainPanel
	AppendImage/T/W=$("NRS_MainPanel#Logo") root:packages:NIKA_RSoXS:NRS_Vars:NRS_Logo
	ModifyImage/W=$("NRS_MainPanel#Logo") NRS_Logo ctab= {*,*,Grays,0}
	ModifyGraph/W=$("NRS_MainPanel#Logo") margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph/W=$("NRS_MainPanel#Logo") mirror=0,nticks=0,standoff=0,axthick=0,btlen=3
	ModifyGraph/W=$("NRS_MainPanel#Logo") tkLblRot(left)=90
	SetAxis/W=$("NRS_MainPanel#Logo")/A/R left
	
end

Function NRS_KillLogo()

	KillWIndow/Z NRS_MainPanel#Logo
	Killwaves/Z root:packages:NIKA_RSoXS:NRS_Vars:NRS_logo

end

//Counts the number of .fits files matching the given scan ID. This loop might be faster than opening files and reading?
//Testing that out right now, looks fine

Function NRS_CountFits(ScanName,FitsName_List)
	String ScanName //Name of file without the -Al ending
	String FitsName_List //Full list of fits files in the folder
	
	string fitsname //Name of individual fits file as it scans
	
	Variable Index = 0 //Index to cycle through the fits files
	Variable NumFits = 0 //Number of fits files that match the given string
	
	Do //Cycle through the fits list
		fitsname = StringFromList(index,FitsName_List,";")
		if(strlen(FitsName) < 1) //No more fits files?
			break //breaks loop
		endif
		if(stringmatch(fitsname,ScanName+"-*.fits"))
			//May want to save something here
			
			
			NumFits += 1
		endif
		index += 1
	while (1)
	
	return NumFits
end

//Returns a list of keywords able to be searched through the given keyword searchable string 'header'
//Designed with the .fits header in mind...Probabily can be written to be used elsewhere. Inpendent assignment of time as well
Function/S NRS_UpdateHeaderKeys(header)

	String header
	
	//Strings to assigned in this function
	SVAR imagetime = root:Packages:NIKA_RSoXS:NRS_Vars:Imagetime
	SVAR imagekeys = root:Packages:NIKA_RSoXS:NRS_Vars:ImageKeys
	SVAR imagekeypick = root:Packages:NIKA_RSoXS:NRS_Vars:ImageKeyPick
	SVAR imagevalue = root:Packages:NIKA_RSoXS:NRS_Vars:ImageValue

	//Parse the date and time from the header
	variable hour,minute,second
	sscanf stringbykey("DATE",header,"=",";"),"'%*d%*[-]%*d%*[-]%*d%*[T]%d%*[:]%d%*[:]%d'",hour,minute,second //Changed "DATETIME" to "DAte" and added key info
	imagetime =  num2str(hour)+":" + num2str(minute)

	//Build keyword string
	imagekeys=""
	string s1 //temp keyword to append to imagekeys
	variable i //to loop through keywords in the header
	for(i=0;i<itemsinlist(header,";");i+=1)
		//String DEbugg= stringfromlist(i,header,";")
		splitstring /e="^([^:]*)=[^:]*$" stringfromlist(i,header,";"),s1 //Chanaged  a ":" to a "=" to account for new header format
		imagekeys+=s1+";"
	endfor

	return imagekeys
//	svar imagekeypick
	//PopupMenu Key,mode=1,popvalue=imagekeypick, win=DataReduction//,value=imagekeys
//	string /g imagevalue = stringbykey(imagekeypick,header)

end


//I am about 99% sure that when importing a .fits file saturated images will haave a negative value in the imported CCD data.
//This Function will read the 2D map and check for negative values. If so it will make a quick mask and append it to the image over the dead pixels
//The mask will remain on the displayed image. As you load data this will be checked really quick and updated. Otherwise it will just be full of NaNs.
Function NRS_CheckSaturatedPixel(map)
	wave map	
	
	NVAR Sat = root:packages:NIKA_RSoXS:NRS_Vars:Saturated_Pixel //LED that turns on if saturated
	
	String CurrentFolder = GetDatafolder(1)
	SetDataFolder root:Packages:Nika_RSoXS:NRS_Vars:
	Wave Saturate_Mask
	wavestats/Z/Q map
	
	if ( V_min < 0 ) //Check to see if any values were saturated
		duplicate/o Map Saturate_Mask
		Saturate_Mask = Saturate_Mask > 0 ? NaN : Saturate_Mask //removes all good data points
		Saturate_Mask = Saturate_Mask < 0 ? 1 : Saturate_Mask //Highlights bad data points
		Make/N=1/o Num_SaturatedPixels
		Sat = 1
		
	else
		Sat =0
		Saturate_Mask = NaN
	endif
	
	SetDataFolder $CurrentFolder
end

//Function NRS_AppendSaturatedPixelMask()
//
//	wave Saturate_Mask = root:Packages:Nika_RSoXS:NRS_Vars:Saturate_Mask
//	
//	Appendimage /W=NRS_MainPanel#Disp Saturate_Mask
//	ModifyImage /W=NRS_MainPanel#Disp Saturate_Mask explicit=1,eval={1,65535,0,0},eval={255,-1,-1,-1}
//	
//end
	
function NRS_convertpathtonika([main,mask,dark,beamcenter])
	variable mask,dark,beamcenter,main
	PathInfo Path_1101fits //CHanged from Path_1101Panel
	doupdate
	if(main)
		svar SampleNameMatchStr = root:Packages:Convert2Dto1D:SampleNameMatchStr
		SampleNameMatchStr = ""
		newpath /O/Q/Z Convert2Dto1DDataPath S_path
		SVAR MainPathInfoStr=root:Packages:Convert2Dto1D:MainPathInfoStr
		MainPathInfoStr=S_path
		NI1A_UpdateDataListBox()	
	endif
	if(mask)
		NI1M_CreateMask()
		newpath /O/Q/Z Convert2Dto1DMaskPath S_path
		popupmenu CCDFileExtension win=NI1M_ImageROIPanel, popmatch="fits"
		SVAR CCDFileExtension=root:Packages:Convert2Dto1D:CCDFileExtension
		CCDFileExtension = "fits"
		NI1M_UpdateMaskListBox()
	endif
	if(dark)
		newpath /O/Q/Z Convert2Dto1DEmptyDarkPath S_path
//		popupmenu SelectBlank2DDataType win=NI1A_Convert2Dto1DPanel, popmatch="*fits"
		nVAR usedarkfield=root:Packages:Convert2Dto1D:UseDarkField
		usedarkfield=1
//		SVAR BlankFileExtension=root:Packages:Convert2Dto1D:BlankFileExtension
//		BlankFileExtension = ".fits"
//		SVAR DataFileExtension=root:Packages:Convert2Dto1D:DataFileExtension
//		DataFileExtension = ".fits"
		svar EmptyDarkNameMatchStr = root:Packages:Convert2Dto1D:EmptyDarkNameMatchStr
		EmptyDarkNameMatchStr = ""
		NI1A_UpdateEmptyDarkListBox()	
	endif
	//Sets up the Beam Center and Calibration Panel for use in RSoXS
	//Just reroutes some panel options to other locations
	if(beamcenter)
		NI1_CreateBmCntrFile()
		newpath /O/Q/Z Convert2Dto1DBmCntrPath S_path
		popupmenu BmCntrFileType win=NI1_CreateBmCntrFieldPanel, popmatch="fits"
		PopupMenu BmCalibrantName value="User;PS 100nm Spheres;PS 300nm Spheres;Ag behenate", proc=NRS_popupmenuproc
		SVAR BmCntrFileType=root:Packages:Convert2Dto1D:BmCntrFileType
		BmCntrFileType = "fits"
		SVAR BCPathInfoStr=root:Packages:Convert2Dto1D:BCPathInfoStr
		BCPathInfoStr=S_Path
		NI1BC_UpdateBmCntrListBox()
	endif
end
	

function /t getfilename([list,selwave])
	variable list
	wave selwave
	wave/t files = root:Packages:Nika_RSoXS:NRS_Vars:files
	if(waveexists(selwave))
		string filelistsel=""
		variable i
		for(i=0;i<numpnts(files);i+=1)
			if(selwave[i])
				filelistsel = AddListItem(files[i], filelistsel)
			endif
		endfor
		return filelistsel
	elseif(paramisdefault(list)||list==0)
		controlInfo /W=DataReduction fitslist
		variable row = V_Value
		return files[row]
	else
		string filelist=""
		variable j
		for(j=0;j<numpnts(files);j+=1)
			filelist = AddListItem(files[j], filelist)
		endfor
		return filelist
	endif
end
	
	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////
//This function loads the appropriate mask selected by the user into NIKA. This is identical to 'NI1A_LoadMask()' with the controlInfo looking at a different panel
//All functions below were created by Jan Ilavsky of APS
Function NRS_LoadMask()

	string OldDf=GetDataFOlder(1)
	setDataFOlder root:Packages:Convert2Dto1D
	Wave/T  ListOf2DMaskData=root:Packages:Convert2Dto1D:ListOf2DMaskData
	SVAR CurrentMaskFileName=root:Packages:Convert2Dto1D:CurrentMaskFileName

	controlInfo /W=NRS_MainPanel#NRS_AdvancedOptions MaskListBoxSelection //changed the panel that you look at here.
	variable selection = V_Value
	if(selection<0)
		setDataFolder OldDf
		abort
	endif
	SVAR FileNameToLoad
	FileNameToLoad=ListOf2DMaskData[selection]
	if(stringmatch(FileNameToLoad[strlen(FileNameToLoad)-4,inf],".tif"))
		NI1A_UniversalLoader("Convert2Dto1DMaskPath",FileNameToLoad,"tiff","M_ROIMask")
	else
		NI1_MaskHDFLoader("Convert2Dto1DMaskPath",FileNameToLoad,".hdf","M_ROIMask")
		//NI1A_UniversalLoader
	endif

	CurrentMaskFileName = FileNameToLoad
	wave M_ROIMask
	Redimension/B/U M_ROIMask
	M_ROIMask=M_ROIMask>0.5 ? 1 : 0
	setDataFolder oldDf
end


//This is from NIKA to append sectors to the displayed image...I changed the window that you append stuff too...That is all.
// Version 1.0 - Fixed a bug where the image in the preview was inverted...but not the sectors. It will now appropriately display the sectors.
Function NRS_DrawSectorsIn2DGraph()

	IN2G_PrintDebugStatement(IrenaDebugLevel, 5,"")
	string oldDf=GetDataFOlder(1)
	setDataFolder root:Packages:Convert2Dto1D

 	DoWindow /W=NRS_MainPanel disp
	if(V_Flag)
	     setDrawLayer/W=NRS_MainPanel#disp ProgFront //changed
		NVAR ycenter=root:Packages:Convert2Dto1D:BeamCenterY
		NVAR xcenter=root:Packages:Convert2Dto1D:BeamCenterX
		NVAR DoSectorAverages=root:Packages:Convert2Dto1D:DoSectorAverages
		NVAR UseSectors = root:Packages:Convert2Dto1D:UseSectors
		NVAR UseLineProfile=root:Packages:Convert2Dto1D:UseLineProfile
		NVAR NumberOfSectors=root:Packages:Convert2Dto1D:NumberOfSectors
		NVAR SectorsStartAngle=root:Packages:Convert2Dto1D:SectorsStartAngle
		NVAR SectorsHalfWidth=root:Packages:Convert2Dto1D:SectorsHalfWidth
		
		Wave CCDImageToConvert=root:Packages:Nika_RSoXS:NRS_Vars:Data_disp //changed
		
		NVAR SectorsStepInAngle=root:Packages:Convert2Dto1D:SectorsStepInAngle
		variable i, tempEndX, tempEndY, sectorCenterAngle, tempLength
		variable temp1, temp2, temp3, temp4
		
		if(DoSectorAverages && UseSectors)
			For(i=0;i<NumberOfSectors;i+=1)
				//calculate coordinates for lines...
				sectorCenterAngle = SectorsStartAngle+90 + i*(SectorsStepInAngle)
				if(sectorCenterAngle>=90 && sectorCenterAngle<180)
					temp1 = DimSize(CCDImageToConvert, 0 )-xcenter
					temp2 = ycenter
				elseif(sectorCenterAngle>=180 && sectorCenterAngle<270)
					temp1 = xcenter
					temp2= ycenter
				elseif(sectorCenterAngle>=270 && sectorCenterAngle<360)
					temp1 = xcenter
					temp2= DimSize(CCDImageToConvert, 1)-ycenter
				elseif(sectorCenterAngle>=360 && sectorCenterAngle<450)
					temp1 = DimSize(CCDImageToConvert, 0 )-xcenter
					temp2= DimSize(CCDImageToConvert, 1)-ycenter
				endif
				tempLength = sqrt((temp1 * sin(pi/180*sectorCenterAngle))^2+ (temp2 * cos(pi/180*sectorCenterAngle))^2)
				//center line
				tempEndX= (xcenter + (tempLength)*sin(pi/180*(sectorCenterAngle)))
				tempEndY=(ycenter + (tempLength)*-cos(pi/180*(sectorCenterAngle))) //Flipped the sign to follow the flipping image
				string AxList= AxisList("NRS_MainPanel#disp" )
				if(stringMatch(axlist,"*top*"))
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=top,ycoord=left,save
				else
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=bottom,ycoord=left,save
				endif
				SetDrawEnv/W=NRS_MainPanel#disp linefgc= (8704,8704,8704),dash= 7  
				SetDrawEnv /W=NRS_MainPanel#disp linethick=2
				Drawline/W=NRS_MainPanel#disp xcenter, ycenter, tempEndX, tempEndY
				//side lines
				tempEndX= (xcenter + (tempLength)*sin(pi/180*(sectorCenterAngle-SectorsHalfWidth)))
				tempEndY=(ycenter + (tempLength)*-cos(pi/180*(sectorCenterAngle-SectorsHalfWidth))) //Flipped the sign to follow the flipping image
				if(stringMatch(axlist,"*top*"))
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=top,ycoord=left,save
				else
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=bottom,ycoord=left,save
				endif
				SetDrawEnv/W=NRS_MainPanel#disp linefgc= (65280,65280,0)
				SetDrawEnv /W=NRS_MainPanel#disp dash= 2,linethick= 1.00
				Drawline/W=NRS_MainPanel#disp xcenter, ycenter, tempEndX, tempEndY
				tempEndX=(xcenter + (tempLength)*sin(pi/180*(sectorCenterAngle+SectorsHalfWidth)))
				tempEndY=(ycenter + (tempLength)*-cos(pi/180*(sectorCenterAngle+SectorsHalfWidth))) //Flipped the sign to follow the flipping image
				if(stringMatch(axlist,"*top*"))
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=top,ycoord=left,save
				else
					setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=bottom,ycoord=left,save
				endif
				SetDrawEnv/W=NRS_MainPanel#disp linefgc= (65280,65280,0)
				SetDrawEnv /W=NRS_MainPanel#disp dash= 2,linethick= 1.00
				Drawline/W=NRS_MainPanel#disp xcenter, ycenter, tempEndX, tempEndY
			 endfor
	      endif
	endif
	setDataFolder OldDf
End


//Builds a cross on the beam center...Otherwise this code comes from NIKA
Function NRS_DrawCenterIn2DGraph()
	string oldDf=GetDataFOlder(1)
	setDataFolder root:Packages:Convert2Dto1D
 	DoWindow /W=NRS_Mainpanel disp
 	Variable plot = V_Flag
 	ControlInfo/W=NRS_MainPanel NRS_BeamCenterPreview

	if(Plot && V_Value)
	     setDrawLayer/W=NRS_MainPanel#disp Overlay
	     NVAR ycenter=root:Packages:Convert2Dto1D:BeamCenterY
	     NVAR xcenter=root:Packages:Convert2Dto1D:BeamCenterX
			if(stringMatch(AxisList("NRS_MainPanel#Disp"),"*top*"))
				setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=top,ycoord=left,save
			else
				setdrawenv/W=NRS_MainPanel#disp fillpat=0,xcoord=bottom,ycoord=left,save
			endif
			
			SetDrawEnv/W=NRS_MainPanel#disp linefgc=(65535, 0,0 )
			SetDrawEnv/W=NRS_MainPanel#disp linethick=2
			DrawLine/W=NRS_MainPanel#disp xcenter-20, ycenter, xcenter+20, ycenter
			
			SetDrawEnv/W=NRS_MainPanel#disp linefgc=(65535, 0,0 )
			SetDrawEnv/W=NRS_MainPanel#disp linethick=2
			DrawLine/W=NRS_MainPanel#disp xcenter, ycenter-20, xcenter, ycenter+20
	endif 	
end

//This funciton will check to see if an I0 has been loaded.
//If nothing is found it will create a dummy I0 for NIKA
Function NRS_CheckForI0()

	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:NIKA_RSOXS
	NVAR SampleExposure = root:packages:NIKA_RSoXS:NRS_Vars:Preview_EPUPol
	String ListofI0  = wavelist("CorrectionFactor_*",";","")
	
	variable i //for loops
	
	for(i=0;i<itemsinlist(ListofI0,";");i++)
		//check for P100
		if(stringmatch(stringfromlist(i,listofI0,";"),"*Pol"+num2str(SampleExposure)+"*"))
			//I0 exists so we move on
			return 0
		endif
	endfor
					
	String TempCF = "CorrectionFactor_Pol"+num2str(SampleExposure)
	String TempBE = "Beamline_Energy_Pol"+num2str(SampleExposure)
	Make/O/N=2 $TempCF = {1,1}	
	Make/O/N=2 $TempBE = {1,10000}
	print "make a fake I0 here"
end
//Similarly makes a dark matching the exposure time of the selected file.

Function NRS_CheckforDark()

	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:packages:Convert2Dto1D
	
	NVAR SampleExposure = root:packages:NIKA_RSoXS:NRS_Vars:Preview_Exposure
	Wave/Z DataMap = root:Packages:Convert2Dto1D:CCDImageToConvert
//	print SampleExposure
	
	String ListofDarks  = wavelist("DarkFieldData*",";","")
//	print ListofDarks
	Variable i
	for(i=0;i<itemsinlist(ListofDarks,";");i+=1)
		//Check for a dark		
		//Replace all decimal points with a p
		String dummystring = replacestring(".",Stringfromlist(i,ListofDarks,";"),"p")
		if(stringmatch(DummyString,"*"+Num2str(SampleExposure)))
			//Dark exists so move on
			return 0
		endif
	endfor
	//No dark exists
	String DarkImage = ReplaceString(".","DarkFieldData_"+num2str(SampleExposure),"p")
	duplicate/o DataMap $DarkImage
	Wave Dark_Dummy = $DarkImage
	Dark_Dummy = 0
end
