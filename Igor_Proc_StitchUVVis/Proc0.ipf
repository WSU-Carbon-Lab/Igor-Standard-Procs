#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//This routine will stitch UV-Vis spectra obtained from thw two light sources in our lab

//1. Load the UV spectrum
//2. Load the Vis spectrum
//3. Select stitch point/wavelength
//4. Find intensity for both spectra at selected stitch point
//5. Calculate scale factor based on intensity of both spectra at each point
//6. Multiply spectra by scale factor

Function/WAVE stitchUVVisSpec(wl,xw1,yw1,yw2)

	Variable wl	//Wavelength that we want to stitch things at
	Wave xw1 //Wavelength wave
	Wave yw1,yw2//Absorbance wave from UV source and Vis source respectively
	
	
	FindLevel/Q xw1,wl
	Variable p1 = round(V_LevelX)
	
	Variable int1 = yw1[p1]
	Variable int2 = yw2[p1]
	Variable scaleFactor = int1/int2
	
	print p1,int1,int2
	
	//Duplicate/O yw1,stitchSpec
	Duplicate/O yw2,stitchSpec
	Variable n = numpnts(yw1)
	
	//stitchSpec[p1,n-1] = yw2 * scaleFactor
	stitchSpec[0,p1-1] = yw1 / scaleFactor
	return stitchSpec
End 

Function ImportQE(sampleName,uvName,visName,pathName,[sp])	//Imports QEPro absorbance data
	
	String sampleName //Name for data folder
	String uvName	//Name for UV Spec
	String visName	//Name for Vis Spec
	String pathName	//Path to data. Just enter ""
	Variable sp //Wavelength to stitch data at. Default value is 400nm
	
	if(ParamIsDefault(sp))
		sp = 440
	endif
	
	String dataFolder = GetDataFolder(1)
	String foldername = "root:" +RemoveEnding(uvName,".txt")
	String dataName2 = RemoveEnding(uvName,".txt")
	
	NewDataFolder/O/S $sampleName
	String columnInfoStr = " "
	
	columnInfoStr += "C=1,F=0,N=Wavelength;"
	columnInfoStr += "C=1,F=0,N="+ uvName +";"
	
	LoadWave/G/B=columnInfoStr/D/W/N/O/Q/P=$pathName uvName
	
	columnInfoStr = " "
	columnInfoStr += "C=1,F=0,N=Wavelength;"
	columnInfoStr += "C=1,F=0,N="+ visName +";"
	LoadWave/G/B=columnInfoStr/D/W/N/O/Q/P=$pathName visName
	
	String xwName = "Wavelength"
	Wave wl = $xwName
	Wave uvSpec = $uvName
	Wave visSpec = $visName
	
	//Stitch the spectra
	Wave  stitchSpec = stitchUVVisSpec(sp,wl,uvSpec,visSpec)
	
	//Plot the spectra
	DoWindow UV_Vis_Plot
	if(!V_Flag)
		Display/N=UV_Vis_Plot/w=(0,0,350,400)/K=1 uvSpec,visSpec,stitchSpec vs wl
		Label left "Absorbance [a.u]"
		Label bottom "Wavelength[nm]"
		ModifyGraph mirror=1,minor=1,fStyle=1
		Legend/C/N=text0/J/F=0/M/A=RT "\\s("+uvName+") UV\r\\s("+visName+") VIS\r\\s(stitchSpec) Stitch"
		ModifyGraph rgb($visName)=(0,0,65535),lstyle(stitchSpec)=3,rgb(stitchSpec)=(0,0,0)
		SetAxis left 0,*
	endif
	
	SetDataFolder $dataFolder
End