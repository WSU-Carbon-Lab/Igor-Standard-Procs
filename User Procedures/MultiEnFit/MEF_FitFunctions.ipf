#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function MEF_CalcInt(pop,En)
	Variable pop, En //Population and Energy for the fit
	
	//Folder Stuff
	String CurrentFolder = GetDataFolder(1)
	SetDataFolder root:Packages:MEF_Vars
	
	//Start with initializing the model wave
	String ModelQs = "Qmodel_Set"+num2str(En)
	Wave ModeLQ = $ModelQs
	if(!waveexists(ModelQ))
		ABort "Qmodel not found"
	endif
	//INitialize the modelInt
	Duplicate/o ModelQ, $("IntensityModel_Set"+num2str(En)+"_pop"+num2str(pop))
	Wave ModelInt = $("IntensityModel_Set"+num2str(En)+"_pop"+num2str(pop))
	ModelInt = 0
	
	//Figure out what fit we are running
	SVAR Pop_FF = $("root:Packages:MEF_Vars:FormFactor_pop"+num2str(pop))
	//Cycle though the options and do the right function
	
	if(!cmpstr(Pop_FF,"Power_Law")) //Background is added later
		//A * x^pow + B * x
		NVAR A = $("root:Packages:MEF_Vars:PorodScale_En"+num2str(En)+"_pop"+num2str(pop))
		NVAR B = $("root:Packages:MEF_Vars:Porod_line_pop"+num2str(pop))
		NVAR pow = $("root:Packages:MEF_Vars:PorodPower_En"+num2str(En)+"_pop"+num2str(pop))
		
	//	if(En==1)
	//		ModelInt = A*ModelQ^pow //+ B*ModelQ	
	//	elseif(En==2)
	//		ModelInt = B*ModelQ^pow
	//	endif
		
		ModelInt = A*ModelQ^pow
		
	endif
	
	
	SetDataFolder $CurrentFolder
	
	
end