Function/WAVE scaleExptoDFTStep(allExpEnergy,stepSum,anchorStep1,anchorStep2,mol,expSpecName,method)
	
	Wave allExpEnergy,stepSum
	Variable anchorStep1,anchorStep2,method
	String mol,expSpecName
	
	Variable tol = 1e-3,expScale
	//Make the concatenated NEXAFS wave
	Wave allExpSpec = makeLongSpecWave(expSpecName)
	
	//Find the bare atom waves
	Wave mu_energy = findBAenergy()
	Wave mu = findBAMA(mol)
	WaveStats/Q stepSum
	
	Duplicate/O allExpSpec,allExpSpec2
	//Find the experimental energies used for scaling   
	Variable expLo = findWaveValAtEergy(allExpEnergy,allExpSpec,anchorStep1)
	Variable expHi = findWaveValAtEergy(allExpEnergy,allExpSpec,anchorStep2)
	
	//Find the bare atom energies used for scaling   
	Variable bareStepLo = findWaveValAtEergy(mu_energy,mu,anchorStep1)
	Variable bareStepHi = findWaveValAtEergy(mu_energy,mu,anchorStep2)
	
	//Choose 1-Point vs 2-Point Scaling
	if(method == 1)
		expScale = expLo/bareStepLo //Try anchoring at just preedge
		print expLo,bareStepLo,expScale
	elseif(method == 2)
		expScale = (expHi-expLo) / (bareStepHi-bareStepLo)//Anchor at post and preedge	
	endif
	Variable eta = 1
	allExpSpec2 = allExpSpec / (expScale*eta)
	Variable expLo_2 = findWaveValAtEergy(allExpEnergy,allExpSpec2,anchorStep1)
	Variable dif = expLo_2 - bareStepLo
	if(abs(dif) <= tol) 
		print dif//This  should be zero!
	else
		print "Difference between experiment and step is larger than tolerance. Check."
		if(expLo_2 > bareStepLo)
			allExpSpec2 -= abs(dif)
		elseif(expLo_2 < bareStepLo)
			allExpSpec2 += abs(dif)
		endif
		expLo_2 = findWaveValAtEergy(allExpEnergy,allExpSpec2,anchorStep1)
		dif = expLo_2 - bareStepLo
		print dif
	endif
	return allExpSpec2
End

Function/WAVE makeLongSpecWave(expSpecName)
	
	String expSpecName
	String listOfExpSpec   = SortList(WaveList(expSpecName   + "*",";",""),";",24)
	Make/O/D/N=0 allExpSpec = 0
	Concatenate/O/NP listOfExpSpec,allExpSpec
	
	return allExpSpec
End

Function/WAVE findBAenergy()
	
	String currentFolder=GetdataFolder(1)
	SetDataFolder root:Packages:NXA
	WAVE/Z mu_energy
	String newNameEn = currentFolder + "mu_energy"
	Duplicate/O mu_energy,$newNameEn
	SetDataFolder currentFolder
	
	return mu_energy
End

Function/WAVE findBAMA(mol)
	String mol
	
	String currentFolder=GetdataFolder(1)
	SetDataFolder root:Packages:NXA
	WAVE/Z mu=$(mol+"_mu")
	String newNameMu = currentFolder + "mu"
	Duplicate/O mu,$newNameMu
	
	If( !WaveExists(mu) )
		Abort "Couldn't find the indicated mass absorption wave! Aborting!"
	endif
	
	SetDataFolder currentFolder
	
	return mu
End

Function findWaveValAtEergy(xw,yw,val)
	
	Wave xw,yw
	Variable val
	
	Variable p1 = round(BinarySearchInterp(xw,val))
	Variable dftLo = yw[p1]
	
	return dftLo
End