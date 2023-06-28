#pragma rtGlobals=1		// Use modern global access method.
#include "CollinsProcs_Igor8_NSLSII"


//Combines two autostitch commands to 
Function AutostitchBL(sets,qpts,outstr,[qpwr,thickness,anchor])
	String Sets,qpts,outstr // Same from Autostitch, "String of sets, qpts to stich and output name
	Variable qpwr,thickness,anchor
	
	qpwr= ParamIsdefault(qpwr) ? 2 : qpwr
	thickness= ParamIsdefault(thickness) ? 1 : thickness
	anchor= ParamIsdefault(anchor) ? 0 : anchor

	
	
	Autostitch(Sets,qpts,outstr,quiet=1,sectsize=360,qpwr=qpwr,thickness=thickness,anchor=anchor)
	Autostitch(Sets,qpts,outstr,ani=1,d=1,qpwr=qpwr,thickness=thickness,anchor=anchor)
	
end



//Makes a duplicate image plot of the current CCDImageToConvertFig that won't change with new data processing
Function NIKA2Dduplicate([n,f,overwrite])
	String n //optional name of data.  Default=NIKA2D## using "UniqueName"
	String f //optional folder name to put the copied data
	Variable overwrite //overwrite if wave already exists
	String CurrentFolder=GetDataFolder(1)
	If( ParamIsDefault(f) )
		f=CurrentFolder
	endif
	If( ParamIsDefault(overwrite) )
		overwrite=0
	endif
	If( ParamIsDefault(n) )
		n="NIKA2D"
	elseif( CheckName(n,1) !=0 && overwrite==0 )
		n=UniqueName(n,1,0)
	endif
	DoWindow/F CCDImagetoConvertFig
	if( V_Flag != 1 )
		Print "Couldn't find Window!!! Abort!"
		return 0
	endif
	NewDataFolder/O/S $f
	String wn=f+":"+n
	SetDataFolder root:Packages:Convert2Dto1D
	Duplicate/o CCDImageToConvert_dis $wn
	DoWindow/F CCDImagetoConvertFig
	DoIgorMenu "Edit" "Duplicate"
	String winNm=UniqueName(n+"_win",6,0)
	DoWindow/C $winNm
	ReplaceWave/W=$winNm image=CCDImageToConvert_dis, $wn
	SetDataFolder $CurrentFolder
End


//Processes scattering profiles at various energies and two linear x-ray polarizations that have been processed by NIKA
Function/S ProcEscan(name, sets, qpts, [tagQ, tagE,Bkg,XRFeff,qPwr,thick,anchor,d,q,xFilter])
	String name, sets, qpts //name of scan, list of match strings to identify the datasets, qpt list for stitching
	Variable TagQ, TagE //extract the Q/Energy profile at the tagged Energy/Q
	String Bkg //folder containing abs. spectra of film and SiN to be used in background correction
	Variable XRFeff //XRF efficiency assuming a K-edge photon has been absorbed (for Carbon = 2.05e-4/sr)
	Variable qPwr, thick //intensity factors
	Variable anchor //optionally scale the stitched intensity traces not to the lowest qrange
	Variable d //display the data
	Variable q //don't print warning messages to history
	String xFilter //extra filter to rid of spurious processed data from NIKA (kind of a reverse match)
	clearTimers()
	Variable timer=StartMSTimer
	//defaults
	qPwr= ParamIsDefault(qPwr) ? 2 : qPwr
	Bkg= SelectString( ParamIsDefault(Bkg), Bkg, "")
	thick= ParamIsDefault(thick) ? 1 : thick
	d= ParamIsDefault(d) ? 1 : d
	q= ParamIsdefault(q) ? 1 : q
	xFilter= SelectString( ParamIsDefault(xFilter), xfilter, "")

	String CurrentFolder=GetDataFolder(1), dFldr="root:Escans:"+name, mFldr=dFldr+":Meta", procFldr=dFldr+":Proc", pFldr=dFldr+":Prof",wnote
	String sFldr=dFldr+":Stitch"
	NewDataFolder/o/s root:Escans
	NewDataFolder/o/s $dFldr
	NewDataFolder/o/s $procFldr
	NewDataFolder/o/s $pFldr
	NewDataFolder/o/s $mFldr
	
	qpts=AssembleScanInfo(sets,qpts,mFldr,10,xFilter=xfilter)  //assembles info on all the profiles for this scan processed by NIKA
	
	//Ferron Changes here, working on allowing stitching together different energies ranges within an identical Q range
	String EnergiesList = AssembleMetaWave("Energies",mFLDR)
	Wave 	Energies0,Pols0, Energy = UniqueVals_Set(EnergiesList,0.05,":Energy") //Energy is now a list of all possible energies within the data sets given
	Duplicate/o Energy $(dFLDR + ":Energy")
	//Returns identical results as of 1:54pm on 11/26/2018***********
	
	//WAVE Energies0, Pols0, Energy=UniqueVals(Energies0,0.05,wn=":Energy") //Set up unique energies (within 0.1 eV) //Ferron commented for line 90
	//END CHANGES HERE.
	
	SetDataFolder $procFldr
	
	//prepare backgrounds for subtraction
	WAVE RawAbs=$("root:NEXAFS:"+Bkg+":"+Bkg+"_absSub"), BkgE=$("root:NEXAFS:"+Bkg+":Energy")
	WAVE RawSiNAbs=$("root:NEXAFS:"+Bkg+":SiN_abs") //SiN abs to correct for reduced intensities
	Duplicate/d/o Energy, AbsEff, SiNtrans
	AbsEff=0; SiNtrans=1
	
	If( !ParamIsDefault(XRFeff) && WaveExists(RawAbs) && WaveExists(BkgE) && WaveExists(RawSiNabs) )
		AbsEff=1-exp(-interp(Energy,BkgE, RawAbs)) // convert to abs efficiency = 1-Transmission = (1-exp^-OD)
		SiNtrans=exp(-interp(Energy,BkgE,RawSiNAbs)) //convert to transmission = I/I0 = exp^(-OD)
	endif
	
	make/o/n=(0,0) HD_Mat //initialize all data matricies
//	sprintf wnote, "ScanNum:%g;XRFspec:%s;XRFamp:%g;qPwr:%g;thick:%g;", snum, Bkg, fact, qPwr, thick
//	Note/K/NOCR HD_mat, wnote
	Duplicate/d/o HD_Mat, HDu_Mat, VD_Mat, VDu_Mat, HaniMat, HaniUmat, VaniMat, VaniUmat //initialize 2D data Matricies
	Duplicate/d/o HD_Mat, HH_mat, HHu_Mat, HV_mat, HVu_mat, VH_mat, VHu_mat, VV_mat, VVu_mat
	
	Variable i, j, k, pol, maxQ, p0, pf, nsets=itemsInList(sets)
	Variable/C TSIcmplx, rng
	String bn, tg, pstr
	wavestats pols0
	For(i=0; i<numpnts(energy); i+=1 )		//cycle through energies
		For( j=0; j<2; j+=1 )				// cycle through polarizations
			IF( j==1 && V_max ==100 )
				continue 
			endif
			bn=pFldr+":"+name+SelectString(j,"_HP","_VP")
			pol= j==0 ? 100 : 190; pstr=SelectString(j,"H","V")
			k= nsets>1 ? RXS_StitchQranges(name,qpts,energy[i],pol,q,anchor) : 0 //stitches all q-ranges for all sectors at this En & pol
			RXS_CombineSectors(bn,energy[i],pol,mFldr,sFldr,itemsInlist(qpts)+1) //combine 8 sectors into 3 (Vertical, Horizontal, and Diagonal)
			For( k=-1; k<2; k+=1 )  //do the background subtraction
				tg=SelectString(k,"H","V","D")
				WAVE qw=$(bn+"_"+tg+"x"), dw=$(bn+"_"+tg+"y"), uw=$(bn+"_"+tg+"u") 
				dw/=SiNTrans[i] //account for the reduced incident intensity due to limited SiN window transmission
				uw/=SiNTrans[i]
				IF( energy[i]>283 ) 
					dw -= XRFeff*AbsEff[i]/(4*pi) //subtracting off XRF background
				endif
				If( i==0 && j==0 ) //inialize the 2D intensity matrix via the first energy (may need to change since q-range changes with energy)
					Duplicate/d/o qw, $("qOut"+tg), $("qOut"+tg+"i"), $("intProf"+tg), $("intProfU"+tg), $("AniProf"+tg), $("AniProfU"+tg)
					WAVE qOut=$("qOut"+tg), qOuti=$("qOut"+tg+"i") //rescale q-values to encompas q-range all energies will possibly see
					maxQ=qOut[numpnts(qOut)-1] * waveMax(energy)/wavemin(energy)
					setscale/i x, log(qOut[0]), log(maxQ), qOut
					qOut=10^x
					imageScale(qOuti) //incase user wants to display these matrices later
				endif
				WAVE qOut=$("qOut"+tg), intProf=$("intProf"+tg), intProfU=$("intProfU"+tg), AniProf=$("AniProf"+tg), AniProfU=$("AniProfU"+tg)
				IntProf= Range(qw, qOut[p]) ? interp(qOut[p],qw,dw) : nan
				IntProfU= Range(qw, qOut[p]) ? interp(qOut[p],qw,uw) : nan
				WAVE intMat=$(pstr+tg+"_Mat"), intUmat=$(pstr+tg+"u_Mat")
				Concatenate {intProf}, intMat
				Concatenate {intProfU}, intUmat
			endfor //ends with the diagonal waves loaded which are used for the average profiles (max q-range, min CCD edge-effects)
			CalcAni(bn,useD=1) //calculates the anisotropy from the three primary sectors
			WAVE qw=$(bn+"_aniX"), dw=$(bn+"_aniY"), uw=$(bn+"_aniU")
			AniProf= Range(qw,qOut[p]) ? interp(qOut,qw,dw) : nan
			aniProfU= Range(qw,qOut[p]) ? interp(qOut,qw,uw) : nan
			WAVE aniMat=$(pstr+"aniMat"), aniUmat=$(pstr+"aniUmat")
			Concatenate {AniProf},aniMat
			Concatenate {AniProfU}, aniUmat
		endfor
	endfor
	
	//access the combined data matricies
	WAVE HD_Mat, HDu_Mat, VD_Mat, VDu_Mat, HaniMat, HaniUmat, VaniMat, VaniUmat, qOutH, qOutV, qOutD
	qOutH*=10; qOutV*=10; qOutD*=10 //convert from A to nm
	
	//most important waves to parent folder
	SetDataFolder $dFldr
	Duplicate/d/o HD_Mat, Intensity, IntensityU, Anisotropy, AnisotropyU
	Duplicate/d/o Energy, eImg, TSI, TSIu, TSIani, TSIaniU //output Intensity waves
	Duplicate/d/o qOut, qImg, qVals
	ImageScale(qImg) //make X & Y waves against which to plot the matrix data
	ImageScale(eImg)
	Make/o/n=0 nanWave=nan
	Variable/G ProfEn, SpecQ
	SetDataFolder $procFldr
	
	//combine results from the two polarizations & calc TSI
	Intensity= real(wAvg( {HD_Mat,VD_Mat},{HDu_mat,VDu_mat} ))/thick  //uses weighted averaging based on uncertainty
	WAVE uncert=combineUncert( HD_mat, VD_mat, HDu_mat, VDu_mat); IntensityU=uncert/thick
//	IntensityU= (HDu_mat+VDu_mat)/thick //Alternative if V & H are very different (occurs if anisotropy not sinusoidal)
	Anisotropy= real(wAvg( {HaniMat,-VaniMat},{HaniUmat,VaniUmat} ))
	WAVE uncert=combineUncert( HaniMat,VaniMat,HaniUmat,VaniUmat,f=-1); AnisotropyU=uncert
	RXS_CalcTSI(name,thick) //calculates the TSI and TSI-based anisotropy waves
		
	//Add Q-scaling to data
	Intensity*= qOut[p]^qPwr * (qPwr==0 && thick!=1 ? 1e7 : 1) //convert thickness from nm into cm if there is no Q-scaling
	IntensityU*=qOut[p]^qPwr * (qPwr==0 && thick!=1 ? 1e7 : 1) //convert thickness from nm into cm if there is no Q-scaling
	Intensity= Intensity<0 ? nan : intensity //rid any negative intensity data (due to over subtraction of backgrounds)
	
	//calculates the profile and spectrum slices of the data
	DoWindow $(name+"_data"); i=10; j=0
	If( V_Flag!=0 ) //handle preexisting cursors
		bn=CsrInfo(A,name+"_Data")
		i = NumberByKey("POINT",bn)
		j = NumberByKey("YPOINT",bn)
	endif
	k= ParamIsDefault(tagQ) ? RXS_TagQ(name,qVals[i]) : RXS_TagQ(name,tagQ)
	k= ParamIsDefault(tagE) ? RXS_TagE(name,energy[j]) : RXS_TagE(name,tagE)
	If( d )  //displays the data
		DoWindow $(name+"_Layout")
		If( V_Flag==0 && d==1 )
			NewLayout/N=$(name+"_Layout") as name+" Layout"
			PrintSettings/I margins={0.5,0.5,0.5,0.5}
		endif
		RXS_DisplayMat(name,csrX=binarySearch(qVals, (tagQ==0 ? qVals[10] : tagQ) ),csrY=binarySearch(energy,(tagE==0 ? energy[0] : tagE)) )
		RXS_DisplayProfile(name)
		RXS_DisplaySpectra(name)
		RXS_DisplayTSI(name)
	endif
	RXS_UpdateLabels(name,qpwr,!ParamIsDefault(thick))
	SetDataFolder $CurrentFolder
	string returnStr
	sprintf returnStr, "Took %1.2f seconds." , stopMStimer(timer)*1e-6
	return returnStr
End	


//Stitches data from multple Q-ranges (CCD images) given an energy and polarization
Function RXS_StitchQranges(dn,qpts,en,pol,q,anchor)
	String dn, qpts
	Variable en,pol //current energy and polarization
	Variable q //don't report stitch warnings
	Variable anchor //which set should the intensity be anchored to?
	String CurrentFolder=GetDataFolder(1), dFldr="root:Escans:"+dn, mFldr=dFldr+":Meta", procFldr=dFldr+":Proc", pFldr=dFldr+":Prof",pStr
	String sFldr=dFldr+":Stitch", fStr
	NewDataFolder/o $sFldr
	SetDataFolder $mFldr
	Variable pair, i, nSectors, qpt, R, uR, qRange
	WAVE/T Scans0, Notes //load all the metadata wave lists made by "AssembleScanInfo"
	WAVE Energies0, Pols0, Angles0, Dwells0
	SetDataFolder $procFldr
	WAVE qOutD //Below: determine conservative lowest Q-val for case of missing energy in lowest q-range dataset.
	Variable q0= waveExists(qOutD) ? CommonRange("qOutH;qOutV;qOutD") : 0.0005  //special case: missing first energy of first q-range
			
	//Cycle through the stitch pairs
	For( pair=1; pair<=itemsInList(qpts); pair+=1 )
		SetDataFolder $mFldr //load a pair of meta data
		pStr=num2str(pair)
		WAVE/Z/T ScansN=$("Scans"+pStr)
		WAVE/Z EnergiesN=$("Energies"+pStr), PolsN=$("Pols"+pStr), AnglesN=$("Angles"+pStr), DwellsN=$("Dwells"+pStr)
		WAVE/Z/C ratio=$("StitchRatios"+num2str(pair))
		Make/o/n=8 sectors={0,45,90,135,180,225,270,315}
//		Make/o/n=4 sectors={45,135,225,315}
//		WAVE/Z sectors=RXS_compileSectors(dn,en,pol) //returns wave with all appropriate sectors at en
		Duplicate/d/o sectors, rTerms, uTerms, idx0, idxN
		rTerms=nan; uTerms=nan; idx0=nan; idxN=nan
		qpt=Str2Num(stringfromList(pair-1,qpts))
		//cycle through the sectors & calculate the stitch ratio
		For( i=0; i<numpnts(sectors); i+=1 )
			idx0[i]=QueryMeta(sectors[i],En,pol,mfldr) //store lookup table of matches
			idxN[i]=QueryMeta(sectors[i],En,pol,mfldr,tg=pStr)
			If( pair==1 )
				fStr="root:SAS:'"+scans0[idx0[i]]
				WAVE/Z y1=$(fstr+"':'r_"+scans0[idx0[i]]+"'"), x1=$(fstr+"':'q_"+scans0[idx0[i]]+"'"), u1=$(fstr+"':'s_"+scans0[idx0[i]]+"'")
			else //already stitched two together and created new profile waves
				WAVE/Z y1=$(sFldr+":y_"+num2str(sectors[i])), x1=$(sFldr+":x_"+num2str(sectors[i])), u1=$(sFldr+":u_"+num2str(sectors[i]))
			endif
			fStr="root:SAS:'"+scansN[idxN[i]] //Moved into If statement to check if QueryMeta Failed
			//WAVE/Z  y2=$(fstr+"':'r_"+scansN[idxN[i]]+"'"), x2=$(fstr+"':'q_"+scansN[idxN[i]]+"'"), u2=$(fstr+"':'s_"+scansN[idxN[i]]+"'") //Second set found appropriate sector/energy/thing
			//WAVE/Z  y2=$(fstr+"':'r_"+scans0[idx0[i]]+"'"), x2=$(fstr+"':'q_"+scans0[idx0[i]]+"'"), u2=$(fstr+"':'s_"+scans0[idx0[i]]+"'")

//
			if(!QueryMetaFailCheck(idxN)) //Check to see if one of the QueryMetas fails based on the initial parameters (skipped energies or sectors)
				fStr="root:SAS:'"+scans0[idx0[i]]
				WAVE/Z y2=$(fstr+"':'r_"+scans0[idx0[i]]+"'"), x2=$(fstr+"':'q_"+scans0[idx0[i]]+"'"), u2=$(fstr+"':'s_"+scans0[idx0[i]]+"'") //Second set does not contain something
			else
				fStr="root:SAS:'"+scansN[idxN[i]]
				WAVE/Z  y2=$(fstr+"':'r_"+scansN[idxN[i]]+"'"), x2=$(fstr+"':'q_"+scansN[idxN[i]]+"'"), u2=$(fstr+"':'s_"+scansN[idxN[i]]+"'") //Second set found appropriate sector/energy/thing
			endif
//			
			IF( !WaveExists(y1) || !waveExists(y2) )
				continue //if sectors don't line up
			elseIf( !Range(x1,qpt,tol=3) )
				If( !q )
					Printf "Stitch Warning: %g out of range for %s [%1.4f,%1.4f].\r", qpt, NameOfWave(x1) ,x1[0],x1[inf]
				endif
				continue
			elseif( !Range(x2,qpt,tol=3) )
				If( !q )
					Printf "Stitch Warning: %g out of range for %s [%1.4f,%1.4f]\r", qpt, NameOfWave(x2),x2[0],x2[inf]
				endif
				continue
			endif
			rTerms[i]=interp( qpt, x1, y1) / interp(qpt, x2, y2) //gives up to 8 votes on what the stitch ratio should be
			uTerms[i]=rTerms[i]*sqrt( (interp(qpt,x1,u1)/interp(qpt,x1,y1))^2 + (interp(qpt,x2,u2)/interp(qpt,x2,y2))^2 )
		Endfor
		wavestats/Q rTerms
		If( V_npnts ) //at least one sector to stitch
			Ratio[idxN[0]]=wAvg(rTerms,uTerms) //store ratio in the first sector
		else
			printf "StitchRatio Error for set %g at E = %1.1f eV and P = %d\r", pair, en, pol
		endif
		//cycle back through the sectors stitching the pair together based on the known stitch ratio
		For( i=0; i<numpnts(sectors); i+=1 )
			R=real(Ratio[idxN[0]]); uR=imag(Ratio[idxN[0]])
			If( pair==1 )
				fStr="root:SAS:'"+scans0[idx0[i]]
				WAVE y1=$(fstr+"':'r_"+scans0[idx0[i]]+"'"), x1=$(fstr+"':'q_"+scans0[idx0[i]]+"'"), u1=$(fstr+"':'s_"+scans0[idx0[i]]+"'")
			else //already stitched two together and created new profile waves
				WAVE y1=$(sFldr+":y_"+num2str(sectors[i])), x1=$(sFldr+":x_"+num2str(sectors[i]))
				WAVE u1=$(sFldr+":u_"+num2str(sectors[i])), ua1=$(sFldr+":ua_"+num2str(sectors[i]))
			endif
			fStr="root:SAS:'"+scansN[idxN[i]]
			WAVE  y2=$(fstr+"':'r_"+scansN[idxN[i]]+"'"), x2=$(fstr+"':'q_"+scansN[idxN[i]]+"'"), u2=$(fstr+"':'s_"+scansN[idxN[i]]+"'")
			If( WaveExists(y1) && pair==1 )
				Duplicate/d/o y1, $(sFldr+":y_"+num2str(sectors[i]))
				Duplicate/d/o x1, $(sFldr+":x_"+num2str(sectors[i]))
				Duplicate/d/o u1, $(sFldr+":u_"+num2str(sectors[i])), $(sFldr+":ua_"+num2str(sectors[i]))
			elseif( !WaveExists(y1) ) //will automatically be pair=1
				make/o/d/n=2 $(sFldr+":y_"+num2str(sectors[i])), $(sFldr+":x_"+num2str(sectors[i])), $(sFldr+":u_"+num2str(sectors[i])), $(sFldr+":ua_"+num2str(sectors[i]))
				WAVE/Z y1=$(sFldr+":y_"+num2str(sectors[i])), x1=$(sFldr+":x_"+num2str(sectors[i]))
				WAVE/Z u1=$(sFldr+":u_"+num2str(sectors[i])), ua1=$(sFldr+":ua_"+num2str(sectors[i]))
				y1=nan; x1={q0,qpt}; u1=nan;	ua1=nan //initialize with 2-point nan data if missing the first set.
			elseif( !WaveExists(y2) )
				make/o/d/n=1 dmmy, dmmyX
				dmmy=nan; dmmyX=qpt
				WAVE y2=dmmy, x2=dmmyX, u2=dmmy, ua2=dmmy
			endif
			WAVE yOut=$(sFldr+":y_"+num2str(sectors[i])), xOut=$(sFldr+":x_"+num2str(sectors[i]))
			WAVE uOut=$(sFldr+":u_"+num2str(sectors[i])), uaOut=$(sFldr+":ua_"+num2str(sectors[i]))

			Variable pt1= BinarySearch(x1,qpt), pt2= BinarySearch(x2,qpt)
			pt1= pt1<0 ? numpnts(x1) : pt1
			pt2= pt2<0 ? 0 : pt2
			Variable nPts= pt1 + numpnts(x2) - pt2
			Redimension/n=(npts) yOut, xOut, uOut, uaOut
			//Fill from the second wave
			xOut[pt1,npts-1]=x2[p-pt1+pt2]
			If( anchor>=pair )
				yOut/=R
				uOut=yOut * sqrt( (uOut/yOut)^2 + (uR/R)^2 )
				uaOut=uaOut/R
				yOut[pt1,npts-1]=y2[p-pt1+pt2]
				uOut[pt1,npts-1]=u2[p-pt1+pt2]//yOut * sqrt( (u2[p-pt1+pt2]/y2[p-pt1+pt2])^2 + (uR/R)^2 )
				uaOut[pt1,npts-1]=u2[p-pt1+pt2]
			else
				yOut[pt1,npts-1]=y2[p-pt1+pt2]*R
				uOut[pt1,npts-1]=yOut * sqrt( (u2[p-pt1+pt2]/y2[p-pt1+pt2])^2 + (uR/R)^2 )
				uaOut[pt1,npts-1]=u2[p-pt1+pt2]*R
			endif
		endfor
	endfor
	SetDataFolder $CurrentFolder
end

//returns a wave with all the appropriate sectors given a data pair, energy and polarization and location of the metadata
Function/WAVE RXS_compileSectors(dn,en,pol)
	String dn//data name
	Variable en,pol //energy and polarization
	make/d/o/n=0 AllAngles
	variable i
	string bn= "root:Escans:"+dn+":Meta:Angles"
	For( i=0; ; i+=1 )
		WAVE Angles=$(bn+num2str(i)), Energies=$(bn+num2str(i)), Pols=$(bn+num2str(i))
		If( !waveExists(angles) )
			break
		endif
		Concatenate {Angles}, AllAngles
	Endfor
	WAVE Sectors=UniqueVals(AllAngles,1,wn="Sectors")
	return Sectors
end


//combine sectors (hozontal, vertical and 45) if multiple exist, and rename
//can handle one or two H & V sectors and up 1-4 45 deg sectors
Function RXS_CombineSectors(bname, en,pol,mf,sf,nRanges)//(bName, th0, n,[en,pol,mf,outf])
	String bname
	Variable en, pol //get the waves directly from the SAS folder rather than prev. processed data in current folder
	string mf, sf //mdatadata folder and sector folder (folders where stitched sectors are recorded)
	Variable nRanges
	String Currentfolder=GetDataFolder(1)

	//Combine vertical sectors based on availablility and user selection
	If( nRanges<2 )  //didn't stitch anything yet so all data is still in SAS folder
		WAVE/Z x90=$SASref(en,pol,90,1,fldr=mf), y90=$SASref(en,pol,90,0,fldr=mf), u90=$SASref(en,pol,90,-1,fldr=mf)
		WAVE/Z x270=$SASref(en,pol,270,1,fldr=mf), y270=$SASref(en,pol,270,0,fldr=mf), u270=$SASref(en,pol,270,-1,fldr=mf)
	else
		WAVE/Z x90=$(sf+":x_90"), y90=$(sf+":y_90"), u90=$(sf+":u_90")
		WAVE/Z x270=$(sf+":x_270"), y270=$(sf+":y_270"), u270=$(sf+":u_270")
	endif
	If( !WaveExists(x90) && !WaveExists(x270) )
		print "Couldn't find any horizontal profiles! ("+bName+"_90_x OR "+bname+"_270_x)"
		return -1
	elseif( !WaveExists(x270) )
		Duplicate/d/o x90 $(bName+"_Vx"); Duplicate/d/o y90 $(bName+"_Vy"); Duplicate/d/o u90 $(bName+"_Vu")
	elseif( !WaveExists(x90) )
		Duplicate/d/o x270 $(bName+"_Vx"); Duplicate/d/o y270 $(bName+"_Vy"); Duplicate/d/o u270 $(bName+"_Vu")
	else
		CombineProfiles(y90,x90,u90,y270,x270,u270)
		WAVE AvX, AvY, AvU
		Duplicate/d/o AvX $(bName+"_Vx"); Duplicate/d/o AvY $(bName+"_Vy"); Duplicate/d/o Avu $(bName+"_Vu")
	endif
	
	//Combine horizontal sectors based on availablility and user selection
	If( nRanges<2 )
		WAVE/Z x0=$SASref(en,pol,0,1,fldr=mf), y0=$SASref(en,pol,0,0,fldr=mf), u0=$SASref(en,pol,0,-1,fldr=mf)
		WAVE/Z x180=$SASref(en,pol,180,1,fldr=mf), y180=$SASref(en,pol,180,0,fldr=mf), u180=$SASref(en,pol,180,-1,fldr=mf)
	else
		WAVE/Z x0=$(sf+":x_0"), y0=$(sf+":y_0"), u0=$(sf+":u_0")
		WAVE/Z x180=$(sf+":x_180"), y180=$(sf+":y_180"), u180=$(sf+":u_180")
	endif
	If( !WaveExists(x0) && !WaveExists(x180) )
		print "Couldn't find any vertical profiles! ("+bName+"_0_x OR "+bname+"_180_x)"
		return -1
	elseif( !WaveExists(x180) )
		Duplicate/d/o x0 $(bName+"_Hx"); Duplicate/d/o y0 $(bName+"_Hy"); Duplicate/d/o u0 $(bName+"_Hu")
	elseif( !WaveExists(x0) )
		Duplicate/d/o x180 $(bName+"_Hx"); Duplicate/d/o y180 $(bName+"_Hy"); Duplicate/d/o u180 $(bName+"_Hu")
	else
		CombineProfiles(y0,x0,u0,y180,x180,u180)
		Duplicate/d/o AvX $(bName+"_Hx"); Duplicate/d/o AvY $(bName+"_Hy"); Duplicate/d/o AvU $(bName+"_Hu")
	endif
	
	//Combine diagonal sectors based on availablility and user selection
	If( nRanges<2 )
		WAVE/Z x45=$SASref(en,pol,45,1,fldr=mf), y45=$SASref(en,pol,45,0,fldr=mf), u45=$SASref(en,pol,45,-1,fldr=mf)
		WAVE/Z x135=$SASref(en,pol,135,1,fldr=mf), y135=$SASref(en,pol,135,0,fldr=mf), u135=$SASref(en,pol,135,-1,fldr=mf)
		WAVE/Z x225=$SASref(en,pol,225,1,fldr=mf), y225=$SASref(en,pol,225,0,fldr=mf), u225=$SASref(en,pol,225,-1,fldr=mf)
		WAVE/Z x315=$SASref(en,pol,315,1,fldr=mf), y315=$SASref(en,pol,315,0,fldr=mf), u315=$SASref(en,pol,315,-1,fldr=mf)
	else
		WAVE/Z x45=$(sf+":x_45"), y45=$(sf+":y_45"), u45=$(sf+":u_45")
		WAVE/Z x135=$(sf+":x_135"), y135=$(sf+":y_135"), u135=$(sf+":u_135")
		WAVE/Z x225=$(sf+":x_225"), y225=$(sf+":y_225"), u225=$(sf+":u_225")
		WAVE/Z x315=$(sf+":x_315"), y315=$(sf+":y_315"), u315=$(sf+":u_315")
	endif
	If( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x225) && !WaveExists(x315))
		print "Couldn't find any 45 deg profiles! ("+bName+"_45_x, OR 135, 225, 315)"
		return -1
	elseif( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x225) ) //only one
		Duplicate/d/o x315 $(bName+"_Dx"); Duplicate/d/o y315 $(bName+"_Dy"); Duplicate/d/o u315 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x225) && !WaveExists(x315) ) //only one
		Duplicate/d/o x135 $(bName+"_Dx"); Duplicate/d/o y135 $(bName+"_Dy"); Duplicate/d/o u135 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x315) ) //only one
		Duplicate/d/o x225 $(bName+"_Dx"); Duplicate/d/o y225 $(bName+"_Dy"); Duplicate/d/o u225 $(bName+"_Du"); return 0
	elseif( !WaveExists(x135) && !WaveExists(x225) && !WaveExists(x315) )//only one
		Duplicate/d/o x45 $(bName+"_Dx"); Duplicate/d/o y45 $(bName+"_Dy"); Duplicate/d/o u45 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x135) ) //only two
		CombineProfiles(y225, x225,u225, y315,x315,u315)
	elseif( !WaveExists(x45) && !WaveExists(x225) ) //only two
		CombineProfiles(y135, x135, u135, y315,x315,u315)
	elseif( !WaveExists(x45) && !WaveExists(x315) ) //only two
		CombineProfiles(y135, x135, u135, y225,x225,u225)
	elseif( !WaveExists(x135) && !WaveExists(x225) ) //only two
		CombineProfiles(y45, x45, u45,y315,x315,u315)
	elseif( !WaveExists(x135) && !WaveExists(x315) ) //only two
		CombineProfiles(y45, x45, u45,y225,x225,u225)
	elseif( !WaveExists(x225) && !WaveExists(x315) ) //only two
		CombineProfiles(y45, x45, u45,y135,x135,u135)
	elseif( !WaveExists(x45) ) //only three
		CombineProfiles(y135,x135,u135,y225,x225,u225)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU, y315,x315,u315)
	elseif( !WaveExists(x135) ) //only three
		CombineProfiles(y45,x45,u45,y225,x225,u225)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y315,x315,u315)
	elseif( !WaveExists(x225) ) //only three
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y315,x315,u315)
		elseif( !WaveExists(x315) ) //only three
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y225,x225,u225)
	else //all 4 are present
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX1; Duplicate/d/o AvY tmpY1; Duplicate/d/o AvU tmpU1
		CombineProfiles(y225,x225,u225,y315,x315,u315)
		Duplicate/d/o AvX tmpX2; Duplicate/d/o AvY tmpY2; Duplicate/d/o AvU tmpU2
		CombineProfiles(tmpY1,tmpX1,tmpU1,tmpY2,tmpX2,tmpU2)
	endif
	Duplicate/d/o AvX $(bName+"_Dx"); Duplicate/d/o AvY $(bName+"_Dy"); Duplicate/d/o AvU $(bName+"_Du")
	
	SetDataFolder $CurrentFolder
End


//Calculates the TSI(E) and TSIani(E) spectra from the RXS data
Function RXS_CalcTSI(dn,thick)
	string dn
	Variable thick
	String CurrentFolder=GetDataFolder(1), df="root:Escans:"+dn
	SetdataFolder $df
	WAVE Energy, TSI, TSIu, TSIani, TSIaniU
	SetDataFolder :proc
	WAVE HD_Mat, HDu_Mat, VD_Mat, VDu_Mat, HaniMat, HaniUmat, VaniMat, VaniUmat
	WAVE HH_mat, HHu_Mat, HV_mat, HVu_mat, VH_mat, VHu_mat, VV_mat, VVu_mat
	Duplicate/d/o HD_Mat, LorInt, LorIntU, TotLorInt, TotLorIntU
	Duplicate/d/o Energy, HH_TSI, HHu_TSI, HV_TSI, HVu_TSI, VH_TSI, VHu_TSI, VV_TSI, VVu_TSI //for calc TSI Anisotropy
	Duplicate/d/o Energy, HD_TSI, HDu_TSI, VD_TSI, VDu_TSI, AniH, AniHu, AniV, AniVu //for cacl TSI Anisotropy
	Duplicate/d/o Energy, TSIH, TSIuH, TSIV, TSIuV
	WAVE qOutH, qOutV, qOutD
	Variable j, k
	Variable/C Lmt
	string pstr, tg

	//calc TSI
	For( j=0; j<2; j+=1 )
		pstr=SelectString(j,"H","V")
		WAVE intMat=$(pstr+"D_Mat"), intUmat=$(pstr+"Du_Mat"), qOut=$("qOutD")
		WAVE TSIspec=$("TSI"+pstr), TSIspecU=$("TSIu"+pstr)
		Duplicate/d/o intMat, LorInt, LorIntU, TotLorInt, TotLorIntU
		LorInt=intMat*qOut[p]^2 // part of the Jacobian for integrating sphereically symmetric scattering
		LorIntU = ( intUmat*qOut[p]^3 )^2
		Lmt=BAC_Integrate2Dnans(LorInt,qOut,TotLorInt)
		BAC_Integrate2Dnans(LorIntU,qOut,TotLorIntU) //uncertainty in the sum is the sum of the uncertainties
		totLorIntU /= qOut[numpnts(qOut)-1]-qOut[0] //renormalize the uncertainties over the integration range
		TSIspec=TotLorInt[imag(Lmt)][p]/(2*pi^2)/Thick  //TSI is the high-limit Q & add rest of Jacobian & opt. thickness
		TSIspecU= sqrt(TotLorIntU[imag(Lmt)][p] )/(2*pi^2)/thick //done in quadrature
	endfor
	TSI= real(wAvg( {TSIH,TSIV}, {TSIuH,TSIuV} ))  //weighted average of the TSI's from the two polarizations
	combineUncert( TSIH,TSIV,TSIuH,TSIuV )
	WAVE uncert
	TSIu= max( uncert, TSI/sqrt( 2^16/energy*10) ) //limit uncertainties to the maximum CCD single-pixel photon-counting statistics (~2300)
	
	//calc TSIani
	String wvStr="HH_mat;HV_mat;HD_mat", xWvStr="qOutH;qOutV;qOutD;"
	Variable/C cRange=nanRange(wvStr,xWvStr) //range for the TSI-based anisotropy data
	Variable p0, pF
	For( j=0; j<2; j+=1 )  //cycle through polarizations
		pstr=SelectString(j,"H","V")
		For( k=-1; k<2; k+=1 ) //cycle through sectors
			tg=SelectString(k,"H","V","D")
			WAVE intMat=$(pstr+tg+"_Mat"), intUmat=$(pstr+tg+"u_Mat"), qOut=$("qOut"+tg)
			WAVE TSIspec=$(pstr+tg+"_TSI"), TSIspecU=$(pstr+tg+"u_TSI")
			Duplicate/d/o intMat, LorInt, LorIntU, TotLorInt, TotLorIntU
			LorInt= intMat*qOut[p]^2
			LorIntU= ( intUmat*qOut[p]^3 )^2
			BAC_Integrate2Dnans(LorInt,qOut,TotLorInt)
			BAC_Integrate2Dnans(LorIntU,qOut,TotLorIntU) //uncertainty in the sum is the sum of the uncertainties
			p0=BinarySearch(qOut, real(cRange)); pf=BinarySearch(qOut, imag(cRange))
			TotLorIntU /= qOut[pf]-qOut[p0] //renormalize uncertainties over integration range
			TSIspec = TotLorInt[pf][p] - TotLorInt[p0][p] //Save the integral over the common range
			TSIspecU = sqrt(TotLorIntU[pf][p] - TotLorIntU[p0][p])
		endfor
		WAVE H_TSI=$(pstr+"H_TSI"), Hu_TSI=$(pstr+"Hu_TSI"), V_TSI=$(pstr+"V_TSI"), Vu_TSI=$(pstr+"Vu_TSI")
		WAVE AniPol=$("Ani"+pstr), AniUpol=$("Ani"+pstr+"u")
		AniPol= (H_TSI-V_TSI)/(H_TSI+V_TSI)*100 //calculation of anisotropy
		AniUpol= sqrt( ( 2*H_TSI*Vu_TSI )^2 + ( 2*V_TSI*Hu_TSI )^2 ) / (H_TSI+V_TSI)^2 * 100
	endfor
	TSIani= real(wAvg( {-AniH,AniV},{AniHu,AniVu} )) //combining anisotropy measurement from two polarizations
	combineUncert( AniH,AniV,AniHu,AniVu,f=-1)
	TSIaniU=uncert
	
	SetDataFolder $CurrentFolder
end

//Finds the comon range for a list of 2D waves where there are nans in the beginning/end of the x-direction (rows)
Function/C nanRange(wList,xwList)
	String wList, xwList
	Variable i, neg=-inf, pos=inf, np, nq
	For( i=0; i<itemsInList(wList); i+=1 )
		WAVE w=$stringfromlist(i,wList), xw=$StringFromList(i,xwList)
		Duplicate/d/o w, nanLoc
		np=dimSize(w,0); nq=dimSize(w,1)
		NanLoc= NumType(w)==2 ? p : nan  //fill matrix with p-values where nans were originally (rest fill with nans)
		wavestats/Q nanLoc
		If( V_npnts<1 )
			continue
		endif
		ImageStats/G={0,(np-1)/2,0,nq-1} nanLoc  //search first half of the matrix
		neg= xw[V_max+1]>neg ? xw[V_max+1] : neg //x-value of the first row without nans (in the first half of matrix)
		ImageStats/G={(np-1)/2,np-1,0,nq-1} nanLoc  //search last half of the matrix
		pos= xw[V_min-1]<pos ? xw[V_min-1] : pos //x-value of the first row with nans (in the last half of matrix)
	endfor
	return cmplx(neg,pos)
end


//Integrates a 2D wave in the x-direction (rows) ignoring any nans at the beginning or end of the x-dimension
// such that the x-integration (row) range is the same for all y-values (columns)
Function/C BAC_Integrate2Dnans(inw,xw,outw)
	WAVE inW, xw, outw
	Variable i, j, np=DimSize(inW,0), nq=DimSize(inW,1), p0, pf
	Duplicate/d/o inW, nanLoc, inCrop, outCrop
	Duplicate/d/o xw, newX
	nanLoc= NumType(inW)==2 ? p : nan  //fill matrix with p-values where nans were originally (rest fill with nans)
	wavestats/Q nanLoc
	If( V_npnts>0 )
		ImageStats/G={0,(np-1)/2,0,nq-1} nanLoc  //search first half of the matrix
		p0=V_max+1 //p-value of the first row without nans (in the first half of matrix)
		ImageStats/G={(np-1)/2,np-1,0,nq-1} nanLoc  //search last half of the matrix
		pf=V_min //p-value of the first row with nans (in the last half of matrix)
		DeletePoints/M=0 pf, np-pf, inCrop, outCrop, newX
		Deletepoints/M=0 0, p0, inCrop, outCrop, newX
		outCrop=0; outw=0
		Integrate/T/Dim=0 inCrop /X=newX /D=outCrop
		outw[p0,pf-1]=outCrop[p-p0][q]
		return cmplx(p0,pf-1)
	else
		Integrate/T/Dim=0 inw /X=xw /D=outw
		return cmplx(0,numpnts(xw)-1)  //return the range of good data
	endif
end


//Extracts a tagged scattering and anisotropy profile from processed data given an energy
Function RXS_TagE(dn,En)
	String dn //data name
	Variable en
	String CurrentFolder=GetDataFolder(1), dFldr="root:Escans:"+dn
	SetdataFolder $dFldr
	WAVE Intensity, IntensityU, Anisotropy, AnisotropyU, Energy, qVals
	Duplicate/d/o qVals IntProf, IntProfU, AniProf, AniProfU
	Variable pt=BinarySearch(energy,en)
	IntProf=Intensity[p][pt]
	IntProfU=IntensityU[p][pt]
	AniProf=Anisotropy[p][pt]
	AniProfU=AnisotropyU[p][pt]
end


//Extracts a tagged scattering & anisotropy spectrum from processed data given a Q-value
Function RXS_TagQ(dn,Q)
	String dn //data name
	Variable Q
	String CurrentFolder=GetDataFolder(1), dFldr="root:Escans:"+dn
	SetdataFolder $dFldr
	WAVE Intensity, IntensityU, Anisotropy, AnisotropyU, Energy, qVals
	Duplicate/d/o Energy IntSpec, IntSpecU, AniSpec, AniSpecU
	Variable pt=BinarySearch(qVals,q)
	IntSpec=Intensity[pt][p]
	IntSpecU=IntensityU[pt][p]
	AniSpec=Anisotropy[pt][p]
	AniSpecU=AnisotropyU[pt][p]
end

//Extracts a slice of RXS data based on a cursor position in the 2D data graph
Function RXS_csrHook(s)
	STRUCT WMWinHookStruct &s
	String CurrentFolder=GetDataFolder(1), csrinfoStr
	Variable Xpt, Ypt
	If( (s.eventCode==7 ) && s.pointNumber*0==0 ) //Cursor moved to valid point
		String fn=GetWavesDataFolder(ImageNameToWaveRef(s.winName,s.traceName),1)
		setDataFolder $fn
		WAVE Intensity, IntensityU, Anisotropy, AnisotropyU, Energy, qVals
		NVAR profEn, specQ
		Duplicate/d/o Energy IntSpec, IntSpecU, AniSpec, AniSpecU
		specQ=qVals[s.pointNumber]
		IntSpec=Intensity[s.pointNumber][p]
		IntSpecU=IntensityU[s.pointNumber][p]
		AniSpec=Anisotropy[s.pointNumber][p]
		AniSpecU=AnisotropyU[s.pointNumber][p]
		Duplicate/d/o qVals IntProf, IntProfU, AniProf, AniProfU
		profEn=Energy[s.YpointNumber]
		IntProf=Intensity[p][s.YpointNumber]
		IntProfU=IntensityU[p][s.YpointNumber]
		AniProf=Anisotropy[p][s.YpointNumber]
		AniProfU=AnisotropyU[p][s.YpointNumber]
		doupdate
	endif
	SetDataFolder $CurrentFolder
end


//Displays raw 2D Intensity and Anisotropy data
Function RXS_DisplayMat(dn,[csrX,csrY]) : Graph
	String dn
	Variable csrX, csrY
	String CurrentFolder= GetDataFolder(1), df, wn
	df= "root:Escans:"+dn
	wn= dn+"_Data"
	DoWindow/F $wn
	IF( V_Flag )
		return 0
	endif
	SetDataFolder $df
	NVAR ProfEn, SpecQ
	WAVE nanWave, Intensity, Anisotropy, qVals, energy, qImg, eImg
	ProfEn=energy[csrY]; SpecQ=qVals[csrX]
	Display /W=(10,40,270,340)/B=middle/N=$wn nanWave vs qVals as replaceString("_",wn," ")
	AppendImage/W=$wn Intensity vs {qImg,eImg}
	ModifyImage/W=$wn Intensity ctab= {*,*,Terrain,0}, log= 1
	AppendImage/W=$wn/L=ani Anisotropy vs {qImg,eImg}
	ModifyImage/W=$wn Anisotropy ctab= {-20,20,RedWhiteBlue,0}
	ModifyGraph/W=$wn margin(left)=36,margin(bottom)=29,margin(top)=29,margin(right)=72
	ModifyGraph/W=$wn grid(left)=2,grid(bottom)=2,grid(ani)=2, log(bottom)=1,log(middle)=1
	ModifyGraph/W=$wn tick(left)=2,tick(bottom)=2,tick(ani)=2,tick(middle)=1
	ModifyGraph/W=$wn mirror=1,mirror(middle)=0,noLabel(middle)=2, mirror(bottom)=0
	ModifyGraph/W=$wn standoff=0, lblPosMode(ani)=1, lblPos(left)=36,lblPos(bottom)=27
	ModifyGraph/W=$wn freePos(ani)=0, freePos(middle)={0.5,kwFraction}
	ModifyGraph/W=$wn axisEnab(left)={0.5,1}, axisEnab(ani)={0,0.5}
	Label/W=$wn left "Photon Energy [eV]"
	Label/W=$wn bottom "Q [nm\\S-1\\M]"
	Label/W=$wn ani "Photon Energy [eV]"
	Cursor/W=$wn/P/I/S=2/H=1/N A Intensity csrX, csrY
	ColorScale/W=$wn/C/N=intScale/Z=1/A=RC/X=0.00/Y=20.00/E image=Intensity, heightPct=40
	ColorScale/W=$wn/C/N=intScale tickLen=3, log=1, lblMargin=15, logLTrip=0.1
	AppendText/W=$wn/N=intScale "Intensity*Q\\S2\\M [sr\\S-1\\Mnm\\S-3\\M]"
	ColorScale/W=$wn/C/N=aniScale/A=RC/X=0.00/Y=-20.00/E image=Anisotropy, heightPct=40
	AppendText/W=$wn/N=aniScale "Anisotropy [%]"
	SetWindow $wn,hook(cursorHook)=RXS_csrHook
	Add_2piOverQticks2top()
	BAC_appendToLayout(dn+"_Layout",wn,w=252,h=288)
	SetDataFolder $CurrentFolder
End

//Updates the various axis labels with the correct units based on the desired Q-scaling and thickness normalization
Function RXS_UpdateLabels(dn,qPwr,th)
	String dn
	Variable qPwr, th //th is not thickness but truth of whether user supplied a thickness
	String axisLbl
	If( qPwr>0 )
		sprintf axisLbl, "Intensity*Q\\S%d\\M [sr\\S-1\\Mnm\\S%d\\M]", qpwr, -qpwr-th
	else
		if( th )
			axisLbl="Intensity [sr\\S-1\\Mcm\\S-1\\M]"
		else
			axisLbl="Intensity [sr\\S-1\\M]"
		endif
	endif
	Dowindow $(dn+"_Data")
	If( V_Flag )
		ReplaceText/W=$(dn+"_Data")/N=intScale axisLbl
	endif
	Dowindow $(dn+"_Profiles")
	If( V_Flag )
		Label/Z/W=$(dn+"_Profiles") left axisLbl
	endif
	Dowindow $(dn+"_Spectra")
	IF( V_Flag )
		Label/Z/W=$(dn+"_Spectra") left axisLbl
	endif
	Dowindow $(dn+"_TSI")
	If( V_Flag )
		if( th )
			Label/Z/W=$(dn+"_TSI") left "TSI [nm\S-4\M]"
		else
			Label/Z/W=$(dn+"_TSI") left "TSI [nm\S-3\M]"
		endif
	endif
end


Function RXS_DisplayProfile(dn) : Graph
	String dn
	String CurrentFolder= GetDataFolder(1), df, wn
	df= "root:Escans:"+dn
	wn= dn+"_Profiles"
	DoWindow/F $wn
	IF( V_Flag )
		return 0
	endif
	SetDataFolder $df
	WAVE IntProf, IntProfU, AniProf, AniProfU, qVals, Energy, nanWave
	Display /W=(6,43.25,220.5,329.75)/N=$wn IntProf vs qVals as replaceString("_",wn," ")
	AutoPositionWindow
	AppendToGraph/W=$wn/L=ani AniProf vs qVals
	AppendToGraph/W=$wn/B=IntBtm nanWave vs qVals
	AppendToGraph/W=$wn/T=aniTop nanWave vs qVals
	ModifyGraph/W=$wn margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph/W=$wn mode=4, marker=19, msize=1, mskip=5, grid=2
	ModifyGraph/W=$wn log(left)=1,log(bottom)=1,log(IntBtm)=1,log(aniTop)=1
	ModifyGraph/W=$wn tick=2, zero(ani)=1, mirror=1, minor=1, noLabel(IntBtm)=2,noLabel(aniTop)=2, standoff=0
	ModifyGraph/W=$wn logHTrip(left)=100, logLTrip(left)=0.01
	ModifyGraph/W=$wn lblPosMode(left)=1,lblPosMode(ani)=1, lblPos(left)=38,lblPos(bottom)=28
	ModifyGraph/W=$wn freePos(ani)=0, freePos(IntBtm)={0.51,kwFraction}, freePos(aniTop)={0.51,kwFraction}
	ModifyGraph/W=$wn axisEnab(left)={0.51,1}, axisEnab(ani)={0,0.49}
	Label/W=$wn left "Intensity*Q\\S2\\M [sr\\S-1\\Mnm\\S-3\\M]"
	Label/W=$wn bottom "Q [nm\\S-1\\M]"
	Label/W=$wn ani "Anisotropy [%]"
	SetAxis/W=$wn ani -20,30
	ErrorBars/W=$wn/Y=4 IntProf Y,wave=(IntProfU,IntProfU)
	ErrorBars/W=$wn/Y=4 AniProf Y,wave=(AniProfU,AniProfU)
	TextBox/C/N=EnTag "Energy: \\{root:Escans:"+dn+":ProfEn} eV"
	BAC_appendToLayout(dn+"_Layout",wn,w=216,h=288)
	SetDataFolder $CurrentFolder
End


Function RXS_DisplaySpectra(dn) : Graph
	String dn
	String CurrentFolder= GetDataFolder(1), df, wn
	df= "root:Escans:"+dn
	wn= dn+"_Spectra"
	DoWindow/F $wn
	IF( V_Flag )
		return 0
	endif
	SetDataFolder $df
	WAVE IntSpec, IntSpecU, AniSpec, AniSpecU, qVals, Energy, nanWave
	Display /W=(6,43.25,220.5,329.75)/N=$wn IntSpec vs Energy as replaceString("_",wn," ")
	AutoPositionWindow
	AppendToGraph/W=$wn/L=ani AniSpec vs Energy
	AppendToGraph/W=$wn/B=IntBtm nanWave vs Energy
	AppendToGraph/W=$wn/T=aniTop nanWave vs Energy
	ModifyGraph/W=$wn margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph/W=$wn mode=0, grid=2
	ModifyGraph/W=$wn log(left)=1, tick=2, zero(ani)=1, mirror=1, minor=1, noLabel(IntBtm)=2,noLabel(aniTop)=2, standoff=0
	ModifyGraph/W=$wn logHTrip(left)=100, logLTrip(left)=0.01
	ModifyGraph/W=$wn lblPosMode(left)=1,lblPosMode(ani)=1, lblPos(left)=38,lblPos(bottom)=28
	ModifyGraph/W=$wn freePos(ani)=0, freePos(IntBtm)={0.51,kwFraction}, freePos(aniTop)={0.51,kwFraction}
	ModifyGraph/W=$wn axisEnab(left)={0.51,1}, axisEnab(ani)={0,0.49}
	Label/W=$wn left "Intensity*Q\\S2\\M [sr\\S-1\\Mnm\\S-3\\M]"
	Label/W=$wn bottom "Photon Energy [eV]"
	Label/W=$wn ani "Anisotropy [%]"
	SetAxis/W=$wn ani -20,30
	ErrorBars/W=$wn/Y=4 IntSpec Y,wave=(IntSpecU,IntSpecU)
	ErrorBars/W=$wn/Y=4 AniSpec Y,wave=(AniSpecU,AniSpecU)
	TextBox/C/N=QTag "\\{\"Q: %1.2e\",root:Escans:"+dn+":specQ} nm\\S-1\\M"
	BAC_appendToLayout(dn+"_Layout",wn,w=216,h=288)
	SetDataFolder $CurrentFolder
End


Function RXS_DisplayTSI(dn) : Graph
	String dn
	String CurrentFolder= GetDataFolder(1), df, wn
	df= "root:Escans:"+dn
	wn= dn+"_TSI"
	DoWindow/F $wn
	IF( V_Flag )
		return 0
	endif
	SetDataFolder $df
	WAVE TSI, TSIu, TSIani, TSIaniU, qVals, Energy, nanWave
	Display /W=(6,43.25,220.5,329.75)/N=$wn TSI vs Energy as replaceString("_",wn," ")
	AutoPositionWindow
	AppendToGraph/L=ani TSIani vs Energy
	AppendToGraph/B=IntBtm nanWave vs Energy
	AppendToGraph/T=aniTop nanWave vs Energy
	ModifyGraph/W=$wn margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph/W=$wn mode=0, grid=2
	ModifyGraph/W=$wn log(left)=1, tick=2, zero(ani)=1, mirror=1, minor=1, noLabel(IntBtm)=2,noLabel(aniTop)=2, standoff=0
	ModifyGraph/W=$wn logHTrip(left)=100, logLTrip(left)=0.01
	ModifyGraph/W=$wn lblPosMode(left)=1,lblPosMode(ani)=1, lblPos(left)=38,lblPos(bottom)=28
	ModifyGraph/W=$wn freePos(ani)=0, freePos(IntBtm)={0.51,kwFraction}, freePos(aniTop)={0.51,kwFraction}
	ModifyGraph/W=$wn axisEnab(left)={0.51,1}, axisEnab(ani)={0,0.49}
	Label left "TSI [nm\\S-4\\M]"
	Label/W=$wn bottom "Photon Energy [eV]"
	Label ani "TSI Anisotropy [%]"
	SetAxis ani -20,30
	ErrorBars/Y=4 TSI Y,wave=(TSIu,TSIu)
	ErrorBars/Y=4 TSIani Y,wave=(TSIaniU,TSIaniU)
	BAC_appendToLayout(dn+"_Layout",wn,w=216,h=288)
	SetDataFolder $CurrentFolder
End

Function Ediff(Intensities,scalePt)
	WAVE intensities
	Variable scalePt
	Duplicate/d/o intensities Diff
	diff=(intensities*intensities[scalePt][0]/intensities[scalePt][q]-intensities[p][0])/intensities[p][0]*100
end


//correct combined profiles (Horizontal, Vertical & Diagonal) for SiN atten, and XRF background
//Function SubBkg(bn,pFldr,idx,A,XRF,SiNtrans)  
//	String bn, pFldr
//	Variable idx //idx of the XRF SiNTrans waves to do background correction
//	Variable A //amplitude of the XRF intensity
//	WAVE XRF, SiNtrans
//	Variable i
//	For( i=-1; i<2; i+=1 )
//		WAVE dw=$(pFldr+bn+SelectString(i,"_Hy","_Vy","_Dy")), uw=$(pFldr+bn+SelectString(i,"_Hu","_Vu","_Du")) 
//		dw/=SiNTrans[i]
//		uw/=SiNTrans[i]
//		dw-=A*XRF[i]
//	endfor
//end
//

// function that stitches all appropriate (energy & pol) datasets from two different sets of data and calculates the anisotropic signal
// Uses meta data in wave note captured by NIKA for determining the profile's energy, polarization, dwell and angle
// It assumes that the user has processed 5 10deg sectors (optionally 'sectSize') in 45deg intervals starting at 90 or 'startAngle'
Function AutoStitch(sets,qpts,OutStr, [d, quiet, sm, nSectors, sectSize, th0,pStyle,aStyle, NoStitch, dupFilter, xFilter,ani, a2nm, qPwr,grp,aniExt,noU,crop,meta,nonsine,constPix,pSummary,anchor,Energies,thickness,w])
	String sets // a list of match strings (Usually just the scan number) for the datasets that you want to stitch
	String qpts // a list of Q-values where the stitch will take place (should have one less stitch point than datasets)
	String outstr // BaseName for the new dataset
	String pStyle, aStyle //profile plot style or anisotropy plot style (don't include parenthesis)
	String grp //Group by either "Pols" or "Energies" (maybe "Angles" in the future)
	Variable dupFilter // filter for duplicate profiles -1=shortest dwell; -2=longest dwell; -3=first name (alpha); -4=last name (alpha); #>0=stitch together at #
	String xFilter // Remove any profile names that match this string or list string (remove via multiple filters)
	Variable ani // calc the anisotropy from the data (need the 5 sectors for this)
	Variable d, quiet, NoStitch // just plots the raw data and doesn't stitch (tries to calculate anisotropy from the data)
	Variable sm, a2nm, qPwr //standard data processing/adjustments
	Variable nSectors, sectSize, th0
	Variable noU //don't plot the uncertainties (can become very combersome when used in the autoloader)
	Variable crop //number of points to delete from the high-q end of the final stitched profile
	Variable meta // input 1 if it already exsists.  This way the user can alter the scan lists once they are created
	Variable aniExt //Q value at which the anisotropy is extened via calculation of just the H-sectors or V-sectors from each pol (whichever pair is left at highest q)
	Variable nonsine //anisotropy is not sinusoidal with Az, so don't use the diagonal sectors to calculate the anisotropy (just use the H & V sectors)
	Variable constPix //stitch points occur at the same point on the detector (qpts adjusted per energy anchoring to first energy). Use when overlap is minimal
	Variable pSummary //plot summary at end consisting of average of all diagonals as the (average profile) and the average anisotropy for both pols)
	Variable anchor // which set to be the absolute scattering intensity achor (default is zero - the lowest-Q set)
	String Energies // List of energies (within 0.1 eV) that autostitch will process (ignoring any others - can be important if an entire scattering spectrum has been acquired)
	Variable thickness //film thickness in nm for this sample.  AutoStitch will normalize the scattering intensity to this thickness.
	Variable w // Picks to do a weighted average (1) or not (0) for the final average of averaged intensities
		//NOTE Q-point must still have the 45 deg sectors for each polarization so as to normalize the overall intensity between the two pols correctly

	//set defaults
	th0= ParamIsDefault(th0) ? 0 : th0
	nSectors= ParamIsDefault(nSectors) ? (360-th0)/45 : nSectors
	sectSize= ParamIsDefault(sectSize) ? 10 : sectSize
	dupFilter= ParamIsDefault(dupFilter) ? -2 : dupFilter
	xFilter= SelectString(ParamIsDefault(xFilter), xfilter, "")
	pStyle= SelectString(ParamIsDefault(pStyle), pStyle, "profStyle")
	aStyle= SelectString(ParamIsDefault(aStyle), aStyle, "aniStyle")
//	grp= SelectString(ParamIsDefault(grp), grp, "Pols")
	a2nm= ParamIsDefault(a2nm) ? 1 : a2nm //default is to do the conversion
	qPwr = ParamIsDefault(qPwr) ? 2 : qPwr //default is to do Q^2 scaling
	quiet= ParamIsDefault(quiet) ? 1 : quiet //usually doesn't print messages to history (quiet)
	nonsine= ParamIsDefault(nonsine) ? 0 : nonsine //usually assume a sinusoidal anisotropy
	constPix= ParamIsDefault(constPix) ? 0 : constPix //usually don't adjust qpt
	pSummary= ParamIsDefault(pSummary) ? (ani==1 ? 1 : 0 ) : pSummary
	anchor= ParamIsDefault(anchor) ? 0 : (anchor>itemsInList(sets)-1 ? itemsInList(sets)-1 : anchor ) //don't let anchor be greater than the number of sets
	w = ParamIsDefault(w) ? 1 : w
			
	
	//Special case of stitching Circular Profiles
	nSectors= sectSize==360 ? 1 : nSectors
	ani= sectSize==360 ? 0 : ani
	th0= sectSize==360 ? 45 : th0

	Variable i, j, k, pair, V_fndMatch, V_calcedSR, qpt, stitchangle, OneSet=ItemsInList(sets)<2 ? 1 : 0 //if user just supplies one dataset, just plot
	Variable hc=12398.4 //Planck's constant*Speed of Light in eV*Angstrom
	String CurrentFolder=GetDataFolder(1)
	String TopFldr="root:StitchedProfiles:"
	String DataFldr=topFldr+Outstr, MetaFldr=DataFldr+":Metadata", ProfFldr=DataFldr+":Profiles", AniFldr=DataFldr+":Anisotropy"
	String PrTblName="", profStr="", wNote, tmpStr, avName, denomName
	NewDataFolder/s/o $(removeEnding(topFldr,":")) //parent folder for all profiles calculated with AutoStitch
	NewDataFolder/s/o $DataFldr //make new data folder to put all the results
	NewDataFolder/s/o $ProfFldr //make new data folder to put all the results
	NewDataFolder/s/o $AniFldr //make new data folder to put all the results
	NewDataFolder/s/o $MetaFldr //make new data folder to put all the results
	

	//makes all the wave lists containing the filtered scan metadata
	If( !meta ) //If meta=1, the sets are already assembled and shouldn't be reassembled (likely a user tweak in them)
		AssembleScanInfo(sets, qpts, MetaFldr, sectSize,dupFilter=dupFilter,xFilter=xFilter )
	endif

	WAVE/T Scans0, Notes //load all the metadata wave lists made by "AssembleScanInfo"
	WAVE Energies0, Pols0, Angles0, Dwells0
	Make/o/n=0 ProcessedEnergies, Pols //for plotting
	Variable nE=NumUnique(Energies0,0.1), nP=NumUnique(Pols0,1), nS=NumUnique(Angles0,1)
	Notes="" //clear notes
	IF( nS!=nSectors )
		printf "Warning: Only found %d sectors instead of %d expected.\r", nS, nSectors
		nSectors=nS
	endif
	If( numpnts(Scans0) < 1 )
		print "Error: Didn't find any scans matching the first matchstring!"
		SetDataFolder $CurrentFolder
		return 0
	endif
	
	// Cycle through the stitch pairs
	For( pair=1; pair<ItemsInList(sets); pair+=1 )
		
		//load a pair of dataset info
		SetDataFolder $MetaFldr
		WAVE/Z/T ScansN=$("Scans"+num2str(pair))
		WAVE/Z EnergiesN=$("Energies"+num2str(pair)), PolsN=$("Pols"+num2str(pair))
		WAVE/Z AnglesN=$("Angles"+num2str(pair)), DwellsN=$("Dwells"+num2str(pair))
		WAVE/Z/C ratio=$("StitchRatios"+num2str(pair))
		If( numpnts(ScansN) < 1 )
			print "Error: Didn't find any scans matching the Nth matchstring! (N="+num2str(pair+1)+")"
			SetDataFolder $CurrentFolder
			return 0
		endif
		
		//Cycle through the energies & polarizations
		For( i=0; i<numpnts(Scans0); i+=k )
			If( !ParamIsDefault(Energies) && !GoodEnergy(Energies,Energies0[i]) ) // Only process &/or display the energies indicated by the user
				continue //if user wants to only process/plot subset of measured energies 
			endif
			qpt=Str2Num( StringFromList(pair-1, qpts) )
			If( constPix ) //adjust current stitch point to be at const scattering angle (or CCD pixel) instead of q anchored to the first energy
				StitchAngle=asin( qpt*hc/(Energies0[0]*4*pi) )
				qpt= 4*pi*Energies0[i]/hc*sin(stitchAngle)
//				printf "E:%1.2f, P:%d, A:%1.2f, Q:%1.4f\r", Energies0[i], Pols0[i], StitchAngle*180/pi, qpt
			endif
			
			//Cycle through the sectors
			V_calcedSR=0
			For( k=0; k<nSectors; k+=1 )
				//Find the matching first sector in the scan B list
				V_fndMatch=0
				SetDataFolder $MetaFldr
				j=FindScanMatch(pair,i+k,sets,outstr,V_FndMatch, quiet) //returns the index in the Nth scan list that matched the ith scan in scan0 & sets V_fndMatch
				SetDataFolder $ProfFldr
		
				//Stitch each sector in the scan
				sprintf profStr, "%s_%d_%1.1f_%s", OutStr, Pols0[i+k], Energies0[i+k], SelectString( sectSize==360, num2str(round(Angles0[i+k])), "C")
				//Determine the average ratio for all sectors (do only once for each scan). For now, I'm assuming all the angles are there.
				If( V_CalcedSR==0 && V_fndMatch )
					ratio[i]= CalcStitchRatio( SelectString(pair>1,Scans0[i+k],profStr), ScansN[j], nSectors, sectSize, th0, qpt, pair>1, quiet )
					V_CalcedSR=1
				endif
				If( real(ratio[i])==-1 ) //error in finding the ratio
					Notes[i+k]="Stitch Ratio Error"
//					ratio=1  //reset ratio so if there's no match next time, It doesn't kill the Notes
					continue
				endif
			
				If( !NoStitch && V_FndMatch)
					Stitch( SelectString(pair>1,Scans0[i+k],profStr), ScansN[j],profStr,qpt,ratio[i],q=quiet,appnd=pair>1,rescale=(anchor>=pair) )
					WAVE yOut=$(profStr+"_y"), xOut=$(profStr+"_x"), uOut=$(profStr+"_u"), uaOut=$(profStr+"_uAni")
					Duplicate/d/o yOut $(profStr+"_yTmp") //store current stitch profile temporarily for next round in case there are more datasets to stitch
					Duplicate/d/o xOut $(profStr+"_xTmp")
					Duplicate/d/o uOut $(profStr+"_uTmp")
					Duplicate/d/o uaOut $(profStr+"_uaTmp")
					wNote=Note(yOut)
					sprintf tmpStr, "s%dName:%s;End%d:%s;Start%d:%s;A%d:%s;uA%d:%s;", pair, ScansN[i+k], pair-1, StringByKey("End0",wNote), pair, StringByKey("Start1",wNote), pair, StringbyKey("a0",wNote), pair, StringByKey("uA0",wNote)
					Notes[i+k]+= SelectString(pair>1, wNote, tmpStr)
				elseif( NoStitch ) // user just want's the raw data processed and displayed w/out stitching
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'r_"+Scans0[i+k]+"'") $(profStr+"_y0")
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'q_"+Scans0[i+k]+"'") $(profStr+"_x0")
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'s_"+Scans0[i+k]+"'") $(profStr+"_u0")
					duplicate/d/o $("root:SAS:'"+ScansN[i+k]+"':'r_"+ScansN[i+k]+"'") $(profStr+"_y"+num2str(pair) )
					duplicate/d/o $("root:SAS:'"+ScansN[i+k]+"':'q_"+ScansN[i+k]+"'") $(profStr+"_x"+num2str(pair) )
					duplicate/d/o $("root:SAS:'"+ScansN[i+k]+"':'s_"+ScansN[i+k]+"'") $(profStr+"_u"+num2str(pair) )
				elseif( !V_FndMatch && pair==1 )
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'r_"+Scans0[i+k]+"'") $(profStr+"_yTmp")
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'q_"+Scans0[i+k]+"'") $(profStr+"_xTmp")
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'s_"+Scans0[i+k]+"'") $(profStr+"_uTmp")
					duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'s_"+Scans0[i+k]+"'") $(profStr+"_uaTmp")
					WAVE yOut=$(profStr+"_yTmp")
					sprintf tmpStr, "s%dName:%s;", pair-1, profStr
					Notes[i+k]+= tmpStr
				Endif
			endfor // loop over sectors
		endfor // loop over scans
	endfor // loop over pairs
	
	//apply various corrections, notes, display the data, and calculate anisotropy if requested
	If( d==1 )
		DoWindow $Outstr
		If( V_Flag==0 ) //make new layout for all the graphs
			NewLayout/N=$Outstr as OutStr
			PrintSettings/I margins={0.5,0.5,0.5,0.5}
		elseif( WinType(OutStr)!=3 ) //it isn't a layout
			DoWindow/F $Outstr
			print "Name Conflict: A window already exists with this name."
			SetDataFolder $CurrentFolder
			abort
		endif
	endif
//	print OutStr, StopMStimer(0), StopMStimer(1), StopMStimer(2), StopMStimer(3), StopMStimer(4), StopMStimer(5), StopMStimer(6), StopMStimer(7), StopMStimer(8), StopMStimer(9)
	For( i=0; i<numpnts(Scans0); i+=k )
		If( !ParamIsDefault(Energies) && !GoodEnergy(Energies,Energies0[i]) ) // Only process &/or display the energies indicated by the user
			continue //if user wants to only process/plot subset of measured energies 
		endif
		InsertPoints inf, 1, ProcessedEnergies,Pols
		ProcessedEnergies[numpnts(ProcessedEnergies)-1]=Energies0[i]
		Pols[numpnts(Pols)-1]=Pols0[i]
		For( k=0; k<nSectors; k+=1 )
			If( StringMatch(Notes[i+k], "Stitch Ratio Error") )
				continue
			endif
			SetDataFolder $profFldr
			sprintf profStr, "%s_%d_%1.1f_%s", OutStr, Pols0[i+k], Energies0[i+k], SelectString( sectSize==360, num2str(round(Angles0[i+k])), "C")
			if( OneSet ) // No stitching
				duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'r_"+Scans0[i+k]+"'") $(profStr+"_y")
				duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'q_"+Scans0[i+k]+"'") $(profStr+"_x")
				duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'s_"+Scans0[i+k]+"'") $(profStr+"_u")
				duplicate/d/o $("root:SAS:'"+Scans0[i+k]+"':'s_"+Scans0[i+k]+"'") $(profStr+"_uAni")
			else //remove temp waves
				Duplicate/d/o $(profStr+"_yTmp") $(profStr+"_y")
				Duplicate/d/o $(profStr+"_xTmp") $(profStr+"_x")
				Duplicate/d/o $(profStr+"_uTmp") $(profStr+"_u")
				Duplicate/d/o $(profStr+"_uaTmp") $(profStr+"_uAni")
				KillWaves/Z $(profStr+"_yTmp"), $(profStr+"_xTmp"), $(profStr+"_uTmp"), $(profStr+"_uaTmp")
			endif
			WAVE yOut=$(profStr+"_y"), xOut=$(profStr+"_x"), uOut=$(profStr+"_u"), uaOut=$(profStr+"_uAni")
			//assuming no uncert in x values (not really true, but they SHOULD be much smaller)
			//normalize to thickness if user gave input - nm is qPwr >0 otherwise in cm
			
			//I CHANGED THE NEXT PORTION
			
			yOut *= xOut^qPwr / (thickness==0 ? 1 : thickness* (qPwr==0 ? 1e-7 : 1 ) )
			uOut *= xOut^qPwr / (thickness==0 ? 1 : thickness* (qPwr==0 ? 1e-7 : 1 )  )
			
			//THE TWO LINES ABOVE
			
			
			//uncertainties w/out Ratio uncertaines added
			uaOut *= xOut^qPwr / (thickness==0 ? 1 : thickness* (qPwr==0 ? 1e-7 : 1 ) )
			xOut *= a2nm ? 10 : 1
			
			
			
			yOut = yOut<0 ? 0 : yOut //don't allow negative values
			
			
			
//			Duplicate/d/o uOut $(profStr+"_uAni")  //add the ratio uncertainties back into the calculation for anisotropy uncertainty
			If( sm>0 )
				smooth /E=3 sm, yOut
			endif
			Variable/C Integral=IntData(yOut,xOut,uOut)
			Deletepoints numpnts(xOut)-1, crop, yOut, xout, uOut, uaOut //crop the data
			NOTE/K/NOCR yOut, Notes[i+k]+"Integral:"+num2str(real(integral))+";IntU:"+num2str(imag(integral))
			NOTE/K/NOCR xOut, Notes[i+k]
			NOTE/K/NOCR uOut, Notes[i+k]
			avName=CombinePols(profStr+"_", w=w) //avName is profstr with the pol info stripped
			
			nE=NumUnique(ProcessedEnergies,0.1)
			If( d )
//				If( (ParamIsDefault(grp) && OnePol(Pols0)) || (!ParamIsDefault(grp) && StringMatch("Pols",grp)) ) //groups plot by polarization - only plots if user specified this or just on pol was acquired
					PlotProf( yOut, xOut, uOut, Outstr, "P"+num2str(Pols0[i+k]), "", 0, pStyle, th0, nSectors, nE, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
				if( !ParamIsDefault(grp) && StringMatch("Energies",grp) ) //plots energies into separate plots
					PlotProf( yOut, xOut, uOut, Outstr, "En"+ReplaceString(".",num2str(Energies0[i+k]),"p"), "", 0, pStyle, th0, nSectors, nP, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
				endif
				If( StrLen(avName) && Pols0[i+k]==190 ) //doubly sure that it only plots it once.
					WAVE avY=$(avName+"y"), avX=$(avName+"x"), avU=$(avName+"u")
					Integral=IntData(avY, avX, avU)
					NOTE/NOCR avY, "AvInt:"+num2str(integral)+";"
					//PlotProf( avY, avX, avU, Outstr, "Avg", "", 0, pStyle, 45, 1, nE, noU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness,pols=pols)
				endif
			endif
		endFor //sectors
		
		If( ani==1 && !StringMatch(Notes[i], "Stitch Ratio Error") )  //calculate the anisotropy from the sector plots
			sprintf profStr, "%s_%d_%1.1f", OutStr, Pols0[i], Energies0[i]
			CombineSectors(profSTr, th0, nSectors)//combine the appropriate sectors into just Horizontal, Vertical & Diagnoal sectors (used by CalcAni)
			CalcAni(profStr, fldr=aniFldr,th0=th0, nSectors=nSectors, useD=!nonsine)
			SetDataFolder $aniFldr
			string avAniName=CombinePols(profStr+"_ani",f=-1,w=!nonsine) //avName is profstr with the pol info stripped
			SetDataFolder $profFldr
			String avDname=CombinePols(profStr+"_D",w=!nonsine) //diagnoal average for all polarizations (used as averaged profile)
//			denomName= MakeAniDenominators(profstr)
			If( aniExt>0 )
				If( constPix ) //adjust current stitch point to be at const scattering angle (or CCD pixel) instead of q anchored to the first energy
					StitchAngle=asin( aniExt*hc/(Energies0[0]*4*pi) )
					qpt= 4*pi*Energies0[i]/hc*sin(stitchAngle) * 10 //get it in invAngstroms
//					printf "E:%1.2f, P:%d, A:%1.2f, Q:%1.4f\r", Energies0[i], Pols0[i], StitchAngle*180/pi, qpt
				endif
				ExtendAni(profStr,qpt, aniFldr, quiet)
			endif
			If( d )
				WAVE aniY=$(aniFldr+":'"+profStr+"_aniY'"), aniX=$(aniFldr+":'"+profStr+"_aniX'"), aniU=$(aniFldr+":'"+profStr+"_aniU'")
				WAVE Hy=$(profStr+"_Hy"), Hx=$(profStr+"_Hx"), Hu=$(profStr+"_Hu")
				WAVE Vy=$(profStr+"_Vy"), Vx=$(profStr+"_Vx"), Vu=$(profStr+"_Vu")
				WAVE Dy=$(profStr+"_Dy"), Dx=$(profStr+"_Dx"), Du=$(profStr+"_Du")
				PlotProf( aniY, aniX, aniU, Outstr, "AniPol", "", 0, aStyle, 45, 1, nE, noU,hide=pSummary,d=d,pols=pols)
				PlotProf( Hy, Hx, Hu, Outstr, "Profs", "", 0, pStyle, th0, 3, nE, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
				PlotProf( Dy, Dx, Du, Outstr, "Profs", "", 0, pStyle, th0, 3, nE, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
				PlotProf( Vy, Vx, Vu, Outstr, "Profs", "", 0, pStyle, th0, 3, nE, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
				//Don't need this if plotting the summary
//				If( StrLen(avAniName) && Pols0[i]==190 ) //doubly sure that it only plots it once.
//					WAVE avY=$(aniFldr+":'"+avAniName+"y'"), avX=$(aniFldr+":'"+avAniName+"x'"), avU=$(aniFldr+":'"+avAniName+"u'")
//					PlotProf( avY, avX, avU, Outstr, "Ani", "", 0, aStyle, 45, 1, nE, NoU,hide=pSummary,d=d,pols=pols)
//					WAVE avY=$(profFldr+":'"+avDname+"y'"), avX=$(profFldr+":'"+avDname+"x'"), avU=$(profFldr+":'"+avDname+"u'")
//					PlotProf( avY, avX, avU, Outstr, "Diagonal", "", 0, pStyle, 45, 1, nE, NoU,hide=pSummary,d=d,qPwr=qPwr,thickness=thickness)
//				endif
			endif
		endif
	Endfor //scans
	If( pSummary && ani && d )
		PlotSummary(outStr,qpwr,thickness>0)
	endif
	SetDataFolder $CurrentFolder
End


//Checks to see if there is only one polarization in the metadata: return 1 if only one pol, 0 if more than one
Function OnePol(pols)
	WAVE pols
	Variable pol0=pols[0], i
	For( i=0; i<numpnts(pols); i+=1 )
		If( pols[i] != pol0 )
			return 0 //found a different polarization
		endif
	endfor
	return 1 //all polarizations equaled the first one
end


//Searches the energy list to see if current energy is on the list. Returns 1 if it is, 0 if it's not on the list
Function GoodEnergy(Energies,CurrentEnergy)
	String Energies
	Variable CurrentEnergy
	Variable k
	For( k=0; k<itemsInList(Energies) ; k+=1 )
		If( abs( CurrentEnergy - Str2Num(StringFromList(k,Energies))) < 0.1 )
			return 1
		endif
	Endfor
	return 0
end

//returns number of matching energies within tolerance
Function NumMatches(eList,eWave,tol)
	string eList
	WAVE eWave
	variable tol
	Variable i, j, ectr=0
	For( i=0; i<numpnts(eWave); i+=1 )
		For( j=0; j<itemsInList(eList); j+=1 )
			If( abs( eWave[i] - Str2Num(StringFromList(j,elist))) < tol )
				ectr+=1
			endif
		endfor
	endfor
	return ectr
end

//searches for a match profile in the paired list. Returns index of the pair and sets V_FndMatch=0 if can't find one
Function FindScanMatch(pair, i, sets, dName, V_FndMatch, quiet)
	Variable pair, i, quiet
	Variable &V_FndMatch
	String sets, dName
	WAVE/Z/T Scans0, ScansN=$("Scans"+num2str(pair))
	WAVE/Z Energies0, Pols0, EnergiesN=$("Energies"+num2str(pair)), PolsN=$("Pols"+num2str(pair))
	WAVE/Z Angles0, Dwells0, AnglesN=$("Angles"+num2str(pair)), DwellsN=$("Dwells"+num2str(pair))
	Variable j, k
	For( j=0; j<numpnts(ScansN); j+=1 )
		If( abs(Pols0[i]-PolsN[j])<1 && abs(Energies0[i]-EnergiesN[j])<0.1 && abs(Angles0[i]-AnglesN[j])<1 )
			V_fndMatch+=1
			break
		endif
	EndFor
	For( k=j+1; k<numpnts(ScansN); k+=1 ) //look for any more matches
		If( abs(Pols0[i]-PolsN[k])<1 && abs(Energies0[i]-EnergiesN[k])<0.1 && abs(Angles0[i]-AnglesN[k])<1 )
			V_fndMatch+=1
			break
		endif
	EndFor
	
	//Handle Errors in matching the  scans
	If( V_FndMatch==0 && quiet==0 ) //No match found
		Printf "Match Error! No Match for P=%d, E=%1.1f, A=%d in set: %s\r", Pols0[i], Energies0[i], Angles0[i], StringFromList(pair,sets)
		EditMetaData(dName) //Generate diagnostics
	elseif( V_FndMatch>1 && quiet==0 )	//Check for and filter any duplicates still remaining by asking user to decide
		Printf "Match Error! Multiple Matchs for, P=%d, E=%1.1f, A=%d to set: %s\r", Pols0[i], Energies0[i], Angles0[i], StringFromList(pair,sets)
		EditMetaData(dName) //Generate diagnostics
		//at this point it will just use the first match it finds and print this error
		// in future could have user apply another filter and/or select the one they want
	endif
	return j
end

Function/WAVE MakeEnergies(Pols, energies)
	WAVE pols, energies
	FindLevel/P/Q Pols, 189
	If( !V_Flag ) //found crossing
		Make/o/n=(ceil(V_levelX)) Energy
	else
		Make/o/n=(numpnts(energies)) Energy
	endif
	Energy=energies[p]
	return Energy
end

Function QueryMeta(Angle, Energy, Pol, fldr,[tg])
	Variable angle, energy, pol
	string fldr,tg
	WAVE aW=$(fldr+":Angles0"), eW=$(fldr+":Energies0"), pW=$(fldr+":Pols0")
	If( !ParamIsDefault(tg) )
		WAVE aW=$(fldr+":Angles"+tg), eW=$(fldr+":Energies"+tg), pW=$(fldr+":Pols"+tg)
	endif		
	Variable i
	For( i=0; i<numpnts(aW); i+=1 )
		If( aW[i]==angle && eW[i]==energy && pW[i]==Pol )
			return i //return index of the match
		endif
	endfor
	return -1 //didn't find a match
end

//BETA CODE!!! Check to see if QueryMeta returns all '-1' for a specific match...
//Returns -1 for FAILURE (QueryMeta fails its search)
//Returns 1 for PASS

Function QueryMetaFailCheck(w)
	wave w
	
	variable i
	for(i=0;i<numpnts(w);i+=1)
		if(!(w[i] == -1))
			return 1
		endif
	endfor
	return -1
	
end



//makes a table with all the metadata for the indicated dataset
Function EditMetaData(dName,[type])
	string dName
	variable type // 1=Escan (just to look in a different folder for the metadata)
	String CurrentFolder=GetDataFolder(1), fldr="root:StitchedProfiles"+dName+"_data:Metadata", TblName="MetaData_"+dName
	If( type==1 )
		fldr="root:Escans:"+dName+":Meta"
	else
		fldr=GetFldrNameASmeta(dName)
	endif
	If( !DataFolderExists(fldr) )
		print "Can't find the data!"
		return 0
	endif
	SetDataFolder $fldr
	variable i
	DoWindow/F $TblName
	If( V_Flag )
		SetDataFolder $CurrentFolder
		return 0
	endif
	edit/n=$Tblname
	For( i=0; ; i+=1 )
		WAVE/T/Z Scans=$("Scans"+num2str(i))
		WAVE/Z Energies=$("Energies"+num2str(i)), Pols=$("Pols"+num2str(i)), Angles=$("Angles"+num2str(i)), Dwells=$("Dwells"+num2str(i))
		WAVE/Z/C Ratios=$("StitchRatios"+num2str(i))
		If( !WaveExists(Scans) )
			break
		endif
		If( WaveExists(Ratios) )
			AppendtoTable/W=$Tblname Ratios
			ModifyTable/W=$Tblname width($NameOfWave(Ratios))=30,sigDigits($NameOfWave(Ratios))=4
		endif
		Appendtotable/W=$Tblname Scans, Pols, Energies, Angles, Dwells
		ModifyTable/W=$Tblname size=8, width(Point)=22, width($NameOfWave(scans))=100,width($NameOfWave(Pols))=30,width($NameOfWave(Energies))=30
		ModifyTable/W=$Tblname width($NameOfWave(Angles))=30,width($NameOfWave(Dwells))=30
	endfor
	SetDataFolder $CurrentFolder
end

//checks a bunch of possible locations for the folder with the knowledge of my past naming conventions (for bkwrds compatibility)
Function/S GetFldrNameASmeta(dn)
	String dn
	string fldr="root:"+dn
	If( !DataFolderExists(fldr) )
		fldr="root:"+dn+"_data"
	endif
	If( !DataFolderExists(fldr) )
		fldr="root:StitchedProfiles:"+dn
	endif
	If( !DataFolderExists(fldr) )
		fldr="root:StitchedProfiles:"+dn+"_data"
	endif
	fldr+=":Metadata"
	return fldr
end
	

//Combines (averages) profiles from Horizontal and Vertical polarizations.  Assumes the polarziations are in the name.
Function/S CombinePols(inName,[f,w])
	String inName
	Variable f //factor to multiply by the 190 polarization prior to averaging
	Variable w //use weighted average
	f= ParamIsDefault(f) ? 1 : f
	w= ParamIsDefault(w) ? 1 : w
	String w100n=replacestring("_190_",inName,"_100_")
	String w190n=replacestring("_100_",inName,"_190_")
	w100n=removeEnding(w100n,"y"); w100n=removeEnding(w100n,"x")
	w190n=removeEnding(w190n,"y"); w190n=removeEnding(w190n,"x")
	String avName=replacestring("_100_",w100n,"_")
	WAVE/Z w100y=$(w100n+"y"), w100x=$(w100n+"x"), w100u=$(w100n+"u")
	WAVE/Z w190y=$(w190n+"y"), w190x=$(w190n+"x"), w190u=$(w190n+"u")
	If( WaveExists(w100y) && WaveExists(w190y) )
		CombineProfiles(w100y,w100x,w100u,w190y,w190x, w190u,f=f, w=w)
		Duplicate/d/o AvY, $(avName+"y"); Duplicate/d/o AvX, $(avName+"x"); Duplicate/d/o AvU, $(avName+"u")
		return avName
	elseif( WaveExists(w100y) )
		Duplicate/d/o w100y, $(avName+"y"); Duplicate/d/o w100x, $(avName+"x"); Duplicate/d/o w100u, $(avName+"u")
		return avName
	elseif( WaveExists(w190y) )
		Duplicate/d/o w190y, $(avName+"y"); Duplicate/d/o w190x, $(avName+"x"); Duplicate/d/o w190u, $(avName+"u")
		return avName
	else
		return ""
	endif
end



//Macro asdf
//					WAVE avY=$(denomName+"_dh"), avX=$(profStr+"_hx"), avU=$(profStr+"_hu")
//					avName=replacestring(".",num2str(round(Energies0[i]*10)/10),"p")
//					PlotProf( avY, avX, avU, Outstr,"Denom_"+avName,"",0,pStyle,45, 1, 10,1)
//					WAVE avY=$(denomName+"_dv"), avX=$(profStr+"_vx"), avU=$(profStr+"_vu")
//					PlotProf( avY, avX, avU, Outstr,"Denom_"+avName,"",0,pStyle,45, 1, 10,1)
//					WAVE avY=$(denomName+"_d100"), avX=$(profStr+"_45x"), avU=$(profStr+"_45u")
//					PlotProf( avY, avX, avU, Outstr,"Denom_"+avName,"",0,pStyle,45, 1, 10,1)
//					WAVE avY=$(denomName+"_d190"), avX=$(profStr+"_45x"), avU=$(profStr+"_45u")
//					PlotProf( avY, avX, avU, Outstr,"Denom_"+avName,"",0,pStyle,45, 1, 10,1)
//					WAVE Cy=$(denomName+"_Cy"), Cx=$(denomName+"_Cx"), Cu=$(denomName+"_Cu")
//					PlotProf( Cy, Cx, Cu, Outstr,"Ani_denom","",0,aStyle,45, 1, 10,1)
//EndMacro
//

//Some academic exercise that failed
Function/S MakeAniDenominators(inName)
	String inName
	String w100n=replacestring("_190_",inName,"_100_")
	String w190n=replacestring("_100_",inName,"_190_")
	WAVE/Z w100hy=$(w100n+"_Hy"), w100hx=$(w100n+"_Hx"), w100hu=$(w100n+"_Hu")
	WAVE/Z w100vy=$(w100n+"_Vy"), w100vx=$(w100n+"_Vx"), w100vu=$(w100n+"_Vu")
	WAVE/Z w100dy=$(w100n+"_Dy"), w100dx=$(w100n+"_Dx"), w100du=$(w100n+"_Du")
	WAVE/Z w190hy=$(w190n+"_Hy"), w190hx=$(w190n+"_Hx"), w190hu=$(w190n+"_Hu")
	WAVE/Z w190vy=$(w190n+"_Vy"), w190vx=$(w190n+"_Vx"), w190vu=$(w190n+"_Vu")
	WAVE/Z w190dy=$(w190n+"_Dy"), w190dx=$(w190n+"_Dx"), w190du=$(w190n+"_Du")
	String avName=replacestring("_100_",w100n,"_")
	If( WaveExists(w100hy) && WaveExists(w190hy) )
		Duplicate/d/o w100hy $(avName+"_dh")
		Duplicate/d/o w100vy $(avName+"_dv")
		Duplicate/d/o w100dy $(avName+"_d100")
		Duplicate/d/o w190dy $(avName+"_d190")
		WAVE dh=$(avName+"_dh"), dv=$(avName+"_dv"), d100=$(avName+"_d100"), d190=$(avName+"_d190")
		CombineProfiles(w100hy,w100hx,w100hu,w190hy,w190hx,w190hu,f=-1)
		Duplicate/d/o AvY, $(avName+"_Ay"); Duplicate/d/o AvX, $(avName+"_Ax"); Duplicate/d/o AvU $(avName+"_Au")
		WAVE/Z Ay=$(avName+"_Ay"), Ax=$(avName+"_Ax"), Au=$(avName+"_Au")
		CombineProfiles(w100vy,w100vx,w100vu,w190vy,w190vx,w190vu,f=-1)
		Duplicate/d/o AvY, $(avName+"_By"); Duplicate/d/o AvX, $(avName+"_Bx"); Duplicate/d/o AvU $(avName+"_Bu")
		WAVE/Z By=$(avName+"_By"), Bx=$(avName+"_Bx"), Bu=$(avName+"_Bu")
		CombineProfiles(Ay,Ax,Au,By,Bx,Bu)
		Duplicate/d/o AvY, $(avName+"_Cy"); Duplicate/d/o AvX, $(avName+"_Cx"); Duplicate/d/o AvU $(avName+"_Cu")
		WAVE/Z Cy=$(avName+"_Cy"), Cx=$(avName+"_Cx"), Cu=$(avName+"_Cu")
		dh=w100hy+w190hy*0.91
		dv=w100vy+w190vy*0.91
		d100=2*(w100dy)
		d190=2*(w190dy)
		Cy/=(interp(Cx,w100dx,w100dy)+interp(Cx,w190dx,w190dy))*.01
		return avName
	Endif
end	
//

//USED IN CONJUCTION WITH ASSEMBLESCANINFO //Returns a string filled with waves given a specific header...
//For example, will grab all Energies scans into a list "Energy0;Energy1;Energy2" etc.


Function/S AssembleMetaWave(Prefix,MetaFldr)

	String Prefix //The specific prefix that you are looking to assemble.
	String MetaFLDR //The folder that holds all the meta data
	//Loop variable
	Variable i

	//Given the meta folder we grab a list of all waves within the folder
	String MetaWaveList = ReplaceString(",",Stringbykey("WAVES",DataFolderDir(2,$MetaFldr),":",";"),";") //Sets up a semicolon separated list of waves in the meta folder
	String PrefixWaveList = "" //List of all waves given the appropriate Prefix.
	for(i=0;i<itemsinlist(MetaWaveList,";");i+=1)
		String ListItem = Stringfromlist(i,MetaWaveList,";")
		
		
		if(StringMatch(Listitem,Prefix+"*"))
			PrefixWaveList += MetaFldr+":"+ListItem + ";"
		endif
	endfor

//	print metaWavelist
//	print PrefixWaveList
	return PrefixWaveList

end



//Function AssembleScanInfo(sets, qpts, fldr, nSectors, sectSize,startAngle,dupFilter,xFilter )
Function/S AssembleScanInfo(sets, qpts, fldr, sectSize,[dupFilter,xFilter,nSectors] )
	String sets, qpts //tells if user want's to splice scans into same Q-range (b/c of det saturation, or forgot an energy or polarization, etc.)
	String fldr
	String xFilter
//	Variable nSectors, sectSize,startAngle,dupFilter
	Variable sectSize,dupFilter, nSectors
	Variable i, j, k, range, angle, scanAngle, qpt, insertion
	nSectors= ParamIsDefault(nSectors) ? 8 : nSectors
	dupFilter= ParamIsdefault(dupFilter) ? -2 : dupFilter
	xFilter= SelectString(ParamIsDefault(xFilter), xFilter, "")
	String CurrentFolder=GetDataFolder(1), sName //set name which is usually just a number but can be '2b' if several sets are inserted into one Q-range
	String rangeList="" //list mirroring "sets" that describes which Q-range each set should go in
	String newQpts="" //new list of q-points stripped of potential blank spots due to insertions
	SetDataFolder root:SAS
	String FldrList=ReplaceString(",",StringByKey("FOLDERS",DataFolderdir(1)),";")	//Full list of profiles processed by NIKA

	For( i=0, range=0; i<ItemsInList(sets); i+=1, range+=1 ) //range is the q-range hierarchy, i is the index in the set list
		SetDataFolder $fldr
		If( i>0 )
			qpt=Str2Num( StringFromList(i-1, qpts) )
			Insertion= Numtype(qpt)==2 && i<ItemsInList(sets)-1  ? Insertion+1 : 0 //if  qpt=nan, then set flag to insert into prev scan
			range -= Insertion ? 1 : 0 //demote the range since this one is just to be inserted into previous range
			If( insertion==1 ) //first  insertion into this range, duplicate the earlier dataset's metadata into range + "a"
				WAVE/T Scans=$("Scans"+num2str(range))
				WAVE Energies=$("Energies"+num2str(range)), Pols=$("Pols"+num2str(range)), Angles=$("Angles"+num2str(range)), Dwells=$("Dwells"+num2str(range))
				Duplicate/t/o Scans, $("Scans"+num2str(range)+"a")
				Duplicate/d/o Energies, $("Energies"+num2str(range)+"a")
				Duplicate/d/o Pols, $("Pols"+num2str(range)+"a")
				Duplicate/d/o Angles, $("Angles"+num2str(range)+"a")
				Duplicate/d/o Dwells, $("Dwells"+num2str(range)+"a")		
			Endif
			sName=num2str(range)+selectString( insertion, "", num2char(97+insertion) )
		else
			sName="0"
		endif
		rangeList+=num2str(range)+";"
		Make/t/o/n=0 $("Scans"+sName) //make the wave lists for each set
		Make/o/n=0 $("Energies"+sName), $("Pols"+sName), $("Angles"+sName), $("Dwells"+sName)
		WAVE/T Scans=$("Scans"+sName)
		WAVE Energies=$("Energies"+sName), Pols=$("Pols"+sName), Angles=$("Angles"+sName), Dwells=$("Dwells"+sName)
		If( range>0 )  //only make stitchratio waves for range pairs
			Make/o/c/n=0 $("StitchRatios"+sName)
			WAVE/C Ratios=$("StitchRatios"+sName)
		endif

		SetDataFolder root:SAS
		String profStr, tmpList, tmpFilter
		String ScanListA=ListMatch(FldrList,"*"+StringFromList(i,sets)+"*"),ScanListB="",eList="", pList1="" //filter using user's match string
		//filter for sector size replacing decimal with "p".  If sectSize=360, then look for circular average
		ScanListA= ListMatch(ScanListA,  SelectString(sectSize==360, "*_"+replacestring(".", num2str(sectSize), "p"), "*_C" ) )

		If( StrLen(xFilter) )	//filter for xfilter
			For( j=0; j<ItemsInList(xFilter); j+=1 )
				tmpFilter=StringFromList(j,xFilter)
				tmpList=ListMatch(ScanListA,"*"+tmpFilter+"*")
				ScanListA=RemoveFromList(tmpList,ScanListA)
			endfor
		endif

		For( k=0; k<nSectors; k+=1 ) //filter for correct angles
			angle=k*45
			For( j=0; j<ItemsInList(ScanListA); j+=1 ) //find any in the list at this angle
				scanAngle= str2Num( ParseFilePath(0, StringFromList(j, ScanListA), "_",1,1) ) //extract the profile's angle from the name
				If( abs( scanAngle-angle) < 0.1 )
					ScanListB+= StringFromList(j,ScanListA)+";" //make a new list with just the correct scan angles
				endif
			endfor
		Endfor
		If( sectSize==360 )  //only do if not looking for Circ-Average
			ScanListB=ScanListA
		endif

		Redimension/n=(ItemsInList(ScanListB)) Scans, Energies, Pols, Angles, Dwells //make waves instead of lists for sorting later
		Scans= StringFromList(p,ScanListB)
		
		For( j=0; j<numpnts(Scans); j+=1 )//Get meta information on each scan on energy and polarization from the wave note
			WAVE w=$(":'"+Scans[j]+"':'r_"+Scans[j])+"'"
			
			//QUICK ADDITION TO IGOR 7 - T Ferron
			//String EnergiesString = Stringbykey("Beamline Energy",Note(w),"=",";")
			//String PolsString = Stringbykey("EPU polarization",Note(w),"=",";")
			//String DwellsString = Stringbykey("Exposure",Note(w),"=",";")	
			
			//Quick changes for compatibility with data from Eliot's beamline (NSLSII) -Obaid 
			//String EnergiesString = Stringbykey("en_energy_setpoint",Note(w),":",";")
			//String PolsString = Stringbykey("en_polarization_setpoint",Note(w),":",";")
			//String DwellsString = Stringbykey("SampleMeasurementTime",Note(w),"=",";")
			
			String EnergiesString = Stringbykey("X-ray_energy",Note(w),":",";")
			String PolsString = Stringbykey("X-ray_polarization_setpoint",Note(w),":",";")
			String DwellsString = Stringbykey("Shutter_Opening_time",Note(w),":",";")
			//No units/spaces to remove (for NSLSII data) so the follwoing three lines of code are commentized  -Obaid
			//EnergiesString = RemoveUnits(EnergiesString)
			//PolsString = RemoveUnits(PolsString)
			//DwellsString = RemoveUnits(DwellsString)
				
					
			Energies[j] = round(20*str2num(EnergiesString))/20 //to nearest 0.05eV
			Pols[j] = str2num(PolsString)
			Dwells[j] = Str2num(DwellsString)
		
	//OLD CODE		
//			Energies[j]= round(20*NumberByKey("BeamlineEnergy",Note(w),"="))/20 //to nearest 0.05eV
//			Pols[j]= NumberByKey("EPUpolarization",Note(w),"=")
//			Dwells[j]= NumberByKey("Exposure",Note(w),"=")

//Return to old code
			Angles[j]= sectSize==360 ? 360 : Str2Num( ParseFilePath(0, Scans[j], "_", 1, 1) )
			
			//EDIT 5/23/2016 Thomas Ferron  -- correct for Horseshit Polarization 200 nonsense
//			if (Pols[j]==200)
//				Pols[j]=100
//			elseif(Pols[j]==380)
//				Pols[j] = 190
//			endif
//			
//		Endfor
		if (Pols[j]==0)
				Pols[j]=100
			elseif(Pols[j]==90)
				Pols[j] = 190
			endif
			
		Endfor
		//Sort by Polarization, then energies within that, then angles within that, and finally dwells within that
//		If( i==0 )
//			Energies+=0.1
//		endif
		Sort {Pols,Energies,Angles,Dwells}, Scans, Energies, Pols, Angles, Dwells

		//Check for and filter any duplicate angles still remaining by asking user to decide what to do
		For( j=0; j<numpnts(Scans)-1; j+=1 )
			If( abs(Energies[j]-Energies[j+1]) < 0.1 && abs(Angles[j]-Angles[j+1]) < 1 )
				ApplyUserFilter(j, dupFilter, Pols, Energies, Angles, Dwells, Scans)
			endif
		Endfor
		If( i==0 )
			SetDataFolder $fldr
			Duplicate/t/o Scans Notes; Notes=""
		endif
		If( range>0 )
			Redimension/n=(numpnts(scans)) Ratios
			Ratios=0
		endif
	Endfor
	Variable nRanges=Range, jRange, en, pol, sect, idx
	SetDataFolder $fldr
	For( i=0; i<nRanges; i+=1 )
		WAVE/T Scans=$("Scans"+num2str(i))
		WAVE Energies=$("Energies"+num2str(i)), Pols=$("Pols"+num2str(i)), Angles=$("Angles"+num2str(i)), Dwells=$("Dwells"+num2str(i))
		WAVE/C/Z Ratios=$("Ratios"+num2str(i))
		range=str2Num(stringFromList(i,rangeList))
		insertion=0
		For( j=i+1; j<itemsInlist(rangeList); j+=1 )
			jRange=str2Num(stringFromList(j,rangeList))
			If( range==jRange )  //these two sets are in the same range
				Insertion= j-i
				sName=num2str(Jrange)+selectString( insertion, "", num2char(97+insertion) )
				WAVE/T ScansJ=$("Scans"+sName)
				WAVE EnergiesJ=$("Energies"+sName), PolsJ=$("Pols"+sName), AnglesJ=$("Angles"+sName), DwellsJ=$("Dwells"+sName)
				For( k=0; k<Numpnts(ScansJ); k+=1 ) //cycle through inserted set and insert the scans into the original set
					en=EnergiesJ[k]
					pol=PolsJ[k]
					sect=AnglesJ[k]
					idx=QueryMeta(sect, En, Pol, fldr,tg=num2str(i)) //find the matching profile
					If( idx>=0 ) //found a match
						Scans[idx]=scansJ[k]
						Energies[idx]=EnergiesJ[k]
						Pols[idx]=PolsJ[k]
						Angles[idx]=AnglesJ[k]
						Dwells[idx]=DwellsJ[k]
					else //no match found, this is a unique dataset
						InsertPoints inf, 1, Energies, Pols, Angles, Dwells, Scans
						If( WaveExists(Ratios) )
							InsertPoints inf, 1, Ratios
						endif
					endif
				endfor
			else
				i+=insertion //if there were insertions, skip these for the next range
				break
			endif
		endfor
	Endfor
	For( i=0; i<itemsInList(qPts); i+=1 )
		qpt=Str2Num(StringFromList(i,qpts) )
		If( NumType(qpt)!=2 )
			newQpts+=num2str(qpt)+";"
		endif
	endfor
	SetDataFolder $CurrentFolder
	return newQpts
end


//returns the name of the dataset for AssembleSetInfo
Function/S MakeSetName(range,insertion)
	Variable range, insertion
	If( insertion==1 ) //first  insertion into this range, duplicate the earlier dataset's metadata into range + "a"
		WAVE/T Scans=$("Scans"+num2str(range))
		WAVE Energies=$("Energies"+num2str(range)), Pols=$("Pols"+num2str(range)), Angles=$("Angles"+num2str(range)), Dwells=$("Dwells"+num2str(range))
		Duplicate/t/o Scans, $("Scans"+num2str(range)+"a")
		Duplicate/d/o Energies, $("Energies"+num2str(range)+"a")
		Duplicate/d/o Pols, $("Pols"+num2str(range)+"a")
		Duplicate/d/o Angles, $("Angles"+num2str(range)+"a")
		Duplicate/d/o Dwells, $("Dwells"+num2str(range)+"a")		
	Endif
	return num2str(range)+selectString( insertion, "", num2char(97+insertion) )
end

Function SortMetaData(set, [type])
	Variable set
	String type
	WAVE/T Scans=$("Scans"+num2str(set))
	WAVE Energies=$("Energies"+num2str(set)), Pols=$("Pols"+num2str(set)), Angles=$("Angles"+num2str(set)), Dwells=$("Dwells"+num2str(set))
	IF( ParamIsDefault(type) )
		Sort {Pols,Energies,Angles,Dwells}, Scans, Energies, Pols, Angles, Dwells
	elseif( stringmatch("Energy",type) )
		Sort {Pols,Angles,Dwells,Energies}, Scans, Energies, Pols, Angles, Dwells
	endif
end


//Caculates the average ratio between the two scans for all sectors based on the given wave names
Function/C CalcStitchRatio(wn1, wn2, nSectors,sectSize,startAngle, qpt,appnd, quiet)
	String wn1, wn2
	Variable nSectors, sectSize, startAngle, qpt //q-value where the stich takes place
	Variable appnd //appending raw data wn2 to PROCESSED data in current folder under wn1
	Variable quiet
	Variable i, avgRatio=0, failed=0, uncert=0
	Make/o/n=(nSectors) rTerms=nan, uTerms=nan
	String bn1, bn2
	bn2=ParseFilePath(1,wn2,"_",1,1)//remove the angle and sector size info from the data name
	bn1=ParseFilePath(1,wn1,"_",1,!appnd)  //if appending doesn't have sector size at the end
	For( i=0; i<nSectors; i+=1 )
		If( appnd )
			wn1= SelectString( sectSize==360, (bn1+num2str(startAngle+i*45)), wn1 ) //reg vs Circ profiles
			WAVE/Z y1=$( wn1+"_yTmp" ), x1=$( wn1+"_xTmp"), u1=$(wn1+"_uTmp")
		else
			wn1= SelectString( sectSize==360, (bn1+num2str(startAngle+i*45)+"_"+ReplaceString(".", num2str(sectSize), "p")), wn1 )
			WAVE/Z y1=$("root:SAS:'"+wn1+"':'r_"+wn1+"'"), x1=$("root:SAS:'"+wn1+"':'q_"+wn1+"'"), u1=$("root:SAS:'"+wn1+"':'s_"+wn1+"'")
		endif
		wn2= SelectString( sectSize==360, (bn2+num2str(startAngle+i*45)+"_"+ReplaceString(".", num2str(sectSize), "p")), wn2 )
		WAVE/Z y2=$("root:SAS:'"+wn2+"':'r_"+wn2+"'"), x2=$("root:SAS:'"+wn2+"':'q_"+wn2+"'"), u2=$("root:SAS:'"+wn2+"':'s_"+wn2+"'")
		
		If( !WaveExists(y1) || !WaveExists(x1) || !WaveExists(y2) || !WaveExists(x2) )
			failed+=1 //this error will already be reported
			continue
		elseif( !Range(x1,qpt,tol=3) )
			If( !quiet )
				Printf "Stitch Ratio Warning: %g out of range for %s [%1.4f,%1.4f].\r", qpt, wn1,x1[0],x1[inf]
			endif
			failed+=1
			continue
		elseif( !Range(x2,qpt,tol=3) )
			If( !quiet )
				Printf "Stitch Ratio Warning: %g out of range for %s [%1.4f,%1.4f]\r", qpt, wn2,x2[0],x2[inf]
			endif
			failed+=1
			continue
		endif
		rTerms[i]=interp( qpt, x1, y1) / interp(qpt, x2, y2)
		uTerms[i]=rTerms[i]*sqrt( (interp(qpt,x1,u1)/interp(qpt,x1,y1))^2 + (interp(qpt,x2,u2)/interp(qpt,x2,y2))^2 )
		avgRatio+=interp( qpt, x1, y1) / interp(qpt, x2, y2)
		uncert+= (interp(qpt,x1,u1)/interp(qpt,x1,y1))^2 + (interp(qpt,x2,u2)/interp(qpt,x2,y2))^2
	Endfor
	If( failed<nSectors )
		return wAvg(rTerms,uTerms)
		avgRatio/=nSectors-failed
		If( nSectors-failed > 2 ) //use population stats for uncert if the pop is at least 3
			wavestats/Q rTerms
			Uncert= min( sqrt(uncert)*avgRatio, V_sdev )/sqrt(nSectors-failed)
		else // if there are not enough to do population stats, use error propagation
			Uncert=sqrt(uncert)*avgRatio/sqrt(nSectors-failed)
		endif
	else
		printf "Stitch Ratio Error: for %s and %s!\r", bn1, bn2
		return -1
	endif
	return cmplx(avgRatio,Uncert)
end


//Stitches two specific waves, saves the waves in Outstr+_proc folder, and displays the data if desired
Function Stitch(wn1,wn2, OutStr, qpt, a, [d,qPwr,a2nm,sm,q,appnd,rescale])
	String wn1, wn2, OutStr
	Variable qpt //q-value where the stitch will take place
	Variable/C a // muliplication factor for second dataset to match the first
	Variable d //1 = new dispaly, 2 = append to top graph
	Variable qPwr //multiplies the data by q^qPwr - DEFAULT = 2
	Variable a2nm //transforms the q from angrstroms to nm - DEFAULT = 1 (Yes, transform it)
	Variable sm //smooths the profile data - DEFAULT = no smooth
	Variable q //don't print output to history
	Variable appnd //appending raw data wn2 to stitched data in current folder under wn1
	Variable rescale //rescale stiched data to second data in pair rather than the first
	
	String noteStr
	Variable i, j, R=real(a), uR=imag(a)
	// Load data: first wave loaded depends on if we're appeding to an already stitched dataset
	WAVE y1=$(Selectstring( appnd, "root:SAS:'"+wn1+"':'r_"+wn1+"'", wn1+"_yTmp" ))
	WAVE x1=$(SelectString( appnd, "root:SAS:'"+wn1+"':'q_"+wn1+"'", wn1+"_xTmp" ))
	WAVE u1=$(SelectString( appnd, "root:SAS:'"+wn1+"':'s_"+wn1+"'", wn1+"_uTmp" ))
	WAVE ua1=$(SelectString( appnd, "root:SAS:'"+wn1+"':'s_"+wn1+"'", wn1+"_uaTmp" ))
	WAVE y2=$("root:SAS:'"+wn2+"':'r_"+wn2+"'"), x2=$("root:SAS:'"+wn2+"':'q_"+wn2+"'"), u2=$("root:SAS:'"+wn2+"':'s_"+wn2+"'")
	If( numpnts(x2)<1 )  //if there is a wave but there's no data in it, don't stitch
		return 0
	endif
	// Stitch the data
	Variable pt1= BinarySearch(x1,qpt), pt2= BinarySearch(x2,qpt)
	IF( pt1<0 )
		Printf "Stitch Warning: Q=%g out of range for '%s'[%g,%g]\r", qpt, NameOfWave(y1),x1[0],x1[inf]
		pt1=numpnts(x1) //just use all the points
	elseif( pt2<0 )
		Printf "Stitch Warning: Q=%g out of range for '%s'[%g,%g]\r", qpt, NameOfWave(y2),x2[0],x2[inf]
		pt2=0 //just start at the begining
	endif
	Variable nPts= pt1 + numpnts(x2) - pt2
	Duplicate/d/o y1 $(Outstr+"_y")
	Duplicate/d/o x1 $(Outstr+"_x")
	Duplicate/d/o u1 $(Outstr+"_u")
	Duplicate/d/o ua1 $(Outstr+"_uAni") // don't include Ratio uncert. for anisotropy uncert. calculation
	WAVE yOut=$(Outstr+"_y"), xOut=$(Outstr+"_x"), uOut=$(Outstr+"_u"), uaOut=$(Outstr+"_uAni")
	Redimension/n=(npts) yOut, xOut, uOut, uaOut
	//Save parameters in the wave note
	sprintf noteStr, "s0Name:%s;s1Name:%s;End0:%d;Start1:%d;A0:%g;uA0:%g;", wn1, wn2, pt1, pt2, real(a), imag(a)
	NOTE/K/NOCR yOut, noteStr
	NOTE/K/NOCR xOut, noteStr
	NOTE/K/NOCR uOut, noteStr
	//Fill from the second wave
	xOut[pt1,numpnts(xOut)]=x2[p-pt1+pt2]
	If( rescale )
		yOut/=R
		uOut=yOut * sqrt( (u1/y1)^2 + (uR/R)^2 )
		uaOut=u1/R
		yOut[pt1,numpnts(xOut)]=y2[p-pt1+pt2]
		uOut[pt1,numpnts(xOut)]=u2[p-pt1+pt2]//yOut * sqrt( (u2[p-pt1+pt2]/y2[p-pt1+pt2])^2 + (uR/R)^2 )
		uaOut[pt1,numpnts(xOut)]=u2[p-pt1+pt2]
	else
		yOut[pt1,numpnts(xOut)]=y2[p-pt1+pt2]*R
		uOut[pt1,numpnts(xOut)]=yOut * sqrt( (u2[p-pt1+pt2]/y2[p-pt1+pt2])^2 + (uR/R)^2 )
		uaOut[pt1,numpnts(xOut)]=u2[p-pt1+pt2]*R
	endif
	If( !q )
		printf "Stitched '%s' to '%s'. Saved as '%s' with ratio=%s\r", wn1, wn2, OutStr+"_y", FormatNumber(real(a),imag(a))
	endif
end	


Function PlotProf( yw, xw, uw, s1, s2, s3, new, style, Cphase, num, resetC,NoU,[hide, d, qpwr thickness,pols])
	WAVE yw, xw, uw
	String s1, s2, s3 //name of the display
	String style //user-applied graph style
	Variable new //force new display
	Variable NoU //add/take off the uncertainty bars
	Variable Cphase, num, resetC // controls the color scheme
	Variable hide //hides the window once it is created
	Variable d //standard display variable (only cares if it d=2)
	Variable qPwr // intensities are multiplied by Q^qPwr
	Variable thickness // nonzero if normalized to thickness
	WAVE pols //polarizations for the plotted data

	Variable i
	
	String wName=SearchWindowName(s1,s2,s3) //search for existing windows
	String styleMacro
	
	String lgndtxt, str
	If( (new==1 || strLen(wName)==0) && d!=2 ) //make a new window
		If( !strLen(s1) )
			return 0 //error with s1 name
		elseif( !strLen(s2) )
			wName=uniquename( s1+"_", 6, 0)
		elseif( !strLen(s3) ) 
			wName=uniquename( s1+"_"+s2+"_", 6, 0)
		else
			wName=uniquename( s1+"_"+s2+"_"+s3+"_", 6, 0)		
		endif
		
		variable tNW=startMStimer //New Window Timer
		display/W=(5,40,300,280)/N=$wname as wname

		Dowindow/C $wname
		doWindow/hide=(hide) $wname
		sprintf lgndTxt, "\\f01%s %s %s\\f00", s1, s2, s3
		TextBox/W=$wname/C/N=lgnd lgndTxt
		
		BAC_AppendToLayout(s1,wname, w=230, h=216)
	elseif( d==2 )
		wName=WinName(0,1)
	endif
	
	DoWindow $wname
	
	If( !Plotted(yw,wname) ) //don't plot again if already plotted
		appendtograph/W=$wname yw vs xw
		ModifyGraph/W=$wname/Z grid=2,tick=2,mirror=1,standoff=0,muloffset={1,1} //very basic styling
		ModifyGraph/W=$wname/Z mode=4,marker=19,msize=1,mskip=5
		ModifyGraph/W=$wname/Z margin(left)=40,margin(bottom)=36,margin(top)=29,margin(right)=14
		ModifyGraph logHTrip(left)=100,logLTrip(left)=0.01
		Str=RemoveFromList(s1,NameOfWave(yw),"_")
		Str=ParseFilePath(1,Str,"_",1,0)
		Str=RemoveEnding(Str,"_")
		Str=ReplaceString("_", Str, " ")
		sprintf lgndtxt, "\\s('%s') %s", NameOfWave(yw), Str
		AppendText/W=$wname/N=lgnd lgndtxt
		GroupStyle(phase=Cphase, n=num, wn=wname, resetC=resetC)
		If( strLen(style) )
			If( StringMatch( style, "profStyle") )
				sprintf styleMacro "%s(wn=\"%s\", qPwr=%d, th=%d)", removeEnding(style,"()"), wname, qPwr, thickness>0
			elseif( StringMatch( style, "aniStyle") )
				sprintf styleMacro "%s(wn=\"%s\",pols=%s)", removeEnding(style,"()"), wname, GetWavesDataFolder(pols,2)
			else
				sprintf styleMacro "%s(wn=\"%s\")", removeEnding(style,"()"), wname
			endif			
			Execute styleMacro
		endif
	endif
	
	If( NoU ) //add or turn off uncertainty bars
		ErrorBars/W=$wname $NameOfWave(yw) OFF
	else
		ErrorBars/W=$wname/T=0 $NameOfWave(yw) Y, wave=(uw,uw)
	endif
End

//Updates intensity axis depending on the qPwr and thickness correction
Function/S UpdateIntAxis(qpwr,th,wn)
	Variable qpwr, th
	String wn
	String intLabel
	IF( qpwr>0 ) //update Intensity axis
		sprintf intLabel, "Intensity*Q\\S%d\\M [sr\\S-1\\Mnm\\S%d\\M]", qpwr, -qpwr-th
		Label/W=$wn/Z left intLabel
		If( th )
			SetAxis/W=$wn left 1e-10,1e-07
		else
			SetAxis/W=$wn left 1e-08,1e-05
		endif
	else
		If( th )
			intLabel="Intensity [sr\\S-1\\Mcm\\S-1\\M]"
			Label/W=$wn/Z left "Intensity [sr\\S-1\\Mcm\\S-1\\M]"
			//assumes thickness is ~100 nm or 1e-5 cm
			SetAxis/W=$wn left 1e+1, 1e+6
		else
			intLabel="Intensity [sr\\S-1\\M]"
			Label/W=$wn/Z left "Intensity [sr\\S-1\\M]"
			SetAxis/W=$wn left 1e-4, 1e+1
		endif
	endif
	return intLabel
end


//Plots the circularly averaged profiles for each energy in an upper plot and the anisotropy in a lower plot
Function PlotSummary(dname,qpwr,th)
	String dName
	Variable qpwr, th //ways the intensity may have been normalized (I*Q^qPwr/thickness) - th is a 0 or 1 (not actual film thickness!)
	String CurrentFolder=GetDataFolder(1), wn=dname+"_Summary", avgN, aniN, lgndTxt, Str
	String dFldr=removeending(GetFldrNameASmeta(dName),":metadata")
	SetDataFolder $(dFldr+":metadata")
	WAVE energies0
	Variable i, nE=numUnique(Energies0,0.1), Ectr=0
	string intLabel
	SortMetaData(0,type="Energy")
	DoWindow/F $wn
	If( V_Flag ) //don't make new one if it exsists already
		UpdateIntAxis(qpwr,th,wn)
		return 0
	endif
	Display /W=(5,40,300,450)/N=$wn as wn
	sprintf lgndTxt, "\\JC\\f01%s\\f00", dname
	TextBox/C/N=lgnd lgndTxt
	AppendText/N=lgnd "\\JL"
//	SetDataFolder $("root:"+dName+"_data:Profiles")
	SetDataFolder $(dFldr+":Profiles")
	For( i=0; i<ne; i+=1 )
		sprintf avgN, "%s_%1.1f_C", dName, energies0[i]
//		sprintf avgN, "%s_%1.1f_45", ReplaceString("Sect",dName,""), energies0[i]
		WAVE/Z avgY=$(avgN+"_y"), avgX=$(avgN+"_x"), avgU=$(avgN+"_u")
		IF( WaveExists(avgY) && WaveExists(avgX) )
			appendtograph/W=$wn avgy vs avgx
			ErrorBars/W=$wn/T=0 $NameOfWave(avgY) Y, wave=(avgU, avgU)
			sprintf lgndtxt, "\\s('%s') %1.1f", NameOfWave(avgY), Energies0[i]
			AppendText/W=$wn/N=lgnd lgndtxt
			Ectr+=1
		else
			print "Summary Error: Couldn't find "+avgN
		endif
	Endfor
	SetDataFolder $(dFldr+":Anisotropy")
	For( i=0; i<ne; i+=1 ) //have to do it this way to get the order right
		sprintf aniN, "%s_%1.1f_ani", dName, energies0[i]
		WAVE/Z aniY=$(aniN+"y"), aniX=$(aniN+"x"), avgU=$(aniN+"u")
		If( WaveExists(aniY) && WaveExists(aniX) )
			appendtograph/W=$wn/L=ani aniy vs anix
			ErrorBars/W=$wn/T=0 $NameOfWave(aniY) Y, wave=(avgU, avgU)
		else
			print "Summary Error: Couldn't find "+aniN
		endif
	Endfor
	Add_2piOverQticks2top()
	ModifyGraph/W=$wn margin(left)=40,margin(bottom)=36,margin(top)=29,margin(right)=14
	ModifyGraph/W=$wn lSize=1.5, grid=2, log(left)=1,log(bottom)=1,log(MT_bottom)=1,tick=2
	ModifyGraph/W=$wn mirror(left)=1,mirror(ani)=1, minor(ani)=1, standoff=0
	ModifyGraph/W=$wn lblMargin(MT_bottom)=6, lblPosMode(left)=1,lblPosMode(MT_bottom)=1,lblPosMode(ani)=1
	ModifyGraph/W=$wn lblPos(left)=39, lblLatPos(MT_bottom)=-40, lblLatPos(bottom)=-40, freePos(MT_bottom)=0, freePos(ani)=0
	ModifyGraph/W=$wn axisEnab(left)={0.51,1}, axisEnab(ani)={0,0.49}, zero(ani)=1
	ModifyGraph/W=$wn logHTrip(left)=100,logLTrip(left)=0.01
	Label/W=$wn bottom "Q [nm\\S-1\\M]"
	Label/W=$wn MT_bottom "2\\F'Symbol'p\\F'Arial'/q [nm]"
	Label/W=$wn ani "Anisotropy [%]"
	UpdateIntAxis(qpwr,th,wn)
	SetAxis/W=$wn bottom 0.005,1
	SetAxis/W=$wn ani -20,60
	GroupStyle(phase=45, n=1, wn=wn, resetC=Ectr)
	DoWindow $dName
	If( V_Flag )
		Str= LayoutInfo(dName, num2str(i) )
		Str= StringbyKey("NAME",str)
		if( StrLen(str)==0 ) //window is not in the layout
			AppendLayoutObject/W=$dName graph $wn
			ModifyLayout/W=$dName frame($wn)=0, mag=1, width($wn)=252, height($wn)=324
		endif
	endif
	BAC_AppendToLayout(dName,wn, w=230, h=324)
	SetDataFolder $CurrentFolder
end

//appends a graph to a specified layout
Function BAC_AppendToLayout(Lname,Wname,[w,h])
	String Lname, Wname
	Variable w,h
	w=ParamIsDefault(w) ? 252 : w
	h=ParamIsDefault(h) ? 216 : h
	DoWindow $Lname //Look for layout
	String str
	Variable i
	If( V_Flag ) //Layout exists
		For( i=0; ; i+=1 )
			Str= LayoutInfo(Lname, num2str(i) )
			Str= StringbyKey("NAME",str)
			If( StringMatch(Str,Wname) ) //window is already in the layout
				break
			elseif( StrLen(str)==0 ) //window is not in the layout
				AppendLayoutObject/W=$Lname graph $wname
				ModifyLayout/W=$Lname frame($wname)=0, mag=1, width($wname)=w, height($wname)=h
				break
			endif
		endfor
	endif
end


//Searches current window names based on the strings provided and retuns the name if it finds it.  If not, returns an empty string.
Function/S SearchWindowName(s1, s2, s3)
	String s1, s2, s3
	String List
	Variable i=0
	If( !strLen(s1) )
		return ""
	elseif( !strLen(s2) )
		List=winList(s1+"_*",";","")
	elseif( !strLen(s3) )
		List=winList(s1+"_"+s2+"_*",";","")
	else
		List=winList(s1+"_"+s2+"_"+s3+"_*",";","")
	endif
	Do
		If( Str2Num(ParseFilePath(0,StringFromList(i,List),"_",1,0))*0!=0 ) //rest of win name is not a number
			List=RemoveListItem(i,List)
		else
			i+=1
		endif
	While( i<itemsInList(List) )
	Return StringFromList(0,SortList(List,";",1))
end


//Combines (averages) the anisotropy profiles from Horizontal and Vertical polarizations to eliminate parasitic scatter from the signal.
Function CombineAniProfs(bName,[d])
	String bName
	variable d
	string P100List=waveList(bName+"_100_*_ani",";",""), w100n, w190n, aName, wName
	variable i
	For( i=0; i<ItemsInList(P100List); i+=1 )
		w100n=StringFromList(i,P100List)
		w190n=replacestring("_100_",w100n,"_190_")
		aName=replacestring("_100_",w100n,"_")
		WAVE/Z w100=$w100n, w100q=$(w100n+"q"), w100u=$(w100n+"u")
		WAVE/Z w190=$w190n, w190q=$(w190n+"q"), w190u=$(w190n+"u")
		If( WaveExists(w100) && WaveExists(w190) )
			CombineProfiles(w100,w100q,w100u,w190,w100u,w190q,f=-1)
			Duplicate/d/o AvY, $aName; Duplicate/d/o AvX, $(aName+"Q"); Duplicate/d/o AvU, $(aName+"u")
			WAVE ani=$aName, aniQ=$(aName+"Q"), aniU=$(aName+"u")
			If( i==0 && d==1 )
				wName=uniquename( bName+"_Ani_", 6, 0)
				display/N=$wname as wname
			endif
			If( d==1 || d==2 )
				appendtograph/W=$wName ani vs aniQ
			endif
			Add_2piOverQticks2top()
			ModifyGraph/W=$wname grid=2,tick=2,mirror(left)=1,standoff=0
			Legend/C/N=lgnd
		endif
	Endfor
end


//Filters out duplicates (same energies, angles, & polarizations) based on the dwell or name
Function ApplyUserFilter(index, mode, Pols, Energies, Angles, Dwells, Scans)
	WAVE Pols, Energies, Angles, Dwells
	WAVE/T Scans
	Variable index
	Variable mode// -1=shortest dwell; -2=longest dwell; -3=first name (alpha); -4=last name (alpha); #>0=stitch together at #
	Variable En=Energies[index], th=Angles[index], pol=Pols[index], index2=index
	Do  //find the end of the set of scans at the same energy & angle and polarization (but multiple dwells/scan names)
		index2+=1
	While ( abs(en-Energies[index2])<0.1 && abs(th-Angles[index2])<1 && abs(pol-Pols[index2])<1 && index2<numpnts(energies)-1 )
	//duplicate the subset of the scans that have the same energy
	Make/o/n=(index2-index) tmpPol, tmpEn, tmpAng, tmpDwell
	Make/t/o/n=(index2-index) tmpScn
	tmpPol=Pols[index+p]
	tmpEn=Energies[index+p]
	tmpAng=Angles[index+p]
	tmpDwell=Dwells[index+p]
	tmpScn=Scans[index+p]

	//Eventually this part will be in if statements baed on MODE, but for now get rid of short dwells
	Sort {tmpDwell} tmpScn, tmpPol, tmpEn, tmpAng, tmpDwell
	variable dwell=tmpDwell[index2-index], i, j
	For( i=index2-index-1; i>=0; i-=1 )  //delete the short dwells from the sublist
		If( tmpDwell[i]!=dwell )
			DeletePoints i, 1, tmpScn, tmpPol, tmpEn, tmpAng, tmpDwell
		endif
	Endfor
	Sort {tmpAng}, tmpScn, tmpPol, tmpEn, tmpAng, tmpDwell //resort for angle
	DeletePoints index, index2-index, Pols, Energies, Angles, Dwells, Scans  //delete the subset from the main list
	variable newPts=numpnts(tmpDwell)
	InsertPoints index, newPts, Pols, Energies, Angles, Dwells, Scans
	For( i=index, j=0; i<index+newPts; i+=1, j+=1 )  //put the data back into the main list
		Pols[i]=tmpPol[j]
		Energies[i]=tmpEn[j]
		Angles[i]=tmpAng[j]
		Dwells[i]=tmpDwell[j]
		Scans[i]=tmpScn[j]
	Endfor
	return index //eventually might stitch the two together here first and then skip these.
end


//Default style for displaying scattering profiles (log-log plot, axis labels, etc.)
Function ProfStyle([wn,qpwr,th]) : GraphStyle
	string wn
	Variable qpwr, th  //possible intensity normalizations that were applied to the data reflected in Left axis lable and range
	wn= SelectString( ParamIsDefault(wn), wn, WinName(0,1) )
	Add_2piOverQticks2top()
	ModifyGraph/W=$wn/Z grid=2, log=1, tick=2, mirror(left)=1, mirror(bottom)=0, lblMargin(MT_bottom)=6, standoff=0
	ModifyGraph/W=$wn/Z lblPosMode(MT_bottom)=1, lblLatPos(MT_bottom)=-40, freePos(MT_bottom)=0
	ModifyGraph/W=$wn/Z mirror(left)=1, mirror(bottom)=0, mirror(MT_bottom)=0
	UpdateIntAxis(qpwr,th,wn)
	Label/W=$wn/Z bottom "Q [nm\\S-1\\M]"
	Label/W=$wn/Z MT_bottom "2\\F'Symbol'p\\F'Arial'/q [nm]"
	SetAxis/W=$wn/Z bottom 0.005,1
	ModifyGraph/W=$wn/Z lSize=1.5
End


//Default style for displaying anisotropy (log-lin plot, axis labels, etc.)
Function aniStyle([wn,pols]) : GraphStyle
	string wn
	WAVE pols
	wn= SelectString( ParamIsDefault(wn), wn, WinName(0,1) )
	Add_2piOverQticks2top()
	ModifyGraph/W=$wn/Z margin(left)=36,margin(bottom)=29,margin(top)=22,margin(right)=14,lSize=1.5
	variable i
	For( i=0; i<numpnts(pols); i+=1 )
		ModifyGraph/W=$wn/Z muloffset[i]={1,pols[i]==100 ? 1 : -1}
	endfor
	ModifyGraph/W=$wn/Z mirror(left)=1, mirror(bottom)=0, mirror(MT_bottom)=0, grid=2, standoff=0
	ModifyGraph/W=$wn/Z log(bottom)=1,log(MT_bottom)=1,log(left)=0
	ModifyGraph/W=$wn/Z tick=2, zero(left)=1,mirror(left)=1,mirror(bottom)=2,minor(left)=1,minor(bottom)=1
	ModifyGraph/W=$wn/Z lblMargin(MT_bottom)=6
	ModifyGraph/W=$wn/Z lblPosMode(MT_bottom)=1, lblLatPos(MT_bottom)=-40, freePos(MT_bottom)=0
	Label/W=$wn/Z left "Anisotropy [%]"
	Label/W=$wn/Z bottom "Q [nm\\S-1\\M]"
	Label/W=$wn/Z MT_bottom "2\\F'Symbol'p\\F'Arial'/q [nm]"
	SetAxis/W=$wn/Z left -20,60
	SetAxis/W=$wn/Z bottom 0.005,1
End


//returns whether user wants to use this theta.
Function useTh(th, th0, num)
	Variable th, th0, num
	Variable i
	for( i=0; i<num; i+=1 )
		If( abs(th-(th0+i*45))<1 )
			return 1
		endif
	endfor
	return 0
end


//returns a string representing the complete path to the wave representing NIKA processed data in the SAS folder
Function/S SASref(en,pol,sect,type,[fldr])
	Variable en,pol,sect,type //energy, polarization, sector, and type of the desired wave - type: -1=uncert., 0=data, +1=q-values
	String fldr //folder containing the Metadata
	fldr= Selectstring( ParamIsDefault(fldr), fldr, GetDataFolder(1))
	Variable idx=QueryMeta(sect, en, Pol, fldr)
	If( idx==-1 )
		return ""
	endif
	WAVE/T scans=$(fldr+":scans0")
	string refS="root:SAS:'"+scans[idx]+SelectString(type,"':'s_","':'r_","':'q_")+scans[idx]+"'"
	return refS
end
	
	
//combine sectors (hozontal, vertical and 45) if multiple exist, and rename
//can handle one or two H & V sectors and up 1-4 45 deg sectors
Function CombineSectors(bName, th0, n,[en,pol,mf,outf])
	String bname
	Variable th0, n //first angle & number of sectors, respectively
	Variable en, pol //get the waves directly from the SAS folder rather than prev. processed data in current folder
	string mf, outf //mf=metadata folder, outf=output folder
	String Currentfolder=GetDataFolder(1)
	mf= Selectstring( ParamIsDefault(mf), mf, CurrentFolder)
	outf= Selectstring( ParamIsDefault(outf), outf, CurrentFolder)
	SetDataFolder $outf

	//Combine vertical sectors based on availablility and user selection
	If( !ParamIsDefault(en) && !ParamIsDefault(pol)  )
		WAVE/Z x90=$SASref(en,pol,90,1,fldr=mf), y90=$SASref(en,pol,90,0,fldr=mf), u90=$SASref(en,pol,90,-1,fldr=mf)
		WAVE/Z x270=$SASref(en,pol,270,1,fldr=mf), y270=$SASref(en,pol,270,0,fldr=mf), u270=$SASref(en,pol,270,-1,fldr=mf)
	else
		WAVE/Z x90=$(Selectstring(useTh(90,th0,n), "", bName+"_90_x")), x270=$(Selectstring(useTh(270,th0,n), "", bName+"_270_x"))
		WAVE/Z y90=$(bName+"_90_y"), y270=$(bName+"_270_y"), u90=$(bName+"_90_uAni"), u270=$(bName+"_270_uAni")
	endif
	If( !WaveExists(x90) && !WaveExists(x270) )
		print "Couldn't find any horizontal profiles! ("+bName+"_90_x OR "+bname+"_270_x)"
		return -1
	elseif( !WaveExists(x270) )
		Duplicate/d/o x90 $(bName+"_Vx"); Duplicate/d/o y90 $(bName+"_Vy"); Duplicate/d/o u90 $(bName+"_Vu")
	elseif( !WaveExists(x90) )
		Duplicate/d/o x270 $(bName+"_Vx"); Duplicate/d/o y270 $(bName+"_Vy"); Duplicate/d/o u270 $(bName+"_Vu")
	else
		CombineProfiles(y90,x90,u90,y270,x270,u270)
		Duplicate/d/o AvX $(bName+"_Vx"); Duplicate/d/o AvY $(bName+"_Vy"); Duplicate/d/o Avu $(bName+"_Vu")
	endif
	
	//Combine horizontal sectors based on availablility and user selection
	If( !ParamIsDefault(en) && !ParamIsDefault(pol)  )
		WAVE/Z x0=$SASref(en,pol,0,1,fldr=mf), y0=$SASref(en,pol,0,0,fldr=mf), u0=$SASref(en,pol,0,-1,fldr=mf)
		WAVE/Z x180=$SASref(en,pol,180,1,fldr=mf), y180=$SASref(en,pol,180,0,fldr=mf), u180=$SASref(en,pol,180,-1,fldr=mf)
	else
		WAVE/Z x0=$(Selectstring(useTh(0,th0,n), "", bName+"_0_x")), x180=$(Selectstring(useTh(180,th0,n), "", bName+"_180_x"))
		WAVE/Z y0=$(bName+"_0_y"), y180=$(bName+"_180_y"),u0=$(bName+"_0_uAni"), u180=$(bName+"_180_uAni")
	endif
	If( !WaveExists(x0) && !WaveExists(x180) )
		print "Couldn't find any vertical profiles! ("+bName+"_0_x OR "+bname+"_180_x)"
		return -1
	elseif( !WaveExists(x180) )
		Duplicate/d/o x0 $(bName+"_Hx"); Duplicate/d/o y0 $(bName+"_Hy"); Duplicate/d/o u0 $(bName+"_Hu")
	elseif( !WaveExists(x0) )
		Duplicate/d/o x180 $(bName+"_Hx"); Duplicate/d/o y180 $(bName+"_Hy"); Duplicate/d/o u180 $(bName+"_Hu")
	else
		CombineProfiles(y0,x0,u0,y180,x180,u180)
		Duplicate/d/o AvX $(bName+"_Hx"); Duplicate/d/o AvY $(bName+"_Hy"); Duplicate/d/o AvU $(bName+"_Hu")
	endif
	
	//Combine diagonal sectors based on availablility and user selection
	If( !ParamIsDefault(en) && !ParamIsDefault(pol)  )
		WAVE/Z x45=$SASref(en,pol,45,1,fldr=mf), y45=$SASref(en,pol,45,0,fldr=mf), u45=$SASref(en,pol,45,-1,fldr=mf)
		WAVE/Z x135=$SASref(en,pol,135,1,fldr=mf), y135=$SASref(en,pol,135,0,fldr=mf), u135=$SASref(en,pol,135,-1,fldr=mf)
		WAVE/Z x225=$SASref(en,pol,225,1,fldr=mf), y225=$SASref(en,pol,225,0,fldr=mf), u225=$SASref(en,pol,225,-1,fldr=mf)
		WAVE/Z x315=$SASref(en,pol,315,1,fldr=mf), y315=$SASref(en,pol,315,0,fldr=mf), u315=$SASref(en,pol,315,-1,fldr=mf)
	else
		WAVE/Z x45=$(Selectstring(useTh(45,th0,n), "", bName+"_45_x")), x135=$(Selectstring(useTh(135,th0,n), "", bName+"_135_x"))
		WAVE/Z x225=$(Selectstring(useTh(225,th0,n), "", bName+"_225_x")), x315=$(Selectstring(useTh(315,th0,n), "", bName+"_315_x"))
		WAVE/Z y45=$(bName+"_45_y"), y135=$(bName+"_135_y"), y225=$(bName+"_225_y"), y315=$(bName+"_315_y")
		WAVE/Z u45=$(bName+"_45_uAni"), u135=$(bName+"_135_uAni"), u225=$(bName+"_225_uAni"), u315=$(bName+"_315_uAni")
	endif
	If( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x225) && !WaveExists(x315))
		print "Couldn't find any 45 deg profiles! ("+bName+"_45_x, OR 135, 225, 315)"
		return -1
	elseif( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x225) ) //only one
		Duplicate/d/o x315 $(bName+"_Dx"); Duplicate/d/o y315 $(bName+"_Dy"); Duplicate/d/o u315 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x225) && !WaveExists(x315) ) //only one
		Duplicate/d/o x135 $(bName+"_Dx"); Duplicate/d/o y135 $(bName+"_Dy"); Duplicate/d/o u135 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x135) && !WaveExists(x315) ) //only one
		Duplicate/d/o x225 $(bName+"_Dx"); Duplicate/d/o y225 $(bName+"_Dy"); Duplicate/d/o u225 $(bName+"_Du"); return 0
	elseif( !WaveExists(x135) && !WaveExists(x225) && !WaveExists(x315) )//only one
		Duplicate/d/o x45 $(bName+"_Dx"); Duplicate/d/o y145 $(bName+"_Dy"); Duplicate/d/o u145 $(bName+"_Du"); return 0
	elseif( !WaveExists(x45) && !WaveExists(x135) ) //only two
		CombineProfiles(y225, x225,u225, y315,x315,u315)
	elseif( !WaveExists(x45) && !WaveExists(x225) ) //only two
		CombineProfiles(y135, x135, u135, y315,x315,u315)
	elseif( !WaveExists(x45) && !WaveExists(x315) ) //only two
		CombineProfiles(y135, x135, u135, y225,x225,u225)
	elseif( !WaveExists(x135) && !WaveExists(x225) ) //only two
		CombineProfiles(y45, x45, u45,y315,x315,u315)
	elseif( !WaveExists(x135) && !WaveExists(x315) ) //only two
		CombineProfiles(y45, x45, u45,y225,x225,u225)
	elseif( !WaveExists(x225) && !WaveExists(x315) ) //only two
		CombineProfiles(y45, x45, u45,y135,x135,u135)
	elseif( !WaveExists(x45) ) //only three
		CombineProfiles(y135,x135,u135,y225,x225,u225)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU, y315,x315,u315)
	elseif( !WaveExists(x135) ) //only three
		CombineProfiles(y45,x45,u45,y225,x225,u225)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y315,x315,u315)
	elseif( !WaveExists(x225) ) //only three
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y315,x315,u315)
		elseif( !WaveExists(x315) ) //only three
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX; Duplicate/d/o AvY tmpY; Duplicate/d/o AvU tmpU
		CombineProfiles(tmpY, tmpX, tmpU,y225,x225,u225)
	else //all 4 are present
		CombineProfiles(y45,x45,u45,y135,x135,u135)
		Duplicate/d/o AvX tmpX1; Duplicate/d/o AvY tmpY1; Duplicate/d/o AvU tmpU1
		CombineProfiles(y225,x225,u225,y315,x315,u315)
		Duplicate/d/o AvX tmpX2; Duplicate/d/o AvY tmpY2; Duplicate/d/o AvU tmpU2
		CombineProfiles(tmpY1,tmpX1,tmpU1,tmpY2,tmpX2,tmpU2)
	endif
	Duplicate/d/o AvX $(bName+"_Dx"); Duplicate/d/o AvY $(bName+"_Dy"); Duplicate/d/o AvU $(bName+"_Du")
	
	SetDataFolder $CurrentFolder
End

//Averages two supposedly identical sets of data that have different x-ranges
//averages the overlaping range and takes data from datasets having an extended range
//x-data is indexed in a log scale
Function CombineProfiles(y1, x1, u1, y2, x2, u2,[f,w])
	wave y1, x1, u1, y2, x2, u2
	variable f //factor to multiply by y2
	Variable w //weighted average
	f= ParamIsDefault(f) ? 1 : f
	w= ParamIsDefault(w) ? 1 : w
	Variable LoX=min(x1[0],x2[0]), HiX=max(x1[numpnts(x1)],x2[numpnts(x2)]), npts=max(numpnts(x1),numpnts(x2))
	Make/o/n=(npts) AvX, AvY, AvU
	setscale/i x, log(LoX), log(HiX), AvX
	AvX=10^x; //spaces the new x-data evenly in a log scale
	If( w )
		AvY= Range(x1,AvX[p]) && Range(x2,AvX[p]) ? real(wAvg( {interp(AvX,x1,y1),f*interp(AvX,x2,y2)},{interp(AvX,x1,u1),interp(AvX,x2,u2)} )) : ( Range(x1,AvX[p]) ? interp(AvX,x1,y1) : f*interp(AvX,x2,y2) )
		AvU= Range(x1,AvX[p]) && Range(x2,AvX[p]) ? imag(wAvg( {interp(AvX,x1,y1),f*interp(AvX,x2,y2)},{interp(AvX,x1,u1),interp(AvX,x2,u2)} )) : ( Range(x1,AvX[p]) ? interp(AvX,x1,u1) : f*interp(AvX,x2,u2) )
	else
		AvY= Range(x1,AvX[p]) && Range(x2,AvX[p]) ? (interp(AvX,x1,y1)+f*interp(AvX,x2,y2))/2: nan
		AvU= Range(x1,AvX[p]) && Range(x2,AvX[p]) ? Max(interp(AvX,x1,u1),interp(AvX,x2,u2)) : nan
	endif
	string n1=note(y1), n2=note(y2), combNote
	Note/k/NOCR AvY, n1
	Note/k/NOCR AvX, n1
//	combNote="!!!!!Source1:"+NameOfWave(y1)+";"+n1+";"
////	sprintf combNote, "!!!!!Source1:%s;%s;", NameOfWave(y1),n1
//	Note/k/NOCR AvY, combNote
//	Note/k/NOCR AvX, combNote
//	combNote="!!!!!Source2:"+NameOfWave(y2)+";"+n2+";"
////	sprintf combNote, "!!!!!Source2:%s;%s;",NameOfWave(y2),n2
//	Note/NOCR AvY, combNote
//	Note/NOCR AvX, combNote
end

//Returns the truth of whether the value "val" is in the range of the wave "w" values.
//Assumes "w" values are ordered from most neg to most pos.
Function Range(w,val,[tol])
	wave w
	variable val, tol //tolerances of how far out of range the value can be while believing the interpolations
	tol= ParamIsDefault(tol) ? 0.1 : tol // default10% of the first and last spacing (pretty conservative)
	variable tol0=abs(w[1]-w[0])*tol, tolEnd=abs(w[numpnts(w)-1]-w[numpnts(w)-2])*tol
	return (w[0]-val)<tol0 && (val-w[numpnts(w)])<tolEnd ? 1 : 0
end

//Calculates the Anisotropy signal (Diff over Sum) for a given set of profiles
//looks in current folder for 90, 135, 180, 225, and 270 degree data that's already stitched & averages the duplicates (e.g. 90 & 270)
Function CalcAni(bName, [fldr,d,style,wname,pWinName,test,th0,nSectors,useD])
	String bName, fldr, style, wname, pWinName
	Variable d, test, th0,nSectors, useD
	th0= ParamIsDefault(th0) ? 0 : th0
	nSectors= ParamIsDefault(nSectors) ? 8 : nSectors
	If( ParamIsDefault(wname) )
		wname=WinName(0,1)
	endif
	String CurrentFolder=GetDataFolder(1)
	If( ParamIsDefault(fldr) )
		fldr=CurrentFolder+"Ani"
	endif

	//combine the appropriate waves from the list
//	Variable Success=CombineSectors(bName, th0, nSectors)
//	if( Success<0 )
//		return -1
//	endif
	WAVE Hy=$(bName+"_Hy"), Hx=$(bName+"_Hx"), Hu=$(bName+"_Hu")
	WAVE Vy=$(bName+"_Vy"), Vx=$(bName+"_Vx"), Vu=$(bName+"_Vu")
	WAVE Dy=$(bName+"_Dy"), Dx=$(bName+"_Dx"), Du=$(bName+"_Du")
	Variable int//get integral for the sectors
	Int=IntData(Hy,Hx,Hu); NOTE/NOCR Hy, "Int:"+num2str(int)+";"
	Int=IntData(Vy,Vx,Vu); NOTE/NOCR Vy, "Int:"+num2str(int)+";"
	Int=IntData(Dy,Dx,Du); NOTE/NOCR Dy, "Int:"+num2str(int)+";"
	Variable loH=max(Hx[0], Dx[0]), loV=max(Vx[0],Dx[0]) //Low Q
	Variable hiH=min(Hx[numpnts(Hx)],Dx[numpnts(Dx)]), hiV=min(Vx[numpnts(Vx)],Dx[numpnts(Dx)]) //high Q
	Variable nptsH=min(numpnts(Hx),numpnts(Dx)), nptsV=min(numpnts(Vx),numpnts(Dx))
	Variable lo=max(Hx[0],Vx[0]), hi=min(Hx[numpnts(Hx)],Vx[numpnts(Vx)])
	Variable npts=min(numpnts(Hx),numpnts(Vx))

	NewDataFolder/O/S $fldr
	If( !useD )
		Make/o/n=(npts) $(bName+"_aniY"), $(bName+"_aniX"), $(bName+"_aniU")
		WAVE ani=$(bName+"_aniY"), qw=$(bName+"_aniX"), aniU=$(bName+"_aniU")
		setscale/I x, log(lo), log(Hi), qw
		qw=10^x
		Ani= (interp(qw,Vx,Vy) - interp(qw,Hx,Hy)) / (interp(qw,Vx,Vy) + interp(qw,Hx,Hy)) *100
		AniU= ( 2*interp(qw,Hx,Hy)*interp(qw,Vx,Vu) )^2 
		AniU+= ( 2*interp(qw,Vx,Vy)*interp(qw,Hx,Hu) )^2 
		AniU= sqrt( AniU ) /(interp(qw,Hx,Hy)+interp(qw,Vx,Vy))^2 *100
	elseif (useD) 
		Make/o/n=(nptsH) $(bName+"_ani_Hy"), $(bName+"_ani_Hx"), $(bName+"_ani_Hu") //ani from diff btwn H & diagonal sectors
		WAVE aniH=$(bName+"_ani_Hy"), aniHx=$(bName+"_ani_Hx"), aniHu=$(bName+"_ani_Hu")
		setscale/I x, log(loH), log(HiH), aniHx
		aniHx=10^x
		aniH= (interp(aniHx,Dx,Dy) - interp(aniHx,Hx,Hy)) / interp(aniHx,Dx,Dy) *100
		aniHu= sqrt( (interp(aniHx,Hx,Hy)*interp(aniHx,Dx,Du))^2 + (interp(aniHx,Dx,Dy)*interp(aniHx,Hx,Hu))^2 )/ interp(aniHx,Dx,Dy)^2 *100
		
		Make/o/n=(nptsV) $(bName+"_ani_Vy"), $(bName+"_ani_Vx"), $(bName+"_ani_Vu") //ani from diff btwn V & diagonal sectors
		WAVE aniV=$(bName+"_ani_Vy"), aniVx=$(bName+"_ani_Vx"), aniVu=$(bName+"_ani_Vu")
		setscale/I x, log(loV), log(HiV), aniVx
		aniVx=10^x
		aniV= ( interp(aniVx,Vx,Vy) - interp(aniVx,Dx,Dy) ) / interp(aniVx,Dx,Dy) * 100
		aniVu= sqrt( (interp(aniVx,Vx,Vy)*interp(aniVx,Dx,Du))^2 + (interp(aniVx,Dx,Dy)*interp(aniVx,Vx,Vu))^2 ) / interp(aniVx,Dx,Dy)^2 * 100
		
		CombineProfiles(aniH,aniHx,aniHu, aniV,aniVx,aniVu)
		Duplicate/d/o avY $(bName+"_aniY"); Duplicate/d/o avX $(bName+"_aniX"); Duplicate/d/o avU $(bName+"_aniU")
		WAVE ani=$(bName+"_aniY"), qw=$(bName+"_aniX"), aniU=$(bName+"_aniU")
	endif
	String wnote=note(Hu)
//	sprintf wnote, "SourceH:%s;SourceD:%s;SourceV:%s",NameOfWave(Hy),NameOfWave(Dy),NameOfWave(Vy)
//	sprintf wnote, "SourceH:%s;SourceV:%s",NameOfWave(Hy),NameOfWave(Vy)
	Note/K/NOCR ani, wnote
	Note/K/NOCR qw, wnote
	Note/K/NOCR aniU, wnote
	
	//plot results if necessary
	If( d==1 )
		wName=uniquename( bName+"_AniPol_", 6, 0)
		display/N=$wname as wname
		pWinName= uniqueName( bName+"_Profs_",6,0)
		display/N=$pWinName as pWinName
	endif
	If( d==1 || d==2 )
		DoWindow/F $wname
		If( !plotted(ani,wname) )
			appendtograph/W=$wname ani vs qw //don't append if already there!!
		endif
		Add_2piOverQticks2top()
		ModifyGraph/W=$wname grid=2,log(bottom)=1,tick=2,mirror(left)=1,standoff=0
		Label/W=$wname left "Anisotropy [%]"
		Label/W=$wname bottom "Q [nm\S-1\M]"
		If( !ParamIsDefault(style) )
			Execute style+"()"
		endif
		DoWindow/F $pwinName
		If( !plotted(Hy,pwinName) )//don't append if already there!!
			appendtograph/W=$pWinName Hy vs Hx
			appendtograph/W=$pWinName Dy vs Dx
			appendtograph/W=$pWinName Vy vs Vx
		endif
		Add_2piOverQticks2top()
		GroupStyle(phase=0,n=3,wn=pWinName)
	elseif( d==3 )
		DoWindow/F $wname
		appendtograph/W=$wname/r ani vs qw
		ModifyGraph/W=$wname log(bottom)=1
		Label/W=$wname right "Anisotropy [%]"
		Label/W=$wname bottom "Q [nm\S-1\M]"
		If( !ParamIsDefault(style) )
			Execute style+"()"
		endif
	endif
	
	SetDataFolder $CurrentFolder
End


//extends the anisotropy data based on one sector from two orthogonal polarizations
Function ExtendAni(profStr,Qpt,fldr,quiet,[mode])
	String profStr, fldr
	Variable qpt, quiet
	Variable mode //  -1=have two sectors D and another (H or V); 0=only have H sector; 1=only have V sector;
	mode= ParamIsDefault(mode) ? 0 : mode
	String CurrentFolder=GetDataFolder(1)
	String w100n=replaceString("_190_",profStr,"_100_")
	String w190n=replaceString("_100_",profStr,"_190_")
	String aniN=replaceString("_100_",w100n,"_") + "_ani"
	WAVE/Z H100y=$(w100n+"_Hy"), H100x=$(w100n+"_Hx"), H100u=$(w100n+"_Hu")
	WAVE/Z D100y=$(w100n+"_Dy"), D100x=$(w100n+"_Dx"), D100u=$(w100n+"_Du")
	WAVE/Z V100y=$(w100n+"_Vy"), V100x=$(w100n+"_Vx"), V100u=$(w100n+"_Vu")
	WAVE/Z H190y=$(w190n+"_Hy"), H190x=$(w190n+"_Hx"), H190u=$(w190n+"_Hu")
	WAVE/Z D190y=$(w190n+"_Dy"), D190x=$(w190n+"_Dx"), D190u=$(w190n+"_Du")
	WAVE/Z V190y=$(w190n+"_Vy"), V190x=$(w190n+"_Vx"), V190u=$(w190n+"_Vu")
	IF( !WaveExists(H100y) || !WaveExists(V100y) || !WaveExists(H190y) || !WaveExists(V190y) )
		return 0
	elseif( qpt*0!=0 ) //qpt=nan b/c constPix conversion selected & user typed ExtAni in nm^-1 instead of A^-1
		print "AniExtension Error: ExtAni likely entered in nm^-1 instead of A^-1!"
		return 0
	endif
	Variable ratio, Uratio, i, j
	String noteStr, mStr=SelectString(mode, "_D","_H","_V") //based on the mode select which data we have to generate the correction factor "Ratio"
	mStr="_D"  //changed my mind: base the ratio on the Diagonal profile at the qPt designated by the ExtAni input since this is more robust btwn polarizaitons
	WAVE/Z R100y=$(w100n+mStr+"y"), R100x=$(w100n+mStr+"x"), R100u=$(w100n+mStr+"u")
	WAVE/Z R190y=$(w190n+mStr+"y"), R190x=$(w190n+mStr+"x"), R190u=$(w190n+mStr+"u")
	If( Range(R100x,qpt) && Range(R190x,qpt) )
		ratio= interp(qpt, R100x, R100y)/interp(qpt,R190x,R190y)
		Uratio= ratio* sqrt( (interp(qpt,R100x,R100u)/interp(qpt,R100x,R100y))^2 + (interp(qpt,R190x,R190u)/interp(qpt,R190x,R190y))^2 )
	else
		printf "Ani Extension Error: %g out of range (R100[%1.3f,%1.3f], R190[%1.3f,%1.3f])\r", qpt, R100x[0],R100x[numpnts(R100x)],R190x[0],R190x[numpnts(R190x)]

	endif
	SetDataFolder $fldr
	//recalculates the entire Anisotropy Q-spectrum based on just these two sectors, will patch into the original anisotropy data at Qpt later
	If( mode==0 || H100x[numpnts(H100x)]>V100x[numpnts(V100x)] ) //Horizontal sectors go to highest Q (Or only have H-sectors)
		Duplicate/d/o H100x $(aniN+"Hy"), $(aniN+"Hx"), $(aniN+"Hu"), denom
		WAVE aniSy=$(aniN+"Hy"), aniSx=$(aniN+"Hx"), aniSu=$(aniN+"Hu"), denom
		denom= ( H100y + interp(aniSx, H190x, H190y)*ratio ) //denominator
		aniSy= ( interp(aniSx, H190x, H190y)*ratio - H100y ) / denom * 100
//		aniSu= 2*ratio*H100y*interp(aniSx,H190x,H190y) * sqrt( (Uratio/ratio)^2 + (H100u/H100y)^2 + (interp(aniSx,H190x,H190u)/interp(aniSx,H190x,H190y))^2 ) / denom^2 * 100
		aniSu= 2*ratio*H100y*interp(aniSx,H190x,H190y) * ( (Uratio/ratio) + (H100u/H100y) + (interp(aniSx,H190x,H190u)/interp(aniSx,H190x,H190y)) ) / denom^2 * 100
	else
		Duplicate/d/o V100x $(aniN+"Vy"), $(aniN+"Vx"), $(aniN+"Vu"), denom //Vertical sectors go to highest Q
		WAVE aniSy=$(aniN+"Vy"), aniSx=$(aniN+"Vx"), aniSu=$(aniN+"Vu"), denom
		denom=  ( V100y + interp(aniSx, V190x, V190y)*ratio )  //denominator
		aniSy= ( interp(aniSx, V190x, V190y)*ratio - V100y ) / denom * 100
		aniSu= 2*ratio*V100y*interp(aniSx,V190x,V190y) * sqrt( (Uratio/ratio)^2 + (V100u/V100y)^2 + (interp(aniSx,V190x,V190u)/interp(aniSx,V190x,V190y))^2 ) / denom	^2 * 100
	endif

	// Stitch the new data to the averaged anisotropy data
	WAVE/Z aniY=$(aniN+"y"), aniX=$(aniN+"x"), aniU=$(aniN+"u")
	Variable pt1= BinarySearch(anix,qpt), pt2= BinarySearch(aniSx,qpt)
	IF( pt1<0 )
		Printf "ExtAni Stitch Warning: Q=%g out of range for '%s'[%g,%g]\r", qpt, NameOfWave(aniY),anix[0],anix[inf]
		pt1=numpnts(x1) //just use all the points
	elseif( pt2<0 )
		Printf "ExtAni Stitch Warning: Q=%g out of range for '%s'[%g,%g]\r", qpt, NameOfWave(aniSy),aniSx[0],aniSx[inf]
		pt2=0 //just start at the begining
	endif
	Variable nPts= pt1 + numpnts(aniSx) - pt2
	Redimension/n=(nPts) aniY, aniX, aniU
	//Save parameters in the wave note
	sprintf noteStr, "AniExtQ:%g;AniExtRatio:%1.3f;",qpt,ratio
	NOTE/NOCR aniY, noteStr
	NOTE/NOCR aniX, noteStr
	NOTE/NOCR aniU, noteStr
	//Fill from the new anisotropy wave
//	print aniX[pt1], aniY[pt1], aniSx[pt2], aniSy[pt2]
	aniX[pt1,numpnts(aniX)]=aniSx[p-pt1+pt2]
	aniY[pt1,numpnts(aniX)]=aniSy[p-pt1+pt2]
	aniU[pt1,numpnts(aniX)]=aniSu[p-pt1+pt2]
//	print Uratio/ratio, h100u[pt2]/H100y[pt2], interp(aniSx[pt2],H190x,H190u)/interp(aniSx[pt2],H190x,H190y)
	SetDataFolder $CurrentFolder
end


Function Add_2piOverQticks2top()
	WAVE/Z tickWave=root:Packages:userTicks:logticks_2piOverQ
	If( !WaveExists(tickWave) )
		logTicks_2piOverQ()
	endif
	NewFreeAxis/O/T MT_bottom
	ModifyFreeAxis/Z MT_bottom,master= bottom,hook= TransformMirrorAxisHook
	ModifyGraph userticks(MT_bottom)={root:Packages:userTicks:logticks_2piOverQ_vals,root:Packages:userTicks:logticks_2piOverQ_Labels}
	ModifyGraph lblPosMode(MT_bottom)=1//,lblMargin(MT_bottom)=6,lblLatPos(MT_bottom)=-35
	ModifyGraph freePos(MT_bottom)=0, tick(MT_bottom)=2, log(MT_bottom)=1
	Label MT_bottom "2\\F'Symbol'p\\F'Arial'/Q [nm]"
end

Function logTicks_2piOverQ()
	String CurrentFolder=GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S root:Packages:userTicks
	If( WaveExists(logticks_2piOverQ_Vals) )
		SetDataFolder $CurrentFolder
		return 0
	endif
	Make/n=0/o/d logticks_2piOverQ_Vals
	Make/n=(0,2)/o/t logticks_2piOverQ_Labels
	WAVE vals=logticks_2piOverQ_Vals
	WAVE/T Labels=logticks_2piOverQ_Labels
	variable i, j
	For( i=1, j=0; j<4; i*=10, j+=1 )
		insertpoints inf, 9, vals, Labels
		vals[j*9, j*10+8]=2*pi/( (p+1-j*9)*i )
		Labels[j*9][0]=num2str(i)
		Labels[j*9][1]="Major"
		Labels[j*9+1,j*9+8][1]="Minor"
	Endfor
	SetDimLabel 1,1,'Tick Type',Labels
	SetDataFolder $CurrentFolder
End

Function MPFstyle() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	String tList=TraceNameList("",";",1)
	String ftName=ListMatch(tList,"fit_*")
	If( !strLen(ftname) )
		ftName=ListMatch(tList,"'fit_*") //add single quote
	endif
	String dName=stringfromlist(1,ftName,"_")
	ModifyGraph/Z margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
	ModifyGraph/Z rgb[1]=(65280,0,0),rgb[2]=(1,4,52428),rgb[3]=(65280,0,0),rgb[4]=(2,39321,1)
	ModifyGraph/Z tick=2
	ModifyGraph/Z zero(Res_left)=1
	ModifyGraph/Z mirror=1
	ModifyGraph/Z nticks(left)=4,nticks(Peaks_Left)=3,nticks(Res_left)=2
	ModifyGraph/Z minor=1
	ModifyGraph/Z standoff=0
	ModifyGraph/Z notation(left)=1,notation(Peaks_Left)=1,notation(Res_left)=1
	ModifyGraph/Z lblPosMode(Peaks_Left)=1,lblPosMode(Res_left)=1
	ModifyGraph/Z lblPos(left)=48
	ModifyGraph/Z freePos(Peaks_Left)={0,kwFraction}
	ModifyGraph/Z freePos(Res_left)={0,kwFraction}
	ModifyGraph/Z axisEnab(left)={0.25,0.75}
	ModifyGraph/Z axisEnab(Peaks_Left)={0,0.2}
	ModifyGraph/Z axisEnab(Res_left)={0.8,1}
	Label/Z left "Intensity \\u"
	Label/Z bottom "Q [nm\\S-1\\M]"
	Label/Z Peaks_Left "Peaks \\u"
	Label/Z Res_left "Res \\u"
	TextBox/C/N=text0/B=2/A=LT dName
EndMacro


//Extracts the average scattering and anisotropic scattering ratio at a given Q value from a processed energy scan
Function ProcEscanOld(name, snum, qVal, fact,[tagE])
	String name //name of scan
	Variable snum //scan number to process
	Variable qVal //which q-value to process at
	Variable fact //factor to multiply by the film absorbance to subtract off the XRF backgrounds
	Variable TagE //extract the profile at the tagged energy

//	WAVE RawAbs=root:NEXAFS:PSPMMAfresh:PSPMMAsub, RawAbsEn=root:NEXAFS:PSPMMAfresh:energy //original subtractions
//	WAVE RawAbs=root:NEXAFS:PSPMMAfresh:PSPMMAfresh_absSub, RawAbsEn=root:NEXAFS:PSPMMAfresh:energy //actual BL11 NEXAFS
	WAVE RawAbs=root:Spectra:Scan034_fit0, RawAbsEn=root:Spectra:Energy0 //my model by hand (frozen at 28 wt.% PMMA etc.)
	WAVE RawSiNTrans=root:NEXAFS:PSPMMAfresh:SiNtrans, SiNen=root:NEXAFS:PSPMMAfresh:Energy

	String CurrentFolder=GetDataFolder(1), dFldr="root:Escans:"+name, SectFldr="root:Escans:"+name+":Meta", cFldr="root:Escans:"+name+":MetaC"
	NewDataFolder/o/s root:Escans
	NewDataFolder/o/s $dFldr
	NewDataFolder/o/s $SectFldr
	AssembleScanInfo(num2str(snum),"1","root:Escans:"+name+":Meta",10)  //assembles info for the sectors
	WAVE AnglesSect=Angles0, DwellsSect=Dwells0, EnergiesSect=Energies0, PolsSect=Pols0, NotesSect=Notes0
	WAVE/T ScansSect=Scans0
	NewDataFolder/o/s $cFldr
	AssembleScanInfo(num2str(snum),"1","root:Escans:"+name+":MetaC",360)  //assembles info for the average profiles
	WAVE DwellsC=dwells0, EnergiesC=energies0, PolsC=pols0, NotesC=notes0
	WAVE/T ScansC=scans0
	SetDataFolder $dFldr
	WAVE Energy=MakeEnergies(PolsC,EnergiesC) //Set up energies
	Duplicate/d/o Energy, Absorption, SiNtrans
	Absorption=interp(Energy,RawAbsEn, RawAbs)
	SiNtrans=interp(Energy,SiNen,RawSiNTrans)
	Duplicate/d/o Energy, ew2Di, IntH, IntHu, IntV, IntVu, IntAvg, IntAvgU, TSIH, TSIHu, TSIV, TSIVu, TSIdiff TSIavg, TSIavgU //output Intensity waves
	Duplicate/d/o Energy ,AniH, AniHu, AniV, AniVu, AniAvg, AniAvgU  //output Anisotropy waves
	Variable i, idx, uDenom, init=0 //2D data matrix initialized 
	Variable/C TSIcmplx
	For(i=0; i<numpnts(energy); i+=1 )
		// Average Intensities H-pol
		idx=QueryMeta(360,energy[i],100,cFldr)
		If( idx<0 )
			print "Error! HC Query Failed for E="+num2str(energy[i])
			break
		endif
		WAVE qwH=$("root:SAS:'"+ScansC[idx]+"':'q_"+ScansC[idx]+"'")
		WAVE dwH=$("root:SAS:'"+ScansC[idx]+"':'r_"+ScansC[idx]+"'")
		WAVE uwH=$("root:SAS:'"+ScansC[idx]+"':'s_"+ScansC[idx]+"'")
		If( !init )
			Duplicate/d/o qwH, qw2D, IntProf, uProf, qw2Di
			If( WaveExists(Hint2D) )
				Redimension/n=(0,0) Hint2D, Hint2Du, Vint2D, Vint2D
			Endif
			ImageScale(qw2Di) //make X & Y waves to plto the 2D data against in an image
			ImageScale(ew2Di)
			init=1
		endif
		Duplicate/d/o dwH, subProf, subProfU
		subProf=dwH[p]*energy[i]/siNTrans[i]
		subProfU=uwH[p]*energy[i]/SiNTrans[i]
		If( energy[i]>=283 ) //don't subtract any XRF below the edge (should be none)
			subProf-=fact*Absorption[i] //subtract off the XRF background from the data which should be proportional to Abs
		endif
		//making 2D wave
		IntProf=interp(qw2D,qwH,subProf)
		uProf=interp(qw2D,qwH,subProfU)
		Concatenate {IntProf}, Hint2D
		Concatenate {uProf}, Hint2Du
		//just energy dependence
		IntH[i]=interp(qVal,qwH,subProf)
		IntHu[i]=interp(qVal,qwH,subProfU)
		If( abs(energy[i]-tagE)<0.1 )
			Duplicate/d/o subProf TagIntH
			Duplicate/d/o qwH, TagQH
			TagQH*=10 //from A^-1 into nm^-1
		endif
		//Invariant
		TSIcmplx=IntData(subProf,qwH,subProfU,Xpwr=2)*100/(2*pi^2) //convert from 1/A^2 to 1/nm^2 & multiply the Jacobian
		TSIH[i]=real(TSIcmplx)
		TSIHu[i]=imag(TSIcmplx)
		
		// Average Intensities V-pol
		idx=QueryMeta(360,energy[i],190,cFldr)
		If( idx<0 )
			print "Error! VC Query Failed for E="+num2str(energy[i])
			break
		endif
		WAVE qwV=$("root:SAS:'"+ScansC[idx]+"':'q_"+ScansC[idx]+"'")
		WAVE dwV=$("root:SAS:'"+ScansC[idx]+"':'r_"+ScansC[idx]+"'")
		WAVE uwV=$("root:SAS:'"+ScansC[idx]+"':'s_"+ScansC[idx]+"'")
		Duplicate/d/o dwV, subProf, subProfU
		subProf=dwV[p]*energy[i]/SiNTrans[i]
		subProfU=uwV[p]*energy[i]/SiNTrans[i]
		If( energy[i]>=283 ) //don't subtract any XRF below the edge (should be none)
			subProf-=fact*Absorption[i]
		endif
		//making 2D wave
		IntProf=interp(qw2D,qwV,subProf)
		uProf=interp(qw2D,qwV,subProfU)
		Concatenate {IntProf}, Vint2D
		Concatenate {uProf}, Vint2Du
		//Just energy dependence
		IntV[i]=interp(qVal,qwV,subProf)
		IntVu[i]=interp(qVal,qwV,subProfU)
		//Invariant
		TSIcmplx=IntData(subProf,qwV,subProfU,Xpwr=2)*100/(2*pi^2)
		TSIV[i]=real(TSIcmplx)
		TSIVu[i]=imag(TSIcmplx)
		If( abs(energy[i]-tagE)<0.1 )
			Duplicate/d/o dwV TagIntV
			Duplicate/d/o qwV TagQV
			TagQV*=10 //from A^-1 into nm^-1
		endif
		
		// Anisotropy H-pol
		// Horizontal Scattering (p-polarized scattering)
		idx=queryMeta(0,energy[i],100,sectFldr) //horizontal to the right
		WAVE q0=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d0=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u0=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		idx=queryMeta(180,energy[i],100,sectFldr) //horizontal to the left
		WAVE q1=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d1=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u1=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		CombineProfiles(d0, q0, u0, d1, q1, u1) //combine the two directions
		WAVE AvY, AvX, AvU
		Duplicate/d/o AvY, dHH
		Duplicate/d/o AvX, qHH
		Duplicate/d/o AvU, uHH
		// Vertical Scattering (s-polarized scattering)
		idx=queryMeta(90,energy[i],100,sectFldr) //vertical up
		WAVE q0=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d0=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u0=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		idx=queryMeta(270,energy[i],100,sectFldr) //vertical down
		WAVE q1=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d1=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u1=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		CombineProfiles(d0, q0, u0, d1, q1, u1) //combine the two directions
		WAVE AvY, AvX, AvU
		Duplicate/d/o AvY, dVH
		Duplicate/d/o AvX, qVH
		Duplicate/d/o AvU, uVH
		
		AniH[i]=(interp(qVal,qHH,dHH)-interp(qVal,qVH,dVH)) / (interp(qVal,qHH,dHH)+interp(qVal,qVH,dVH))*100
		uDenom= (interp(qVal,qHH,dHH) + interp(qVal,qVH,dVH) )^2
		AniHu[i]= sqrt( (2*interp(qVal,qVH,dVH)*interp(qVal,qHH,uHH))^2 + (2*interp(qVal,qHH,dHH)*interp(qVal,qVH,uVH))^2 )*100/uDenom
		
		// Anisotropy V-pol
		// Horizontal Scattering (s-polarized scattering)
		idx=queryMeta(0,energy[i],190,sectFldr) //horizontal to the right
		WAVE q0=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d0=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u0=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		idx=queryMeta(180,energy[i],190,sectFldr) //horizontal to the left
		WAVE q1=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d1=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u1=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		CombineProfiles(d0, q0, u0, d1, q1, u1) //combine the two directions
		WAVE AvY, AvX, AvU
		Duplicate/d/o AvY, dHV
		Duplicate/d/o AvX, qHV
		Duplicate/d/o AvU, uHV
		// Vertical Scattering (p-polarized scattering)
		idx=queryMeta(90,energy[i],190,sectFldr) //vertical up
		WAVE q0=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d0=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u0=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		idx=queryMeta(270,energy[i],190,sectFldr) //vertical down
		WAVE q1=$("root:SAS:'"+ScansSect[idx]+"':'q_"+ScansSect[idx]+"'")
		WAVE d1=$("root:SAS:'"+ScansSect[idx]+"':'r_"+ScansSect[idx]+"'")
		WAVE u1=$("root:SAS:'"+ScansSect[idx]+"':'s_"+ScansSect[idx]+"'")
		CombineProfiles(d0, q0, u0, d1, q1, u1) //combine the two directions
		WAVE AvY, AvX, AvU
		Duplicate/d/o AvY, dVV
		Duplicate/d/o AvX, qVV
		Duplicate/d/o AvU, uVV
		
		AniV[i]=(interp(qVal,qVV,dVV)-interp(qVal,qHV,dHV)) / (interp(qVal,qVV,dVV)+interp(qVal,qHV,dHV))*100
		uDenom= (interp(qVal,qVV,dVV) + interp(qVal,qHV,dHV) )^2
		AniVu[i]= sqrt( (2*interp(qVal,qHV,dHV)*interp(qVal,qVV,uVV))^2 + (2*interp(qVal,qVV,dVV)*interp(qVal,qHV,uHV))^2 )*100/uDenom
		
		AniAvg[i]= (AniH[i]+AniV[i])/2
		AniAvgU[i]= sqrt( (AniHu[i]/AniH[i])^2 + (AniVu[i]/AniV[i])^2 )
	endfor

	IntAvg=(IntH+IntV)/2
	IntAvgU=IntHu+IntVu
	TSIavg=(TSIH+TSIV)/2
	TSIdiff=abs(TSIH-TSIV)
	TSIavgU= max( TSIavg * sqrt( (TSIHu/TSIH)^2 + (TSIVu/TSIV)^2 ), TSIdiff/2 )
	qw2D*=10; qw2Di*=10
	Duplicate/d/o Hint2D AvgInt2D, AvgInt2Du, DiffInt2D, RatioInt2D, AvgLor2D
	Duplicate/d/o IntProf, DiffProf
	Duplicate/d/o energy HiQdiff
	Duplicate/d/o TagIntH TagIntAvg
	Duplicate/d/o TagQH TagQ
	TagIntAvg=(TagIntH+interp(TagQH, TagQV,TagIntV))/2
	Duplicate/d/o TagIntAvg TagIntAvg2
	TagIntAvg2*=TagQ^2
	Variable nQpts=numpnts(IntProf), pts2avg=20
	RatioInt2D=(Hint2D/Vint2D-1)*100
	AvgInt2D=(Hint2D+Vint2D)/2// - fact*Absorption[q]
	AvgLor2D=AvgInt2D[p][q]*qw2D[p]^2
	AvgInt2Du=Hint2Du+Vint2Du
	For( i=0; i<Dimsize(Avgint2D,1); i+=1 )
		IntProf=AvgLor2D[p][i]
//		IntProf=(Avgint2D[p][i] - Fact*Absorption[i])*qw2D[p]^2
		Differentiate IntProf /X=qw2D /D=DiffProf
		DiffInt2D[][i]=DiffProf[p]
		Wavestats/Q/R=[nQpts-pts2avg-5,nQpts-5] DiffProf
		HiQdiff[i]=V_avg
//		Duplicate/d/o HiQdiff HiQdiff_raw
	endfor
End	

//Quick script to remove the units at the end of each wavenote key. This is new in Igor 7 due to new way to import Data

Function/T RemoveUnits(s)

	String s
	String EndString
	
	Variable spaceloc = strsearch(s," ",0)
	
	EndString = s[0,spaceloc-1]
	return EndString
	
end