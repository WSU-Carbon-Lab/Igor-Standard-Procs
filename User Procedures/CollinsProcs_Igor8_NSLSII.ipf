#pragma rtGlobals=1		// Use modern global access method.



///This is a compilable version of CollinsProc.ipf into Igor 7...NOT ALL FUNCTIONS HAVE BEEN TESTED.

//Takes a number and it's uncertainty and returns a string in the shorthand format "3.54(3)e+6"
Function/S FormatNumber(Val,Uncert)
	Variable Val, Uncert
	If( Val*0!=0 || Uncert*0!=0 )
		return "nan(nan)"
	elseIf( Uncert==0 )
		return num2str(Val)+"(0)"
	endif
	Variable sigDigits=ceil(log(abs(Val)))-floor(log(Uncert)) //number digits in value based on size of uncert
	String fstr, str
	If( sigDigits<0 && (uncert>950 || uncert<0.0095 ) ) //report range in which zero is a good value
		sprintf str, "%1.0e", uncert  //stores the exponential info
		return "0("+TruncExp(Uncert,1)+")e"+num2str(str2num(str[2,inf]) )
	elseif( sigDigits<0 && uncert>0.95 )
		return "0("+TruncDigits(Uncert,1)+")"
	elseif( sigDigits<0 )
		return "0("+num2str(Uncert)+")"
	elseif( abs(Val)>950 || abs(Val)<0.0095 ) //report exponential format
		sprintf str, "%1.0e", abs(Val)  //stores the exponential info
		return TruncExp(Val,sigDigits)+"("+TruncExp(Uncert,1)+")e"+num2str(str2num(str[2,inf]) )
	elseif( uncert>0.95 ) //uncertainty rounds to one or more
		val= val<1 ? round(val) : val
		return TruncDigits(Val,sigDigits)+"("+TruncDigits(Uncert,1)+")"
	else
		Variable nDigits= abs(Val)>Uncert ? ceil(log(abs(val))) : floor(log(uncert)) //order of mag of the value
		sprintf fstr, "%%1.%df",sigDigits-nDigits;	sprintf str, fstr, val
		return str+"("+TruncExp(Uncert,1)+")"
	endif
End


//takes # and returns it as a string rounded to the # of digits given assuming # is larger than desired
//digits e.g. 427 -> 430 or 400
Function/S TruncDigits(Val,Digits)
	Variable Val, Digits
	Variable nDigitsOver1=floor(log(val)) //# digits in front of decimal
	If( Digits<=1 ) //deal with rounding up an order of magnitude
		String str1=num2str(Val)
		String str2=TruncExp(Val,Digits)
		If( (stringmatch(str1[0],"9") || stringmatch(str1[2],"9")) && stringmatch(str2[0],"1") )
			nDigitsOver1+=1
		endif
	endif
	Variable num=str2num(TruncExp(Val,Digits))*10^nDigitsOver1
	String str
	sprintf str, "%d", num
	return str
end


//returns as a string only the leading digits of a number in exponential format without the exponent
//rounded to that digit
Function/S TruncExp(Val, digits)
	Variable Val, digits
	String fstr, str
	sprintf fstr, "%%1.%de", digits-1
	sprintf str, fstr, Val
	If( digits==1 ) //deal with rounding up an order of magnitude
		String str1=num2str(Val)
		If( (stringmatch(str1[0],"9") || stringmatch(str1[2],"9")) && stringmatch(str[0],"1") )
			digits+=1
			sprintf fstr, "%%1.%de", digits-1
			sprintf str, fstr, str2num(str)
		endif
	endif
	If( digits==1 )
		return str[0]
	else
		return str[0,digits]
	endif
end


//This function alters the input wave so that it can be used to scale an image x or y axis
Function imageScale(w)
	WAVE w
	Variable i, d
	insertpoints 0, 1, w //insert the exta point at the beginning
	w[0]=w[1]-(w[2]-w[1])/2 //new first pt is 1/2 pixel over
	for( i=1; i<numpnts(w)-1; i+=1 )
		d=(w[i+1]-w[i])/2
		w[i]+=d
		If( i==numpnts(w)-2 )
			w[i+1]+=d
		endif
	endfor
End


//returns the center of a window or the screen
Function/C BC_wPos(wName)
	String wName
	Dowindow $wName
	If( stringMatch(wName, "Screen") )
		GetWindow kwFrameInner wsize
	elseif( V_Flag==1 )
		GetWindow $wName wsize
	endif
	// Convert points into pixels
	V_left *= ScreenResolution/72
	V_right *= ScreenResolution/72
	V_top *= ScreenResolution/72
	V_bottom *= ScreenResolution/72
	return cmplx( (V_right+V_left)/2, (V_top+V_bottom)/2 )
end

//returns the number of unique values in a wave within the tolerance
Function NumUnique(w,tol)
	wave w
	Variable tol
	variable i
	String List="", val
	For( i=0; i<Numpnts(w); i+=1 )
		val= num2str(round(w[i]/tol)*tol)
		If( WhichListItem( val, List) == -1)
			List= AddListItem(val,List)
		endif
	EndFor
	return ItemsInList(List)
end	

//outputs a wave containing only the unique values from an input wave within the given tolerance
Function/WAVE UniqueVals(w,tol,[wn])
	wave w
	Variable tol
	String wn //optional output name
	String wDF=GetWavesDataFolder($NameOfWave(w),1)+SelectString(ParamIsDefault(wn),wn,NameOfWave(w)+"_unique")
	Duplicate/d/o w $(wDF)
	WAVE uEn=$(wDF)
	Variable i, j
	For( i=0; i<Numpnts(w); i+=1 )
		uEn= abs(uEn[p]-uEn[i])<tol ? ( p==i ? uEn : nan ) : uEn //mark all duplicates as a nan
	Endfor
	sort uEn, uEn //get nans to the end of the wave
	wavestats/Q uEn
	deletepoints V_npnts, V_numNans, uEn //delete all nans
	return uEn
end

//Combines multiple waves and outputs only the unique vals within a given tolerance

Function/Wave UniqueVals_Set(ws,tol,wn)
	String ws //String of waves that you want to check
	Variable tol
	string wn //Not an optional name anymore
	
	Make/n=0/O $wn
	Wave uVal = $wn //w
	Variable i,j
	//Variable NumVals =0 //Number of values over the multiple input data waves
	//Quick loop to find number of Total elements combined in all the given waves
	For(j=0;j<itemsinlist(ws,";");j+=1)
		wave tempwave = $stringfromlist(j,ws,";")
		InsertPoints 0,numpnts(tempwave), uVal
		uVal[0,numpnts(Tempwave)-1] = Tempwave
	Endfor
	//New loop to remove the expected duplicates
	For( i=0; i<Numpnts(uVal); i+=1 )
		uVal= abs(uVal[p]-uVal[i])<tol ? ( p==i ? uVal : nan ) : uVal //mark all duplicates as a nan
	Endfor
	sort uVal, uVal //get nans to the end of the wave
	wavestats/Q uVal
	deletepoints V_npnts, V_numNans, uVal //delete all nans
	return uVal
	
End

Function TestAboveFunc(ws)
	String ws 
	
	//Make/o Energies
	wave Energy = UniqueVals_Set(ws,0.05,"Energy")
end

//Alters the display options to make the sector graphs look good for anisotropy (assumes 5 sectors per energy)
Function GroupStyle([phase, n, CO, off,wn,resetC]) : GraphStyle
	Variable phase //start 'angle' or shade: 0=dark, 45=mid-going-light, 90=light shade, 135=mid-going-dark (DEFAULT: 0)
	Variable n //number of traces in group (DEFAULT: 3)
	Variable off  //multiplicative offset of each group (DEFAULT: 1)
	Variable resetC //resets the color series after this many groups (DEFAULT: 9)
	WAVE/T CO //Color order for the groups
	String wn //window name (DEFAULT: TopGraph)
	phase= ParamIsDefault(phase) ? 0 : phase*pi/180 //convert to radians
	n= ParamIsDefault(n) ? 3 : n
	off= ParamIsDefault(off) ? 1 : off
	off= off==0 ? 1 : off
	resetC= ParamIsDefault(resetC) ? 9 : resetC
	wn= SelectString( ParamIsDefault(wn), wn, WinName(0,1))
	String CurrentFolder=GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S root:Packages:GroupStyle
	If( !WaveExists(CO) )
		make/o/t/n=9 CO
		CO={"Black", "Red","Blue","Green","Orange","Purple","Yellow","Magenta","Teal"}
	endif
	Make/o/n=(3,3) Black, Red, Blue, Green, Orange, Purple, Yellow, Magenta, Teal
	Black={{0,0,0},{17408,17408,17408},{30464,30464,30464}}
	Red={{39168,0,0},{65280,0,0},{65280,32768,32768}}
	Blue={{0,0,39168},{0,0,65280},{32768,32768,65280}}
	Green={{0,26112,0},{0,52224,0},{16384,65280,16384}}
	Orange={{39168,26112,0},{65280,43520,0},{65280,54528,32768}}
	Purple={{19712,0,39168},{29440,0,58880},{44032,29440,58880}}
	Yellow={{26112,26112,0},{52224,52224,0},{65280,65280,16384}}
	Magenta={{39168,0,31232},{65280,0,52224},{65280,32768,58880}}
	Teal={{0,26112,26112},{0,52224,52224},{16384,65280,65280}}
	Variable i, j, k, c, nTraces=ItemsInList( TraceNameList(wn, ";", 1) )
	For( i=0, c=0; i<nTraces/n; i+=1, c+=1 ) //cycle through groups
		For( j=0; j<n; j+=1) //cycle through traces in a group
			k=abs(round(2*sin(j*pi/4+phase)))
			c= c==resetC ? 0 : c
			ModifyGraph/W=$wn/Z rgb[i*n+j]=(Color(CO,c,k,0),Color(CO,c,k,1),Color(CO,c,k,2))
		Endfor
	Endfor
		
	//Offsets for the groups.
	If( off!=1 )
		For( i=1; i<9; i+=1 )
			For( j=0; j<n; j+=1 )
				ModifyGraph/W=$wn/Z muloffset[i*n+j]={1,off^i}
			endfor
		endfor
	endif
	SetDataFolder $CurrentFolder
End


//returns an rgb value of a color based on a list of color hues and shades
Function Color(CO,hue, shade rgb)
	wave/T CO //color order
	Variable hue, shade, rgb
	WAVE color=$(CO[hue])
	return color[rgb][shade]
end

//returns 1 if wave with same name is already plotted in the window otherwise returns 0
Function Plotted(w,wn)
	wave w //wave that might be plotted
	string wn //window name
	If( strLen(wn)==0 ) //get name of top graph
		wn=StringFromList(0,WinList("*",";",""))
	endif
	return WhichListItem(NameOfWave(w), ReplaceString("'", TraceNameList(wn,";",1), "" ) )==-1 ? 0 : 1
end

//computes the weighted average of the values in 'yw' weighted by the standard error values in 'uw' and calculates the propagated errors along with it
Function/C wAvg(yw,uw,[p1,p2])
	WAVE yw, uw
	Variable p1,p2 //optional beginning & ending parameters
	Variable npts=numpnts(yw), i, avg=0, denom=0, uncert, i0=0
	If( !ParamIsDefault(p1) & !ParamIsDefault(p2) )
		npts=p2-p1+1
		i0=p1
	Endif
	For( i=i0; i<i0+npts; i+=1 )
		If( yw[i]*0==0 && uw[i]*0==0 ) //skip nans
			avg+= yw[i]*uw[i]^-2
			denom+= uw[i]^-2
		endif
	endfor
	avg= avg/denom
	uncert= sqrt( 1/denom )
	return cmplx( avg, uncert )
end

//Finds a peak based on the derivative within the x-range provided
Function BAC_FindPk(ywIn, xwIn, x1,x2)
	WAVE ywIn, xwIn
	Variable x1, x2
	NewDataFolder/S/O :ProcVars
	Duplicate/o ywIn, yw
	Duplicate/o xwIn, xw
	Sort xw, yw, xw //new waves with sorted x values
	variable p1=BinarySearch(xw,x1)
	variable p2=BinarySearch(xw,x2)
	Duplicate/R=[p1,p2]/o yw, dyw, ddyw //just look at data range
	Duplicate/R=[p1,p2]/o xw, dxw
	setscale/p x, 0, 1, dyw //make internal x-values equal to scaling
	setscale/p x, 0, 1, ddyw
	smooth 3, dyw
	differentiate dyw /X=dxw
	differentiate dyw /X=dxw /D=ddyw
	Wavestats/Q ddyw //position of sharpest peak within a pixel
	variable PosIndex=V_MaxLoc
	FindLevel/Q/P/R=[PosIndex,0], ddyw, 0 //lower range of where peak could be
	variable lo=V_LevelX
	FindLevel/Q/P/R=[PosIndex,9999], ddyw, 0 //uppper range of where peak could be
	variable hi=V_LevelX
	FindLevel/Q/P/R=[lo,hi], dyw, 0 //Actual index of peak (sub pixel res using linear interp)
	SetDataFolder ::
//	Killwaves/Z yw, xw, dxw, dyw, ddyw
	return dxw[V_LevelX]
end

//make a 2D dataset w higher res by a factor A
Function InterpImage(w,A)
	Wave w
	Variable a
	Duplicate/d/o w iw
	Variable nx=dimSize(w,0), ny=dimSize(w,1)
	Variable x0=dimOffset(w,0), y0=dimOffset(w,1)
	Variable dx=dimDelta(w,0), dy=dimDelta(w,1)
	Redimension/n=(a*nx-(a-1), a*ny-(a-1)) iw
	SetScale/P x, x0, dx/a, iw
	SetScale/P y, y0, dy/a, iw
	iw=interp2D(w,x,y)
	String wn=nameofwave(w)
	duplicate/d/o iw $(wn+"_interp")
end

//Makes a color table suitable to highlight the servere corelations in a correlation matrix
Function MakeCorrelationColorTable()
	Make/d/o/n=1001 Saturation
	Saturation=exp(0.01*x)
	Saturation/=wavemax(saturation)
	Saturation*=65535
	Saturation=65535-Saturation
	Make/d/o/n=(2001,3) corrColorTab
	SetScale/i x, -1, 1, corrColorTab
	CorrColorTab[0,1000][0]=65535
	CorrColorTab[0,1000][1]=Saturation[1000-p]
	CorrColorTab[0,1000][2]=Saturation[1000-p]
	
	CorrColorTab[1000,2000][0]=Saturation[p-1000]
	CorrColorTab[1000,2000][1]=Saturation[p-1000]
	CorrColorTab[1000,2000][2]=65535
end	

//Interpolates a complex wave
Function/C cInterp(x, xw, yw)
	Variable x
	WAVE xw
	WAVE/C yw
	duplicate/d/o xw yReal, yImag
	yReal=real(yw)
	yImag=imag(yw)
	return cmplx( interp(x, xw, yReal), interp(x, xw, yImag) )
end

//integrate a profile with uncertainties only in the y-data
//returns a complex # where real part is the integral, and the imaginary part is the uncertainty
Function/C IntData(yw,xw,uw,[Xpwr])
	WAVE yw, xw, uw
	Variable Xpwr //optionally multiply xw into yw with a power of Xpwr (used for integrating scattering profiles)
	Duplicate/o yw, intW, intWu
	intW*=xw^Xpwr
	intWu=( uw*xw^Xpwr )^2 //quadrature assumes stochastic uncertainties (rather than systematic)
	return cmplx(AreaXY(xw,intW), sqrt(AreaXY(xw,intWu)))
end

Function GraphIntegrals([wn])
	String wn
	wn= SelectString( ParamIsDefault(wn), wn, winName(0,1) )
	String tList=TraceNameList(wn,";",1)
	Variable i
	Variable integral
	For( i=0; i<itemsInList(tList); i+=1 )
		WAVE yw= TraceNameToWaveRef(wn,StringFromlist(i,tList))
		WAVE xw= XWaveRefFromTrace(wn,StringFromlist(i,tList))
		IF( WaveExists(xw) )
			integral= Intdata(yw,xw,xw)
		else
			integral=Area(yw)
		endif
		printf "%s:\tArea = %1.2e\r",NameOfWave(yw),integral
	endfor
end

//converts an integer to a letter based on the 26-letter alphabet. After 'z' it goes to double letters 'aa' etc.
Function/S Int2Let(num)
	Variable num
	If( num<1 )
		return ""
	elseif( num>=1 && num<27 )  // num2char function does floor(num)
		return num2char(96+num)
	elseif( num>=27 && num<53 )
		return num2char(96+num-26)+num2char(96+num-26)
	else
		return ""
	endif
end

//clears all the MStimers in IGOR
function clearTimers()
	Variable i, dmy
	For( i=0; i<10; i+=1 )
		dmy=stopMStimer(i)
	endfor
end

//returns most neg and most positive values commonly contained in the wave list
Function/C CommonRange(wList)
	String wList
	Variable neg=-inf, pos=inf, i
	For( i=0; i<itemsinList(wList); i+=1 )
		WAVE w=$StringFromList(i,wList)
		wavestats/Q w
		neg= V_min>neg ? V_min : neg
		pos= V_Max<pos ? V_max : pos
	endfor
	return cmplx(neg,pos)
end

//Used for combining uncertainties for data measured twice with their own uncert.
//Uses weighted averaging uncertainties, but also adds in uncert. from differences between the two measurements
//Doesn't use differences btwn measurements if they are systematcially higher or lower than each other
//NOTE: both waves need the same dimensions
Function/WAVE combineUncert(d1,d2,u1,u2,[f])
	WAVE d1,d2,u1,u2 //data and uncertainty waves
	variable f //optional factor to multiply by d2
	f= ParamIsDefault(f) ? 1 : f
	String CurrentFolder=GetDataFolder(1)
	If( !StringMatch("Proc",GetDataFolder(0)) )
		NewDataFolder/o/s :Proc
	endif
	duplicate/d/o d1, diff, uncert
	diff=d1-f*d2
	If( waveMax(diff)>0 && waveMin(diff)<0 )
		//uncert=max( imag(wAvg( {d1,f*d2},{u1,f*u2} )), abs(diff)/2 ) Removed 11/10/2015 DUe to oversubtaction
		uncert=imag(wAvg( {d1,f*d2},{u1,f*u2} ))
	else
		uncert=imag(wAvg( {d1,f*d2},{u1,f*u2} ))
	endif
	SetDataFolder $CurrentFolder
	return uncert
end

//tiles a square array(w) TILES times
Function TileWave(w,tiles)
	wave w
	Variable tiles
	Variable num=DimSize(w,0)
	Make/o/n=(tiles*num, tiles*num) tiled
	variable i, j
	For( i=0; i<tiles; i+=1)
		For( j=0; j<tiles; j+= 1)
			tiled[i*num,(i+1)*num-1][j*num,(j+1)*num-1]=w[p-i*num][q-j*num]
		endfor
	endfor
end

//returns a 3-point wave containing the red, green, & blue values for a color fractionally located at "fract" on the color table "tableStr"
//used for displaying a series of traces in one plot whose color continuously changes
Function/WAVE ColorFromTable(fract,tableStr)
	Variable fract //relative (fractional) position in the color table to extract
	String tableStr
	WAVE tab=$(tableStr+"_tab")
	If( !WaveExists(tab) )
		ColorTab2Wave $tableStr
		WAVE M_colors
		Duplicate M_colors, $(tableStr+"_tab")
		WAVE tab=$(tableStr+"_tab")
	endif
	Variable npts=dimSize(tab,0)
	Make/o/n=3 rgb=tab[fract*(npts-1)][p]
	return rgb
end

