#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function Calc_BCP_SF(Iw,Qw,Rg,f1,vol,chi)

	Wave IW //Intensity wave of the Structure factor S(q)
	Wave Qw // The q-wave that you want to calculate the scattering for
	
	Variable Rg //The full polymer radius of Gyration
	Variable f1 //The volume fraction of block '1'
	Variable vol //Average volume of the polymer in [nm^3] 
	
	Variable chi //For now, this is equal to chi/vu 
	
//CHeck the volume fraction to make sure they are physical
	if(f1 > 100 || f1 < 0)
		print "Selected volume fraction must be between 0 and 100"
		return 1
	endif
	f1 = f1 > 1 ? f1/100 : f1
	Variable f2 = 1-f1
	
	//Make a form-factor wave
	duplicate/o Iw, Iw_FF ; Wave Iw_FF
	
	Calc_BCP_FormFactor(Iw_FF,Qw,Rg,f1)
	
	Iw = Iw_ff - 2*chi*vol
	Iw = 1/Iw	
	
end


//Calc the scattering from a disordered BCP
//Does not include contrast or volume information
Function Calc_BCP_FormFactor(Iw,Qw,Rg,f1)

	Wave IW //Intensity wave of the form factor F(q)
	Wave Qw // The q-wave that you want to calculate the scattering for
	
	Variable Rg //The full polymer radius of Gyration
	Variable f1 //The volume fraction of block '1'
	
//CHeck the volume fraction to make sure they are physical
	if(f1 > 100 || f1 < 0)
		print "Selected volume fraction must be between 0 and 100"
		return 1
	endif
	f1 = f1 > 1 ? f1/100 : f1
	Variable f2 = 1-f1

	//Calculate the Radius of Gyration of the individual components
	//Roe Page 225 (Equation 6.51)
	
	Variable Rg1 = Sqrt(f1*Rg^2)
	Variable Rg2 = Sqrt(f2*Rg^2)
		
	//Calculate the form factor from page 225
	//Make waves for Debye Functions
	Duplicate/o IW, Dx, Dx1, Dx2 //Debye with Rg, Rg1, Rg2
	Wave Dx ; wave Dx1 ; wave Dx2 //Grab the waves for later use
	
	//Fill in the 3 Debye functions
	Calc_Debye_Scattering(Dx,Qw,Rg)
	Calc_Debye_Scattering(Dx1,Qw,Rg1)
	Calc_Debye_Scattering(Dx2,Qw,Rg2)
	
	Iw = Dx
	Iw /= ( (f1^2*Dx1*f2^2*Dx2) - (1/4)*( Dx - (f1^2*Dx1) - (f2^2*Dx2) )^2 )
	
	
end


//Calc the scattering of a Debye Function...No folder structure within function.
//DOES NOT INCLUDE CONTRAST OR VOLUME INFORMATION....Strictly the integral of the polymer chain
//This will be a fully contained function
//Scattering and Q waves must be created before calculation
Function Calc_Debye_Scattering(Iw,qw,Rg)

	Wave Iw // The Wave you want to calculate the intensity
	Wave Qw // The q-wave that you want to calculate the scattering for
	Variable Rg //The polymer Radius of gyration
	
	//Function from Page 164 of Roe Scattering
	//D(x) = 2(e^-x + x - 1) / x^2 -> x = q^2*Rg^2
	
	if(Dimsize(Iw,0) != Dimsize(Qw,0))
		Print "Intensity wave does not match q-wave, please fix and try again"
		return 0	
	endif
	
	Duplicate/O Qw xw ; Wave xw // Make and Assign the xw
	xw = Qw^2 * Rg^2 // Calc x
	
	IW = 2*(Exp(-xw) + xw - 1)
	IW /= xw^2

end


Function Make_Qw_Logspace(npnts,start,stop)

	Variable npnts //number of points you want
	Variable start //start point
	Variable stop  //end point
		
	Variable Div = (log(stop) - log(start))/(npnts-1) //Exponential differentiation

	Make/N=(npnts)/O Qw ; Wave Qw //Make the wave with N=npnts
	
	Qw = start * 10^(div*p) //Calculate the log-spaced values
	
	
end