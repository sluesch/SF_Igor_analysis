#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3			// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include <Reduce Matrix Size>

function ctrans_avg(wave wav, int refit,int dotcondcentering, string kenner_out,[string condfit_prefix, variable minx, variable maxx, int average])
	// wav is the wave containing original CT data
	// refit tells whether to do new fits to each CT line
	// dotcondcentering tells whether to use conductance data to center the CT data
	// kenner_out is the prefix to replace dat for this analysis
	// kenner_out and condfit_prefix can not contain a number otherwise getfirstnu will not work
	print dimsize(wav,0)
		print dimsize(wav,1)

	variable refnum, ms
	//	option to limit fit to indexes [minx,maxx]

	if (paramisdefault(minx))
		minx=pnt2x(wav,0)
	endif

	if (paramisdefault(maxx))
		maxx=pnt2x(wav,dimsize(wav,0))

	endif

	if (paramIsDefault(average))
		average=1
	endif

	//		stopalltimers() // run this line if timer returns nan

	refnum=startmstimer

	display
	string datasetname=nameofWave(wav) // typically datXXXcscurrent or similar
	string kenner=getsuffix(datasetname) //  cscurrent in the above case
	int wavenum=getfirstnum(datasetname) // XXX in the above case

	// these are the new wave names to be made
	string avg = kenner_out + num2str(wavenum) + "cleaned_avg"
	string centered=kenner_out+num2str(wavenum)+"centered"
	string cleaned=kenner_out+num2str(wavenum)+"cleaned"
	string fit_params_name = kenner_out+num2str(wavenum)+"fit_params"
	wave fit_params = $fit_params_name


	variable N=10; // how many sdevs are acceptable?
//	if (dimsize(wav,0)>1e4)
////notch_filters(wav, Hzs="60;180;300",  Qs="20;150;250"); 
//
//string nfname=datasetname+"_nf"; 
//wave temp_wave
////resampleWave($nfname,300 );
//duplicate/o temp_wave $datasetname
//endif

closeallGraphs(); display




	wave W_coef
	wave badthetasx
	wave badgammasx


	string quickavg=avg_wav($datasetname) // averages datasetname and returns the name of the averaged wave

	if (refit==1)
		if (average==1) // sometimes we do not want to average
			get_initial_params($quickavg);// print W_coef
			fit_transition($quickavg,minx,maxx);// print W_coef
		endif

		get_fit_params($datasetname,fit_params_name,minx,maxx) //
	endif

	if (dotcondcentering==0)
	find_plot_thetas(wavenum,N,fit_params_name)
	doupdate

		if (average==1) // sometimes we do not want to average
			duplicate/o/r=[][3] $fit_params_name mids
			centering($datasetname,centered,mids) // centred plot and average plot
			cleaning($centered,badthetasx)
		endif

	elseif(dotcondcentering==1)
		string condfit_params_name=condfit_prefix+num2str(wavenum)+"fit_params"
		print condfit_params_name
		wave condfit_params = $condfit_params_name
		find_plot_gammas(condfit_params_name,N)
		duplicate/o/r=[][2] condfit_params mids

		centering($datasetname,centered,mids)
		cleaning($centered,badgammasx)
		plot_badgammas($centered)

	endif

	if (average==1) // sometimes we do not want to average
		avg_wav($cleaned) // quick average plot
		get_initial_params($quickavg); //print W_coef
		fit_transition($avg,minx,maxx)
		prepfigs(wavenum,N,kenner,kenner_out,minx,maxx)
	endif
	ms=stopmstimer(refnum)
	print ms/1e6
end





//what does this mean in Igor pro: [p][q] > flag ? p : NaN
//In Igor Pro, the expression "[p][q] > flag ? p : NaN" is a conditional statement that checks if the value of the two-dimensional array element located at [p][q] is greater than the value of the variable "flag".
//If the condition is true, the statement returns the value of "p". If the condition is false, the statement returns "NaN", which stands for "Not a Number" and is used to represent undefined or unrepresentable numerical values.

function /wave get_initial_params(sweep)

	// for a given sweep returns a guess of initial parameters for the fit function: Charge transiton

	wave sweep
	duplicate /o sweep x_array
	x_array = x

	variable amp = wavemax(sweep) - wavemin(sweep) //might be worthwile looking for a maximum/minimum with differentiation
	//variable amp = 0.001
	variable const = mean(sweep)
	variable theta = 50

	duplicate /o sweep sweepsmooth
	Smooth/S=4 201, sweepsmooth ;DelayUpdate

	differentiate sweepsmooth
	extract/INDX sweepsmooth, extractedwave, sweepsmooth == wavemin(sweepsmooth)
	variable mid = x_array[extractedwave[0]]

	//extract/INDX sweepsmooth, extractedwave, sweepsmooth == 0 //new
	//variable amp = sweep[extractedwave[0]] - sweep[extractedwave[1]] // new


	variable lin = 0.001  // differentiated value of flat area?

	Make /D/N=6/O W_coef
	W_coef[0] = {amp,const,theta,mid,lin,0,0}; print W_coef

	killwaves extractedwave, sweepsmooth
	return W_coef

end




function /wave fit_transition(current_array,minx,maxx)
	// fits the current_array, If condition is 0 it will get initial params, If 1:
	// define a variable named W_coef_guess = {} with the correct number of arguments

	wave current_array
	variable minx,maxx


	wave W_coef
	wave/t T_Constraints 

	 // W_coef[0]= {-0.03095,1.29277,7.82084,-334.495,4.46489e-05,0}
	   W_coef[0]= {-0.00422886,0.216177,225.6,32.8403,-2.42636e-07,0}
	   T_Constraints[0]= {"K0 > -0.1","K0 < 0","K1 > 0","K1 < 0.5","K2 > 0","K3 > -500","K3 < 500"}



	removefromgraph /all; appendtoGraph current_array
	//FuncFit/q /TBOX=768 CT_faster W_coef current_array(minx,maxx) /D
	FuncFit/q/H="000000"/TBOX=768 Chargetransition W_coef current_array(minx,maxx) /D /C=T_Constraints 
	makecolorful(); 
	doupdate
	//sleep/s 1
end




function /wave get_fit_params(wave wavenm, string fit_params_name,variable minx, variable maxx)
	// returns wave with the name wave "dat"+ wavenum +"fit_params" eg. dat3320fit_params

	//If condition is 0 it will get initial params, If 1:
	// define a variable named W_coef_guess = {} with the correct number of arguments


	variable i
	string w2d=nameofwave(wavenm)
	int wavenum=getfirstnum(w2d)
	int nc
	int nr
	wave temp_wave,slice
	wave W_coef
	wave W_sigma

	duplicate/o wavenm, temp_wave


	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points); 
	//nc=10
	
	make /N= (nc , 12) /o $fit_params_name
	wave fit_params = $fit_params_name
//reduce_and_chop(wavenm, 5000)
print dimsize(wavenm,0); print dimsize(temp_wave,0) 
//resampleWave(wavenm,300 );
print dimsize(wavenm,0); print dimsize(temp_wave,0) 
	for (i=0; i < nc ; i+=1)

rowslice(temp_wave,i)
		fit_transition(slice,minx,maxx)
		fit_params[1 * i][,5] = W_coef[q]
		fit_params[1 * i][6,] = W_sigma[q-6]         //I genuinely cant believe this worked
		// i dont think the q-5 does anything, should double check
	endfor

	return fit_params

end


function find_plot_thetas(int wavenum,variable N,string fit_params_name)

	//If condition is 0 it will get initial params, If 1:
	// define a variable named W_coef_guess = {} with the correct number of arguments

//	string fit_params_name =kenner_out+num2str(wavenum)+"fit_params"
	variable thetamean
	variable thetastd
	variable i
	int nr
	//variable N //how many sdevs?



	wave fit_params = $fit_params_name
	nr = dimsize(fit_params,0)

	duplicate /O/R =[0,nr][2] fit_params thetas


	thetamean = mean(thetas)
	thetastd = sqrt(variance(thetas))

	make /o/n =(nr) meanwave
	make /o/n =(nr) stdwave
	make /o/n =(nr) stdwave2
	make /o/n = 0 goodthetas
	make /o/n = 0 goodthetasx
	make /o/n = 0 badthetas
	make /o/n = 0 badthetasx


	meanwave = thetamean
	stdwave = thetamean - N * thetastd
	stdwave2 = thetamean + N * thetastd


	//display thetas, meanwave, stdwave, stdwave2


	for (i=0; i < nr ; i+=1)

		if (abs(thetas[i] - thetamean) < (N * thetastd))

			insertPoints /v = (thetas[i]) nr, 1, goodthetas // value of theta
			insertpoints /v = (i) nr, 1, goodthetasx        // the repeat

		else

			insertPoints /v = (thetas[i]) nr, 1, badthetas // value of theta
			insertpoints /v = (i) nr, 1, badthetasx        // repeat

		endif

	endfor

execute("theta_graph()")

end


function plot_badthetas(wave wavenm)

	int i
	int nr
	wave badthetasx
	string w2d=nameofwave(wavenm)

	duplicate /o wavenm, wavenmcopy
	nr = dimsize(badthetasx,0)

	display
if (nr>0)
	for(i=0; i < nr; i +=1)
		appendtograph wavenmcopy[][badthetasx[i]]

	endfor

	QuickColorSpectrum2()

	ModifyGraph fSize=24
	ModifyGraph gFont="Gill Sans Light"
	//    ModifyGraph width={Aspect,1.62},height=300
	Label bottom "voltage"
	Label left "current"
	TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 bad thetas of " +w2d
endif
end



function /wave cleaning(wave center, wave badthetasx)
	string w2d=nameofwave(center)
	int wavenum=getfirstnum(w2d)
	string cleaned=getprefix(w2d)+num2str(wavenum)+"cleaned"
	duplicate/o center $cleaned


	// removing lines with bad thetas;

	variable i, idx
	int nc
	int nr
	nr = dimsize(badthetasx,0) //number of rows
	i=0
	if (nr>0)
		do
			idx=badthetasx[i]-i //when deleting, I need the -i because if deleting in the loop the indeces of center change continously as points are deleted
			DeletePoints/M=1 idx,1, $cleaned
			DeletePoints/M=1 idx,1, center

			i=i+1
		while (i<nr)
	endif
end

function find_thetas_smaller(variable  val)
wave thetas
Duplicate/o thetas,reduced;
reduced = thetas[p] > val ? p : NaN
redimension/n=-1 reduced
WaveTransform zapnans reduced
end

function find_thetas_larger(variable  val)
wave thetas
Duplicate/o thetas,reduced;
reduced = thetas[p] > val ? p : NaN
redimension/n=-1 reduced
WaveTransform zapnans reduced
end

function prepfigs(wavenum,N,kenner, kenner_out, minx, maxx)
	variable wavenum,N
	string kenner, kenner_out
	variable minx, maxx
	string datasetname ="dat"+num2str(wavenum)+kenner // this was the original dataset name
	string avg = kenner_out + num2str(wavenum) + "cleaned_avg" // this is the averaged wave produced by avg_wave($cleaned)
	string centered=kenner_out+num2str(wavenum)+"centered" // this is the centered 2D wave
	string cleaned=kenner_out+num2str(wavenum)+"cleaned" // this is the centered 2D wave after removing outliers ("cleaning")
	string fit_params_name = kenner_out+num2str(wavenum)+"fit_params" // this is the fit parameters 
	string quickavg = datasetname+"_avg" // this is the wave produced by avg_wave($datasetname)
	wave W_coef


	/////////////////// quick avg fig  //////////////////////////////////////

	string fit_name = "fit_"+quickavg

	display $quickavg
	fit_transition($quickavg,minx,maxx)
	Label bottom "gate V"
	Label left "csurrent"
	ModifyGraph fSize=24
	ModifyGraph gFont="Gill Sans Light"
	ModifyGraph mode($fit_name)=0,lsize($fit_name)=1,rgb($fit_name)=(65535,0,0)
	ModifyGraph mode($quickavg)=2,lsize($quickavg)=2,rgb($quickavg)=(0,0,0)
	legend
	Legend/C/N=text0/J/A=LB/X=59.50/Y=53.03


	/////////////////// thetas  //////////////////////////////////////


	//find_plot_thetas(wavenum,N,fit_params_name)
	//plot_badthetas($datasetname) // thetas vs repeat plot and bad theta sweep plot
	plot2d_heatmap($datasetname)
	plot2d_heatmap($cleaned)
	plot2d_heatmap($centered)



	/////////////////// plot avg fit  //////////////////////////////////////

	string fit = "fit_"+kenner_out+num2str(wavenum)+"cleaned_avg" //new array
	print fit

	display $avg; W_coef[3]=0
	fit_transition($avg,minx,maxx)

	Label bottom "gate V"
	Label left "csurrent"
	ModifyGraph fSize=24
	ModifyGraph gFont="Gill Sans Light"
	ModifyGraph mode($fit)=0,lsize($fit)=1,rgb($fit)=(65535,0,0)
	ModifyGraph mode($avg)=2,lsize($avg)=2,rgb($avg)=(0,0,0)
	legend

	TileWindows/O=1/C/P
	Legend/C/N=text0/J/A=LB/X=59.50/Y=53.03


	//MultiGraphLayout(WinList("*", ";", "WIN:1"), 3, 20, "AllGraphLayout");


end

Function Chargetransition(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = Amp*tanh((x - Mid)/(2*Theta)) + Linear*x + Const+0*padder
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = Amp
	//CurveFitDialog/ w[1] = Const
	//CurveFitDialog/ w[2] = Theta
	//CurveFitDialog/ w[3] = Mid
	//CurveFitDialog/ w[4] = Linear
	//CurveFitDialog/ w[5] = padder

	return w[0]*tanh((x - w[3])/(2*w[2])) + w[4]*x + w[1]+x^2*w[5]
End

Function CT_faster(w,ys,xs) : FitFunc
	Wave w, xs, ys
	ys= w[0]*tanh((xs - w[3])/(-2*w[2])) + w[4]*xs + w[1]+xs^2*w[5]
End

Function CT2(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = Amp*tanh((x - Mid)/(2*theta)) + Linear*x + Const+Quad*x^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = Amp
	//CurveFitDialog/ w[1] = Const
	//CurveFitDialog/ w[2] = Theta
	//CurveFitDialog/ w[3] = Mid
	//CurveFitDialog/ w[4] = Linear
	//CurveFitDialog/ w[5] = Quad


	return w[0]*tanh(-(x - w[3])/(2*w[2])) + w[4]*x + w[1]+w[5]*x^2
End

Function calc_occ(variable x,variable xo, variable theta) 
	return 0.5*tanh((x - xo)/(2*theta))
end



function dotcond_centering(wave waved, string kenner_out)
	string w2d=nameofwave(waved)
	int wavenum=getfirstnum(w2d)
	string centered=kenner_out+num2str(wavenum)+"centered"
	string fit_params_name = "cond"+num2str(wavenum)+"fit_params"
	wave fit_params = $fit_params_name
	wave new2dwave=$centered
	copyscales waved new2dwave
	new2dwave=interp2d(waved,(x+fit_params[q][2]),(y)) // column 3 is the center fit parameter
End

function/wave sqw_analysis(wave wav, int delay, int wavelen)

// this function separates hot (plus/minus) and cold(plus/minus) and returns  two waves for hot and cold //part of CT
	variable nr, nc
	nr=dimsize(wav,0); print nr
	nc=dimsize(wav,1); print nc
	variable i=0
	variable N
	N=nr/wavelen/4;

	Make/o/N=(nc,(N)) cold1, cold2, hot1, hot2
	wave slice, slice_new

	do
		rowslice(wav,i)
		Redimension/N=(wavelen,4,N) slice
		DeletePoints/M=0 0,delay, slice
		reducematrixSize(slice,0,-1,1,0,-1,4,1,"slice_new")

		cold1[i][]=slice_new[0][0][q]
		cold2[i][]=slice_new[0][2][q]
		hot1[i][]=slice_new[0][1][q]
		hot2[i][]=slice_new[0][3][q]


		i=i+1
	while(i<nc)

	duplicate/o cold1, cold; cold=(cold1+cold2)/2
	duplicate/o hot1, hot; hot=(hot1+hot2)/2

	matrixtranspose hot
	matrixtranspose cold

	CopyScales wav, cold, hot
	
	duplicate/o hot, nument
	nument=cold-hot;

end


function center_demod(int filenum, int delay, int wavelen)
string wname="dat"+num2str(filenum)+"cscurrent_2d_nf_entr";
sqw_analysis($wname,delay,wavelen)//--> nument
wave W_coef, cold, hot  
W_coef[0]= {0.0531997,0.880123,10.688,-12.024,7.28489e-05,7.50215e-08}
W_coef={-0.0039044,0.2688,244.77,-191.16,-9.6209e-07,0}

//ctrans_avg(cold,1,0, "cst", average=0,minx=-1000,maxx=1000)
wave mids
//duplicate/o/r=[][3] ct0fit_params mids
string wname1="dat"+num2str(filenum)+"cscurrentx_2d";

wave demod
demodulate(filenum, 2, "cscurrent_2d_nf")

wave nument
centering(demod,"entropy",mids) // centred plot and average plot
centering(nument,"numentropy",mids) // centred plot and average plot
find_thetas_larger(350) // makes wave called reduced
wave entropy, entropy_avg, numentropy, numentropy_avg, reduced
cleaning(entropy, reduced); 
cleaning(nument, reduced); 

avg_wav(entropy); 
avg_wav(numentropy)
entropy_avg=entropy_avg*2;
//wavetransform/o zapnans entropy_avg
//wavetransform/o zapnans numentropy_avg

//Integrate entropy_avg/D=entropy_avg_INT;
//Integrate numentropy_avg/D=numentropy_avg_INT;

variable factor
factor=calc_scaling( cold, hot,  mids)
//entropy_avg_int=entropy_avg_INT*factor;
//numentropy_avg_int=numentropy_avg_INT*factor;

//execute("intent_graph()")

end


Window hot_cold() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(711,55,1470,497) hot_centr_avg,cold_centr_avg,fit_cold_centr_avg,fit_hot_centr_avg
	ModifyGraph lSize=2
	ModifyGraph lStyle(fit_cold_centr_avg)=7,lStyle(fit_hot_centr_avg)=7
	ModifyGraph rgb(cold_centr_avg)=(0,0,65535),rgb(fit_cold_centr_avg)=(26214,26214,26214)
	ModifyGraph rgb(fit_hot_centr_avg)=(26214,26214,26214)
	TextBox/C/N=CF_cold_centr_avg/X=6.14/Y=6.89 "Coefficient values ± one standard deviation\r\tAmp   \t= -0.049555 ± 6.43e-05\r\tConst \t= 0.88481 ± 2.42e-05"
	AppendText "\tTheta \t= 7.8717 ± 0.0313\r\tMid   \t= 0.15659 ± 0.0306\r\tLinear\t= 9.1999e-05 ± 5.6e-07"
	TextBox/C/N=CF_hot_centr_avg/X=6.46/Y=35.81 "Coefficient values ± one standard deviation\r\tAmp   \t= -0.049523 ± 7.25e-05\r\tConst \t= 0.88481 ± 2.47e-05"
	AppendText "\tTheta \t= 10.005 ± 0.038\r\tMid   \t= -1.4357 ± 0.0354\r\tLinear\t= 9.1736e-05 ± 6.13e-07"
EndMacro


Window theta_graph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(35,53,710,497) meanwave,stdwave,stdwave2
	AppendToGraph goodthetas vs goodthetasx
	AppendToGraph badthetas vs badthetasx
	ModifyGraph gFont="Gill Sans Light"
	ModifyGraph mode(goodthetas)=3,mode(badthetas)=3
	ModifyGraph lSize(goodthetas)=2
	ModifyGraph lStyle(meanwave)=3,lStyle(stdwave)=3,lStyle(stdwave2)=3
	ModifyGraph rgb(meanwave)=(17476,17476,17476),rgb(stdwave)=(52428,1,1),rgb(stdwave2)=(52428,1,1)
	ModifyGraph rgb(goodthetas)=(2,39321,1)
	ModifyGraph fSize=24
	Label left "theta values"
	Label bottom "repeat"
	Legend/C/N=text0/J "\\s(meanwave) mean\r\\s(stdwave) 2*std\r\\s(goodthetas) good\r\\s(badthetas) outliers"
	TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 thetas of dat0"
EndMacro

Window intent_graph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(1471,53,2142,499) entropy_avg,numentropy_avg
	AppendToGraph/R entropy_avg_INT,numentropy_avg_INT
	ModifyGraph lSize(entropy_avg)=2,lSize(numentropy_avg)=2,lSize(entropy_avg_INT)=2
	ModifyGraph lSize(numentropy_avg_INT)=2
	ModifyGraph lStyle(numentropy_avg)=7,lStyle(numentropy_avg_INT)=7
	ModifyGraph rgb(entropy_avg_INT)=(4369,4369,4369),rgb(numentropy_avg_INT)=(4369,4369,4369)
	ModifyGraph zero(right)=15
	Legend/C/N=text1/J/X=67.51/Y=8.48 "\\s(entropy_avg) entropy_avg\r\\s(numentropy_avg) numentropy_avg\r\\s(entropy_avg_INT) entropy_avg_INT"
	AppendText "\\s(numentropy_avg_INT) numentropy_avg_INT"
	SetDrawLayer UserFront
	SetDrawEnv xcoord= axrel,ycoord= right,linethick= 2,linefgc= (65535,0,26214),dash= 7
	DrawLine 0,0.693147,1,0.693147
	SetDrawEnv xcoord= prel,ycoord= right,linethick= 2,linefgc= (1,4,52428)
	DrawLine 0,1.09861,1,1.09861
EndMacro


function calc_scaling(wave cold,wave hot, wave mids)

//first we need to center cold and hot wave
centering(cold,"cold_centr",mids) // centred plot and average plot
centering(hot,"hot_centr",mids) // centred plot and average plot
wave badgammasx
cleaning(cold, badgammasx); 
cleaning(hot, badgammasx); 


wave cold_centr,hot_centr, cold_centr_avg, hot_centr_avg, W_coef
avg_wav(cold_centr);
avg_wav(hot_centr);
 execute("hot_cold()")
//DeletePoints 5,1, W_coef
W_coef={-0.00311856,0.216734,186.321,399.402,-7.41618e-07,0}
wave/t T_Constraints
T_Constraints[0]= {"K0 > -0.1","K0 < 0","K1 > 0","K1 < 0.5","K2 > 0","K3 > -500","K3 < 500"}


    
//    FuncFit/TBOX=768 Chargetransition W_coef cold_centr_avg /D 
    FuncFit/q/H="000000"/TBOX=768 Chargetransition W_coef cold_centr_avg /D /C=T_Constraints 
    duplicate/o W_coef, cold_r
    
    FuncFit/q/H="000000"/TBOX=768 Chargetransition W_coef hot_centr_avg /D /C=T_Constraints 
    duplicate/o W_coef, hot_r
    
 variable   Go= (cold_r[0]+hot_r[0]); print Go
 variable   dT=((hot_r[2]-cold_r[2])); print dT
 variable factor=-1/Go/dT; print factor
 return factor
end

