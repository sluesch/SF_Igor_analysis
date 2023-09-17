#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

Menu "Graph"
	"Close All Graphs/9", CloseAllGraphs()
End
Function CloseAllGraphs()
	String name
	do
		name = WinName(0,1) // name of the front graph
		if (strlen(name) == 0)
			break // all done
		endif
		DoWindow/K $name // Close the graph
	while(1)
End



Function AddLegend(wav,param)
    wave wav        
    string param 
    string graphName
    
        graphName = WinName(0, 1)   // Top graph
    
    
    String list = TraceNameList(graphName, ";", 1)
    String legendText = ""
    Variable numItems = ItemsInList(list)
    Variable i
    for(i=0; i<numItems; i+=1)
        String item = StringFromList(i, list)+"--"+param+num2str(wav[i])
//        if (CmpStr(item,"wave1") == 0)
//            continue            // Skip this trace
//        endif
        String itemText
        sprintf itemText, "\\s(%s) %s", item, item
        if (i > 0)
            legendText += "\r"      // Add CR
        endif
        legendText += itemText
    endfor
    Legend/K/N=text0
    Legend/C/N=MyLegend/W=$graphName legendText
    Legend/C/N=text0/J/A=MT/E=0
End

	string wvname;wvname=stringfromlist(0,imagenamelist("",";"));TextBox/C/N=text1/F=0/A=MT/E wvname


function displayplot2D(start, endnum,whichdat,[delta,xnum, shiftx, shifty])
	variable start, endnum
	string whichdat
	variable delta, xnum, shiftx, shifty
	if(paramisdefault(delta))
		delta=1
	endif
	if(delta==0)
		abort
	endif

	if(paramisdefault(shiftx))
		shiftx=0
	endif
	if(paramisdefault(shifty))
		shifty=0
	endif
		
	variable i=0, totoffx=0, totoffy=0
	string st
	//udh5()
	Display /W=(35,53,960,830)
	i=start
	do
		st="dat"+num2str(i)+whichdat
		appendtograph $st
		wavestats /q $st
		totoffx=shiftx*mod((i-start)/delta,xnum)
		totoffy=shifty*floor((i-start)/delta/xnum)-v_avg
		ModifyGraph offset($st)={totoffx,totoffy}
		i+=delta
	while (i<=endnum)
	makecolorful()
	legend
	Legend/C/N=text0/J/A=RC/E

end

function displayplot(start, endnum,whichdat,prefix[delta,shiftx, shifty])
	variable start, endnum
	string whichdat,prefix
	variable delta, shiftx, shifty
	if(paramisdefault(delta))
		delta=1
	endif
	if(delta==0)
		abort
	endif

	if(paramisdefault(shiftx))
		shiftx=0
	endif
	if(paramisdefault(shifty))
		shifty=0
	endif
		
	variable i=0, totoffx=0, totoffy=0
	string st
	//udh5()
	Display /W=(35,53,960,830)
	i=start
	do
		st=prefix+num2str(i)+whichdat
		appendtograph $st
		ModifyGraph offset($st)={totoffx,totoffy}
		totoffx=totoffx+shiftx
		totoffy=totoffy+shifty
		i+=delta
	while (i<=endnum)
	makecolorful()
	legend
	Legend/C/N=text0/J/A=RC/E

end

function makecolorful([rev, nlines])
	variable rev, nlines
	variable num=0, index=0,colorindex
	string tracename
	string list=tracenamelist("",";",1)
	colortab2wave rainbow
	wave M_colors
	variable n=dimsize(M_colors,0), group
	do
		tracename=stringfromlist(index, list)
		if(strlen(tracename)==0)
			break
		endif
		index+=1
	while(1)
	num=index-1
	if( !ParamIsDefault(nlines))
		group=index/nlines
	endif
	index=0
	do
		tracename=stringfromlist(index, list)
		if( ParamIsDefault(nlines))
			if( ParamIsDefault(rev))
				colorindex=round(n*index/num)
			else
				colorindex=round(n*(num-index)/num)
			endif
		else
			if( ParamIsDefault(rev))
				colorindex=round(n*ceil((index+1)/nlines)/group)
			else
				colorindex=round(n*(group-ceil((index+1)/nlines))/group)
			endif
		endif
		if(colorindex>99)
			colorindex=99
		endif
		ModifyGraph rgb($tracename)=(M_colors[colorindex][0],M_colors[colorindex][1],M_colors[colorindex][2])
		index+=1
	while(index<=num)

end
Function QuickColorSpectrum2()                            // colors traces with 12 different colors
	String Traces    = TraceNameList("",";",1)               // get all the traces from the graph
	Variable Items   = ItemsInList(Traces)                   // count the traces
	Make/FREE/N=(11,3) colors = {{65280,0,0}, {65280,43520,0}, {0,65280,0}, {0,52224,0}, {0,65280,65280}, {0,43520,65280}, {0,15872,65280}, {65280,16384,55552}, {36864,14592,58880}, {0,0,0},{26112,26112,26112}}
	Variable i
	for (i = 0; i <DimSize(colors,1); i += 1)
		ModifyGraph rgb($StringFromList(i,Traces))=(colors[0][i],colors[1][i],colors[2][i])      // set new color offset
	endfor
End

function plot2d_heatmap(wave wav)

	//plots the repeats against the sweeps for dataset cscurrent_2d

	variable num
	string dataset
	string wvname

	wvname=nameOfWave(wav)

	wave wav = $wvname



	display; //start with empty graph
	appendimage wav //append image of data
	ModifyImage $wvname ctab= {*,*,Turbo,0} //setting color (idk why it prefers the pointer)
	ColorScale /A=RC /E width=20 //puts it on the right centre, /E places it outside

	Label bottom "gate(V)"
	Label left "repeats"

	ModifyGraph fSize=24
	ModifyGraph gFont="Gill Sans Light"
	//    ModifyGraph width={Aspect,1.62},height=300
	//	TextBox/C/N=text1/A=MT/E=2 "raw 2D plot of dat" + num2str(num)

end

function setcolorscale2d(percent)
	variable percent
	variable x1, y1, x2, y2, xs, ys, minz, maxz, i=0, j=0
	string filename
	filename=csrwave(A)
	wave mywave = $filename
	x1=pcsr(A)
	y1=qcsr(A)
	x2=pcsr(B)
	y2=qcsr(B)
	duplicate /o/r=[x1,x2][y1,y2] mywave kjhdfgazs7f833jk
	wavestats/q kjhdfgazs7f833jk
	killwaves kjhdfgazs7f833jk
	ModifyImage '' ctab= {V_min,percent*V_max,PlanetEarth,0}
	//killwaves mywave
end



function spectrum_analyzer(wave data, int cutoff, int newplot)
    // This function performs spectrum analysis on a given wave data
    // data: The input wave for spectrum analysis

    // Built-in powerspectrum function
    variable filenum = getfirstnum(nameOfWave(data))
    variable samp_freq = fd_getmeasfreq(filenum);
    if (newplot==1)
    closeallGraphs()
    Display /W=(35,53,1046,664) 
    endif
    

    // Duplicate the data wave to spectrum
    duplicate/o data spectrum

    // Set the x-axis scale of spectrum
    SetScale/P x 0, 1/samp_freq, "", spectrum

    variable nr = dimsize(spectrum, 0)

    wave slice
    wave powerspec, sum_power

    variable i = 0
    variable N, V_avg

    // Slice the spectrum wave at index i
    rowslice(spectrum, i)
    display slice



// Subtract V_avg from slice and calculate the spectrum
wavestats/q slice;slice=slice-V_avg
calculate_spectrum(slice, linear=1)
duplicate/o powerspec, sum_power
string toplot=nameOfWave(data)+"slice"
duplicate/o slice,$toplot

i=1
do
    // Slice the spectrum wave at index i
    rowslice(spectrum,i);
    wavestats/q slice;slice=slice-V_avg

    // Add powerspec and sum_power to sum+power
    sum_power=powerspec+sum_power
    i=i+1
while(i<dimsize(spectrum,1))

// Set the first element of powerspec to 0
sum_power=sum_power/dimsize(spectrum,1);
sum_power[0,x2pnt(sum_power,cutoff)]=0

string spec_name=nameofWave(data)+"spectrum"
duplicate/o sum_power, $spec_name
Resample/RATE=1 $spec_name
Dowindow/F graph0
appendtograph $spec_name; ModifyGraph log(left)=1
string intspecname=nameOfWave(data)+"intspec" 
Integrate sum_power/D=$intspecname;

appendtoGraph/r/w=graph0 $intspecname
makecolorful();legend
//Dowindow/F graph1
//
//appendtoGraph/r/w=graph1 $toplot
//legend




 
 

// display slice

end

function tune_dot_plots(wave filenums, wave nose, wave plunger,string kenner)

// nose and plunger are the gates that are stepped. 
// modify iff statement to filter out gate values that are of interest
variable i=0
variable idx
variable V_npnts
wavestats/q filenums
wave W_coef,result
closeallgraphs()
do
idx=filenums[i]
if (nose[i]<0 && nose[i]>-260)
//if ((mod(check[i],20)==0)&&(check[i]<-250))

string wav_name="dat"+num2str(idx)+kenner;
print wav_name
display; appendimage $wav_name
ModifyImage '' log=1
ModifyImage '' ctab= {0.042,1,ColdWarm,0}

	Label bottom "RC2"
	Label left "CSS"
	TextBox/C/N=text0 wav_name
	TextBox/C/N=text1 "plunger="+num2str(plunger[i])+" / nose="+num2str(nose[i])
	TextBox/C/N=text1/A=LB/X=9.63/Y=17.14



//diffwave($wav_name)


doupdate
endif
i=i+1

while(i<V_npnts)
//TileWindows/A=(2,4)/O=1 
TileWindows/O=1/C/P

end


