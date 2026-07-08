!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 2    last modified 2013-03-14 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2013
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal unknown 
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front linetype -1 linewidth 1.000
#set timefmt z "%d/%m/%y,%H:%M"
set zdata 
#set timefmt y "%d/%m/%y,%H:%M"
set ydata 
#set timefmt x "%d/%m/%y,%H:%M"
set xdata 
#set timefmt cb "%d/%m/%y,%H:%M"
#set timefmt y2 "%d/%m/%y,%H:%M"
set y2data 
#set timefmt x2 "%d/%m/%y,%H:%M"
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc lt -3 fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0, 0 
set style ellipse size graph 0.05, 0.03, first 0 angle 0 units xy
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set format cb "% g"
set format r "% g"
set angles radians
unset grid
set raxis
set key title ""
set key inside right top vertical Right noreverse enhanced autotitles nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title  offset character 0, 0, 0
#set logscale cb
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 1 #1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis linetype -2 linewidth 1.000
set yzeroaxis linetype -2 linewidth 1.000
set zzeroaxis linetype -2 linewidth 1.000
set x2zeroaxis linetype -2 linewidth 1.000
set y2zeroaxis linetype -2 linewidth 1.000
set ticslevel 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set xtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set xtics autofreq  norangelimit
set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set ytics autofreq  norangelimit
set ztics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set ztics autofreq  norangelimit
set nox2tics
set noy2tics
set cbtics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0 autojustify
set cbtics autofreq  norangelimit
set cbtics ("Mutant" -XXXcbmaxXXX, "-XXXcbmaxMidXXX" -XXXcbmaxMidXXX, "XXXcbmaxMidXXX" XXXcbmaxMidXXX, "WT" XXXcbmaxXXX)
set rtics axis in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set rtics autofreq  norangelimit
set title "" 
set title  offset character 0, 0, 0 font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  offset character 0, 0, 0 font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set x2label "" 
set x2label  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "" 
set ylabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
set yrange [ * : * ] reverse nowriteback
set y2range [ * : * ] reverse nowriteback
set zlabel "" 
set zlabel  offset character 0, 0, 0 font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  offset character 0, 0, 0 font "" textcolor lt -1 rotate by -270
#set cbrange [ 0.001*XXXminFactorXXX : 100.*XXXmaxFactorXXX ] noreverse nowriteback
set cbrange [ -XXXcbmaxXXX : XXXcbmaxXXX ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit noerrorvariables

#set palette defined (XXXcbminXXX "", -1. "gray", -1. "dark-blue", XXXcbmaxXXX*0.05*XXXfactorXXX "blue", XXXcbmaxXXX*0.10*XXXfactorXXX "white", XXXcbmaxXXX*0.30*XXXfactorXXX "red", XXXcbmaxXXX*XXXfactorXXX "dark-red") # Linear

#set cbrange [ XXXcbminXXX : XXXcbmaxXXX ] noreverse nowriteback
#set cbrange [ 0.1 : XXXcbmaxXXX ] noreverse nowriteback

#set palette defined (-1.5 "gray", -1. "gray", -1 "dark-blue", 0.01 "dark-blue", 0.01 "blue", 0.02 "blue", 0.02 "white", 0.035 "white", 0.035 "red", 0.045 "red", 0.045 "dark-red", XXXcbmaxXXX "dark-red") # Linear
#set palette defined (0 "dark-blue", 0.25 "blue", 0.50 "white", 0.75 "red", 1.0 "dark-red") # Log
#set palette defined (0 "dark-blue", 0.15 "blue", 0.30 "white", 0.45 "red", 1.0 "dark-red") # Log cHiC

set palette defined (0 "dark-blue", 0.25 "blue", 0.50 "white", 0.75 "red", 1.0 "dark-red") # Log cHiC
#set palette defined (0 "white", 0.5 "white", 0.5 "dark-red", 1.0 "dark-red") # Log cHiC

#set palette defined (0 "white", XXXscaleXXX "red", 1.0 "dark-red") # Log cHiC
set colorbox

set term post color enhanced


set yrange  [ XXXsizeXXX-0.5 : -0.50000 ] noreverse nowriteback
set y2range [ XXXsizeXXX-0.5 : -0.50000 ] noreverse nowriteback
set xrange  [ -0.50000 : XXXsizeXXX-0.5 ] noreverse nowriteback
set x2range [ -0.50000 : XXXsizeXXX-0.5 ] noreverse nowriteback
#set yrange  [ -0.50000 : XXXsizeXXX-0.5 ] noreverse nowriteback
#set y2range [ -0.50000 : XXXsizeXXX-0.5 ] noreverse nowriteback

#set x2label "Bins of 50kb"
unset ylabel

#set arrow 1 from graph XXXRstartXXX , 0 to graph XXXRstartXXX , XXXRstartXXX nohead filled lt 1 lw 1 lc rgb "black" front
#set arrow 2 from graph 0 , XXXRstartXXX to graph XXXRstopXXX  , XXXRstartXXX nohead filled lt 1 lw 1 lc rgb "black" front
#set arrow 3 from graph XXXRstopXXX  , 0 to graph XXXRstopXXX  , XXXRstartXXX nohead filled lt 1 lw 1 lc rgb "black" front
#set arrow 4 from graph 0 ,  XXXRstopXXX to graph XXXRstopXXX  , XXXRstopXXX  nohead filled lt 1 lw 1 lc rgb "black" front

#set xtics ("XXXRstartXXX" XXXoffsetXXX, "XXXRstopXXX" XXXoffsetXXX+XXXRsizeXXX)
#set ytics ("XXXRstartXXX" XXXoffsetXXX, "XXXRstopXXX" XXXoffsetXXX+XXXRsizeXXX)
set xtics ("53.7" 0, "54.2" 300, "56.7" 600)
set ytics ("53.7" 0, "54.2" 300, "56.7" 600)
unset xtics
unset ytics
set colorbox
set output "contact_map.ps"
#unset xtics
#unset ytics



#plot "< awk 'BEGIN{for(i=0;i<XXXsizeXXX;i++)for(j=0;j<XXXsizeXXX;j++) m[i,j]=0.0}{v=$3;m[$1,$2]=v}END{for(i=0;i<XXXsizeXXX;i++)for(j=0;j<XXXsizeXXX;j++) print i,j,m[i,j]}' contact_map.tab" u 1:2:3 notitle w image pixel
plot 'contact_map.tab' u 1:2:3 notitle w image , \
     "< cat ./scripts/barriers.txt | awk '{print ($2-XXXstartXXX)/XXXresolutionXXX,($2-XXXstartXXX)/XXXresolutionXXX}'" u 1:2 not w p pt 7 ps 0.5 lt 1 lw 3 lc rgb "yellow"

#     "< grep -v small ../scripts/barriers_XXXconditionXXX.txt | grep -w plus_minus" u 1:2 not w p pt 7 ps 0.5 lt 1 lw 3 lc rgb "black", \
#     "< grep -v small ../scripts/barriers_XXXconditionXXX.txt | grep -w plus" u 1:2 not w p pt 7 ps 0.5 lt 1 lw 3 lc rgb "yellow", \
#     "< grep -v small ../scripts/barriers_XXXconditionXXX.txt | grep -w minus" u 1:2 not w p pt 7 ps 0.5 lt 1 lw 3 lc rgb "orange"

#    EOF
