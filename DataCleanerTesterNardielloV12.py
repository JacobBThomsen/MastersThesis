# DataCleaner V11 Nardiello
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
import os #to jumo around dir
import pandas as pd #For supressing e+04 data type when converting to txt
#%%Start Variables for grid plot
GridNum=0
#GridData = [[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8],[1],[2],[3],[4],[5],[6],[7],[8]]
GridData = [1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]
MaxNum=0
MaxData=[1,2,3,4,5,6,7,8]
#To check variables inside function
Check=[]
#%%#Getting the data 
#"C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Nardiello\V11_Nardiello\pathos_tess_lightcurve_tic-00461601177-s018_tess_v1_llc.fits"
#Raw data
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Nardiello\V12_Nardiello') 
a18 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s018_tess_v1_llc.fits' )
a20 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s020_tess_v1_llc.fits' )
a25 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s025_tess_v1_llc.fits' )
a26 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s026_tess_v1_llc.fits' )
a40 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s040_tess_v1_llc.fits' )
a52 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s052_tess_v1_llc.fits' )
a53 = fits.getdata( 'pathos_tess_lightcurve_tic-00461618602-s053_tess_v1_llc.fits' )

per = 6.5043035 #period 

ph = 'best_phot_flux_raw'
ph = 'best_phot_flux_cor'


# Finding the magnitude via the flux
m18 = 25.0 - 2.5*np.log10( a18[ ph ] )
m20 = 25.0 - 2.5*np.log10( a20[ ph ] )
m25 = 25.0 - 2.5*np.log10( a25[ ph ] )
m26 = 25.0 - 2.5*np.log10( a26[ ph ] )
m40 = 25.0 - 2.5*np.log10( a40[ ph ] )
m52 = 25.0 - 2.5*np.log10( a52[ ph ] )
m53 = 25.0 - 2.5*np.log10( a53[ ph ] )
#Finding the phase time of data points
ph18 = np.mod( a18['time'], per )
ph20 = np.mod( a20['time'], per )
ph25 = np.mod( a25['time'], per )
ph26 = np.mod( a26['time'], per )
ph40 = np.mod( a40['time'], per )
ph52 = np.mod( a52['time'], per )
ph53 = np.mod( a53['time'], per )
#cutting data
i18 = np.where( (ph18 > 4.1 ) & (ph18 < 4.4) )[0]
i20 = np.where( (ph20 > 4.1 ) & (ph20 < 4.4) )[0]
i25 = np.where( (ph25 > 4.1 ) & (ph25 < 4.4) )[0]
i26 = np.where( (ph26 > 4.1 ) & (ph26 < 4.4) )[0]
i40 = np.where( (ph40 > 4.1 ) & (ph40 < 4.4) )[0]
i52 = np.where( (ph52 > 4.1 ) & (ph52 < 4.4) )[0]
i53 = np.where( (ph53 > 4.1 ) & (ph53 < 4.4) )[0]
#Change from median
me18=m18 - np.median(m18[i18])
me20=m20 - np.median(m20[i20])
me25=m25 - np.median(m25[i25])
me26=m26 - np.median(m26[i26])
me40=m40 - np.median(m40[i40])
me52=m52 - np.median(m52[i52])
me53=m53 - np.median(m53[i53])
#Data groups
Times = [a18['time'], a20['time'], a25['time'], a26['time'], a40['time'], a52['time'], a53['time']]
PhaseTimes= [ph18, ph20, ph25, ph26, ph40, ph52, ph53]
Mangnitudes = [me18, me20, me25, me26, me40, me52, me53]
Sectors = ['18','20','25','26','40','52','53']
Num =[0,1,2,3,4,5,6]
#%% PLot all data
plt.plot( ph18, me18, '.',  label = 'Sector 18', markersize=2)
plt.plot( ph20, me20, '.' , label = 'Sector 20', markersize=2)
plt.plot( ph25, me25, '.' , label = 'Sector 25', markersize=2)
plt.plot( ph26, me26, '.' , label = 'Sector 26', markersize=2)
plt.plot( ph40, me40, '.' , label = 'Sector 40', markersize=2)
plt.plot( ph52, me52, '.' , label = 'Sector 52', markersize=2)
plt.plot( ph53, me53, '.' , label = 'Sector 53', markersize=2)
plt.ylim( -0.4, 0.8 )
ax = plt.gca() 
ax.invert_yaxis()
plt.title( 'NGC 188 - V12 - Raw Data' )
plt.xlabel( 'Phase time (Days)' )
plt.ylabel( 'Δ Magnitude compared to median' )
plt.legend(fontsize="small", loc="lower right")
#%% Plot single Sector function
def Raw_Plot(i):
    plt.figure()
    plt.plot( Times[i], Mangnitudes[i], '.' , label = 'Sector' )
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V12 - Raw Data Sector:'+Sectors[i] )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Δ Magnitude compared to median' )
    plt.legend(fontsize="small", loc="lower left")
#%% Plot all sectors raw
for i in Num:
    Raw_Plot(i)

#%%Clean data of very high values
for i in Num:
    RollM = np.roll(Mangnitudes[i],1)
    DiffM = Mangnitudes[i]-RollM #Over 0.15 er vist dårligt
    pos = np.where( (0.17 > abs(DiffM) ))[0]
    Mangnitudes[i]=Mangnitudes[i][pos]
    Times[i] = Times[i][pos]

#%% Isolating and analysing the peaks
def Nardiello_data( i ):
    #import Define data
    t = Times[i]
    m = Mangnitudes[i]
    s = Sectors[i]
    
    #Find first eclipse
    if i == 4: #due to its nine peaks
        Times[i]=t[60:3922]
        Mangnitudes[i] = m[60:3922]
        StartEclip=next(x[0] for x in enumerate(Mangnitudes[i]) if x[1] > 0.23) #Position of first eclipse
        StartTime=t[StartEclip] #Time of start eclipse Aperture    
    else: #For all the other normal ones
        StartEclip=next(x[0] for x in enumerate(Mangnitudes[i]) if x[1] > 0.23) #Position of first eclipse
        StartTime=t[StartEclip] #Time of start eclipse Aperture    
    
    TimeStart = StartTime-0.5 #Peak start
    TimeStop = StartTime+0.5 #Peak end
    #cut eclipse
    Startindexin = [i for i in range(len(Mangnitudes[i])) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
    
    MaxVal=max(Mangnitudes[i][Startindexin]) #max value of start eclip
    #The position of the mimimum of the first eclipse
    MaxPos=([index for index, item in enumerate(Mangnitudes[i]) if item == MaxVal]) #Position of max value for Aperture Photometry
    TimeEclipse=t[MaxPos] #Time of eclipse
    
    #Finding all peaks
    inclination = 88.62 #Inclination for V12, tabel 6 Frank paper
    Period = 6.504 #Period of V12
    EndEclipse = t[-1]
    AllPeak=np.arange(TimeEclipse, EndEclipse, (Period/2)).tolist() #list of all peak positins in time
    
    #Loop through all peaks
    for AllPeak in AllPeak:
        # Cut out wanted data into groups (<>)30
        TimeStart = AllPeak-0.3 #Peak start
        TimeStop = AllPeak+0.3 #Peak end
        TimeStartCut = TimeStart - 0.4 #Cutoff for fitting data before peak
        TimeStopCut = TimeStop + 0.4 #Cutoff for fitting data after peak
    
        #cut eclipse 
        ineclipse= [i for i in range(len(t)) if (t[i]>TimeStart and t[i]<TimeStop)] #The eclipse in time
        indexin = [i for i in range(len(Mangnitudes[i])) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
        inMag = Mangnitudes[i][indexin] #The eclipse magnitude values
        
        #Cut fitting data not part of eclipse
        outeclipseBefore= [i for i in range(len(t)) if (t[i]>TimeStartCut and t[i]<TimeStart)] #Before eclipse
        indexBefore = [i for i in range(len(Mangnitudes[i])) if (t[i]>TimeStartCut and t[i]<TimeStart)] #The index of the eclipse
        outeclipseAfter= [i for i in range(len(t)) if (t[i]>TimeStop and t[i]<TimeStopCut)] #after eclipse
        indexAfter = [i for i in range(len(Mangnitudes[i])) if (t[i]>TimeStop and t[i]<TimeStopCut)] #The index of the eclipse
        indexOut = indexBefore+indexAfter #Total Out index
        outeclipse = outeclipseBefore+outeclipseAfter #Total out times
        outMag = Mangnitudes[i][indexOut] #non eclipse magnitude values - PSF Photometry
        
        #Fit qudratic equation to outeclipse data
        if len(t[outeclipse])==0: #To make sure it ignores empty data areas
            def Funk(C): return 0*C
            Quard = Funk
        else:
            QuardFit = np.polyfit(t[outeclipse], outMag, 1) #linear fit 
        
            Quard = np.poly1d(QuardFit) #Equation, eks 0.008607 x + 2.793

        #For nice equations in legend
        Num1=round(QuardFit[0],10)
        Num2=round(QuardFit[1],6)
        #NumP3=round(QuardFitP[2],3) #For Quardratic
        #NumA3=round(QuardFitA[2],3)
        #Normalzing the data
        NormaOut=Quard(t[indexOut]) #Difference between zero and fit
        NormaIN=Quard(t[indexin]) #Difference between zero and fit 
        #Normalized data
        NormaDataOut= outMag-NormaOut #Normalized outside 
        NormaDataIn= inMag-NormaIN #Normalized  outside 
        #Plot everything
        #Raw data
        plt.figure()
        plt.plot(t[ineclipse], inMag,'r.', label = 'Inside Eclipse') 
        plt.plot(t[outeclipse], outMag,'m.', label = 'Outside Eclipse')
        #Fit Equations
        if len(t[outeclipse])==0:
            Quardline = np.linspace(0, 0, 10)
        else:    
            Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
        plt.plot(Quardline, Quard(Quardline), label= str(Num1) + "x + " + str(Num2) )
        #Make nice
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Data and Fits' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Δ Magnitude compared to median' )
        plt.legend(fontsize="small", loc="lower left", bbox_to_anchor=(1, 0.5)) #bbox for out of plot legend
        
        #Plot Only fits
        plt.figure()
        plt.plot(t[outeclipse], outMag,'.', label = 'Outside Eclipse')
        #Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
        plt.plot(Quardline, Quard(Quardline), label= str(Num1) + "x + " + str(Num2) )
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Fits' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Δ Magnitude compared to median' )
        plt.legend(fontsize="small", loc="lower center")
        
        #Plot data used for normalization 
        plt.figure()
        plt.plot(t[indexOut], NormaDataOut,'r.', label = 'Inside Eclipse')
        plt.plot(Quardline,np.array([0] * len(Quardline)),'k')
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Normalized data used in normalization' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Normalized Δ magnitude' )
        plt.legend(fontsize="small", loc="lower center")
        
        #Plot all normalized data
        plt.figure()
        plt.plot(t[indexOut], NormaDataOut,'r.', label = 'Outside Eclipse')
        plt.plot(t[indexin], NormaDataIn,'m.', label = 'Inside Eclipse')
        plt.plot(Quardline,np.array([0] * len(Quardline)),'k')  #Zero line
        plt.plot(Quardline,np.array([0.6] * len(Quardline)),'m')  #PSF max top line
        plt.plot(Quardline,np.array([0.66] * len(Quardline)),'c')  #Aperature max top line
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Normalized data' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Normalized Δ magnitude' )
        plt.legend(fontsize="small", loc="lower left")
        
        #Save data for grid plot all data
        T_out=t[indexOut]
        T_in=t[indexin]
        
        global GridNum #To ensure it works as a variable
        global MaxNum #To ensure it works as a variable
        
        GridData[GridNum]=T_out.tolist() #to.list() to ensure that list in list works 
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataOut.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=T_in.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataIn.tolist()
        GridNum=GridNum+1
        
        #Save data for grid plot max point data
        if len(t[outeclipse])==0:
            MaxData[MaxNum]=0
            MaxNum=MaxNum+1
        else:
            MaxVal=max(NormaDataIn) #max value of Photometry start eclip
            MaxPos=([index for index, item in enumerate(NormaDataIn) if item == MaxVal]) #Position of max value
            Time=T_in[MaxPos] #Time of eclipse 
            TimeEclipse= Time #Their average
            MaxData[MaxNum]=TimeEclipse #to.list() to ensure that list in list works 
            MaxNum=MaxNum+1
        
    # Plot the prison diagram
    #%matplotlib qt
    fig = plt.figure() #%matplotlib qt #To make it pop out
    fig.set_size_inches(9, 7)
    gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
    (ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
    TopText=fig.suptitle('Nardiello - V12: Normalised Δ Magnitude Sector:'+Sectors[i], fontsize=16, weight='bold')
    TopText.set_position([.5, .95])
    LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
    LeftText.set_position([0.05, .5])
    BotText=fig.supxlabel('Time (BJD-2457000)', fontsize=14)
    BotText.set_position([.5, 0.05])
    ax1.plot(GridData[0], GridData[1],'r.', label='_nolegend_')
    ax1.plot(GridData[2], GridData[3],'b.', label='_nolegend_')
    ax1.invert_yaxis()
    Quardline = np.linspace(GridData[0][0], GridData[0][-1], 10)  #GridData[4][0] for S26
    ax1.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax1.set_xlim(GridData[0][0], GridData[0][-1]) #GridData[4][0] for S26
    ax1.set_ylim(0.8, -0.05)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)
    
    
    Placement=round(float(MaxData[0]),3) 
    Placement1 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax1.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax2.plot(GridData[4], GridData[5],'r.', label='_nolegend_')
    ax2.plot(GridData[6], GridData[7],'b.', label='_nolegend_')
    Quardline = np.linspace(GridData[4][0], GridData[4][-1], 10) 
    ax2.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax2.set_xlim(GridData[4][0], GridData[4][-1])
    ax2.set_ylim(0.8, -0.05)
    ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    

    Placement=round(float(MaxData[1]),3) 
    Placement2 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax2.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax3.plot(GridData[8], GridData[9],'r.', label='_nolegend_')
    ax3.plot(GridData[10], GridData[11],'b.', label='_nolegend_')
    ax3.invert_yaxis()
    Quardline = np.linspace(GridData[8][0], GridData[8][-1], 10)   #GridData[22][-1] for S18
    ax3.plot(Quardline,np.array([0.56] * len(Quardline)),'m')  #PSF max top line
    ax3.set_xlim(GridData[8][0], GridData[8][-1])  #GridData[22][-1] for S18
    ax3.set_ylim(0.8, -0.05)
    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax3.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[2]),3) 
    Placement3 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax3.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax4.plot(GridData[12], GridData[13],'r.', label='_nolegend_')
    ax4.plot(GridData[14], GridData[15],'b.', label='_nolegend_')
    if GridData[12]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax4.set_xlim(2925, 2926)
    else: 
        Quardline = np.linspace(GridData[12][0], GridData[12][-1], 10)
        ax4.set_xlim(GridData[12][0], GridData[12][-1])
    ax4.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax4.set_ylim(0.8, -0.05)
    ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[3]),3) 
    Placement4 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax4.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax5.plot(GridData[16], GridData[17],'r.', label='_nolegend_')
    ax5.plot(GridData[18], GridData[19],'b.', label='_nolegend_')
    ax5.invert_yaxis()
    if GridData[16]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax5.set_xlim(1855, 1856)
    else: 
        Quardline = np.linspace(GridData[16][0], GridData[16][-1], 10) #GridData[36][0] for S50
        ax5.set_xlim(GridData[16][0], GridData[16][-1]) #GridData[36][0] for S50
    ax5.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax5.set_ylim(0.8, -0.05)
    ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax5.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[4]),3) 
    Placement5 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax5.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax6.plot(GridData[20], GridData[21],'r.', label='_nolegend_')
    ax6.plot(GridData[22], GridData[23],'b.', label='_nolegend_')
    Quardline = np.linspace(GridData[20][0], GridData[20][-1], 10) 
    ax6.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax6.set_xlim(GridData[20][0], GridData[20][-1])
    ax6.set_ylim(0.8, -0.05)
    ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[5]),3) 
    Placement6 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax6.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax7.plot(GridData[24], GridData[25],'r.', label='_nolegend_')
    ax7.plot(GridData[26], GridData[27],'b.', label='_nolegend_')
    ax7.invert_yaxis()
    Quardline = np.linspace(GridData[24][0], GridData[24][-1], 10)  #GridData[54][-1] for S18
    ax7.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax7.set_xlim(GridData[24][0], GridData[24][-1])    #GridData[54][-1] for S18
    ax7.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    ax7.set_ylim(0.8, -0.05)
    ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax7.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[6]),3) 
    Placement7 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax7.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax8.plot(GridData[28], GridData[29],'r.', label='_nolegend_')
    ax8.plot(GridData[30], GridData[31],'b.', label='_nolegend_')  
    if GridData[28]==[] or isinstance(GridData[28], int)==True  :#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax8.set_xlim(2928, 2929)
    else: 
        Quardline = np.linspace(GridData[28][0], GridData[28][-1], 10)
        ax8.set_xlim(GridData[28][0], GridData[28][-1])
    ax8.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  # max top line

    ax8.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    ax8.set_ylim(0.8, -0.05)
    ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[7]),3) 
    Placement8 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax8.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    for ax in fig.get_axes(): #o make grid plot stuff work
        ax.label_outer()


    #Datafiles i added togehter
    l=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
    #j=[0,1,2,3]
    All_Magnitudes=[]
    All_Time=[]
    for l in l:
        if l == 0 or l == 1:
            All_Time=GridData[0]
            All_Magnitudes=GridData[1]
            continue
        else:
            if (isinstance(GridData[l], int)) == True:
                continue
            else: 
                if l % 2 == 0:
                    All_Time=All_Time+GridData[l]
                else:
                    All_Magnitudes=All_Magnitudes+GridData[l]


    All_Data = np.matrix([All_Magnitudes,All_Time])
    All_Data=np.transpose(All_Data) #Swap rows and collums
    #Placement data
    Placements=[Placement1,Placement2,Placement3,Placement4,Placement5,Placement6,Placement7,Placement8]
    #Make into csv file
    mat = All_Data
    df = pd.DataFrame(data=mat.astype(float))
    df.to_csv('Nardiello_V12_all_peaks_S'+Sectors[i]+'.csv', sep=' ', header=False, float_format='%.10f', index=False)
    
    mat = Placements
    df = pd.DataFrame(mat) #pd.DataFrame(data=mat.astype(float))
    df.to_csv('Nardiello_V12_all_placements_S'+Sectors[i]+'.csv', header=False, float_format='%.10f', index=False)
   
    return All_Magnitudes, All_Time, All_Data, GridData

#%% Prison plot for one

Magni, time, All_Data, GridData = Nardiello_data(6)


#%% Plot all sectors
'''
for i in Num:
    Magni, time, All_Data, GridData = Nardiello_data(i)
'''
 
