# DataCleaner V12 TGLC
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
#%%Start Variables for grid plot
GridNum=0
GridData = [1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]
MaxNum=0
MaxData=[1,2,3,4,5,6,7,8]
#%% Get data from TLGC files function
def get_tglc_phot( filename ):

    hdul         = fits.open( filename )
    hdul.info()

    filt              = list(hdul[1].data['TESS_flags'] == 0) and list(hdul[1].data['TGLC_flags'] == 0)
    # filter out bad datapoints from both TESS FFI flags and TGLC flags

    time           = hdul[1].data['time'][filt]          # Time ... BJD ??

    psf_flux       = hdul[1].data['psf_flux'][filt]      # raw psf flux
    psf_mag        = 25.0 - 2.5*np.log10( psf_flux )
    psf_flux_err   = hdul[1].header['PSF_ERR']           # raw psf flux error
    psf_mag_err    = -1.0857*psf_flux_err/psf_flux       # raw aper flux error

    aper_flux      = hdul[1].data['aperture_flux'][filt] # raw aper flux
    aper_mag       = 25.0 - 2.5*np.log10( aper_flux )
    aper_flux_err  = hdul[1].header['APER_ERR']          # raw aper flux error
    aper_mag_err   = -1.0857*aper_flux_err/aper_flux     # raw aper flux error

    return time, psf_mag, aper_mag 
#%%
def TGLC_data( filename ):
    #import get_tglc_phot
    t, p, a = get_tglc_phot(filename) #Data extraction of time, PSF, and Aperature
    #Check to see if correct file
    plt.figure() #This is all data
    plt.plot(t, p,'r.', label = 'PSF Photometry')
    plt.plot(t, a,'b.', label = 'Aperature Photometry')
    plt.title( 'NGC 188 - V12 - Raw Data' )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="upper left")
    
    #For aperature
    StartEclipA=next(x[0] for x in enumerate(a) if x[1] > 19) #Position of first eclipse
    StartTimeA=t[StartEclipA] #Time of start eclipse Aperture       
    #For PSF
    StartEclipP=next(x[0] for x in enumerate(p) if x[1] > 19) #Position of first eclipse
    StartTimeP=t[StartEclipP] #Time of start eclipse Aperture     
    #Finding first Eclipse max
    StartEclipse= (StartTimeA+StartTimeP)/2 #Their average
    TimeStart = StartEclipse-0.5 #Peak start
    TimeStop = StartEclipse+0.5 #Peak end
    #cut eclipse 
    #Startineclipse= [i for i in range(len(t)) if (t[i]>TimeStart and t[i]<TimeStop)] #The eclipse in time
    Startindexin = [i for i in range(len(p)) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
    
    MaxValA=max(a[Startindexin]) #max value of Aperture Photometry start eclip
    MaxValP=max(p[Startindexin]) #max value of PSF Photometry start eclip
   

    #The position of the mimimum of the first eclipse
    MaxPosA=([index for index, item in enumerate(a) if item == MaxValA]) #Position of max value for Aperture Photometry
    MaxPosP=([index for index, item in enumerate(p) if item == MaxValP]) #Position of max value for PSF Photometry
    TimeA=t[MaxPosA] #Time of eclipse Aperture
    TimeP=t[MaxPosP] #Time of eclipse PSF
    TimeEclipse= (TimeA+TimeP)/2 #Their average

    #Finding all peak
    inclination = 88.62 #Inclination for V12, tabel 6 Frank paper
    Period = 6.504 #Period of V12
    EndEclipse = t[-1]
    AllPeak=np.arange(TimeEclipse, EndEclipse, (Period/2)).tolist() #list of all peak positins in time
    
    #Loop through all peaks
    for AllPeak in AllPeak:
        # Cut out wanted data into groups (<>)30
        TimeStart = AllPeak-0.2 #Peak start
        TimeStop = AllPeak+0.2 #Peak end
        TimeStartCut = TimeStart - 0.4 #Cutoff for fitting data before peak
        TimeStopCut = TimeStop + 0.4 #Cutoff for fitting data after peak
    
        #cut eclipse 
        ineclipse= [i for i in range(len(t)) if (t[i]>TimeStart and t[i]<TimeStop)] #The eclipse in time
        indexin = [i for i in range(len(p)) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
        inMagP = p[indexin] #The eclipse magnitude values - PSF Photometry
        inMagA = a[indexin] #The eclipse magnitude values - Aperture Photometry
        
        #Cut fitting data not part of eclipse
        outeclipseBefore= [i for i in range(len(t)) if (t[i]>TimeStartCut and t[i]<TimeStart)] #Before eclipse
        indexBefore = [i for i in range(len(p)) if (t[i]>TimeStartCut and t[i]<TimeStart)] #The index of the eclipse
        outeclipseAfter= [i for i in range(len(t)) if (t[i]>TimeStop and t[i]<TimeStopCut)] #after eclipse
        indexAfter = [i for i in range(len(p)) if (t[i]>TimeStop and t[i]<TimeStopCut)] #The index of the eclipse
        indexOut = indexBefore+indexAfter #Total Out index
        outeclipse = outeclipseBefore+outeclipseAfter #Total out times
        outMagP = p[indexOut] #non eclipse magnitude values - PSF Photometry
        outMagA = a[indexOut] #non eclipse magnitude values - Aperture Photometry
        
        
        
        #Fit qudratic equation to outeclipse data
        if len(t[outeclipse])==0: #To make sure it ignores empty data areas
            def Funk(C): return 0*C
            QuardP = Funk
            QuardA = Funk
        else:
            QuardFitP = np.polyfit(t[outeclipse], outMagP, 1) #linear fit for PSF 
            QuardFitA = np.polyfit(t[outeclipse], outMagA, 1) #linear fit for Aperature
        
            QuardP = np.poly1d(QuardFitP) #Equation, eks 0.008607 x + 2.793
            QuardA = np.poly1d(QuardFitA)

        #For nice equations in legend
        NumP1=round(QuardFitP[0],10)
        NumP2=round(QuardFitP[1],6)
        #NumP3=round(QuardFitP[2],3) #For Quardratic
        NumA1=round(QuardFitA[0],10)
        NumA2=round(QuardFitA[1],6)
        #NumA3=round(QuardFitA[2],3)
        #Normalzing the data
        NormaPOut=QuardP(t[indexOut]) #Difference between zero and fit PSF outside 
        NormaAOut=QuardA(t[indexOut]) #Difference between zero and fit Aperature outside
        NormaPIn=QuardP(t[indexin]) #Difference between zero and fit PSF inside 
        NormaAIN=QuardA(t[indexin]) #Difference between zero and fit Aperature inside
        #Normalized data
        NormaDataPOut= outMagP-NormaPOut #Normalized PSF outside 
        NormaDataAOut= outMagA-NormaAOut #Normalized Aperature inside 
        NormaDataPIn= inMagP-NormaPIn #Normalized PSF outside 
        NormaDataAIn= inMagA-NormaAIN #Normalized Aperature inside 
        #Plot everything
        #Raw data
        plt.figure()
        plt.plot(t[ineclipse], inMagP,'r.', label = 'PSF Photometry') 
        plt.plot(t[ineclipse], inMagA,'b.',label = 'Aperture Photometry' )
        plt.plot(t[outeclipse], outMagP,'m.', label = 'PSF Photometry')
        plt.plot(t[outeclipse], outMagA,'c.',label = 'Aperture Photometry' )
        #Fit Equations
        if len(t[outeclipse])==0:
            Quardline = np.linspace(0, 0, 10)
        else:    
            Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
        plt.plot(Quardline, QuardP(Quardline), label= str(NumP1) + "x + " + str(NumP2) )
        plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x+ " + str(NumA2) )
        #Make nice
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Data and Fits' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Magnitude' )
        plt.legend(fontsize="small", loc="lower left", bbox_to_anchor=(1, 0.5)) #bbox for out of plot legend
        
        #Plot Only fits
        plt.figure()
        plt.plot(t[outeclipse], outMagP,'.', label = 'PSF Photometry')
        plt.plot(t[outeclipse], outMagA,'.',label = 'Aperture Photometry' )
        #Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
        plt.plot(Quardline, QuardP(Quardline), label= str(NumP1) + "x + " + str(NumP2) )
        plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x + " + str(NumA2) )
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Fits' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Magnitude' )
        plt.legend(fontsize="small", loc="lower center")
        
        #Plot data used for normalization 
        plt.figure()
        plt.plot(t[indexOut], NormaDataPOut,'r.', label = 'PSF Photometry')
        plt.plot(t[indexOut], NormaDataAOut,'b.', label = 'Aperture Photometry')
        plt.plot(Quardline,np.array([0] * len(Quardline)),'k')
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Normalized data used in normalization' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Normalized flux' )
        plt.legend(fontsize="small", loc="lower center")
        
        #Plot all normalized data
        plt.figure()
        plt.plot(t[indexOut], NormaDataPOut,'r.', label = 'PSF Photometry')
        plt.plot(t[indexOut], NormaDataAOut,'b.', label = 'Aperture Photometry')
        plt.plot(t[indexin], NormaDataPIn,'m.', label = 'PSF Photometry')
        plt.plot(t[indexin], NormaDataAIn,'c.', label = 'Aperture Photometry')
        plt.plot(Quardline,np.array([0] * len(Quardline)),'k')  #Zero line
        plt.plot(Quardline,np.array([0.6] * len(Quardline)),'m')  #PSF max top line
        plt.plot(Quardline,np.array([0.66] * len(Quardline)),'c')  #Aperature max top line
        ax = plt.gca() 
        ax.invert_yaxis()
        plt.title( 'NGC 188 - V12 - Normalized data' )
        plt.xlabel( 'Time [BJD-2457000]' )
        plt.ylabel( 'Normalized flux' )
        plt.legend(fontsize="small", loc="lower left")
        
        #Save data for grid plot all data
        T_out=t[indexOut]
        T_in=t[indexin]
        
        global GridNum #To ensure it works as a variable
        global MaxNum #To ensure it works as a variable
        
        GridData[GridNum]=T_out.tolist() #to.list() to ensure that list in list works 
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataPOut.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=T_out.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataAOut.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=T_in.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataPIn.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=T_in.tolist()
        GridNum=GridNum+1
        GridData[GridNum]=NormaDataAIn.tolist()
        GridNum=GridNum+1
        
        #Save data for grid plot max point data
        if len(t[outeclipse])==0:
            MaxData[MaxNum]=0
            MaxNum=MaxNum+1
        else:
            MaxValA=max(NormaDataAIn) #max value of Aperture Photometry start eclip
            MaxValP=max(NormaDataPIn) #max value of PSF Photometry start eclip
            MaxPosA=([index for index, item in enumerate(NormaDataAIn) if item == MaxValA]) #Position of max value for Aperture Photometry
            MaxPosP=([index for index, item in enumerate(NormaDataPIn) if item == MaxValP]) #Position of max value for PSF Photometry
            TimeA=T_in[MaxPosA] #Time of eclipse Aperture
            TimeP=T_in[MaxPosP] #Time of eclipse PSF
            TimeEclipse= (TimeA+TimeP)/2 #Their average
            MaxData[MaxNum]=TimeEclipse #to.list() to ensure that list in list works 
            MaxNum=MaxNum+1
        
    # Plot the prison diagram
    #%matplotlib qt
    fig = plt.figure() #%matplotlib qt #To make it pop out
    fig.set_size_inches(9, 7)
    gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
    (ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
    TopText=fig.suptitle('TLGC - V12: Normalised Δ Magnitude Sector:'+filename[143:-27], fontsize=16, weight='bold')
    TopText.set_position([.5, .95])
    LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
    LeftText.set_position([0.05, .5])
    BotText=fig.supxlabel('Time (BJD-2457000)', fontsize=14)
    BotText.set_position([.5, 0.05])
    ax1.plot(GridData[0], GridData[1],'r.', label='_nolegend_')
    ax1.plot(GridData[2], GridData[3],'b.', label='_nolegend_')
    ax1.plot(GridData[4], GridData[5],'m.', label='_nolegend_')
    ax1.plot(GridData[6], GridData[7],'c.', label='_nolegend_')
    ax1.invert_yaxis()
    Quardline = np.linspace(GridData[0][0], GridData[0][-1], 10)  #GridData[4][0] for s26
    ax1.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax1.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax1.set_xlim(GridData[0][0], GridData[0][-1]) #GridData[4][0] for s26
    ax1.set_ylim(1.05, -0.05)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[0]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax1.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax2.plot(GridData[8], GridData[9],'r.', label='_nolegend_')
    ax2.plot(GridData[10], GridData[11],'b.', label='_nolegend_')
    ax2.plot(GridData[12], GridData[13],'m.', label='_nolegend_')
    ax2.plot(GridData[14], GridData[15],'c.', label='_nolegend_')
    Quardline = np.linspace(GridData[8][0], GridData[8][-1], 10) 
    ax2.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax2.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax2.set_xlim(GridData[8][0], GridData[8][-1])
    ax2.set_ylim(1.05, -0.05)
    ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    

    Placement=round(float(MaxData[1]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax2.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax3.plot(GridData[16], GridData[17],'r.', label='_nolegend_')
    ax3.plot(GridData[18], GridData[19],'b.', label='_nolegend_')
    ax3.plot(GridData[20], GridData[21],'m.', label='_nolegend_')
    ax3.plot(GridData[22], GridData[23],'c.', label='_nolegend_')
    ax3.invert_yaxis()
    Quardline = np.linspace(GridData[16][0], GridData[16][-1], 10)   #GridData[22][-1] for S18
    ax3.plot(Quardline,np.array([0.6] * len(Quardline)),'m')  #PSF max top line
    ax3.plot(Quardline,np.array([0.66] * len(Quardline)),'c')  #Aperature max top line
    ax3.set_xlim(GridData[16][0], GridData[16][-1])  #GridData[22][-1] for S18
    ax3.set_ylim(1.05, -0.05)
    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax3.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[2]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax3.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax4.plot(GridData[24], GridData[25],'r.', label='_nolegend_')
    ax4.plot(GridData[26], GridData[27],'b.', label='_nolegend_')
    ax4.plot(GridData[28], GridData[29],'m.', label='_nolegend_')
    ax4.plot(GridData[30], GridData[31],'c.', label='_nolegend_')
    if GridData[24]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax4.set_xlim(2925, 2926)
    else: 
        Quardline = np.linspace(GridData[24][0], GridData[24][-1], 10)
        ax4.set_xlim(GridData[24][0], GridData[24][-1])
    ax4.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax4.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax4.set_ylim(1.05, -0.05)
    ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[3]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax4.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax5.plot(GridData[32], GridData[33],'r.', label='_nolegend_')
    ax5.plot(GridData[34], GridData[35],'b.', label='_nolegend_')
    ax5.plot(GridData[36], GridData[37],'m.', label='_nolegend_')
    ax5.plot(GridData[38], GridData[39],'c.', label='_nolegend_')
    ax5.invert_yaxis()
    if GridData[36]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax5.set_xlim(1855, 1856)
    else: 
        Quardline = np.linspace(GridData[32][0], GridData[32][-1], 10) #GridData[36][0] for S50
        ax5.set_xlim(GridData[32][0], GridData[32][-1]) #GridData[36][0] for S50
    ax5.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax5.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax5.set_ylim(1.05, -0.05)
    ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax5.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[4]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax5.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax6.plot(GridData[40], GridData[41],'r.', label='_nolegend_')
    ax6.plot(GridData[42], GridData[43],'b.', label='_nolegend_')
    ax6.plot(GridData[44], GridData[45],'m.', label='_nolegend_')
    ax6.plot(GridData[46], GridData[47],'c.', label='_nolegend_')
    Quardline = np.linspace(GridData[40][0], GridData[40][-1], 10) 
    ax6.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax6.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax6.set_xlim(GridData[40][0], GridData[40][-1])
    ax6.set_ylim(1.05, -0.05)
    ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[5]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax6.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax7.plot(GridData[48], GridData[49],'r.', label='_nolegend_')
    ax7.plot(GridData[50], GridData[51],'b.', label='_nolegend_')
    ax7.plot(GridData[52], GridData[53],'m.', label='_nolegend_')
    ax7.plot(GridData[54], GridData[55],'c.', label='_nolegend_')
    ax7.invert_yaxis()
    Quardline = np.linspace(GridData[48][0], GridData[48][-1], 10)  #GridData[54][-1] for S18
    ax7.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax7.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
    ax7.set_xlim(GridData[48][0], GridData[48][-1])    #GridData[54][-1] for S18
    ax7.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    ax7.set_ylim(1.05, -0.05)
    ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax7.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[6]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax7.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax8.plot(GridData[56], GridData[57],'r.', label='_nolegend_')
    ax8.plot(GridData[58], GridData[59],'b.', label='_nolegend_')
    ax8.plot(GridData[60], GridData[61],'m.', label='_nolegend_')
    ax8.plot(GridData[62], GridData[63],'c.', label='_nolegend_')
    
        
    if GridData[24]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax8.set_xlim(2928, 2929)
    else: 
        Quardline = np.linspace(GridData[56][0], GridData[56][-1], 10)
        ax8.set_xlim(GridData[56][0], GridData[56][-1])
    ax8.plot(Quardline,np.array([0.6] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax8.plot(Quardline,np.array([0.66] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line

    ax8.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    ax8.set_ylim(1.05, -0.05)
    ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[7]),3) 
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax8.legend(title=TitelDate, fontsize="x-small", loc="lower right")



    for ax in fig.get_axes():
        ax.label_outer()
#%%
TGLC_data( r"C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\TGLC\V12_TGLC\hlsp_tglc_tess_ffi_gaiaid-573937274335905152-s0019-cam3-ccd2_tess_v1_llc.fits")
#%%
#for filename in os.listdir(r'C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/TGLC/V12_TGLC'):
    #TGLC_data( r'C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/TGLC/V12_TGLC'+filename )