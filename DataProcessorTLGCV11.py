# DataCleaner V11 TGLC
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
import os #to jumo around dir
import pandas as pd #For supressing e+0 data type when converting to txt
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
#%% The functon that does stuff to the data
def TGLC_data( filename ):
    #import get_tglc_phot
    t, p, a = get_tglc_phot( filename ) #Data extraction of time, PSF, and Aperature
    '''   
    #Check to see if correct file
    plt.figure() #This is all data
    plt.plot(t, p,'r.', label = 'PSF Photometry')
    plt.plot(t, a,'b.', label = 'Aperature Photometry')
    plt.title( 'NGC 188 - V11 - Raw Data' )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="upper left")
    '''   
    #Find reference time of primary minimum
    MaxValA=max(a) #max value of Aperture Photometry
    MaxValP=max(p) #max value of PSF Photometry
    MaxPosA=([index for index, item in enumerate(a) if item == MaxValA]) #Position of max value for Aperture Photometry
    MaxPosP=([index for index, item in enumerate(p) if item == MaxValP]) #Position of max value for PSF Photometry
    TimeA=t[MaxPosA] #Time of eclipse Aperture
    TimeP=t[MaxPosP] #Time of eclipse PSF
    TimeEclipse= (TimeA+TimeP)/2 #Their average
    # Cut out wanted data into groups (<>)
    TimeStart = TimeEclipse-0.34 #Peak start
    TimeStop = TimeEclipse+0.34 #Peak end
    TimeStartCut = TimeStart - 0.3 #Cutoff for fitting data before peak
    TimeStopCut = TimeStop + 0.3 #Cutoff for fitting data after peak
    
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
    QuardFitP = np.polyfit(t[outeclipse], outMagP, 2) #Quardratic fit for PSF 
    QuardFitA = np.polyfit(t[outeclipse], outMagA, 2) #Quardratic fit for Aperature
    
    QuardP = np.poly1d(QuardFitP)
    QuardA = np.poly1d(QuardFitA)
    
    #For nice equations in legend
    NumP1=round(QuardFitP[0],10)
    NumP2=round(QuardFitP[1],6)
    NumP3=round(QuardFitP[2],3)
    NumA1=round(QuardFitA[0],10)
    NumA2=round(QuardFitA[1],6)
    NumA3=round(QuardFitA[2],3)
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
    #Line of zeros for other stuff
    Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
   
    #Plot everything

    #Raw data
    plt.figure()
    plt.plot(t[ineclipse], inMagP,'r.', label = 'PSF Photometry') 
    plt.plot(t[ineclipse], inMagA,'b.',label = 'Aperture Photometry' )
    plt.plot(t[outeclipse], outMagP,'.', label = 'PSF Photometry')
    plt.plot(t[outeclipse], outMagA,'.',label = 'Aperture Photometry' )
    #Fit Equations
    Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
    plt.plot(Quardline, QuardP(Quardline), label= str(NumP1) + "x^2 + " + str(NumP2) + "x +" + str(NumP3))
    plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x^2 + " + str(NumA2) + "x +" + str(NumA3))
    #Make nice
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Data and Fits' )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="lower left", bbox_to_anchor=(1, 0.5)) #bbox for out of plot legend

    #Plot Only fits
    plt.figure()
    plt.plot(t[outeclipse], outMagP,'.', label = 'PSF Photometry')
    plt.plot(t[outeclipse], outMagA,'.',label = 'Aperture Photometry' )
    plt.plot(Quardline, QuardP(Quardline), label= str(NumP1) + "x^2 + " + str(NumP2) + "x +" + str(NumP3))
    plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x^2 + " + str(NumA2) + "x +" + str(NumA3))
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Fits' )
    plt.xlabel( 'Time [BJD-2457000d]' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="lower center")
    
    #Plot data used for normalization 
    plt.figure()
    plt.plot(t[indexOut], NormaDataPOut,'r.', label = 'PSF Photometry')
    plt.plot(t[indexOut], NormaDataAOut,'b.', label = 'Aperture Photometry')
    plt.plot(Quardline,np.array([0] * len(Quardline)),'k')
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Normalized data used in normalization' )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Delta Mag' )
    plt.legend(fontsize="small", loc="lower center")
 
    #Plot all normalized data
    t_out=t[indexOut]
    t_in=t[indexin]
    plt.figure()
    plt.plot(t_out, NormaDataPOut,'r.', label = 'PSF Photometry')
    plt.plot(t_out, NormaDataAOut,'b.', label = 'Aperture Photometry')
    plt.plot(t_in, NormaDataPIn,'m.', label = 'PSF Photometry')
    plt.plot(t_in, NormaDataAIn,'c.', label = 'Aperture Photometry')
    plt.plot(Quardline,np.array([0] * len(Quardline)),'k')  #Zero line
    plt.plot(Quardline,np.array([0.257] * len(Quardline)),'m')  #PSF max top line
    plt.plot(Quardline,np.array([0.345] * len(Quardline)),'c')  #Aperature max top line
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Normalized data' )
    plt.xlabel( 'Time [BJD-2457000]' )
    plt.ylabel( 'Delta Mag' )
    plt.legend(fontsize="small", loc="lower left")

    #Save data for grid plot all data
    global GridNum #To ensure it works as a variable
    GridData[GridNum]=t_out.tolist() #to.list() to ensure that list in list works 
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataPOut.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=t_out.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataAOut.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=t_in.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataPIn.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=t_in.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataAIn.tolist()
    GridNum=GridNum+1
    
    #Save data for grid plot max point data
    global MaxNum #To ensure it works as a variable
    MaxData[MaxNum]=TimeEclipse #to.list() to ensure that list in list works 
    MaxNum=MaxNum+1
#%% Test for single files
#TGLC_data( r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\TGLC\V11_TGLC\hlsp_tglc_tess_ffi_gaiaid-573941053907094144-s0020-cam3-ccd4_tess_v1_llc.fits')
#%% Do stuff to all the data files in this dir
for filename in os.listdir(r'C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/TGLC/V11_TGLC'):
    TGLC_data( r'C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/TGLC/V11_TGLC/'+filename )
#%% Plot the prison diagram
fig = plt.figure() #%matplotlib qt #To make it pop out
fig.set_size_inches(9, 7)
gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
TopText=fig.suptitle('TLGC - V11: Normalised Δ Magnitude', fontsize=16, weight='bold')
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
Quardline = np.linspace(GridData[0][0], GridData[0][-1], 10) 
ax1.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax1.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax1.set_xlim(GridData[0][0], GridData[0][-1])
ax1.set_ylim(0.6, -0.02)
ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[0]),3) 
Placement1 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 18'
ax1.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax2.plot(GridData[8], GridData[9],'r.', label='_nolegend_')
ax2.plot(GridData[10], GridData[11],'b.', label='_nolegend_')
ax2.plot(GridData[12], GridData[13],'m.', label='_nolegend_')
ax2.plot(GridData[14], GridData[15],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[8][0], GridData[8][-1], 10) 
ax2.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax2.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax2.set_xlim(GridData[8][0], GridData[8][-1])
ax2.set_ylim(0.6, -0.02)
ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[1]),3) 
Placement2 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 20'
ax2.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax3.plot(GridData[16], GridData[17],'r.', label='_nolegend_')
ax3.plot(GridData[18], GridData[19],'b.', label='_nolegend_')
ax3.plot(GridData[20], GridData[21],'m.', label='_nolegend_')
ax3.plot(GridData[22], GridData[23],'c.', label='_nolegend_')
ax3.invert_yaxis()
Quardline = np.linspace(GridData[16][0], GridData[16][-1], 10) 
ax3.plot(Quardline,np.array([0.257] * len(Quardline)),'m')  #PSF max top line
ax3.plot(Quardline,np.array([0.345] * len(Quardline)),'c')  #Aperature max top line
ax3.set_xlim(GridData[16][0], GridData[16][-1])
ax3.set_ylim(0.6, -0.02)
ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[2]),3)
Placement3 = Placement#PEAK PLACEMENTS 
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 25'
ax3.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax4.plot(GridData[24], GridData[25],'r.', label='_nolegend_')
ax4.plot(GridData[26], GridData[27],'b.', label='_nolegend_')
ax4.plot(GridData[28], GridData[29],'m.', label='_nolegend_')
ax4.plot(GridData[30], GridData[31],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[24][0], GridData[24][-1], 10) 
ax4.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax4.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax4.set_xlim(GridData[24][0], GridData[24][-1])
ax4.set_ylim(0.6, -0.02)
ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[3]),3) 
Placement4 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 26'
ax4.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax5.plot(GridData[32], GridData[33],'r.', label='_nolegend_')
ax5.plot(GridData[34], GridData[35],'b.', label='_nolegend_')
ax5.plot(GridData[36], GridData[37],'m.', label='_nolegend_')
ax5.plot(GridData[38], GridData[39],'c.', label='_nolegend_')
ax5.invert_yaxis()
Quardline = np.linspace(GridData[36][0], GridData[32][-1], 10) 
ax5.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax5.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax5.set_xlim(GridData[36][0], GridData[32][-1])
ax5.set_ylim(0.6, -0.02)
ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[4]),3) 
Placement5 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 40'
ax5.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax6.plot(GridData[40], GridData[41],'r.', label='_nolegend_')
ax6.plot(GridData[42], GridData[43],'b.', label='_nolegend_')
ax6.plot(GridData[44], GridData[45],'m.', label='_nolegend_')
ax6.plot(GridData[46], GridData[47],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[40][0], GridData[40][-1], 10) 
ax6.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax6.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax6.set_xlim(GridData[40][0], GridData[40][-1])
ax6.set_ylim(0.6, -0.02)
ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[5]),3) 
Placement6 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 52'
ax6.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax7.plot(GridData[48], GridData[49],'r.', label='_nolegend_')
ax7.plot(GridData[50], GridData[51],'b.', label='_nolegend_')
ax7.plot(GridData[52], GridData[53],'m.', label='_nolegend_')
ax7.plot(GridData[54], GridData[55],'c.', label='_nolegend_')
ax7.invert_yaxis()
Quardline = np.linspace(GridData[48][0], GridData[48][-1], 10) 
ax7.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax7.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax7.set_xlim(GridData[48][0], GridData[48][-1])
ax7.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
ax7.set_ylim(0.6, -0.02)
ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[6]),3) 
Placement7 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 53'
ax7.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax8.plot(GridData[56], GridData[57],'r.', label='_nolegend_')
ax8.plot(GridData[58], GridData[59],'b.', label='_nolegend_')
ax8.plot(GridData[60], GridData[61],'m.', label='_nolegend_')
ax8.plot(GridData[62], GridData[63],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[56][0], GridData[56][-1], 10) 
ax8.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
ax8.plot(Quardline,np.array([0.345] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax8.set_xlim(GridData[56][0], GridData[56][-1])
ax8.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
ax8.set_ylim(0.6, -0.02)
ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[7]),3)
Placement8 = Placement#PEAK PLACEMENTS 
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 59'
ax8.legend(title=TitelDate, fontsize="x-small", loc="lower right")



for ax in fig.get_axes():
    ax.label_outer()
    
#%%
#Datafiles i added togehter
l=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
#j=[0,1,2,3]
All_Magnitudes_PSF=[]
All_Magnitudes_Aperature=[]
All_Time_Aperature=[]
All_Time_PSF=[]
for l in l:
    if l == 0 or l == 1 or l == 2 or l == 3:
        All_Magnitudes_PSF=GridData[1]
        All_Magnitudes_Aperature=GridData[3]
        All_Time_Aperature=GridData[2]
        All_Time_PSF=GridData[0]
        continue
    else:
        if (isinstance(GridData[l], int)) == True:
            continue
        else: 
            if l % 4 == 0:
                All_Time_PSF=All_Time_PSF+GridData[l]
                continue
            else:
                if l % 2 == 0:
                    All_Time_Aperature=All_Time_Aperature+GridData[l]
                    continue
                else:
                    if (l-1) % 4 == 0:
                        All_Magnitudes_PSF=All_Magnitudes_PSF+GridData[l]
                        continue
                    else:
                        All_Magnitudes_Aperature=All_Magnitudes_Aperature+GridData[l]
                        continue

All_Data = np.matrix([All_Magnitudes_PSF,All_Time_PSF,All_Magnitudes_Aperature ,All_Time_Aperature])
All_Data=np.transpose(All_Data) #Swap rows and collums
#Placement data
Placements=[Placement1,Placement2,Placement3,Placement4,Placement5,Placement6,Placement7,Placement8]
#Make into csv file
mat = All_Data
df = pd.DataFrame(data=mat.astype(float))
df.to_csv('TLGC_V11_all_peaks_PSF_Ap_S'+filename[143:-27]+'.csv', sep=' ', header=False, float_format='%.10f', index=False)

mat = Placements
df = pd.DataFrame(mat) #pd.DataFrame(data=mat.astype(float))
df.to_csv('TLGC_V11_all_placements_S'+filename[143:-27]+'.csv', header=False, float_format='%.10f', index=False)
   

