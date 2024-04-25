# DataCleaner V12 TGLC
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
#%%#Getting the data 
#"C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Nardiello\V11_Nardiello\pathos_tess_lightcurve_tic-00461601177-s018_tess_v1_llc.fits"
#Raw data
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Nardiello\V11_Nardiello') 
a18 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s018_tess_v1_llc.fits' )
a20 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s020_tess_v1_llc.fits' )
a25 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s025_tess_v1_llc.fits' )
a26 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s026_tess_v1_llc.fits' )
a40 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s040_tess_v1_llc.fits' )
a52 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s052_tess_v1_llc.fits' )
a53 = fits.getdata( 'pathos_tess_lightcurve_tic-00461601177-s053_tess_v1_llc.fits' )

per = 35.178 #Period

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
#Finding the time of data points
ph18 = np.mod( a18['time'], per )
ph20 = np.mod( a20['time'], per )
ph25 = np.mod( a25['time'], per )
ph26 = np.mod( a26['time'], per )
ph40 = np.mod( a40['time'], per )
ph52 = np.mod( a52['time'], per )
ph53 = np.mod( a53['time'], per )
#cutting data
i18 = np.where( (ph18 > 14.7 ) & (ph18 < 15.2) )[0]
i20 = np.where( (ph20 > 14.7 ) & (ph20 < 15.2) )[0]
i25 = np.where( (ph25 > 14.7 ) & (ph25 < 15.2) )[0]
i26 = np.where( (ph26 > 14.7 ) & (ph26 < 15.2) )[0]
i40 = np.where( (ph40 > 14.7 ) & (ph40 < 15.2) )[0]
i52 = np.where( (ph52 > 14.7 ) & (ph52 < 15.2) )[0]
i53 = np.where( (ph53 > 14.7 ) & (ph53 < 15.2) )[0]
#Change from median
me18=m18 - np.median(m18[i18])
me20=m20 - np.median(m20[i20])
me25=m25 - np.median(m25[i25])
me26=m26 - np.median(m26[i26])
me40=m40 - np.median(m40[i40])
me52=m52 - np.median(m52[i52])
me53=m53 - np.median(m53[i53])
#Cut sector 18, 20, 53 to make them work 
'''
c20 = np.where( (ph20 > 14.0 ) & (ph20 < 15.5) )[0] #For future, remove that one point instead
c53 = np.where( (ph53 > 14.0 ) & (ph53 < 15.0) )[0] #For future, remove that one point instead
ph20 = ph20[c20]
ph53 = ph53[c53]
me20=me20[c20]
me53=me53[c53]
'''
#Sector 18 still needs cutting
c18 = np.where( (ph18 > 13.0 ) & (ph18 < 16.0) )[0]
ph18 = ph18[c18]
me18=me18[c18]
#Data groups
Times = [ph18, ph20, ph25, ph26, ph40, ph52, ph53]
Mangnitudes = [me18, me20, me25, me26, me40, me52, me53]
#BJD
ph18 = a18['time'][c18]
ph20 = a20['time']
ph25 = a25['time']
ph26 = a26['time']
ph40 = a40['time']
ph52 = a52['time']
ph53 = a53['time']
Times = [ph18, ph20, ph25, ph26, ph40, ph52, ph53]
#%%Clean data of very high values
List=[0,1,2,3,4,5,6]
for i in List:
    RollM = np.roll(Mangnitudes[i],1)
    DiffM = Mangnitudes[i]-RollM #Over 0.15 er vist dårligt
    pos = np.where( (0.05 > abs(DiffM) ))[0]
    Mangnitudes[i]=Mangnitudes[i][pos]
    Times[i] = Times[i][pos]

#%%New data groups
Sectors = ['18','20','25','26','40','52','53']
Num =[0,1,2,3,4,5,6]


#%% PLot all data
plt.plot( ph18, me18, '.',  label = 'Sector 18' )
plt.plot( ph20, me20, '.' , label = 'Sector 20')
plt.plot( ph25, me25, '.' , label = 'Sector 25')
plt.plot( ph26, me26, '.' , label = 'Sector 26')
plt.plot( ph40, me40, '.' , label = 'Sector 40')
plt.plot( ph52, me52, '.' , label = 'Sector 52')
plt.plot( ph53, me53, '.' , label = 'Sector 53')
plt.ylim( -0.2, 0.6 )
ax = plt.gca() 
ax.invert_yaxis()
plt.title( 'NGC 188 - V11 - Raw Data' )
plt.xlabel( 'Phase time (Days)' )
plt.ylabel( 'Δ Magnitude compared to median' )
plt.legend(fontsize="small", loc="lower right")
#%% Plot single Sector function
def Raw_Plot(i):
    plt.figure()
    plt.plot( Times[i], Mangnitudes[i], '.' , label = 'Sector' )
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Raw Data Sector:'+Sectors[i] )
    plt.xlabel( 'Phase time (Days)' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="lower left")
#%% Plot all sectors raw
for i in Num:
    Raw_Plot(i)
    
#%%Clean data of very high values
for i in Num:
    RollM = np.roll(Mangnitudes[i],1)
    DiffM = Mangnitudes[i]-RollM #Over 0.15 er vist dårligt
    pos = np.where( (0.18 > abs(DiffM) ))[0]
    Mangnitudes[i]=Mangnitudes[i][pos]
    Times[i] = Times[i][pos]
#%% Isolating and analysing the peaks
def Nardiello_data( i ):
    #import Define data
    t = Times[i]
    m = Mangnitudes[i]
    s = Sectors[i]
    #Find reference time of primary minimum
    MaxValA=max(m) #max value of Aperture Photometry
    MaxPosA=([index for index, item in enumerate(m) if item == MaxValA]) #Position of max value for Aperture Photometry
    TimeA=t[MaxPosA] #Time of eclipse Aperture
    TimeEclipse= TimeA #Their average
    # Cut out wanted data into groups (<>)
    TimeStart = TimeEclipse-0.34 #Peak start
    TimeStop = TimeEclipse+0.34 #Peak end
    TimeStartCut = TimeStart - 0.43 #Cutoff for fitting data before peak
    TimeStopCut = TimeStop + 0.43 #Cutoff for fitting data after peak
    
    #cut eclipse 
    ineclipse= [i for i in range(len(t)) if (t[i]>TimeStart and t[i]<TimeStop)] #The eclipse in time
    indexin = [i for i in range(len(m)) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
    inMagA = m[indexin] #The eclipse magnitude values - Aperture Photometry
    
    #Cut fitting data not part of eclipse
    outeclipseBefore= [i for i in range(len(t)) if (t[i]>TimeStartCut and t[i]<TimeStart)] #Before eclipse
    indexBefore = [i for i in range(len(m)) if (t[i]>TimeStartCut and t[i]<TimeStart)] #The index of the eclipse
    outeclipseAfter= [i for i in range(len(t)) if (t[i]>TimeStop and t[i]<TimeStopCut)] #after eclipse
    indexAfter = [i for i in range(len(m)) if (t[i]>TimeStop and t[i]<TimeStopCut)] #The index of the eclipse
    indexOut = indexBefore+indexAfter #Total Out index
    outeclipse = outeclipseBefore+outeclipseAfter #Total out times
    outMagA = m[indexOut] #non eclipse magnitude values - magnitude Photometry
    
    #Fit qudratic equation to outeclipse data
    QuardFitA = np.polyfit(t[outeclipse], outMagA, 2) #Quardratic fit for magnitude
    
    QuardA = np.poly1d(QuardFitA)
    
    #For nice equations in legend
    NumA1=round(QuardFitA[0],10)
    NumA2=round(QuardFitA[1],6)
    NumA3=round(QuardFitA[2],3)
    #Normalzing the data
    NormaAOut=QuardA(t[indexOut]) #Difference between zero and fit i outside
    NormaAIN=QuardA(t[indexin]) #Difference between zero and fit i inside
    #Normalized data
    NormaDataAOut= outMagA-NormaAOut #Normalized i inside 
    NormaDataAIn= inMagA-NormaAIN #Normalized i inside 
    #Line of zeros for other stuff
    Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
    
    #Plot everything
    
    #Raw data
    plt.figure()
    plt.plot(t[ineclipse], inMagA,'b.',label = 'Photometry' )
    plt.plot(t[outeclipse], outMagA,'.',label = 'Photometry' )
    #Fit Equations
    Quardline = np.linspace(t[outeclipse][0], t[outeclipse][-1], 1000) 
    plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x^2 + " + str(NumA2) + "x +" + str(NumA3))
    #Make nice
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Data and Fits Sector: '+Sectors[i] )
    plt.xlabel( 'Phase Time (Days)' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="lower left", bbox_to_anchor=(1, 0.5)) #bbox for out of plot legend
    
    #Plot Only fits
    plt.figure()
    plt.plot(t[outeclipse], outMagA,'.',label = 'Photometry' )
    plt.plot(Quardline, QuardA(Quardline), label = str(NumA1) + "x^2 + " + str(NumA2) + "x +" + str(NumA3))
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Fits Sector: '+Sectors[i]  )
    plt.xlabel( 'Phase Time (Days)' )
    plt.ylabel( 'Magnitude' )
    plt.legend(fontsize="small", loc="lower center")
    
    #Plot data used for normalization 
    plt.figure()
    plt.plot(t[indexOut], NormaDataAOut,'b.', label = 'Photometry')
    plt.plot(Quardline,np.array([0] * len(Quardline)),'k')
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Normalized data used in normalization Sector: '+Sectors[i] )
    plt.xlabel( 'Phase Time (Days)' )
    plt.ylabel( 'Delta Mag' )
    plt.legend(fontsize="small", loc="lower center")
    
    #Plot all normalized data
    t_out=t[indexOut]
    t_in=t[indexin]
    plt.figure()
    plt.plot(t_out, NormaDataAOut,'b.', label = 'Photometry')
    plt.plot(t_in, NormaDataAIn,'c.', label = 'Photometry')
    plt.plot(Quardline,np.array([0] * len(Quardline)),'k')  #Zero line
    plt.plot(Quardline,np.array([0.257] * len(Quardline)),'m')  #PSF max top line
    plt.plot(Quardline,np.array([0.345] * len(Quardline)),'c')  #magnitude max top line
    ax = plt.gca() 
    ax.invert_yaxis()
    plt.title( 'NGC 188 - V11 - Normalized data Sector: '+Sectors[i] )
    plt.xlabel( 'Phase Time (Days)' )
    plt.ylabel( 'Delta Mag' )
    plt.legend(fontsize="small", loc="lower left")
    
    #Save data for grid plot all data
    global GridNum #To ensure it works as a variable
    GridData[GridNum]=t_out.tolist() #to.list() to ensure that list in list works 
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataAOut.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=t_in.tolist()
    GridNum=GridNum+1
    GridData[GridNum]=NormaDataAIn.tolist()
    GridNum=GridNum+1
        
    #Save data for grid plot max point data
    global MaxNum #To ensure it works as a variable
    MaxData[MaxNum]=TimeEclipse #to.list() to ensure that list in list works 
    MaxNum=MaxNum+1
#%% Plot all sectors
for i in Num:
    Nardiello_data(i)
#%%
#%% Plot the prison diagram
#%matplotlib qt 
fig = plt.figure() #%matplotlib qt #To make it pop out
fig.set_size_inches(9, 7)
gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
TopText=fig.suptitle('Nardiello - V11: Normalised Δ Magnitude', fontsize=16, weight='bold')
TopText.set_position([.5, .95])
LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
LeftText.set_position([0.05, .5])
BotText=fig.supxlabel('Time', fontsize=14)
BotText.set_position([.5, 0.05])

ax1.plot(GridData[0], GridData[1],'b.', label='_nolegend_')
ax1.plot(GridData[2], GridData[3],'c.', label='_nolegend_')
ax1.invert_yaxis()
Quardline = np.linspace(GridData[0][0], GridData[0][-1], 10) 
ax1.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax1.set_xlim(GridData[0][0], GridData[0][-1])
ax1.set_ylim(0.45, -0.05)
ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[0]),3)
Placement1 = Placement#PEAK PLACEMENTS
Placement1=Placement #For saving the peal positions
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 18'
ax1.legend(title=TitelDate, fontsize="x-small", loc="lower right")


ax2.plot(GridData[4], GridData[5],'b.', label='_nolegend_')
ax2.plot(GridData[6], GridData[7],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[4][0], GridData[4][-1], 10) 
ax2.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax2.set_xlim(GridData[4][0], GridData[4][-1])
ax2.set_ylim(0.45, -0.05)
ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[1]),3) 
Placement2 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 20'
ax2.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax3.plot(GridData[8], GridData[9],'b.', label='_nolegend_')
ax3.plot(GridData[10], GridData[11],'c.', label='_nolegend_')
ax3.invert_yaxis()
Quardline = np.linspace(GridData[8][0], GridData[8][-1], 10) 
ax3.plot(Quardline,np.array([0.30] * len(Quardline)),'c')  #Aperature max top line
ax3.set_xlim(GridData[8][0], GridData[8][-1])
ax3.set_ylim(0.45, -0.05)
ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[2]),3) 
Placement3 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 25'
ax3.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax4.plot(GridData[12], GridData[13],'b.', label='_nolegend_')
ax4.plot(GridData[14], GridData[15],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[12][0], GridData[12][-1], 10) 
ax4.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax4.set_xlim(GridData[12][0], GridData[12][-1])
ax4.set_ylim(0.45, -0.05)
ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[3]),3) 
Placement4 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 26'
ax4.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax5.plot(GridData[16], GridData[17],'b.', label='_nolegend_')
ax5.plot(GridData[18], GridData[19],'c.', label='_nolegend_')
ax5.invert_yaxis()
Quardline = np.linspace(GridData[16][0], GridData[16][-1], 10) 
ax5.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax5.set_xlim(GridData[16][0], GridData[16][-1])
ax5.set_ylim(0.45, -0.05)
ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[4]),3) 
Placement5 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 40'
ax5.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax6.plot(GridData[20], GridData[21],'b.', label='_nolegend_')
ax6.plot(GridData[22], GridData[23],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[20][0], GridData[20][-1], 10) 
ax6.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax6.set_xlim(GridData[20][0], GridData[20][-1])
ax6.set_ylim(0.45, -0.05)
ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[5]),3) 
Placement6 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 52'
ax6.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax7.plot(GridData[24], GridData[25],'b.', label='_nolegend_')
ax7.plot(GridData[26], GridData[27],'c.', label='_nolegend_')
ax7.invert_yaxis()
Quardline = np.linspace(GridData[24][0], GridData[24][-1], 10) 
ax7.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax7.set_xlim(GridData[24][0], GridData[24][-1])
ax7.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
ax7.set_ylim(0.45, -0.05)
ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))

Placement=round(float(MaxData[6]),3) 
Placement7 = Placement#PEAK PLACEMENTS
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 53'
ax7.legend(title=TitelDate, fontsize="x-small", loc="lower right")

ax8.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
ax8.set_xlim(0, 0)
ax8.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
ax8.set_ylim(0.45, -0.05)
ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))

for ax in fig.get_axes():
    ax.label_outer()

#Husk at data kan være i fase dage ikke BJD  
 
#Datafiles i added togehter
l=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
#j=[0,1,2,3

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
Placements=[Placement1,Placement2,Placement3,Placement4,Placement5,Placement6,Placement7]
#Make into csv file
mat = All_Data
df = pd.DataFrame(data=mat.astype(float))
df.to_csv('Nardiello_V11_all_peaks_.csv', sep=' ', header=False, float_format='%.10f', index=False)

mat = Placements
df = pd.DataFrame(mat) #pd.DataFrame(data=mat.astype(float))
df.to_csv('Nardiello_V11_all_placements.csv', header=False, float_format='%.10f', index=False)

    
    
    
    
    
'''
ax8.plot(GridData[56], GridData[57],'r.', label='_nolegend_')
ax8.plot(GridData[58], GridData[59],'b.', label='_nolegend_')
ax8.plot(GridData[60], GridData[61],'m.', label='_nolegend_')
ax8.plot(GridData[62], GridData[63],'c.', label='_nolegend_')
Quardline = np.linspace(GridData[56][0], GridData[56][-1], 10) 
ax8.plot(Quardline,np.array([0.257] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line


Placement=round(float(MaxData[7]),3) 
TitelDate = 'Time of Eclipse: ' + str(Placement) + '     Sector 59'
ax8.legend(title=TitelDate, fontsize="x-small", loc="lower right")
'''

