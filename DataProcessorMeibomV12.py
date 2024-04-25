#Meibom plot data
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
import os #to jump around dir
import pandas as pd #For supressing e+04 data type when converting to txt
#%%Start Variables for grid plot
GridNum=0
GridData = [1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]
MaxNum=0
MaxData=[1,2,3,4,5,6,7,8]
TestNum=0
#%%Get the raw data
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Meibom') 
V_t = np.loadtxt("MEIBOM_data_V12_VBand.txt")[:, 0]#Time
I_t = np.loadtxt("MEIBOM_data_V12_IBand.txt")[:, 0]
V_m = np.loadtxt("MEIBOM_data_V12_VBand.txt")[:, 1]#Magnitudes
I_m = np.loadtxt("MEIBOM_data_V12_IBand.txt")[:, 1]
#%% Treat the data
per = 6.5043035 #period 
V_ph = np.mod( V_t, per )#Finding the phase time of data points
I_ph = np.mod( I_t, per )
V_MD = np.where( (V_ph > 0.28 ) & (V_ph < 3.17) )[0]#Data for median calculation
I_MD = np.where( (I_ph > 0.28 ) & (I_ph < 3.17) )[0]
V_Mm=V_m - np.median(V_m[V_MD])#Change from median
I_Mm=I_m - np.median(I_m[I_MD])
#%%#Data groups
Times = [V_t, I_t]
PhaseTimes= [V_ph, I_ph]
Mangnitudes = [V_Mm, I_Mm]
Sectors = ['V band', 'I band']
Num =[0,1]
#%%De plottes fuldt ud.
plt.figure()
plt.plot( V_t, V_Mm, 'r.',  label = 'I-Band', markersize=2)
plt.plot( I_t, I_Mm, 'b.' , label = 'V-Band', markersize=2)
#plt.ylim( -0.4, 0.8 )
ax = plt.gca()
ax.invert_yaxis()
plt.title( 'NGC 188 - V12 - Raw Data Meibom' )
plt.xlabel( 'Time [HJD-2457000]' )
plt.ylabel( 'Magnitude' )
plt.legend(fontsize="small", loc="lower left")

#%%
def Meibom_data( i ):
    #import Define data
    t = Times[i]
    m = Mangnitudes[i]
    s = Sectors[i]
    
    #Find first eclipse  
    StartEclip=next(x[0] for x in enumerate(m) if x[1] > 0.15) #Position of first eclipse
    StartTime=t[StartEclip] #Time of start eclipse Aperture    
    
    TimeStart = StartTime-0.5 #Peak start
    TimeStop = StartTime+0.5 #Peak end
    #cut eclipse
    Startindexin = [i for i in range(len(m)) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
    
    MaxVal=max(m[Startindexin]) #max value of start eclip
    #The position of the mimimum of the first eclipse
    MaxPos=([index for index, item in enumerate(m) if item == MaxVal]) #Position of max value for Aperture Photometry
    TimeEclipse=t[MaxPos] #Time of eclipse
    
    #Finding all peaks
    inclination = 88.62 #Inclination for V12, tabel 6 Frank paper
    Period = per #Period of V12
    EndEclipse = max(Times[i])
    AllPeak=np.arange(TimeEclipse, EndEclipse, (Period/2)).tolist() #list of all peak positins in time
    #Loop through all peaks
    for AllPeak in AllPeak:
        # Cut out wanted data into groups (<>)30
        TimeStart = AllPeak-0.6 #Peak start
        TimeStop = AllPeak+0.6 #Peak end

    
        #cut eclipse 
        #ineclipse= [i for i in range(len(t)) if (t[i]>TimeStart and t[i]<TimeStop)] #The eclipse in time
        indexin = [i for i in range(len(m)) if (t[i]>TimeStart and t[i]<TimeStop)] #The index of the eclipse
        inMag = m[indexin] #The eclipse magnitude values
        inTime = t[indexin]

        global GridNum #To ensure it works as a variable
        global MaxNum #To ensure it works as a variable
        #global TestNum #To find where mistake is
        
        #TestNum=TestNum+1
        #print(TestNum)
        
        if len(inMag) == 0:
            continue
        else:
            GridData[GridNum]=inTime.tolist() #to.list() to ensure that list in list works 
            GridNum=GridNum+1
            GridData[GridNum]=inMag.tolist()
            GridNum=GridNum+1
            
            MaxVal=max(inMag) #max value of Photometry start eclip
            MaxPos=([index for index, item in enumerate(inMag) if item == MaxVal]) #Position of max value
            Time=inTime[MaxPos] #Time of eclipse 
            TimeEclipse= Time #Their average
            MaxData[MaxNum]=TimeEclipse #to.list() to ensure that list in list works 
            MaxNum=MaxNum+1
            
    # Plotot the prison diagram
    fig = plt.figure() #%matplotlib qt#To make it pop out
    fig.set_size_inches(9, 7)
    gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
    (ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
    TopText=fig.suptitle('Meibom - V12: Normalised Δ Magnitude Type:'+Sectors[i], fontsize=16, weight='bold')
    TopText.set_position([.5, .95])
    LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
    LeftText.set_position([0.05, .5])
    BotText=fig.supxlabel('Time (HJD-2457000)', fontsize=14)
    BotText.set_position([.5, 0.05])
    
    ax1.plot(GridData[0], GridData[1],'r.', label='_nolegend_')
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

    ax2.plot(GridData[2], GridData[3],'r.', label='_nolegend_')
    Quardline = np.linspace(GridData[2][0], GridData[2][-1], 10) 
    ax2.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax2.set_xlim(GridData[2][0], GridData[2][-1])
    ax2.set_ylim(0.8, -0.05)
    ax2.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    

    Placement=round(float(MaxData[1]),3) 
    Placement2 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax2.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax3.plot(GridData[4], GridData[5],'r.', label='_nolegend_')
    ax3.invert_yaxis()
    Quardline = np.linspace(GridData[4][0], GridData[4][-1], 10)   #GridData[22][-1] for S18
    ax3.plot(Quardline,np.array([0.56] * len(Quardline)),'m')  #PSF max top line
    ax3.set_xlim(GridData[4][0], GridData[4][-1])  #GridData[22][-1] for S18
    ax3.set_ylim(0.8, -0.05)
    ax3.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax3.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[2]),3) 
    Placement3 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax3.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax4.plot(GridData[6], GridData[7],'r.', label='_nolegend_')
    if GridData[7]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax4.set_xlim(2925, 2926)
    else: 
        Quardline = np.linspace(GridData[6][0], GridData[6][-1], 10)
        ax4.set_xlim(GridData[6][0], GridData[6][-1])
    ax4.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax4.set_ylim(0.8, -0.05)
    ax4.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[3]),3)
    Placement4 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax4.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax5.plot(GridData[8], GridData[9],'r.', label='_nolegend_')
    ax5.invert_yaxis()
    if GridData[16]==[]:#To make it work whne hole in data
        Quardline = np.linspace(0, 0, 10)
        ax5.set_xlim(1855, 1856)
    else: 
        Quardline = np.linspace(GridData[8][0], GridData[8][-1], 10) #GridData[36][0] for S50
        ax5.set_xlim(GridData[8][0], GridData[8][-1]) #GridData[36][0] for S50
    ax5.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax5.set_ylim(0.8, -0.05)
    ax5.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax5.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[4]),3) 
    Placement5 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax5.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax6.plot(GridData[10], GridData[11],'r.', label='_nolegend_')
    Quardline = np.linspace(GridData[10][0], GridData[10][-1], 10) 
    ax6.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax6.set_xlim(GridData[10][0], GridData[10][-1])
    ax6.set_ylim(0.8, -0.05)
    ax6.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    Placement=round(float(MaxData[5]),3)
    Placement6 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax6.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax7.plot(GridData[12], GridData[13],'r.', label='_nolegend_')
    ax7.invert_yaxis()
    Quardline = np.linspace(GridData[12][0], GridData[12][-1], 10)  #GridData[54][-1] for S18
    ax7.plot(Quardline,np.array([0.56] * len(Quardline)),'m', label='_nolegend_')  #PSF max top line
    ax7.set_xlim(GridData[12][0], GridData[12][-1])    #GridData[54][-1] for S18
    ax7.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    ax7.set_ylim(0.8, -0.05)
    ax7.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = ax7.yaxis.get_major_ticks()
    yticks[-1].set_visible(False)

    Placement=round(float(MaxData[6]),3) 
    Placement7 = Placement#PEAK PLACEMENTS
    TitelDate = 'Time of Eclipse: ' + str(Placement)
    ax7.legend(title=TitelDate, fontsize="x-small", loc="lower right")

    ax8.plot(GridData[14], GridData[15],'r.', label='_nolegend_') 

    Quardline = np.linspace(GridData[14][0], GridData[14][-1], 10)
    ax8.set_xlim(GridData[14][0], GridData[14][-1])
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
    
    #Placements
    Placements=[Placement1,Placement2,Placement3,Placement4,Placement5,Placement6,Placement7,Placement8]
    return Placements
#%%    
#Meibom_data( 0 )
Placements = Meibom_data( 1 )
#%% To txt file
All_Magnitudes = GridData[1]+GridData[3]+GridData[5]+GridData[7]+GridData[9]+GridData[11]+GridData[13]+GridData[15]
All_Time = GridData[0]+GridData[2]+GridData[4]+GridData[6]+GridData[8]+GridData[10]+GridData[12]+GridData[14]
All_Data = np.matrix([All_Magnitudes,All_Time])
All_Data=np.transpose(All_Data) #Swap rows and collums

#%%
#Make into csv file
mat = All_Data
df = pd.DataFrame(data=mat.astype(float))
df.to_csv('Meibom_all_peaks_I.csv', sep=' ', header=False, float_format='%.10f', index=False)

mat = Placements
df = pd.DataFrame(mat) #pd.DataFrame(data=mat.astype(float))
df.to_csv('Meibom_all_placements_I.csv', header=False, float_format='%.10f', index=False)
'''
#%% 
#Make into normeret magnitude csv file
NormP4 = np.where( (np.array(GridData[6]) > (Placements[3]+0.2) ))[0]
median = np.median(np.array(GridData[7])[NormP4])
mat = All_Data
df = pd.DataFrame(data=mat.astype(float))
df.to_csv('Meibom_all_peaks_I.csv', sep=' ', header=False, float_format='%.10f', index=False)
'''