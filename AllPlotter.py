# PLot everything from trested data
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
import os #to jumo around dir
from pandas import * #to read csv files
#%% Datafiles used
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\All_Peaks') 
#TLGC
DataTLGC=['TLGC_V11_all_peaks_PSF_Ap.csv','TLGC_V11_all_placements.csv','TLGC_V12_all_peaks_PSF_Ap_S59.csv','TLGC_V12_all_placements_S59.csv','TLGC_V12_all_peaks_PSF_Ap_S53.csv','TLGC_V12_all_placements_S53.csv','TLGC_V12_all_peaks_PSF_Ap_S52.csv','TLGC_V12_all_placements_S52.csv','TLGC_V12_all_peaks_PSF_Ap_S40.csv','TLGC_V12_all_placements_S40.csv','TLGC_V12_all_peaks_PSF_Ap_S26.csv','TLGC_V12_all_placements_S26.csv','TLGC_V12_all_peaks_PSF_Ap_S25.csv','TLGC_V12_all_placements_S25.csv','TLGC_V12_all_peaks_PSF_Ap_S18.csv','TLGC_V12_all_placements_S18.csv','TLGC_V12_all_peaks_PSF_Ap_S20.csv','TLGC_V12_all_placements_S20.csv']
NameTLGC=['TLCG_V11_Data','TLCG_V11_Peaks','TLCG_V12_S59_Data','TLCG_V12_S59_Peaks','TLCG_V12_S53_Data','TLCG_V12_S53_Peaks','TLCG_V12_S52_Data','TLCG_V12_S52_Peaks','TLCG_V12_S40_Data','TLCG_V12_S40_Peaks','TLCG_V12_S26_Data','TLCG_V12_S26_Peaks','TLCG_V12_S25_Data','TLCG_V12_S25_Peaks','TLCG_V12_S18_Data','TLCG_V12_S18_Peaks','TLCG_V12_S20_Data','TLCG_V12_S20_Peaks']

#Nardiello
DataNardiello=['Nardiello_V11_all_peaks_.csv','Nardiello_V11_all_placements.csv','Nardiello_V12_all_peaks_S20.csv','Nardiello_V12_all_placements_S20.csv','Nardiello_V12_all_peaks_S53.csv','Nardiello_V12_all_placements_S53.csv','Nardiello_V12_all_peaks_S52.csv','Nardiello_V12_all_placements_S52.csv','Nardiello_V12_all_peaks_S40.csv','Nardiello_V12_all_placements_S40.csv','Nardiello_V12_all_peaks_S26.csv','Nardiello_V12_all_placements_S26.csv','Nardiello_V12_all_peaks_S25.csv','Nardiello_V12_all_placements_S25.csv','Nardiello_V12_all_peaks_S18.csv','Nardiello_V12_all_placements_S18.csv']
NameNardiello=['Nardiello_V11_Data','Nardiello_V11_Peaks','Nardiello_V12_S20_Data','Nardiello_V12_S20_Peaks','Nardiello_V12_S53_Data','Nardiello_V12_S53_Peaks','Nardiello_V12_S52_Data','Nardiello_V12_S52_Peaks','Nardiello_V12_S40_Data','Nardiello_V12_S40_Peaks','Nardiello_V12_S26_Data','Nardiello_V12_S26_Peaks','Nardiello_V12_S25_Data','Nardiello_V12_S25_Peaks','Nardiello_V12_S18_Data','Nardiello_V12_S18_Peaks']

CompareNardielloData=['Nard_RAW_V12_all_peaks_S40.csv','Nard_AP1_V12_all_peaks_S40.csv','Nard_PSF_V12_all_peaks_S40.csv']
CompareNardielloName=['Nard_V12_RAW_S40','Nard_V12_AP1_S40','Nard_V12_PSF_S40']
#Meibom
DataMeibom=['Meibom_all_peaks_V.csv', 'Meibom_all_placements_V.csv','Meibom_all_peaks_I.csv', 'Meibom_all_placements_I.csv']

#Earthbased V11
DataEarthBased=['V11_Earthbased_V.csv', 'V11_Earthbased_all_placements_V.csv','V11_Earthbased_I.csv', 'V11_Earthbased_all_placements_I.csv']

#Depth
DepthData=['Nardiello_V12_Depth_S18.csv', 'Nardiello_V12_Depth_S20.csv','Nardiello_V12_Depth_S25.csv', 'Nardiello_V12_Depth_S26.csv', 'Nardiello_V12_Depth_S40.csv', 'Nardiello_V12_Depth_S52.csv', 'Nardiello_V12_Depth_S53.csv']
DepthName=['Nard_V12_S18_Depth','Nard_V12_S20_Depth', 'Nard_V12_S25_Depth', 'Nard_V12_S26_Depth', 'Nard_V12_S40_Depth', 'Nard_V12_S52_Depth', 'Nard_V12_S53_Depth']
#%% Getting TLGC data in dictonary
#Constants
i = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
j = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#Run through
AllData = {}
for i in i:
    if DataTLGC[i][14] == 'e':#filename[143:-27]
        AA = np.loadtxt(DataTLGC[i])
        AllData[NameTLGC[i] + "_PSF"] = AA[:,0]
        AllData[NameTLGC[i] + "_PSF_Time"] = AA[:,1]
        AllData[NameTLGC[i] + "_Ap"] = AA[:,2]
        AllData[NameTLGC[i] + "_Ap_Time"] = AA[:,3]
            
    if DataTLGC[i][14] == 'l':#filename[143:-27]
        AA = np.loadtxt(DataTLGC[i])
        AllData[NameTLGC[i]] = AA
#old string version
'''
AllData = {}
for i in i:
    if DataTLGC[i][14] == 'e':#filename[143:-27]
        with open(DataTLGC[i]) as f:
            AllData[NameTLGC[i] + "_PSF"] = [row.split()[0] for row in f]
        with open(DataTLGC[i]) as f:
            AllData[NameTLGC[i] + "_PSF_Time"] = [row.split()[1] for row in f]
        with open(DataTLGC[i]) as f:
            AllData[NameTLGC[i] + "_Ap"] = [row.split()[2] for row in f]
        with open(DataTLGC[i]) as f:
            AllData[NameTLGC[i] + "_Ap_Time"] = [row.split()[3] for row in f]
            
    if DataTLGC[i][14] == 'l':#filename[143:-27]
        with open(DataTLGC[i]) as f:
            AllData[NameTLGC[i]] = [row.split()[0] for row in f]
#dict_read = dict(map(float,x) for x in reader)
'''
# Getting Nardiello data in dictonary
for j in j:
    if DataNardiello[j][19] == 'e':#filename[143:-27]
        AA = np.loadtxt(DataNardiello[j])
        AllData[NameNardiello[j] + "_Photo"] = AA[:,1]
        AllData[NameNardiello[j] + "_Time"] = AA[:,0]
            
    if DataNardiello[j][19] == 'l':#filename[143:-27]
        AA = np.loadtxt(DataNardiello[j])
        AllData[NameNardiello[j]] = AA
#Comp
AllData[CompareNardielloName[0]] = np.loadtxt(CompareNardielloData[0])
AllData[CompareNardielloName[1]] = np.loadtxt(CompareNardielloData[1])
AllData[CompareNardielloName[2]] = np.loadtxt(CompareNardielloData[2])
# Getting Jordbasert V11 data in dictonary
AA = np.loadtxt(DataEarthBased[0])
AllData["Earthbased_V_Photo"] = AA[:,0]
AllData["Earthbased_V_Time"] = AA[:,1]
AA = np.loadtxt(DataEarthBased[1])
AllData["Earthbased_V_Peaks"] = AA
AA = np.loadtxt(DataEarthBased[2])
AllData["Earthbased_I_Photo"] = AA[:,0]
AllData["Earthbased_I_Time"] = AA[:,1]
AA = np.loadtxt(DataEarthBased[3])
AllData["Earthbased_I_Peaks"] = AA
# Getting Meibom data in dictonary
AA = np.loadtxt(DataMeibom[0])
AllData["Meibom_V_Photo"] = AA[:,0]
AllData["Meibom_V_Time"] = AA[:,1]
AA = np.loadtxt(DataMeibom[1])
AllData["Meibom_V_Peaks"] = AA
AA = np.loadtxt(DataMeibom[2])
AllData["Meibom_I_Photo"] = AA[:,0]
AllData["Meibom_I_Time"] = AA[:,1]
AA = np.loadtxt(DataMeibom[3])
AllData["Meibom_I_Peaks"] = AA
#Getting depth in the database
y=[0,1,2,3,4,5,6]
for y in y:
    AA = np.loadtxt(DepthData[y])
    AllData[DepthName[y]] = AA[:,0]
    AllData[DepthName[y] + "_Time"] = AA[:,1]
#Real peaks are found for meibom due to hole in data and Primary and secundary are claculated
for h in [0,1,2,3,4,5,6,7]:
    P=6.504/2 # half period of v12
    Number=round((AllData['Meibom_V_Peaks'][2]-AllData['Meibom_V_Peaks'][h])/P)
    if Number < 0:
        AllData["Meibom_V_Peaks"][h]=AllData['Meibom_V_Peaks'][2]+(abs(Number)*P)
        AllData["Meibom_I_Peaks"][h]=AllData['Meibom_V_Peaks'][2]+(abs(Number)*P)
    else:
        AllData["Meibom_V_Peaks"][h]=AllData['Meibom_V_Peaks'][2]-(abs(Number)*P)
        AllData["Meibom_I_Peaks"][h]=AllData['Meibom_V_Peaks'][2]-(abs(Number)*P)

def determine_closeness(number, PrimaryRef):
    rounded_number = round((number-PrimaryRef)/6.504)
    difference = abs(((number-PrimaryRef)/6.504) - rounded_number)
    if difference < 0.25:
        return "P - "
    elif difference > 0.75:
        return "S - "
    else:
        return "S - "
#%% PLot Meibom - V12 Jordbaseret
#Start values
fig = plt.figure() #%matplotlib qt#To make it pop out
fig.set_size_inches(9, 7)
gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
TopText=fig.suptitle('NGC 188 - V12: Δ Magnitude, V-band (Meibom et al.)', fontsize=16, weight='bold')
TopText.set_position([.5, .95])
LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
LeftText.set_position([0.05, .5])
BotText=fig.supxlabel('Time before/after eclipse (days)', fontsize=14) #Time (HJD-2457000)
BotText.set_position([.5, 0.05])
#Cutting the relevant datapoints
k =[0,1,2,3,4,5,6,7]
L=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
Data='Meibom_V_Photo'
Time='Meibom_V_Time'
Peaks='Meibom_V_Peaks'

for k in k:
    
    P= 0.2#Peak area plottede
    I1 = np.where( (AllData[Time] > AllData[Peaks][k]-P ) & (AllData[Time] < AllData[Peaks][k]+P) )[0]
    P1 = AllData[Data][I1]
    T1 = AllData[Time][I1]
    #Plotting
    L[k].plot(T1-AllData[Peaks][k], P1,'r.', label='_nolegend_')
    L[k].invert_yaxis()
    Quardline = np.linspace(-P, P, 10)
    L[k].plot(Quardline,np.array([0.595] * len(Quardline)),'k',linestyle='--', label='_nolegend_') 
    L[k].plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
    L[k].set_xlim(-P,P) 
    L[k].set_ylim(0.7, -0.05)
    L[k].yaxis.set_major_locator(MaxNLocator(prune='lower'))
    yticks = L[k].yaxis.get_major_ticks()
    yticks[-1].set_visible(False)
    xticks = L[k].xaxis.get_major_ticks()
    xticks[-1].set_visible(False)
    
    PoS=determine_closeness(AllData[Peaks][k])
    TitelDate = PoS + 'ToE: ' + str("%.2f" % AllData[Peaks][k])
    L[k].legend(title=TitelDate, fontsize="x-small", loc="lower right")
    
    #L[k].set_xticklabels([]) #AllData[Peaks][k]
    
for ax in fig.get_axes(): #o make grid plot stuff work
    ax.label_outer()
#%% V11 Jordbasertet
plt.figure()#V band
plt.plot( AllData['Earthbased_V_Time']-AllData['Earthbased_V_Peaks'], AllData['Earthbased_V_Photo'], 'r.',  label = None, markersize=5)
Quardline = np.linspace(-0.5, 0.5, 10)
plt.plot(Quardline,np.array([0.357] * len(Quardline)),'k',linestyle='--', label='_nolegend_') 
plt.plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
plt.ylim([0.53, -0.05])
plt.xlim([-0.5, 0.5])
ax = plt.gca()
plt.title( 'NGC 188 - V11: Δ Magnitude, V-band (Zhuo et al.)' )
plt.xlabel( 'Time before/after eclipse (days)' )
plt.ylabel( 'Δ Magnitude' )
plt.legend(fontsize="small", loc="lower left")
TitelDate = 'ToE: ' + str("%.2f" % AllData['Earthbased_V_Peaks'])
plt.legend(title=TitelDate, fontsize="small", loc="lower right")

plt.figure()#I band
plt.plot( AllData['Earthbased_I_Time']-AllData['Earthbased_I_Peaks'], AllData['Earthbased_I_Photo'], 'r.',  label = None, markersize=5)
Quardline = np.linspace(-0.5, 0.5, 10)
plt.plot(Quardline,np.array([0.357] * len(Quardline)),'k',linestyle='--', label='_nolegend_') 
plt.plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
plt.ylim([0.53, -0.05])
plt.xlim([-0.5, 0.5])
ax = plt.gca()
plt.title( 'NGC 188 - V11: Δ Magnitude, I-band (Zhuo et al.)' )
plt.xlabel( 'Time before/after eclipse (days)' )
plt.ylabel( 'Δ Magnitude' )
plt.legend(fontsize="small", loc="lower left")
TitelDate = 'ToE: ' + str("%.2f" % AllData['Earthbased_I_Peaks'])
plt.legend(title=TitelDate, fontsize="small", loc="lower right")

 #%%PLot one TLGC
#Defining fig
fig = plt.figure() #%matplotlib qt#To make it pop out
fig.set_size_inches(9, 7)
gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
TopText=fig.suptitle('NGC 188 - V12: Δ Magnitude, S59, TGLC, Aperture Photometry', fontsize=16, weight='bold')
#P=35.178/2 # half period of v11
#TopText=fig.suptitle('NGC 188 - V12: Δ Magnitude, TGLC, PSF Photometry', fontsize=16, weight='bold')
TopText.set_position([.5, .95])
LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
LeftText.set_position([0.05, .5])
BotText=fig.supxlabel('Time before/after eclipse (days)', fontsize=14)
BotText.set_position([.5, 0.05])

#Start values
k =[0,1,2,3,4,5,6,7]
L=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
Data='TLCG_V12_S59_Data_Ap'
Time='TLCG_V12_S59_Data_Ap_Time'
Peaks='TLCG_V12_S59_Peaks'

for k in k:
    
    P= 0.2#Peak area plottede
    #P= 0.5#Peak area plottede
    I1 = np.where( (AllData[Time] > AllData[Peaks][k]-P ) & (AllData[Time] < AllData[Peaks][k]+P) )[0]
    P1 = AllData[Data][I1]
    T1 = AllData[Time][I1]
    if AllData[Peaks][k] > 10: #To print no data for no no data
        #Plotting
        L[k].plot(T1-AllData[Peaks][k], P1,'r.', label='_nolegend_')
        L[k].invert_yaxis()
        Quardline = np.linspace(-P, P, 10)
        L[k].plot(Quardline,np.array([0.595] * len(Quardline)),'k',linestyle='--', label='_nolegend_') #V12
        #L[k].plot(Quardline,np.array([0.357] * len(Quardline)),'k',linestyle='--', label='_nolegend_') #V11
        L[k].plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
        L[k].set_xlim(-P, P) 
        L[k].set_ylim(0.8, -0.05) #V12
        #L[k].set_ylim(0.44, -0.05)#V11
        L[k].yaxis.set_major_locator(MaxNLocator(prune='lower'))
        yticks = L[k].yaxis.get_major_ticks()
        yticks[-1].set_visible(False)
        PoS=determine_closeness(AllData[Peaks][k], 1808.76865399)#For V12
        TitelDate = PoS + 'ToE: ' + str("%.2f" % AllData[Peaks][k])
        L[k].legend(title=TitelDate, fontsize="x-small", loc="lower right")
        
        xticks = L[k].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
    else:
        L[k].set_xlim(-P, P)
        TitelDate = 'No Data'
        L[k].legend(title=TitelDate, fontsize="x-small", loc="lower right")
        xticks = L[k].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
        continue      
    
for ax in fig.get_axes(): #o make grid plot stuff work
    ax.label_outer()
    
#%%PLot one Nardiello
#Defining fig
fig = plt.figure() #%matplotlib qt#To make it pop out
fig.set_size_inches(9, 7)
gs = fig.add_gridspec(4, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) = gs.subplots( sharey='row')
TopText=fig.suptitle('NGC 188 - V12: Δ Magnitude, S25, Nardiello', fontsize=16, weight='bold')
#TopText=fig.suptitle('NGC 188 - V11: Δ Magnitude, Nardiello', fontsize=16, weight='bold')
#P=35.178/2 # half period of v11
TopText.set_position([.5, .95])
LeftText=fig.supylabel('Δ Magnitude', fontsize=14)
LeftText.set_position([0.05, .5])
BotText=fig.supxlabel('Time before/after eclipse (days)', fontsize=14)
BotText.set_position([.5, 0.05])
#Start values
Data='Nardiello_V12_S25_Data_Photo'
Time='Nardiello_V12_S25_Data_Time'
Peaks='Nardiello_V12_S25_Peaks'
if Peaks[12]=='1':#For V11
    k =[0,1,2,3,4,5,6]
else:
    k =[0,1,2,3,4,5,6,7]

L=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
for k in k:
    
    P= 0.2#Peak area plottede V12
    #P= 0.5#Peak area plottede V11
    I1 = np.where( (AllData[Time] > AllData[Peaks][k]-P ) & (AllData[Time] < AllData[Peaks][k]+P) )[0]
    P1 = AllData[Data][I1]
    T1 = AllData[Time][I1]
    if AllData[Peaks][k] > 10: #To print no data for no no data
        #Plotting
        L[k].plot(T1- AllData[Peaks][k], P1,'r.', label='_nolegend_')
        L[k].invert_yaxis()
        Quardline = np.linspace(-P, P, 10)
        L[k].plot(Quardline,np.array([0.595] * len(Quardline)),'k',linestyle='--', label='_nolegend_')  #V12
        #L[k].plot(Quardline,np.array([0.357] * len(Quardline)),'k',linestyle='--', label='_nolegend_') #V11
        L[k].plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
        L[k].set_xlim(-P, P) 
        L[k].set_ylim(0.73, -0.05) #V12
        #L[k].set_ylim(0.36, -0.05)#V11
        L[k].yaxis.set_major_locator(MaxNLocator(prune='lower'))
        yticks = L[k].yaxis.get_major_ticks()
        #yticks[-1].set_visible(False)
        yticks[-1].set_visible(False)
        xticks = L[k].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
        PoS=determine_closeness(AllData[Peaks][k],1808.76865399) #For V12
        TitelDate = PoS + 'ToE: ' + str("%.2f" % AllData[Peaks][k])
        L[k].legend(title=TitelDate, fontsize="x-small", loc="lower right")

    else:
        L[k].set_xlim(-P, P)
        L[k].yaxis.set_major_locator(MaxNLocator(prune='lower'))
        TitelDate = 'No Data'
        L[k].legend(title=TitelDate, fontsize="x-small", loc="lower right")
        xticks = L[k].xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
        continue 
    if Peaks[12]=='1' and k==6:#For V11
        #ax8.plot(Quardline,np.array([0.30] * len(Quardline)),'c', label='_nolegend_')  #Aperature max top line
        ax8.set_xlim(-P, P)
        ax8.set_ylim(0.36, -0.05)
        TitelDate = 'No Data'
        ax8.legend(title=TitelDate, fontsize="x-small", loc="lower right")
        #ax8.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
        ax8.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        xticks = ax8.xaxis.get_major_ticks()
        xticks[-1].set_visible(False)
        

for ax in fig.get_axes(): #o make grid plot stuff work
    ax.label_outer()
#%% Sammenligningsplot mellem Jordbasert of Nard
plt.figure()#I band comparision
plt.plot( AllData['Earthbased_I_Time']-AllData['Earthbased_I_Peaks']-0.025, AllData['Earthbased_I_Photo'], 'r.',  label = 'Zhuo', markersize=5)
#Nardiello comparision part
Data='Nardiello_V11_Data_Photo'
Time='Nardiello_V11_Data_Time'
Peaks='Nardiello_V11_Peaks'
k=5#eclipse used
P=0.4
I1 = np.where( (AllData[Time] > AllData[Peaks][k]-P ) & (AllData[Time] < AllData[Peaks][k]+P) )[0]
P1 = AllData[Data][I1]
T1 = AllData[Time][I1]
plt.plot(T1- AllData[Peaks][k]-0.025, P1,'b.', label='Nardiello', markersize=5)
#Rest of plot
Quardline = np.linspace(-0.5, 0.5, 10)
plt.plot(Quardline,np.array([0.357] * len(Quardline)),'k',linestyle='--', label='_nolegend_') 
plt.plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
plt.ylim([0.4, -0.05])
plt.xlim([-0.5, 0.5])
ax = plt.gca()
plt.title( 'NGC 188 - V11: Δ Magnitude, I-band comparison' )
plt.xlabel( 'Time before/after eclipse (days)' )
plt.ylabel( 'Δ Magnitude' )
plt.legend(fontsize="small", loc="lower left")
TitelDate = 'ToE: ' + str("%.2f" % AllData['Earthbased_I_Peaks'])
plt.legend(title=TitelDate, fontsize="small", loc="lower right")
#%% Nardiello Compare internal Data
plt.figure()# comparision
#Nardiello comparision part#['Nard_V12_RAW_S40','Nard_V12_AP1_S40','Nard_V12_PSF_S40']
Data='Nardiello_V12_S40_Data_Photo'
Time='Nardiello_V12_S40_Data_Time'
Peaks='Nardiello_V12_S40_Peaks'
k=2#eclipse used
P=0.4
I1 = np.where( (AllData[Time] > AllData[Peaks][k]-P ) & (AllData[Time] < AllData[Peaks][k]+P) )[0]
P1 = AllData[Data][I1]
T1 = AllData[Time][I1]
plt.plot(T1- AllData[Peaks][k], P1,'v', label='BEST_PHOT_FLUX_COR', markersize=5)

P1 = AllData['Nard_V12_RAW_S40'][:,0][I1]
T1 = AllData['Nard_V12_RAW_S40'][:,1][I1]
plt.plot(T1- AllData[Peaks][k]+0.05, P1,'*', label='BEST_PHOT_FLUX_RAW', markersize=5)

P1 = AllData['Nard_V12_AP1_S40'][:,0][I1]
T1 = AllData['Nard_V12_AP1_S40'][:,1][I1]
plt.plot(T1- AllData[Peaks][k]+0.10, P1,'.', label='PSF_FLUX_COR', markersize=5)

P1 = AllData['Nard_V12_PSF_S40'][:,0][I1]
T1 = AllData['Nard_V12_PSF_S40'][:,1][I1]
plt.plot(T1- AllData[Peaks][k]+0.15, P1,'+', label='AP1_FLUX_COR', markersize=5)

#Rest of plot
Quardline = np.linspace(-0.5, 0.6, 10)
plt.plot(Quardline,np.array([0.0] * len(Quardline)),'k',linestyle='--', label='_nolegend_')
plt.ylim([0.7, -0.05])
plt.xlim([-0.5, 0.6])
ax = plt.gca()
plt.title( 'NGC 188 - V12: Δ Magnitude, Comparison' )
plt.xlabel( 'Time before/after eclipse (days)' )
plt.ylabel( 'Δ Magnitude' )
plt.legend(fontsize="small", loc="lower left")
TitelDate = 'ToE: ' + str("%.2f" % AllData[Peaks][k])
plt.legend(title=TitelDate, fontsize="small", loc="lower left")
#%%Compare depth of eclipse from all sectors
#Nardiello data plottede
plt.figure(figsize=(8, 6), dpi=80)#
plt.plot([1,2,3,4,5],AllData['Nard_V12_S18_Depth'], '.',  label = 'S18', markersize=5)
plt.plot([6,7,8,9,10,11,12,13],AllData['Nard_V12_S20_Depth'], 'v',  label = 'S20', markersize=5)
plt.plot([14,15,16,17,18,19,20,21],AllData['Nard_V12_S25_Depth'], 's',  label = 'S25', markersize=5)
plt.plot([22,23,24,25,26,27,28],AllData['Nard_V12_S26_Depth'], '*',  label = 'S26', markersize=5)
plt.plot([29,30,31,32,33,34],AllData['Nard_V12_S40_Depth'], '+',  label = 'S40', markersize=5)
plt.plot([35,36,37,38,39],AllData['Nard_V12_S52_Depth'], 'd',  label = 'S52', markersize=5)
plt.plot([40,41,42],AllData['Nard_V12_S53_Depth'], '^',  label = 'S53', markersize=5)
#IS it P or S
All_Time_data = np.concatenate((AllData['Nard_V12_S18_Depth_Time'], AllData['Nard_V12_S20_Depth_Time'],AllData['Nard_V12_S25_Depth_Time'],AllData['Nard_V12_S26_Depth_Time'],AllData['Nard_V12_S40_Depth_Time'],AllData['Nard_V12_S52_Depth_Time'],AllData['Nard_V12_S53_Depth_Time']))
All_Depth_data = np.concatenate((AllData['Nard_V12_S18_Depth'], AllData['Nard_V12_S20_Depth'],AllData['Nard_V12_S25_Depth'],AllData['Nard_V12_S26_Depth'],AllData['Nard_V12_S40_Depth'],AllData['Nard_V12_S52_Depth'],AllData['Nard_V12_S53_Depth']))
P_T=[-1]
S_T=[-1]
P_D=[-1]
S_D=[-1]
for i in range(0,42):
    PoS=determine_closeness(All_Time_data[i]+57000, 50906.0389)#1808.76865399
    if PoS == 'P - ':
        P_T=P_T+[i+1]
        P_D=P_D+[All_Depth_data[i]-0.004]
    if PoS == 'S - ':
        S_T=S_T+[i+1]
        S_D=S_D+[All_Depth_data[i]-0.004]
plt.plot(P_T,P_D,'k', marker = "$P$", lw=0, markersize=4)
plt.plot(S_T,S_D,'k', marker = "$S$", lw=0, markersize=4)

Quardline = np.linspace(0, 43, 10)
plt.plot(Quardline,np.array([0.595] * len(Quardline)),'k',linestyle='--', label='Meibom') 
plt.ylim([0.5, 0.7])
plt.xlim([0, 43])
plt.gca().invert_yaxis()

#plt.xlim([-0.5, 0.5])
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.legend(loc='lower right')
plt.title( 'NGC 188 - V12: Eclipse Depth, Sector Comparison' )
plt.xlabel( 'Eclipses' )
plt.ylabel( 'Δ Magnitude at deepest point' )

