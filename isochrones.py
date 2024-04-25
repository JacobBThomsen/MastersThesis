# isochrones - bachelor
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
#%% - Age
Data = 'gsp035y30marcs_middel.txt'
MassM, MbolM, logTeM, logg, MV, B_V, V_R, V_I, V_J, V_K, J_K = np.loadtxt(Data, unpack=True)
Data = 'gsp035y30marcs_Young.txt'
MassY, MbolY, logTeY, logg, MV, B_V, V_R, V_I, V_J, V_K, J_K = np.loadtxt(Data, unpack=True)
Data = 'gsp035y30marcs_Old.txt'
MassO, MbolO, logTeO, logg, MV, B_V, V_R, V_I, V_J, V_K, J_K = np.loadtxt(Data, unpack=True)

#%% - Delta Fe/H
Data = 'Fe_H_30.txt'
MassL, MbolL, logTeL, logg, MV, B_V, V_R, V_I, V_J, V_K, J_K = np.loadtxt(Data, unpack=True)
Data = 'Fe_H_35.txt'
MassH, MbolH, logTeH, logg, MV, B_V, V_R, V_I, V_J, V_K, J_K = np.loadtxt(Data, unpack=True)
#%% 
fig, ax = plt.subplots(1,figsize=(9,6))
#plt.figure()
# V20p står først
#Kepler
VM = [V20Mp, V20Ms, V18Mp, V18Ms] = 1.0882, 0.8283, 0.9954, 0.9291
VR = [V20Rp, V20Rs, V18Rp, V18Rs] = 1.418, 0.800, 1.0855, 0.9589
VSM = [V20SMp, V20SMs, V18SMp, V18SMs] = 0.004, 0.002, 0.0045, 0.0046 # 1 sigma
VSR = [V20SRp, V20SRs, V18SRp, V18SRs] = 0.014, 0.011, 0.0089, 0.0106 # 1 sigma
V2SM = [V20S2Mp, V20S2Ms, V18S2Mp, V18S2Ms] = 0.007, 0.004, 0.009, 0.0091  # 2 sigma
V2SR = [V20SRp, V20S2Rs, V18S2Rp, V18S2Rs] = 0.042, 0.036, 0.022, 0.028 # 2 sigma
#Jordbaseret 2011
OM = [OV20Mp, OV20Ms, OV18Mp, OV18Ms] = 1.0868, 0.8276, 0.9955,  0.9293 
OR = [OV20Rp, OV20Rs, OV18Rp, OV18Rs] = 1.397, 0.7813, 1.1011, 0.9708
OSM = [OV20SMp, OV20SMs, OV18SMp, OV18SMs] = 0.0039, 0.0022, 0.0033, 0.0032 # 1 sigma
OSR =[OV20SRp, OV20SRs, OV18SRp, OV18SRs] = 0.013, 0.0053, 0.0068, 0.0089 # 1 sigma
#Jordbaseret 2021
JM = [JV20Mp, JV20Ms, JV18Mp, JV18Ms] = 1.088, 0.828, 0.9956,  0.9294 
JR = [JV20Rp, JV20Rs, JV18Rp, JV18Rs] = 1.395, 0.778, 1.0993, 0.9687
JSM = [JV20SMp, JV20SMs, JV18SMp, JV18SMs] = 0.004, 0.002, 0.004, 0.005 # 1 sigma
JSR =[JV20SRp, JV20SRs, JV18SRp, JV18SRs] = 0.004, 0.002, 0.007, 0.008 # 1 sigma
'''
StjKep = np.array([[1.088, 1.418],[0.828, 0.798], [0.9955, 1.0924],[0.9288,  0.9534]])# Erstant all med gennemsnit når tiden er inde
StjOld = np.array([[1.074, 1.399],[0.827, 0.768],[0.9953, 1.1011],[0.9291,   0.9708]])
Sigma = np.array([[0.001, 0.010],[0.0007, 0.007],[0.0044, 0.0091],[0.0047,  0.0104]]) #Errors
TwoSigma = np.array([[0.0023612, 0.04247],[0.00139,  0.03557],[0.00837, 0.01435],[0.009169, 0.02043]])
'''
#%%%
#plt.errorbar(StjKep[:,0], StjKep[:,1], Sigma[:,1], Sigma[:,0], ls='none', capsize=5, color='black')

def error(xy,r1,r2,n):
    Art = matplotlib.patches.Ellipse((xy), 2*r1, 2*r2, angle=0, edgecolor=n, facecolor="none")
    ax.add_patch(Art)
  
    # Kepler 1s
error((VM[0], VR[0]),VSM[0], VSR[0],"dimgray")
error((VM[1], VR[1]),VSM[1], VSR[1],"dimgray")
error((VM[2], VR[2]),VSM[2], VSR[2],"dimgray")
error((VM[3],  VR[3]),VSM[3],  VSR[3],"dimgray")
    # Kepler 2s
error((VM[0], VR[0]),V2SM[0], V2SR[0],"Silver")
error((VM[1], VR[1]),V2SM[1], V2SR[1],"Silver")
error((VM[2], VR[2]),V2SM[2], V2SR[2],"Silver")
error((VM[3],  VR[3]),V2SM[3],  V2SR[3],"Silver")

 # old 1s
error((OM[0], OR[0]),OSM[0], OSR[0],"dimgray")
error((OM[1], OR[1]),OSM[1], OSR[1],"dimgray")
error((OM[2], OR[2]),OSM[2], OSR[2],"dimgray")
error((OM[3],  OR[3]),OSM[3],  OSR[3],"dimgray")
''' 
    #Old data (2021)
error((JM[0], JR[0]),JSM[0], JSR[0],"dimgray")
error((JM[1], JR[1]),JSM[1], JSR[1],"dimgray")
error((JM[2], JR[2]),JSM[2], JSR[2],"dimgray")
error((JM[3],  JR[3]),JSM[3],  JSR[3],"dimgray")
 ''' 
#%%

#Radius regnes
def Radius(Mbol,logTe):
    L = 10**(0.4*(4.75-Mbol))*3.83*10**26 # L i W
    Step =7.126*(10**(-7))  # Stephan *4pi
    T = (10**(logTe))**4 # Teff^4
    R = (((L/(Step*(T)))**(1/2)))/(6.957*10**8) # Radius i R_sol
    return R
#Der plottes
def Plott(Mass,R,x,y):
    Cluster = plt.plot(Mass,R,x,label=y) #hvor er R
    plt.xlim([0.82, 1.11])
    plt.ylim([0.75, 1.5])

    #plt.legend([Cluster, Old, Kep], ["Age", "Old", "Kepler"], loc=0, title='test')
    #plt.legend([Cluster, Old, Kep], ["Age", "Old", "Kepler"])

    
RM = Radius(MbolM,logTeM)
RO = Radius(MbolO,logTeO)
RY = Radius(MbolY,logTeY)

RL = Radius(MbolL,logTeL)
RH = Radius(MbolH,logTeH)
#%%%
Plott(MassO,RO,'r-.','8.8 Gyr')
Plott(MassM,RM,'b-','8.3 Gyr')
Plott(MassY,RY,'g-.','7.8 Gyr')

Old = plt.plot(OM,OR,'y*',label='Brogaard et al. (2011)') 
Kep = plt.plot(VM,VR,'k.',label='Thomsen (Kepler)') 
#OldNew = plt.plot(JM,JR,'g^',label='Thomsen (Jordbaseret)') 

plt.text(1.088-0.035, 1.418,'V20p', size=13)
plt.text(0.831, 0.798+0.03,'V20s', size=13)
plt.text(0.9955-0.01, 1.0924+0.03,'V18p', size=13)
plt.text(0.9288,  0.9534+0.03,'V18s', size=13)

leg = ax.legend(prop={'size': 12});
#%%
Plott(MassL,RL,'b-.','[Fe/H]=0,30')
Plott(MassH,RH,'g-.','[Fe/H]=0,35')
Old = plt.plot(OM,OR,'y*',label='Brogaard et al. (2011)') 
Kep = plt.plot(VM,VR,'k.',label='Thomsen (Kepler)') 
plt.text(1.088-0.035, 1.418,'V20p', size=13)
plt.text(0.831, 0.798+0.03,'V20s', size=13)
plt.text(0.9955-0.01, 1.0924+0.03,'V18p', size=13)
plt.text(0.9288,  0.9534+0.03,'V18s', size=13)
leg = ax.legend(prop={'size': 12});
#%%Making it look nice
#plt.title('Isokroner (NGC 6791)', fontsize = 20)
plt.title('Isokroner (NGC 6791 - 8,2 Gyr)', fontsize = 15) # 
plt.ylabel('R ($R_{\odot}$)',fontsize = 15)
plt.xlabel('M ($M_{\odot}$)',fontsize = 15)


plt.rcParams.update({'font.size': 10})
