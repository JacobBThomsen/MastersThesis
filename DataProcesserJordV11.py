#Jordbasert V11
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
#%%Load data
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\V11Jordbaseret') 
V_t = np.loadtxt("V11_50Bin_20160113_calib_V.txt", usecols=3, unpack=True)#Time
I_t = np.loadtxt("V11_50Bin_20160113_calib_I.txt", usecols=3, unpack=True)
V_m = np.loadtxt("V11_50Bin_20160113_calib_V.txt", usecols=4, unpack=True)#Magnitudes
I_m = np.loadtxt("V11_50Bin_20160113_calib_I.txt", usecols=4, unpack=True)
#%% Plot
plt.plot(V_t,V_m,'r.')
plt.plot(I_t,I_m,'r.')
#%%Find median
V_tM = np.where( (V_t > min(V_t) ) & (V_t < min(V_t)+0.08) )[0]#Data for median calculation
I_tM = np.where( (I_t > min(I_t) ) & (I_t < min(I_t)+0.08) )[0]#Data for median calculation
V_Mm=V_m - np.median(V_m[V_tM])#Change from median
I_Mm=I_m - np.median(I_m[I_tM])#Change from median
#%%PLot new plots
plt.figure()
plt.plot( V_t, V_Mm, 'r.',  label = 'I-Band', markersize=2)
plt.plot( I_t, I_Mm, 'b.' , label = 'V-Band', markersize=2)
#plt.xlim(min(V_t),max(V_t))
ax = plt.gca()
ax.invert_yaxis()
plt.title( 'NGC 188 - V11 - Earth Based' )
plt.xlabel( 'Time [HJD-2457000]' )
plt.ylabel( 'Δ Magnitude' )
plt.legend(fontsize="small", loc="lower left")
#%% Save data
All_Magnitudes = I_Mm
All_Time = I_t
All_Data = np.matrix([All_Magnitudes,All_Time])
All_Data=np.transpose(All_Data) #Swap rows and collums

I_1 = np.where( I_Mm==max(I_Mm) )[0]
Placements = I_t[I_1]
#Make into csv file
mat = All_Data
df = pd.DataFrame(data=mat.astype(float))
df.to_csv('V11_Earthbased_I.csv', sep=' ', header=False, float_format='%.10f', index=False)

mat = Placements
df = pd.DataFrame(mat) #pd.DataFrame(data=mat.astype(float))
df.to_csv('V11_Earthbased_all_placements_I.csv', header=False, float_format='%.10f', index=False)
