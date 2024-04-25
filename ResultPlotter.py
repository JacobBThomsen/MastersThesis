#ResultPLotter
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
#%% Results are loaded
#os.chdir(r'C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/Results') 
#Lightcurve = np.loadtxt('MEIBOMV12_lc.v' )
#LightcurveFit = np.loadtxt('MEIBOMV12_lcfit.v' )
#P_Rv = np.loadtxt('out_Velocity_V12_A.txt' )
#S_Rv = np.loadtxt('out_Velocity_V12_B.txt' )
os.chdir(r'C:\Users\jacob\OneDrive\Dokumenter\Aarhus Universitet\10. Semester\Speciale\Data\Results\Nardiello') 
Lightcurve = np.loadtxt('out_lcV12S25.txt' )
LightcurveFit = np.loadtxt('out_lc_fitV12S25.txt' )

#%% Load rigth colums
Time_lc=Lightcurve[:,0]
Phase_lc=Lightcurve[:,3]
Raw_lc=Lightcurve[:,1]
Model_lc=Lightcurve[:,4]

Phase_Model_lc=LightcurveFit[:,0]
Mag_Model_lc=LightcurveFit[:,1]
#%%Plotting the data
plt.figure()
plt.plot(Time_lc,Model_lc,'r.') #The model data
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: JKTEBOP Nardiello S20 lightcurve')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'Magnitude' )
#%%
plt.figure()
plt.plot(Time_lc,Raw_lc,'r.') #The orginal raw data
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: Nardiello S20 lightcurve')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'Magnitude' )
#%%
plt.figure()
plt.plot(Phase_Model_lc,Mag_Model_lc,'r.') #The model data phase diagram
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: JKTEBOP Nardiello S20 phasediagram')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'Magnitude' )
#%%
plt.figure()
plt.plot(Phase_lc,Raw_lc,'r.') #The orginal raw data phase diagram
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: Nardiello S20 phasediagram')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'Magnitude' )
#%%
sort=np.argsort(Phase_lc)#To make the lines go right
plt.figure()
plt.plot(Phase_lc,Raw_lc,'r.') #The orginal raw data phase diagram
plt.plot(Phase_Model_lc,Mag_Model_lc) #The model data phase diagram
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: JKTEBOP Nardiello S20 lightcurve with Data')
plt.xlabel( 'Phase' )
plt.ylabel( 'Magnitude' )
#%%
plt.figure()
plt.plot(Time_lc,Model_lc,'b') #The model data
plt.plot(Time_lc,Raw_lc,'r.') #The orginal raw data
plt.gca().invert_yaxis()
plt.title( 'NGC 188 - V12: JKTEBOP Nardiello S20 lightcurve with Data')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'Magnitude' )
#%%Residual plot
OC_lc=Lightcurve[:,5]

plt.figure()
plt.plot(Time_lc,OC_lc,'r.') #The model data
Quardline=Quardline = np.linspace(Time_lc[0], Time_lc[-1], 10)
plt.plot(Quardline,np.array([0] * len(Quardline)),'k--')  #Zero line
plt.title( 'NGC 188 - V12: O-C diagram for S20 lightcurve over time')
plt.xlabel( 'Time [BJD-2457000]' )
plt.ylabel( 'O-C' )

plt.figure()
plt.plot(Phase_lc,OC_lc,'r.') #The model data
Quardline=Quardline = np.linspace(0, 1, 10)
plt.plot(Quardline,np.array([0] * len(Quardline)),'k--')  #Zero line
plt.title( 'NGC 188 - V12: O-C diagram for S20 lightcurve over the phase')
plt.xlabel( 'Phase' )
plt.ylabel( 'O-C' )
#%%Calculate uncertianty for file
UnArea = np.where(  (Phase_lc > 0.05) & (Phase_lc < 0.45 ))[0]
UnAreaTwo=np.where( (Phase_lc > 0.55 ) & (Phase_lc < 0.95) )[0]
np.append(UnArea,UnAreaTwo)
STD=np.std(Raw_lc[UnArea])#standard deviation
RMS = np.sqrt(np.mean(Raw_lc[UnArea]**2)) #Root mean square
