#Setting datafiles to the same time
#Packdges
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt #   %matplotlib qt and %matplotlib inline
from sympy import S, symbols, printing #For nice equations in legend
import os #For the looping over files
from matplotlib.ticker import MaxNLocator #For nice ticks
from IPython import get_ipython #For plot pop out
import os #to jumo around dir
import pandas as pd
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

#RV
DataRV=['Velocity_V11_A.txt','Velocity_V11_B.txt','Velocity_V12_A.txt','Velocity_V12_B.txt']
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
#%% Getting Nardiello data in dictonary
for j in j:
    if DataNardiello[j][19] == 'e':#filename[143:-27]
        AA = np.loadtxt(DataNardiello[j])
        AllData[NameNardiello[j] + "_Photo"] = AA[:,0]
        AllData[NameNardiello[j] + "_Time"] = AA[:,1]
            
    if DataNardiello[j][19] == 'l':#filename[143:-27]
        AA = np.loadtxt(DataNardiello[j])
        AllData[NameNardiello[j]] = AA
#Comp
AllData[CompareNardielloName[0]] = np.loadtxt(CompareNardielloData[0])
AllData[CompareNardielloName[1]] = np.loadtxt(CompareNardielloData[1])
AllData[CompareNardielloName[2]] = np.loadtxt(CompareNardielloData[2])
#%% Getting Jordbasert V11 data in dictonary
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
#%% Getting Meibom data in dictonary
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
#%% Getting RV data in dictonary
AA = np.loadtxt(DataRV[3])
AllData["RV_V12_B_Time"] = AA[:,0]
AllData["RV_V12_B_Data"] = AA[:,1]
AllData["RV_V12_B_Un"] = AA[:,2]
AA = np.loadtxt(DataRV[2])
AllData["RV_V12_A_Time"] = AA[:,0]
AllData["RV_V12_A_Data"] = AA[:,1]
AllData["RV_V12_A_Un"] = AA[:,2]
AA = np.loadtxt(DataRV[1])
AllData["RV_V11_B_Time"] = AA[:,0]
AllData["RV_V11_B_Data"] = AA[:,1]
AllData["RV_V11_B_Un"] = AA[:,2]
AA = np.loadtxt(DataRV[0])
AllData["RV_V11_A_Time"] = AA[:,0]
AllData["RV_V11_A_Data"] = AA[:,1]
AllData["RV_V11_A_Un"] = AA[:,2]
#%%Time correction for RV data
AllData["RV_V12_B_Time"]=AllData["RV_V12_B_Time"]-50000
AllData["RV_V12_A_Time"]=AllData["RV_V12_A_Time"]-50000
AllData["RV_V11_B_Time"]=AllData["RV_V11_B_Time"]-50000
AllData["RV_V11_A_Time"]=AllData["RV_V11_A_Time"]-50000
#Create new files
pd.DataFrame({'col1':AllData["RV_V12_B_Time"],'col2':AllData["RV_V12_B_Data"],'col3':AllData["RV_V12_B_Un"]}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Velocity_V12_B_TC.txt', index=False, header=None, sep='\t')
pd.DataFrame({'col1':AllData["RV_V12_A_Time"],'col2':AllData["RV_V12_A_Data"],'col3':AllData["RV_V12_A_Un"]}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Velocity_V12_A_TC.txt', index=False, header=None, sep='\t')
pd.DataFrame({'col1':AllData["RV_V11_B_Time"],'col2':AllData["RV_V11_B_Data"],'col3':AllData["RV_V11_B_Un"]}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Velocity_V11_B_TC.txt', index=False, header=None, sep='\t')
pd.DataFrame({'col1':AllData["RV_V11_A_Time"],'col2':AllData["RV_V11_A_Data"],'col3':AllData["RV_V11_A_Un"]}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Velocity_V11_A_TC.txt', index=False, header=None, sep='\t')
#%% Making Time corrected data for data V12
#Dataset
DataTime='Nardiello_V12_S25_Data_Time'
AllData[DataTime]=AllData[DataTime]+7000
#Uncertianties 
Uncern=[0.006]*len(AllData[DataTime])
#Save time corrected data with uncertianties
pd.DataFrame({'col1':AllData["Nardiello_V12_S25_Data_Time"],'col2':AllData["Nardiello_V12_S25_Data_Photo"],'col3':Uncern}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Nardiello_V12_S25_Data_Time_TC.txt', index=False, header=None, sep='\t')
#%% Making Time corrected data for data V11
#Dataset
DataTime='Nardiello_V11_Data_Time'
AllData[DataTime]=AllData[DataTime]+7000
#Uncertianties
Uncern=[0.006]*len(AllData[DataTime])
#Save time corrected data with uncertianties
pd.DataFrame({'col1':AllData["Nardiello_V11_Data_Time"],'col2':AllData["Nardiello_V11_Data_Photo"],'col3':Uncern}).to_csv('C:/Users/jacob/OneDrive/Dokumenter/Aarhus Universitet/10. Semester/Speciale/Data/All_Peaks/TimeCorrected/Nardiello_V11_Data_Time_TC.txt', index=False, header=None, sep='\t')
