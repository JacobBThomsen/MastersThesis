#FullPlotter - Plotter af alle Beachelor grafer i en
#Pakker
import numpy as np
import matplotlib.pyplot as plt
#%%
#Data_rå
V18_Rv_A = 'V18_Rv_prim.txt'
DateA, RVA, UsikA = np.loadtxt(V18_Rv_A, unpack=True)
V18_Rv_B = 'V18_Rv_sekund.txt'
DateB, RVB, UsikB = np.loadtxt(V18_Rv_B, unpack=True)
#Data_model
V18_Mod1 = 'out3m_V18old_Rc.v'
Phase, Mag, L1, L2, L3, Mag, RvAM, RvAB = np.loadtxt(V18_Mod1, unpack=True)
V18_Mod2 = 'out3l_V18old_Rc.v'
DateMod, MagData, UsikMagData, PhaseData, MagMod, OC = np.loadtxt(V18_Mod2, unpack=True)
V18_RVA = 'V18rva_V18old_Rc.out'
TimeRVA, RVAdata, RVAerror, RVAphase, RVAMod, OCRVA = np.loadtxt(V18_RVA, unpack=True)
V18_RVB = 'V18rvb_V18old_Rc.out'
TimeRVB, RVBdata, RVBerror, RVBphase, RVBMod, OCRVB = np.loadtxt(V18_RVB, unpack=True)
#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)
#%%
# Plot 1 - Lysstyrke og fase

fig = plt.figure()
gs = fig.add_gridspec(1, 2, wspace=0)
axs = gs.subplots(sharey=True)
#fig.suptitle('Lysstyrke (V) for V20 over fase')
fig.text(0.53, 0.28, 'Fase', ha='center')
N = 0.037 #aspect
axs[0].set_xlim([0.47, 0.51]) #V18 axs[0].set_xlim([0.47, 0.53]) #
axs[0].set_aspect(N)
axs[0].plot(PhaseData-1, MagData, 'b.', markersize=5)
axs[0].plot(PhaseData, MagData, 'b.', markersize=5)
axs[0].plot(Phase, Mag, 'r-', markersize=8)
axs[0].plot(Phase-1, Mag, 'r-', markersize=8)
axs[0].set_ylabel('Mag')
axs[0].set_xticks(np.arange(0.47, 0.51, 0.011))  # V18 #axs[0].set_xticks(np.arange(0.47, 0.53, 0.016))  # 

axs[1].set_xlim([-0.02, 0.02]) #V18# axs[1].set_xlim([-0.03, 0.03]) #V20 
axs[1].set_aspect(N)
axs[1].plot(PhaseData, MagData, 'b.', markersize=5)
axs[1].plot(PhaseData-1, MagData, 'b.', markersize=5)
axs[1].plot(Phase-1, Mag, 'r-', markersize=8)
axs[1].plot(Phase, Mag, 'r-', markersize=8)
axs[1].set_xticks(np.arange(-0.01, 0.02, 0.010)) #V18# axs[1].set_xticks(np.arange(-0.018, 0.03, 0.015)) # 

plt.gca().invert_yaxis()
#%%
'''
#Plot 3 all Lysstyrke
fig, ax = plt.subplots()
ax.set_aspect(1)
ax.set_ylabel('R')
ax.set_xlabel('Fase')
ax.plot(PhaseData, MagData, 'b.', markersize=5)
ax.plot(Phase, Mag, 'r-', markersize=8)
plt.gca().invert_yaxis()
'''
#%%
'''
#Plot af OC for magnitude
fig, ax = plt.subplots()
ax.set_aspect(7)
fig.text(0.53, 0.22, 'Fase', ha='center')
ax.plot([0, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))
ax.plot(PhaseData, OC, 'r.', markersize=2)
#plt.errorbar(PhaseData, OC, UsikMagData, ls='none', capsize=0.1, color='red')
#plt.title('O-C Størrelsesklassen af V20')
ax.set_ylabel('O-C (R)')
#ax.set_xlabel('Fase')
'''
#%%
#Plot af OC for magnitude
fig = plt.figure()
gs = fig.add_gridspec(1, 2, wspace=0)
axs = gs.subplots(sharey=True)
#fig.suptitle('Lysstyrke (V) for V20 over fase')
fig.text(0.53, 0.28, 'Fase', ha='center')

n = 1 # aspect
#axs[0].set_ylim([-0.07, 0.05])
axs[0].set_xlim([0.47, 0.51]) #V18# axs[0].set_xlim([0.47, 0.53]) #20# # 
axs[0].set_aspect(n)
axs[0].plot(PhaseData-1,  OC, 'r.', markersize=5)
axs[0].plot(PhaseData,  OC, 'r.', markersize=5)
axs[0].set_ylabel('O-C (Mag)')
axs[0].set_xticks(np.arange(0.47, 0.51, 0.011)) # V18# axs[0].set_xticks(np.arange(0.47, 0.53, 0.016)) #V20 #
axs[0].plot([0, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))

#axs[1].set_ylim([-0.07, 0.05])
axs[1].set_xlim([-0.02, 0.02]) # V18 # axs[1].set_xlim([-0.03, 0.03]) # V20 # 
axs[1].set_aspect(n)
axs[1].plot(PhaseData, OC, 'r.', markersize=5)
axs[1].plot(PhaseData-1, OC, 'r.', markersize=5)
axs[1].set_xticks(np.arange(-0.01, 0.02, 0.010)) #V18# axs[1].set_xticks(np.arange(-0.016, 0.03, 0.015)) # V20 #
axs[1].plot([-1, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))
#%%

#Plot 2  - RV og bane

fig = plt.figure()
gs = fig.add_gridspec(2, hspace=0, height_ratios=[3, 1])
axs = gs.subplots(sharex=True)
OCRVT=np.concatenate((OCRVB, OCRVA))
OCrms=np.nanstd(OCRVT) # O-C usikkerhed
#fig.suptitle('Model og Data af Radial hastigheden af V20 or tilhørende O-C diagram')
axs[0].plot(Phase, RvAM, 'r-', markersize=8)
axs[0].plot(Phase, RvAB, 'b-', markersize=8)
axs[0].plot(RVAphase, RVAdata,'r^', markersize=8)
axs[0].plot(RVBphase, RVBdata, 'b.', markersize=12) 
axs[0].set_ylabel('Radialhastigheden (km/s)')

axs[1].plot([0, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))
axs[1].plot(RVAphase, OCRVA,  'r^', markersize=6)
axs[1].errorbar(RVAphase, OCRVA, OCrms, ls='none', capsize=5, color='red')
axs[1].plot(RVBphase, OCRVB,'b.', markersize=8)
axs[1].errorbar(RVBphase, OCRVB, OCrms, ls='none', capsize=5, color='blue')
plt.ylabel('O-C (km/s)')
plt.xlabel('Fase')

#   %matplotlib qt and %matplotlib inline