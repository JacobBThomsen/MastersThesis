#Plot generator til Bachelor
import numpy as np
import matplotlib.pyplot as plt

#Data_rå
V18_Rv_A = 'V20_Rv_prim.txt'
DateA, RVA, UsikA = np.loadtxt(V18_Rv_A, unpack=True)
V18_Rv_B = 'V20_Rv_sekund.txt'
DateB, RVB, UsikB = np.loadtxt(V18_Rv_B, unpack=True)
#Data_model
V18_Mod1 = 'out3m_Q2.v'
Phase, Mag, L1, L2, L3, Mag, RvAM, RvAB = np.loadtxt(V18_Mod1, unpack=True)
V18_Mod2 = 'out3l_Q2.v'
DateMod, MagData, UsikMagData, PhaseData, MagMod, OC = np.loadtxt(V18_Mod2, unpack=True)
V18_RVA = 'V20rva_Q2.out'
TimeRVA, RVAdata, RVAerror, RVAphase, RVAMod, OCRVA = np.loadtxt(V18_RVA, unpack=True)
V18_RVB = 'V20rvb_Q2.out'
TimeRVB, RVBdata, RVBerror, RVBphase, RVBMod, OCRVB = np.loadtxt(V18_RVB, unpack=True)


#Plot af RV Data over tid
plt.figure()
plt.plot(DateA, RVA, 'r^', markersize=8)
plt.plot(TimeRVA, RVAMod, 'y^', markersize=8)
#plt.errorbar(DateA, RVA, yerr=0.30, ls='none', capsize=5, color='red')
plt.plot(DateB, RVB, 'b.', markersize=11)
plt.plot(TimeRVB, RVBMod, 'c.', markersize=11)
#plt.errorbar(DateB, RVB, yerr=0.30, ls='none', capsize=5, color='blue')
plt.title('Data af RV for V20')
plt.ylabel('RV-A and RV-B')
plt.xlabel('Date')

#Plot af RV model og data over phase
plt.figure()
plt.plot(Phase, RvAM, 'r-', markersize=8)
plt.plot(Phase, RvAB, 'b-', markersize=8)
plt.plot(RVAphase, RVAdata, 'r^', markersize=8)
plt.plot(RVBphase, RVBdata, 'b.', markersize=12)
plt.title('Model og Data af RV for V18')
plt.ylabel('RV')
plt.xlabel('Phase')

#Plot af fase lysstyrke model og data
plt.figure()
plt.plot(PhaseData, MagData, 'b.', markersize=5)
plt.plot(PhaseData-1, MagData, 'b.', markersize=5) #Så den ser pæn ud
#plt.plot(PhaseData, MagMod, 'b.', markersize=8)
#plt.title('Magnitude, meassured and model')
#plt.ylabel('Magnitude')
#plt.xlabel('Phase')

#Plot af fase lysstyrke, model fra m
 #plt.figure()
plt.plot(Phase, Mag, 'r-', markersize=8)
plt.plot(Phase-1, Mag, 'r-', markersize=8) #Så den ser pæn ud
plt.gca().invert_yaxis()
plt.title('The Magnitude over phase, model and data')
plt.ylabel('Magnitude')
plt.xlabel('Phase')

#Plot af tids lysstyrke model og data
plt.figure()
plt.plot(DateMod, MagData, 'b.', markersize=4)
plt.plot(DateMod, MagMod, 'r.', markersize=4)
plt.gca().invert_yaxis()
plt.title('Magnitude, meassured and model')
plt.ylabel('Magnitude')
plt.xlabel('date')

#Plot af OC for magnitude
plt.figure()
plt.plot([0, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))
plt.plot(PhaseData, OC, 'r.', markersize=1)
#plt.errorbar(PhaseData, OC, UsikMagData, ls='none', capsize=0.1, color='red')
plt.title('O-C for the magnitude of V20')
plt.ylabel('O-C (Magnitude)')
plt.xlabel('Phase')

#Plot af OC for RV
plt.figure()
plt.plot([0, 1], [0, 0], color='black', linestyle='--', dashes=(8, 8))
plt.plot(RVAphase, OCRVA, 'r.', markersize=11)
plt.errorbar(RVAphase, OCRVA, RVAerror, ls='none', capsize=5, color='red')
plt.plot(RVBphase, OCRVB, 'b^', markersize=8)
plt.errorbar(RVBphase, OCRVB, RVBerror, ls='none', capsize=5, color='blue')
plt.title('O-C for the RV of V20')
plt.ylabel('O-C (RV)')
plt.xlabel('Phase')

#%matplotlib qt and %matplotlib inline