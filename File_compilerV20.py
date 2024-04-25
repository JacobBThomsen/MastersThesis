#Filecompiler

import shutil



#%%
with open('V20_Qouaters_4.dat','wb') as wfd:
    for f in ['V20_kepler_lc_out_16.dat', 'V20_kepler_lc_out_17.dat', 'V20_kepler_lc_out_15.dat', 'V20_kepler_lc_out_14.dat']:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)
            
#for f in ['V18_kepler_lc_out_1.dat','V18_kepler_lc_out_2.dat','V18_kepler_lc_out_3.dat','V18_kepler_lc_out_6.dat','V18_kepler_lc_out_8.dat','V18_kepler_lc_out_9.dat','V18_kepler_lc_out_10.dat','V18_kepler_lc_out_11.dat','V18_kepler_lc_out_12.dat','V18_kepler_lc_out_13.dat','V18_kepler_lc_out_14.dat','V18_kepler_lc_out_15.dat','V18_kepler_lc_out_16.dat','V18_kepler_lc_out_17.dat']: