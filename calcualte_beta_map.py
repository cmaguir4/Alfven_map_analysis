import scipy.io
# Scientific libraries
import numpy as np
from scipy.optimize import curve_fit
# MatPlotlib
import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy as np
import numpy.ma as ma
from matplotlib import cm
import matplotlib.colors as colors
import math
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
w=scipy.io.readsav('Cartesian_density_map20170902.sav')
s=scipy.io.readsav('Cartesian_Alfven_map20170902.sav')
b=scipy.io.readsav('Cartesian_bfield_map20170902.sav')


###############################################
dens=w['car_den_map']
data_dens=dens['data'][0]

vel=s['cartesian_alfven']

b_field=b['bfield_map']
#prep_dens=np.log10(ma.masked_invalid(data_dens[0]))

n=data_dens
K=1.38064e-16       #Boltzmann Constant 1.3806488(13)×10−16 	erg K−1 (CGS)
T=1e6 # 1MK
p=n*K*T

B=b_field['data'][0]

mu_bar  = 0.6      # used in Mann et al. spherical density model, so used here for consistency
p_mag= B**2 / (2* mu_bar )

beta = p/p_mag

fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(1,2,1)
plt.title('Plasma beta map')
plt.imshow(beta,cmap='Spectral_r',origin='lower',extent=(-6000,6000,-6000,6000),vmin=0,vmax=3)
plt.colorbar(fraction=0.046, pad=0.04)
plt.xlabel('Arcsec')
plt.ylabel('Arcsec')


ax=fig.add_subplot(1,2,2)
plt.imshow(beta,cmap='Spectral_r',origin='lower',extent=(-6000,6000,-6000,6000),vmin=0,vmax=.3)
plt.xlim(-2500,2500)
plt.ylim(-2500,2500)
plt.colorbar(fraction=0.046, pad=0.04)

fig.subplots_adjust(wspace=0.5, hspace=0.0)	
plt.savefig('beta_map.png')
#axv=fig.add_subpl

