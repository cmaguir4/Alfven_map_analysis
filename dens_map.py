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
import math
from matplotlib.colors import Normalize
from scipy.optimize import curve_fit
np.set_printoptions(formatter={'all':lambda x: str(x)})
from scipy.interpolate import interp1d

w=scipy.io.readsav('Dens_Cube-2.sav')
harmon_dens=np.load('/Volumes/Storage/ILOFAR/fundon_dens.npy')
t=2 # 90 degrees

dens=w.img_cube[t]
dens_AR=[]
for j in dens:
#	print(j[0])
	dens_AR.append(j[0])
dens_AR=np.array(dens_AR)	
incr=(5000-960)/len(dens_AR)
dist=np.arange(960,5000,incr )

dif=1000
fa = interp1d( dist,dens_AR)
newh = np.linspace(dist[0], dist[-1], num=dif, endpoint=True)
dens_AR=fa(newh)

incr=(5000-960)/len(dens_AR)
dist=np.arange(960,5000,incr )


###############################################


fig=plt.figure()
ax=fig.add_subplot(1,1,1)
plt.plot(dist,dens_AR,'.')
plt.ylabel('Electron density')
plt.xlabel('Distance along apex (Arcsec)')
plt.ylim(0,0.5e8)
plt.xlim(960,2000)

#plt.hlines(harmon_dens[0]-.2e6,dist[0],dist[-1],color='r' )
heights=[]
for i in harmon_dens:
	
	plt.hlines((np.round(i,2)),dist[0],dist[-1] )
	g =np.full((1, len(dens_AR)),i)
	h= np.isclose(np.round(dens_AR,0),(np.round(i,2)),atol=.35e6)
	ind=(np.where(h==True))
	print(ind)
	plt.plot(dist[ind[0][0]], dens_AR[ind[0][0]], 'ro')
	heights.append(dist[ind[0][0]])
plt.hlines(harmon_dens[0],dist[0],dist[-1],color='b' )

heights_rsun=(np.array(heights)-960)*725e3/696e6
np.save('/Volumes/Storage/Zucca/fund_heights_rsun.npy',heights_rsun)
plt.show()
