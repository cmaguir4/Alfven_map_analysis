from matplotlib import cm
import scipy.io
import math
from matplotlib.colors import Normalize
from sunpy import map
import matplotlib.colors as colors
import glob
from scipy import ndimage
smooth = ndimage.uniform_filter
import numpy as np
from astropy.io import fits
from astropy import units as u
import pdb
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
import sunpy.cm

w=scipy.io.readsav('/Volumes/Storage/LASCO/Cube.sav')
vel_ar=vel['cartesian_alfven']

dif=1000
fd = interp1d( dist,dens_AR)
fv= interp1d( dist,dens_AR)
newh = np.linspace(dist[0], dist[-1], num=dif, endpoint=True)

vel_AR=fa(newh)
dens_AR=fd(newh)

incr=(5000-960)/len(dens_AR)
dist=np.arange(960,5000,incr )



plt.figure()

for i in range(1,len(filenames)-1):
    w=scipy.io.readsav(filenames[i])
    
    ys=np.arange(0,560,20) #960 to 1500 arcsec
    ys=(ys*735e3)/1e6 #Mm above solar surface
    if i>5:
        plt.plot(ys,w.new_val_vel,'--',label=str(legends[i])+' degrees')
    else:
        plt.plot(ys,w.new_val_vel,label=str(legends[i])+' degrees')
    print(len(w.new_val_vel))
    plt.legend()
    plt.ylabel('Alfven Speed km/s')
    
    plt.xlabel('Mm from solar surface')
plt.savefig('alfven_speed_zucca.png')


