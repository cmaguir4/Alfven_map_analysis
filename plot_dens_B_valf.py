# read in mag fiel, electron deensity , alfven speed

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
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec

# calculate line length
R_apex=4040 #arcsed.5000-960
All_H_lines=[]
for i in range(-20,30,10):
    H_line=R_apex/math.cos(np.radians(abs(i)))
    All_H_lines.append(H_line)


for t in range(5):

    # mag fields
    bf=scipy.io.readsav('/Volumes/Storage/LASCO/Cube_bfield_map.sav')
    b_ar=bf.img_cube[t]

    #vel alf
    vel=scipy.io.readsav('/Volumes/Storage/LASCO/Cube.sav')
    vel_ar=vel.img_cube[t]
    
    heights_actual=np.linspace(0,All_H_lines[t],len(vel_ar)) #960 to 4000 arcsec
    heights_actual=np.around(heights_actual,0)


    #dens 
    den=scipy.io.readsav('/Volumes/Storage/Zucca/Dens_Cube-2.sav')
    dens_ar=den.img_cube[t]

    dens_AR=[]
    vel_AR=[]
    b_AR=[]

    for j in range(len(dens_ar)):
        dens_AR.append(dens_ar[j][0])
        b_AR.append(b_ar[j][0] )
        vel_AR.append(vel_ar[j][0])

    dens_AR=np.array(dens_AR)   
    vel_AR=np.array(vel_AR) 
    b_AR=np.array(b_AR) 

    incr=(5000-960)/len(dens_AR)
    dist=np.arange(960,5000,incr )

    dif=1000
    fd = interp1d( heights_actual,dens_AR)
    fv= interp1d( heights_actual,vel_AR)
    fb=interp1d( heights_actual,b_AR)

    newh = np.linspace(heights_actual[0], heights_actual[-1], num=dif, endpoint=True)


    vel_AR=fv(newh)
    dens_AR=fd(newh)
    b_AR=fb(newh)
    #test
    mu_bar=0.6
    m_p=1.67e-24
    new_v_alf=b_AR/np.sqrt(4*np.pi*dens_AR*mu_bar*m_p) # cgs units
    #new_v_alf=(new_v_alf/1e5)*1.2+(120)


    # in r sun 
    newh= newh/960
    fig=plt.figure(figsize=(7,9))
    
    plt.subplot(411)
    plt.title('Angle ' + str(70+(t*10)))
    plt.plot(newh,dens_AR,'--',label='Density',c='k')
    plt.ylabel('Density (cm^-3)')
    plt.xlim(0,2.5)
    plt.grid()

    plt.subplot(412)
    plt.plot(newh,b_AR,'--',label='Magnetic field (G)',c='b')
    plt.ylabel('Magnetic Field')
    plt.xlim(0,2.5)  
    plt.grid()
    
    plt.subplot(413) 
    plt.plot(newh,vel_AR,'--',label='Alfven velocity',c='g')
    plt.ylabel('Alfven Speed (km/s)')
    plt.grid()
    plt.xlim(0,2.5)

    plt.subplot(414) 
    plt.plot(newh,new_v_alf,'--',label='Alfven velocity',c='g')
    plt.ylabel('Calculated Speed (km/s)')
    plt.xlim(0,2.5)

    plt.xlabel('Distance along trace in R sun ')
    plt.grid()
    plt.xlim(0,2.5)
  #  fig.subplots_adjust(wspace=0, hspace=0.0)   
    plt.savefig('zucca_model_output'+str(t)+'.png')

