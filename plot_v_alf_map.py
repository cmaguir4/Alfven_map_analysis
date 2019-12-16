"""
Author:Ciara Maguire
Date: May 2019
Purpose: plot out v alf map from zucca model , takes in .sav file, extracts data array and displays plot
"""
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
import matplotlib.colors as colors
import pdb
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
import sunpy.cm
import sunpy
import matplotlib.patheffects as path_effects
legends=[110,100,90,80,70]
font_deg={'size': 10 ,'color':'white'}

def angle_line(angle):
    y=(x*np.tan(0.0174*angle)) - (x*np.tan(0.0174*angle))[0]
    return y
w=scipy.io.readsav('/Volumes/Storage/Zucca/Cartesian_Alfven_map20170902.sav')
vel_ar=w['cartesian_alfven']

col=np.load('/Volumes/Storage/LASCO/colors.npy')

prep_vel=(ma.masked_invalid(vel_ar['DATA'][0]))

#prep_vel=(ma.masked_where(prep_vel>8e3,prep_vel))

norm=colors.Normalize(vmin=100.,vmax=1500.)
#norm=colors.LogNorm()
fig, ax = plt.subplots(figsize=(8,8))
plt.imshow(prep_vel, aspect = 'equal', norm=norm, cmap='Spectral_r',origin='lower',extent=(-6000,6000,-6000,6000))
plt.colorbar( fraction=0.046, pad=0.04).set_label(label='km s$^{-1}$',size=14)
xi = [i for i in range(0, 600)]
x=[-3000,-2000,-1000,0,1000,2000,3000]
plt.xlim(-3000,3000)
plt.xticks(x)
plt.yticks(x)
plt.xlabel('Arcsec',fontsize=14)
plt.ylabel('Arcsec',fontsize=14)
plt.ylim(-3000,3000)
#plt.title('Alfvén Speed')
x=np.arange(966,3000,1)
t=0
for i in range(-20,30,10):
    y=angle_line(i)
    plt.plot(x,y,linewidth=1.2,color='k')
    plt.plot(x,y,color=col[t],linewidth=.8)

    t=t+1
circle1 = plt.Circle((0, 0), 960, color='white')
circle1_rim=plt.Circle((0, 0), 960, color='k',fill=False)
circle2 = plt.Circle((0, 0), 1920, color='white',fill=False)
circle2_rim = plt.Circle((0, 0), 1920,linewidth=2, color='k',fill=False)
circle3 = plt.Circle((0, 0), 2880, color='white',fill=False)
circle3_rim = plt.Circle((0, 0), 2880, color='k',fill=False,linewidth=2)

ax.text(2500, 680, (str(legends[4])+'$^\circ$'), fontdict=font_deg)
ax.text(2500, 350, (str(legends[3])+'$^\circ$'), fontdict=font_deg)
ax.text(2500, 40, (str(legends[2])+'$^\circ$'), fontdict=font_deg)
ax.text(2400, -220, (str(legends[1])+'$^\circ$'), fontdict=font_deg)
ax.text(2400, -500, (str(legends[0])+'$^\circ$'), fontdict=font_deg)

ax.add_artist(circle1)
ax.add_artist(circle1_rim)

ax.add_artist(circle2_rim)
ax.add_artist(circle2)

ax.add_artist(circle3_rim)
ax.add_artist(circle3)

font_white_a = {'size': 12, 'color':'white','weight':'bold'}


text=ax.text(1500, 1350, '2 R$\odot$', fontdict=font_white_a)
text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                       path_effects.Normal()])
text=ax.text(2600, 1350, '3 R$\odot$', fontdict=font_white_a)
text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                       path_effects.Normal()])
plt.xlim(800,3000)
plt.ylim(-1500,1500)
plt.savefig('figure3.png',dpi=500)#,quality=80)

print('plot made - figure3.png')
