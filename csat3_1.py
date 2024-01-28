# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 12:34:35 2024

@author: Lukas Monrad-Krohn
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


# function to change instrument coordinates to geographic coords
def inst2geo_coords(u_inst, v_inst, v_az):
    ''' function to calculate geographic coordinates from instrument coordinates and V_azimuth
    
        u_inst --> u wind in instrument coord (towards instrument from front, x-axis)
        v_inst --> v wind in instrument coord (towards instrument from right, y-axis in right handed coordinate system)
        v_az   --> azimuth correction angle (v_az = ang(x_inst, mag. north) + magnetic declination - 90°)
        
        returns: u_geo and v_geo'''
    
    RperD = np.pi/180
    u_geo = u_inst * np.cos(v_az*RperD) + v_inst * np.sin(v_az*RperD)
    v_geo = -u_inst * np.cos(v_az*RperD) + v_inst * np.cos(v_az*RperD)
    return u_geo, v_geo

# function to calculate wind direction from u and v wind
def wdir_wspd_from_uv(u, v):
    ''' function to calculate wind direction and speed from u and v wind speeds
    
        u,v  --> the horizontal components of the wind vector (x and y direction)
        
        returns: wind direction (wdir) and wind speed (wspd)
    '''
    DperR = 180/np.pi
    wspd = np.sqrt(u**2+v**2)
    wdir = np.arctan2(-u, -v) * DperR
    
    return wdir, wspd
    
    
### read the original data ####################################################
data_path = 'C:/Users/Lukas Monrad-Krohn/Desktop/master/3_unis24/agf350/fieldwork/data_raw/20240126/cards/'
n_df = 'TOA5_7687.Flux.dat'

df = pd.read_csv(data_path+n_df, skiprows=4, delimiter=',')
df.columns = ["time","RECORD","Ux_mps","Uy_mps","Uz_mps",
                  "T_sonic_degC","Diag_sonic"]
df['time'] = pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S')
df = df.loc[df['time']>pd.datetime(2024, 1, 25, 23, 59, 59)]
df['Ux_mps'] = df.Ux_mps.astype(float)
df['Uy_mps'] = df.Uy_mps.astype(float)
df['Uz_mps'] = df.Uz_mps.astype(float)
df.reset_index(inplace=True)
print(df)


### correction of coordinate systems ##########################################
df['u_geo'], df['v_geo'] = inst2geo_coords(df['Ux_mps'], df['Uy_mps'], -3.160)

### plotting all data timeseries ##############################################
fig, [ax1, ax2, ax3] = plt.subplots(3,1,figsize=(10,5))
fig.suptitle('Sonic at 4m (u and v not corresponing to northerly and easterly)')

ax1.set_ylabel('u [m/s]')
#ax1.scatter(df['time'][::200], df['Ux_mps'][::200], s=1) # instrument coordinates
ax1.scatter(df['time'][::200], df['u_geo'][::200], s=1)

ax2.set_ylabel('v [m/s]')
#ax2.scatter(df['time'][::200], df['Uy_mps'][::200], s=1) # instrument coordinates
ax2.scatter(df['time'][::200], df['v_geo'][::200], s=1)

ax3.set_ylabel('w [m/s]')
ax3.set_xlabel('time [UTC]')
ax3.scatter(df['time'][::200], df['Uz_mps'][::200], s=1) # assume instrument coordinates to be good

plt.savefig('sonic_4m_20200127.png', dpi=300, bbox_inches='tight')
plt.show()

#%%
### averaging for 1 min #######################################################
df['min'] = df['time'].round('min') # currently rounding after normal rules, could also do a floor rounding in list comprehension
df2 = df.groupby(['min']).mean()

df2['wdir'], df2['wspd'] = wdir_wspd_from_uv(df2['u_geo'], df2['v_geo'])
print(df2)
















