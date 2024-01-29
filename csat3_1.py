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
    #wdir = np.arctan2(-u, -v) * DperR
    wdir = np.mod(180+np.rad2deg(np.arctan2(u, v)),360)
    return wdir, wspd
    
    
### read the original data ####################################################
data_path = 'C:/Users/Lukas Monrad-Krohn/Desktop/master/3_unis24/agf350/fieldwork/data_raw/20240127/'
n_df = '20240127_CR3000-SN7687_Sonic_4m.dat'

df = pd.read_csv(data_path+n_df, skiprows=4, delimiter=',')
df.columns = ["time","RECORD","Ux_mps","Uy_mps","Uz_mps",
                  "T_sonic_degC","Diag_sonic"]
df['time'] = pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S', exact=False) #, format='%Y-%m-%d %H:%M:%S.%f'
df = df.loc[df['time']>pd.to_datetime('2024-01-25 23:59:59', format='%Y-%m-%d %H:%M:%S', exact=False)]
df['Ux_mps'] = df.Ux_mps.astype(float)
df['Uy_mps'] = df.Uy_mps.astype(float)
df['Uz_mps'] = df.Uz_mps.astype(float)
df['Diag_sonic'] = df.Diag_sonic.astype(int)
df['T_sonic_degC'] = df.T_sonic_degC.astype(float)
df.reset_index(inplace=True, drop=True)
df.drop(columns='RECORD', inplace=True)
print(df)

### drop many values to use resources better
df = df.loc[df.index < 200000]


### correction of coordinate systems ##########################################
df['u_geo'], df['v_geo'] = inst2geo_coords(df['Ux_mps'], df['Uy_mps'], -3.160)

### plotting all data timeseries ##############################################
fig, [ax1, ax2, ax3] = plt.subplots(3,1,figsize=(10,5))
ax1.set_title('Sonic at 4m (u and v not corresponing to northerly and easterly)')

ax1.set_ylabel('u [m/s]')
#ax1.scatter(df['time'][::200], df['Ux_mps'][::200], s=1) # instrument coordinates
ax1.scatter(df['time'][::200], df['u_geo'][::200], s=1)

ax2.set_ylabel('v [m/s]')
#ax2.scatter(df['time'][::200], df['Uy_mps'][::200], s=1) # instrument coordinates
ax2.scatter(df['time'][::200], df['v_geo'][::200], s=1)

ax3.set_ylabel('w [m/s]')
ax3.set_xlabel('time [UTC]')
ax3.scatter(df['time'][::200], df['Uz_mps'][::200], s=1) # assume instrument coordinates to be good

#plt.savefig('../plots/testing_sonic/sonic_4m_20200127.png', dpi=300, bbox_inches='tight')
plt.show()

#%%
### averaging for 1 min #######################################################
df['min'] = df['time'].round('min') # currently rounding after normal rules, could also do a floor rounding in list comprehension
# drop first and last minute because they are not full
# !!!! not necessary anymore
#df.drop(df[df.time <= df['min'][0]].index, inplace=True)
#df.drop(df[df.time >= df['min'][df.index[-1]]].index, inplace=True)
#df.reset_index(inplace=True)
df2 = df.groupby(['min'], as_index=False).mean()
sizes = np.array(df.groupby(['min']).size())

df2['wdir'], df2['wspd'] = wdir_wspd_from_uv(df2['u_geo'], df2['v_geo'])
#print(df2)

### plotting direction and speed timeseries ##############################################
fig, [ax1, ax2] = plt.subplots(2,1,figsize=(10,5))
ax1.set_title('Sonic-anemometer (CSAT3, 4m, 1min average)')

ax1.set_ylabel('wdir [°]')
ax1.set_ylim([0,360])
ax1.set_yticks([0,90,180,270,360])
ax1.set_xlim([df2['min'].min(), df2['min'].max()])
ax1.grid()
ax1.scatter(df2['min'], df2['wdir'], s=1)

ax2.set_ylabel('U [m/s]')
ax2.set_ylim([0,12.5])
ax2.grid()
ax2.set_xlim([df2['min'].min(), df2['min'].max()])
ax2.scatter(df2['min'], df2['wspd'], s=1)
ax2.set_xlabel('time [UTC]')

#plt.savefig('../plots/testing_sonic/sonic_wdirwspd_4m_20200127.png', dpi=300, bbox_inches='tight')
plt.show()


### calculate tke #############################################################

### first perturbation
# create time array for calculations
#time_per = np.array([[i]*1200 for i in df2.time]).flatten()
def calculate_perturbation(para, sizes):
    count = 0
    per_list = []
    
    for i in range(len(sizes)):
        current_mean = df2[para][i]
        for j in range(sizes[i]):
            perturbation = df[para][count]
            per_list.append(perturbation)
            count +=1
    print('count of perturbations calculated:', count)
    return np.array(per_list)

uprime = calculate_perturbation('Ux_mps', sizes) 
vprime = calculate_perturbation('Uy_mps', sizes)
wprime = calculate_perturbation('Uz_mps', sizes)
#%%
### then the TKE itself
tke = 0.5 * np.sqrt(uprime*2 + vprime*2 + wprime*2)

### plotting TKE ##############################################
fig, ax1 = plt.subplots(1,1,figsize=(10,6))
ax1.set_title('TKE')

ax1.set_ylabel(r'TKE [$m^2/s^2$]')
ax1.scatter(df['time'], tke, s=1)

plt.savefig('../plots/testing_sonic/tke_timeseries_20240127.png', dpi=300, bbox_inches='tight')
plt.show()



### plot of TKE versus wind speed

count = 0
windspeed = []

for i in range(len(sizes)):
    wspd = df2['wspd'][i]
    for j in range(sizes[i]):
        windspeed.append(wspd)
        count += 1

#%%
### plotting TKE ##############################################
fig, ax1 = plt.subplots(1,1,figsize=(10,6))
ax1.set_title('square root of TKE vs. 1 min averaged horizontal wind speed')

ax1.set_ylabel(r'$\sqrt{TKE} [m/s]$')
ax1.set_xlabel('wind speed [m/s]')
ax1.scatter(windspeed, np.sqrt(tke), s=1)

plt.savefig('../plots/testing_sonic/tke_vs_windspeed.png', dpi=300, bbox_inches='tight')
plt.show()

#%%
df['wdir'], df['wspd'] = wdir_wspd_from_uv(df['u_geo'], df['v_geo'])

fig, ax1 = plt.subplots(1,1,figsize=(10,6))
ax1.set_title('square root of TKE vs. horizontal wind speed')

ax1.set_ylabel(r'$\sqrt{TKE} [m/s]$')
ax1.set_xlabel('wind speed [m/s]')
ax1.scatter(df['wspd'], np.sqrt(tke), s=1)

plt.savefig('../plots/testing_sonic/tke_vs_windspeed2.png', dpi=300, bbox_inches='tight')
plt.show()

