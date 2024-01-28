# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 12:34:35 2024

@author: Lukas Monrad-Krohn
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


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
print(df)


fig, [ax1, ax2, ax3] = plt.subplots(3,1,figsize=(10,5))
fig.suptitle('Sonic at 4m (u and v not corresponing to northerly and easterly)')

ax1.set_ylabel('u [m/s]')
ax1.scatter(df['time'][::200], df['Ux_mps'][::200], s=1)

ax2.set_ylabel('v [m/s]')
ax2.scatter(df['time'][::200], df['Uy_mps'][::200], s=1)

ax3.set_ylabel('w [m/s]')
ax3.set_xlabel('time [UTC]')
ax3.scatter(df['time'][::200], df['Uz_mps'][::200], s=1)

plt.savefig('sonic_4m_20200127.png', dpi=300, bbox_inches='tight')
plt.show()