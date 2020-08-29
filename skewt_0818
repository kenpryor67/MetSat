#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 17:12 2020

@author: Ken.Pryor
"""
import matplotlib.pyplot as plt
from skewt import SkewT
from six import StringIO
import numpy as np

fig = plt.figure(figsize=(8, 8))

data_txt = '''
905 32 21
850 28 14
780 19 10
700 11 3
670 9 1
620 5 -2
570 0 -4
500 -6 -8
475 -9 -12
430 -13 -19
400 -17 -24
350 -23 -32
300 -32 -41
250 -43 -51
200 -54 -61
150 -63 -73
135 -66 -76
115 -69 -79
100 -70 -81
'''

sound_data = StringIO(data_txt)
pressure,temperature,dewpoint = np.loadtxt(sound_data, usecols=range(0, 3), unpack=True)

data_std_atm = '''
0.8  92000 
1.5  84559 
2.0  79501 
3.0  70121 
3.5  65780 
4.0  61660 
4.5  57752 
5.5  50539 
6.0  47217 
6.5  44075 
7.0  41105 
8.0  35651 
9.0  30800 
10.5  24540 
11.5  20984 
13.5  15327 
14.5  13100 
15.5  11197 
16.0  10352 
'''

std_atm = StringIO(data_std_atm)
height = np.loadtxt(std_atm, usecols=range(0, 1), unpack=True)
print(len(height))
height_m = height * 1000
print(len(height_m))
pressure_pa = pressure
temperature_c = temperature
dewpoint_c = dewpoint

mydata=dict(zip(('hght','pres','temp','dwpt'),(height_m,pressure_pa,temperature_c,dewpoint_c)))
print(mydata)
S=SkewT.Sounding(soundingdata=mydata)
S.make_skewt_axes(tmin=-40.,tmax=40.,pmin=100.,pmax=1050.)
S.add_profile()
parcel=S.get_parcel(method='sb')
S.lift_parcel(*parcel)
print(parcel)
P_lcl,P_lfc,P_el,CAPE,CIN=S.get_cape(*parcel)
print("CAPE = ",CAPE)
fig.suptitle('NOAA-20 NUCAPS 2109 UTC 18 August 2020')
plt.suptitle('NOAA-20 NUCAPS 2109 UTC 18 August 2020',size=12)
plt.savefig("skewt_0818_sb.png",dpi=250,bbox_inches='tight')
print("Figure saved")
plt.show()

#Compute the Microburst Windspeed Potential Index (MWPI)
Z_UP = 3.5
print("Z_UP = ", Z_UP)
P_UP = pressure[4]
print("P_UP = ", P_UP)
T_UP = temperature[4]
print("T_UP = ", T_UP)
TD_UP = dewpoint[4]
print("TD_UP = ", TD_UP)
Z_LO = 1.5
print("Z_LO = ", Z_LO)
P_LO = pressure[1]
print("P_LO = ", P_LO)
T_LO = temperature[1]
print("T_LO = ", T_LO)
TD_LO = dewpoint[1]
print("TD_LO = ", TD_LO)

def MWPI(Z_UP, Z_LO, T_UP, T_LO, TD_UP, TD_LO, CAPE):
    gamma = (T_LO - T_UP)/(Z_UP - Z_LO)
    DD_UP = T_UP - TD_UP
    print("DD_UP = ", DD_UP)
    DD_LO = T_LO - TD_LO
    print("DD_LO = ", DD_LO)
    DDD = DD_LO - DD_UP
    if DDD < 0:
        DDD = 0
    print("DDD = ", DDD)
    MWPI_IRv1 = (CAPE/100) + gamma + DDD
    MWPI_IRv2 = (CAPE/1000) + (gamma/5) + (DDD/5)
    WGP_IR = (0.4553 * MWPI_IRv1) + 28.769
    WGP_IRv2 = (0.35435365777*(MWPI_IRv2**2)) + (1.29598552473*MWPI_IRv2) + 33.8176788073
    return gamma, MWPI_IRv1, MWPI_IRv2, WGP_IR, WGP_IRv2
    
gamma, MWPI_IRv1, MWPI_IRv2, WGP_IR, WGP_IRv2 = MWPI(Z_UP, Z_LO, T_UP, T_LO, TD_UP, TD_LO, CAPE)

print("Gamma = ", gamma)
print("MWPI_IRv1 = ", MWPI_IRv1)
print("WGP_IR = ", WGP_IR)
print("MWPI_IRv2 = ", MWPI_IRv2)
print("WGP_IRv2 = ", WGP_IRv2)

def Haines_H(T_UP, T_LO, TD_LO):
    Tdiff = T_LO - T_UP
    print("Tdiff = ", Tdiff)
    DD_LO = T_LO - TD_LO
    print("DD_LO = ", DD_LO)
    if Tdiff < 17:
        ST = 1
    elif Tdiff >= 17 and Tdiff <= 21:
        ST = 2
    else:
        ST = 3   
    if DD_LO < 14:
        MT = 1
    elif DD_LO >= 14 and DD_LO <= 20:
        MT = 2
    else:    
        MT = 3
    HI = ST + MT    
    print("ST = ", ST)
    print("MT = ", MT)
    print("HI = ", HI)
    return HI
 
def Haines_M(T_UP, T_LO, TD_LO):
    Tdiff = T_LO - T_UP
    print("Tdiff = ", Tdiff)
    DD_LO = T_LO - TD_LO
    print("DD_LO = ", DD_LO)
    if Tdiff < 5:
        ST = 1
    elif Tdiff >= 5 and Tdiff <= 10:
        ST = 2
    else:
        ST = 3   
    if DD_LO < 5:
        MT = 1
    elif DD_LO >= 5 and DD_LO <= 12:
        MT = 2
    else:    
        MT = 3
    HI = ST + MT    
    print("ST = ", ST)
    print("MT = ", MT)
    print("HI = ", HI)
    return HI
        
def C_Haines(T_UP, T_LO, TD_LO):
    Tdiff = T_LO - T_UP
    print("Tdiff = ", Tdiff)
    DD_LO = T_LO - TD_LO
    print("DD_LO = ", DD_LO)
    if DD_LO >30:
        DD_LO=30
    CA=((T_LO-T_UP)/2)-2
    CB=((DD_LO)/3)-1
    if CB>5:
        CB=5+(CB-5)/2
    CH=CA+CB
    return CH
         
idx_pup_mid = np.where(pressure == 670)
idx_plo_mid = np.where(pressure == 850)
T_UP_mid = temperature_c[idx_pup_mid]
print("T_UP_mid = ", T_UP_mid)
T_LO_mid = temperature_c[idx_plo_mid]
print("T_LO_mid = ", T_LO_mid)
TD_LO_mid = dewpoint[idx_plo_mid]
print("TD_LO_mid = ", TD_LO_mid)
    
HI_M = Haines_M(T_UP_mid, T_LO_mid, TD_LO_mid)
HI_H = Haines_H(T_UP, T_LO, TD_LO)
CH = Haines_M(T_UP_mid, T_LO_mid, TD_LO_mid)
    
print("Haines Index MID = ", HI_M)
print("Haines Index HIGH = ", HI_H)
print("C-Haines Index = ", CH)

fig = plt.figure(figsize=(8, 8))
parcel_2=(850.0, 28.0, 14.0, 'mwpi')
S.make_skewt_axes(tmin=-40.,tmax=40.,pmin=100.,pmax=1050.)
S.add_profile()
S.lift_parcel(*parcel_2)
print(parcel_2)
P_lcl,P_lfc,P_el,CAPE,CIN=S.get_cape(*parcel)
print("CAPE = ",CAPE)
fig.suptitle('NOAA-20 NUCAPS 2109 UTC 18 August 2020')
plt.suptitle('NOAA-20 NUCAPS 2109 UTC 18 August 2020',size=12)
plt.savefig("skewt_mwpi_0818.png",dpi=250,bbox_inches='tight')
print("Figure saved")
plt.show()
