#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 16:01:53 2020

@author: Ken.Pryor
"""
import matplotlib.pyplot as plt
from skewt import SkewT
from six import StringIO
import numpy as np

fig = plt.figure(figsize=(8, 8))

'''
MIRS NPP Profile at 30.2N/93.2W at 0751 UTC 27 August 2020
'''

data_txt = '''
1000  290.64  290.4
972   289.27  289.05
944   288.0   287.77
891   287.03  286.81
839   286.62  283.88
789   286.05  279.69
741   285.27  276.11
695   283.89  273.04
650   281.83  270.7
606   279.39  268.84
545   275.3   265.44
506   272.34  262.5
487   270.88  260.93
450   268.02  257.86
400   263.73  252.65
350   259.3   248.01
293   252.58  241.58
250   246.56  237.2
200   232.17  225.14
'''

sound_data = StringIO(data_txt)
pressure,temperature,dewpoint = np.loadtxt(sound_data, usecols=range(0, 3), unpack=True)

data_std_atm = '''
0.0  101325
0.3  97000
0.6  94000 
1.0  89876
1.5  84559 
2.0  79501
2.6  74000 
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
10.5 24540 
11.5 20984
'''

std_atm = StringIO(data_std_atm)
height = np.loadtxt(std_atm, usecols=range(0, 1), unpack=True)
print(len(height))
height_m = height * 1000
print(len(height_m))
pressure_pa = pressure
temperature_c = temperature - 273.15
dewpoint_c = dewpoint - 273.15

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
fig.suptitle('NPP MIRS 0751 UTC 27 August 2020')
plt.suptitle('NPP MIRS 0751 UTC 27 August 2020',size=12)
plt.savefig("skewt_0827_sb.png",dpi=250,bbox_inches='tight')
print("Figure saved")
plt.show()

#Compute the Microburst Windspeed Potential Index (MWPI)
Z_UP = 5.5
print("Z_UP = ", Z_UP)
P_UP = pressure[11]
print("P_UP = ", P_UP)
T_UP = temperature[11]
print("T_UP = ", T_UP)
TD_UP = dewpoint[11]
print("TD_UP = ", TD_UP)
Z_LO = 3.0
print("Z_LO = ", Z_LO)
P_LO = pressure[7]
print("P_LO = ", P_LO)
T_LO = temperature[7]
print("T_LO = ", T_LO)
TD_LO = dewpoint[7]
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
         
idx_pup_mid = np.where(pressure == 650)
idx_plo_mid = np.where(pressure == 839)
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
parcel_2=(839.0, 13.47, 10.73, 'mwpi')
S.make_skewt_axes(tmin=-40.,tmax=40.,pmin=100.,pmax=1050.)
S.add_profile()
S.lift_parcel(*parcel_2)
print(parcel_2)
P_lcl,P_lfc,P_el,CAPE,CIN=S.get_cape(*parcel)
print("CAPE = ",CAPE)
fig.suptitle('NPP MIRS 0751 UTC 27 August 2020')
plt.suptitle('NPP MIRS 0751 UTC 27 August 2020',size=12)
plt.savefig("skewt_mwpi_0827.png",dpi=250,bbox_inches='tight')
print("Figure saved")
plt.show()

