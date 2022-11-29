import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_excel('V&L_Output_Bats_PassesONLY.xlsx', sheet_name='Sheet1')#, 
excel = df.to_numpy()

F14C_OXII_norm = 1.34066
d13C_OXII_norm = -0.178



AC_full = excel[1:20, 6:20]
BL_full = excel[20:77, 6:20]
Braun_full = excel[77:96, 6:20]
LV1_full = excel[96:115, 6:20]
LV2_full = excel[115:134, 6:20]
LV3_full = excel[134:153, 6:20]
LV4_full = excel[153:172, 6:20]
LV5_full = excel[172:191, 6:20]
LV6_full = excel[191:210, 6:20]
OXA_full = excel[210:266, 6:20]

# sample [ p_i, R_mol, d14/12_mol, d13/12 ]

AC = np.array([ AC_full[:,1] * AC_full[:,12], AC_full[:,4], AC_full[:,5], AC_full[:,9] ]).astype('float')
BL = np.array([ BL_full[:,1] * BL_full[:,12], BL_full[:,4], BL_full[:,5], BL_full[:,9] ]).astype('float')
Braun = np.array([ Braun_full[:,1] * Braun_full[:,12], Braun_full[:,4], Braun_full[:,5], Braun_full[:,9] ]).astype('float')
LV1 = np.array([ LV1_full[:,1] * LV1_full[:,12], LV1_full[:,4], LV1_full[:,5], LV1_full[:,9] ]).astype('float')
LV2 = np.array([ LV2_full[:,1] * LV2_full[:,12], LV2_full[:,4], LV2_full[:,5], LV2_full[:,9] ]).astype('float')
LV3 = np.array([ LV3_full[:,1] * LV3_full[:,12], LV3_full[:,4], LV3_full[:,5], LV3_full[:,9] ]).astype('float')
LV4 = np.array([ LV4_full[:,1] * LV4_full[:,12], LV4_full[:,4], LV4_full[:,5], LV4_full[:,9] ]).astype('float')
LV5 = np.array([ LV5_full[:,1] * LV5_full[:,12], LV5_full[:,4], LV5_full[:,5], LV5_full[:,9] ]).astype('float')
LV6 = np.array([ LV6_full[:,1] * LV6_full[:,12], LV6_full[:,4], LV6_full[:,5], LV6_full[:,9] ]).astype('float')
OXA = np.array([ OXA_full[:,1] * OXA_full[:,12], OXA_full[:,4], OXA_full[:,5], OXA_full[:,9] ]).astype('float')

        
def blank_substraction(R,d):
    R_sub = R - sum(BL[0] * BL[1]) / sum(BL[0])
    d_bl = np.average(np.abs(BL[2]))
    d_sub = np.sqrt(d_bl**2 + d**2)
    return R_sub, d_sub

def fractionation_correction(R,d,d13C):
    R_frac = R * (0.975 / (1+d13C/1000))**2
    d_frac = d * (0.975 / (1+d13C/1000))**2
    return R_frac, d_frac

def standard_normalisation(R,d,p):
    R_std_bl, d_std_bl =  blank_substraction(OXA[1], OXA[2])
    R_std_bl_f, d_std_bl_f = fractionation_correction(R_std_bl, d_std_bl, OXA[3])
    p_std = OXA[0]
    F14C = (sum(p * R) / sum(p)) * (F14C_OXII_norm / ( sum(p_std * R_std_bl_f) / sum(p_std) ))
    dF14C = F14C * np.sqrt( (d/R)**2 + ( (sum(p_std * d_std_bl_f)/sum(p_std)) / (sum(p_std * R_std_bl_f)/sum(p_std)) )**2 )
    return F14C, dF14C

def F14C_calculation(R_mol, d_mol, p, d13C):
    R_mol_bl, d_mol_bl = blank_substraction(R_mol, d_mol)
    R_mol_bl_f, d_mol_bl_f = fractionation_correction(R_mol_bl, d_mol_bl, d13C)
    F14C, dF14C = standard_normalisation(R_mol_bl_f, d_mol_bl_f, p)
    dF14C = sum(dF14C * p) / sum(p)
    return F14C, dF14C

F14C_1, dF14C_1 = F14C_calculation(LV1[1], LV1[2], LV1[0], LV1[3])
F14C_2, dF14C_2 = F14C_calculation(LV2[1], LV2[2], LV2[0], LV2[3])
F14C_3, dF14C_3 = F14C_calculation(LV3[1], LV3[2], LV3[0], LV3[3])
F14C_4, dF14C_4 = F14C_calculation(LV4[1], LV4[2], LV4[0], LV4[3])
F14C_5, dF14C_5 = F14C_calculation(LV5[1], LV5[2], LV5[0], LV5[3])
F14C_6, dF14C_6 = F14C_calculation(LV6[1], LV6[2], LV6[0], LV6[3])
F14C_Braun, dF14C_Braun = F14C_calculation(Braun[1], Braun[2], Braun[0], Braun[3])
F14C_AC, dF14C_AC = F14C_calculation(AC[1], AC[2], AC[0], AC[3])



print(F14C_1, dF14C_1)
print(F14C_2, dF14C_2)
print(F14C_3, dF14C_3)
print(F14C_4, dF14C_4)
print(F14C_5, dF14C_5)
print(F14C_6, dF14C_6)
print(F14C_Braun, dF14C_Braun)
print(F14C_AC, dF14C_AC)










