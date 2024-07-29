import CoolProp
import json
import thermo
import math
import numpy as np
import sys

N = 6000

def phi_exp(Name, t, p):
    try:
        handle = CoolProp.AbstractState("HEOS", Name)
        handle.update(CoolProp.PT_INPUTS, p, t)
        return handle.fugacity_coefficient(0)
    except:
        return math.nan

for i in [sys.argv[1]]:
    substance = {}
    Phi_PR = []
    Phi_cool = []
    handle = CoolProp.AbstractState("HEOS", i)
    
    pc = CoolProp.AbstractState.p_critical(handle)
    
    T = np.linspace(CoolProp.AbstractState.Tmin(handle), CoolProp.AbstractState.Tmax(handle), N)
    P = np.linspace(0.01 * pc, CoolProp.AbstractState.pmax(handle), N)
    
    for pr in P:
        Phi_PR.append([phi_exp(i, t, pr) for t in T])    

    with open('CoolProp ' + i + '.npz', 'wb') as f:
        np.savez_compressed(f, Phi_PR=Phi_PR)
