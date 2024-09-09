import CoolProp
import numpy as np
import sys



def vp(handle, temperature):
    handle.update(CoolProp.QT_INPUTS, 1, temperature)
    return handle.p()

def Hv(handle, t):
    handle.update(CoolProp.QT_INPUTS, 0, t)  # 0 for saturated liquid
    h_liquid = handle.hmolar()  # Get the enthalpy of the liquid (J/kg)
    
    handle.update(CoolProp.QT_INPUTS, 1, t)  # 1 for saturated vapor
    h_vapor = handle.hmolar()  # Get the enthalpy of the vapor (J/kg)
    
    enthalpy_of_vaporization = h_vapor - h_liquid    
    return enthalpy_of_vaporization

def saturation_property(compound, property_name):
    N = 200
    handle = CoolProp.AbstractState("HEOS", compound)

    Tc = handle.T_critical()
    
    Tmin = handle.Tmin()
    Tmax = 0.9999 * Tc
    
    T = np.linspace(Tmin, Tmax, N)
    
    # Calculate the property for each temperature
    property_values = [property_name(handle, t) for t in T]  
    
    # Save the data with the property name in the filename
    filename = f'{property_name.__name__} CoolProp {compound}.npz'
    with open(filename, 'wb') as f:
        np.savez_compressed(f, **{"vp": property_values})

def density(compound):
    Phi_PR = []
    
    N = 6000
    handle = CoolProp.AbstractState("HEOS", compound)

    pc = handle.p_critical()
    
    Tmin = handle.Tmin()
    Tmax = handle.Tmax()
    
    pmin = 0.001*pc
    pmax = handle.pmax()
    
    T = np.linspace(Tmin, Tmax, N)
    P = np.linspace(pmin, pmax, N)
    
    for pr in P:
        Phi_cool = []
        for t in T:
            handle.update(CoolProp.PT_INPUTS, pr, t)
            Phi_cool.append(handle.rhomolar())
        Phi_PR.append(Phi_cool)    
        
    filename = f'Ld CoolProp {compound}.npz'
    with open(filename, 'wb') as f:
        np.savez_compressed(f, **{"Phi_PR": Phi_PR})

    return Phi_PR
    
def C(type,compound):
    Phi_PR = []
    
    N = 6000
    handle = CoolProp.AbstractState("HEOS", compound)

    pc = handle.p_critical()
    
    Tmin = handle.Tmin()
    Tmax = handle.Tmax()
    
    pmin = 0.001*pc
    pmax = handle.pmax()
    
    T = np.linspace(Tmin, Tmax, N)
    P = np.linspace(pmin, pmax, N)
    for pr in P:
        Phi_cool = []
        for t in T:
            handle.update(CoolProp.PT_INPUTS, pr, t)
            if type=="cv":
                Phi_cool.append(handle.cvmolar())
            elif type=="cp":
                Phi_cool.append(handle.cpmolar())
        Phi_PR.append(Phi_cool)    
        
    filename = f'{type} CoolProp {compound}.npz'
    with open(filename, 'wb') as f:
        np.savez_compressed(f, **{"Phi_PR": Phi_PR})

    return Phi_PR

density(sys.argv[1])

saturation_property(sys.argv[1],Hv)
saturation_property(sys.argv[1],vp)
C("cp",sys.argv[1])
