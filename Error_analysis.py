import numpy as np
import json
import CoolProp
import statistics as st
from bisect import bisect
import os

CES = ["ADPCSAFT", "BACKSAFT", "Berthelot", "CKSAFT", "Clausius", "CPA", "cPR", "EPPR78", "GEPCSAFT", 
    "heterogcPCPSAFT", "homogcPCPSAFT", 
    "KU", "LJ SAFT", "ogSAFT", "Patel Teja", "PCPSAFT", "PCSAFT", "pharmaPCSAFT",
    "PR", "PR78", "PSRK", "PTV", "QCPR", "QPCSAFT",
    "RK", "RKPR", "SAFTgammaMie", "SAFTVRMie", "SAFTVRMie15", "SAFTVRQMie",
    "SAFTVRSW", "sCKSAFT", "sCPA", "softSAFT2016", "sPCSAFT", "SRK",
    "structSAFTgammaMie", "tcPR", "tcRK", "TVTPR", "TWUSRK", "UMRPR", "vdW", "VTPR"]
subs = ['n-Nonane', 'MethylLinolenate', 'DimethylCarbonate', 'R21', 'DiethylEther', 'trans-2-Butene', 'R245fa', 'ParaDeuterium', 'OrthoDeuterium', 'Isohexane', 'R365MFC', 'n-Dodecane', 'R410A', 'Deuterium', 'D4', 'R13', 'MD2M', 'n-Hexane', 'Methane', 'Ethane', 'CarbonylSulfide', 'EthylBenzene', 'CarbonMonoxide', 'Isopentane', 'Xenon', 'cis-2-Butene', 'R152A', 'Oxygen', 'EthyleneOxide', 'R1234ze(E)', 'n-Octane', 'R404A', 'R236EA', 'CycloHexane', 'n-Heptane', 'R22', 'R113', 'n-Pentane', 'MethylLinoleate', 'R11', 'SulfurDioxide', 'R23', 'Helium', 'R32', 'R227EA', 'R407C', 'HydrogenSulfide', 'Air', 'R245ca', 'Novec649', 'R143a', 'D5', 'R507A', 'R134a', 'Dichloroethane', 'ParaHydrogen', 'R1233zd(E)', 'Acetone', 'n-Decane', 'HeavyWater', 'MethylPalmitate', 'n-Propane', 'R115', 'R1234yf', 'R236FA', 'Ethylene', 'R116', 'MD4M', 'Benzene', 'Methanol', 'SulfurHexafluoride', 'o-Xylene', 'R125', 'Fluorine', 'R1234ze(Z)', 'CarbonDioxide', 'IsoButane', 'n-Butane', 'NitrousOxide', 'DimethylEther', 'RC318', 'Toluene', 'IsoButene', 'MethylStearate', 'Ammonia', 'Argon', 'R218', 'R41', 'Neon', 'Propyne', 'CycloPropane', 'R12', 'Nitrogen', 'Water', 'MethylOleate', 'R161', 'D6', 'SES36', 'HFE143m', 'n-Undecane', 'R123', 'HydrogenChloride', 'm-Xylene', 'R141b', 'R124', '1-Butene', 'Propylene', 'R14', 'p-Xylene', 'Cyclopentane', 'MDM', 'Hydrogen', 'Neopentane', 'Ethanol', 'OrthoHydrogen', 'R114', 'Krypton', 'MD3M', 'R1243zf', 'MM', 'R142b', 'R40', 'R13I1']


#"CPPCSAFT", "DAPT" , "iPCSAFT",
# Dictionary to store missing files
Missing_files = {}

# Function to check if a file exists and log missing files
def check_file_exists(file_name, equation, component):
    if not os.path.exists(file_name):
        if equation not in Missing_files:
            Missing_files[equation] = {}
        if component not in Missing_files[equation]:
            Missing_files[equation][component] = []
        Missing_files[equation][component].append(file_name)
        return False
    return True

# Function to calculate vapor pressure
def vapor_pressure(temperature, model):
    model.update(CoolProp.QT_INPUTS, 1, temperature)
    return model.p()

# Function to calculate the median of errors
def Median(x):
    try:
        return st.median(x)
    except:
        return []

# Function to merge dictionaries
def merge_dicts(*dicts):
    merged_dict = {}
    for d in dicts:
        merged_dict.update(d)
    return merged_dict

# Medi function definition
def Medi(Equation, Substance, Property):
    Eq = Equation
    comp = Substance
    #loading files
    if Property == "fugacity":
        if not check_file_exists(f'{Eq} {comp} 1.npz', Eq, comp) or not check_file_exists(f'CoolProp {comp}.npz', Eq, comp):
            return {}
        phi_EOS_data = np.load(f'{Eq} {comp} 1.npz')
        phi_exp_data = np.load(f'CoolProp {comp}.npz')
    elif Property == "density":
        if not check_file_exists(f'{Eq} {comp} density.npz', Eq, comp) or not check_file_exists(f'Ld CoolProp {comp}.npz', Eq, comp):
            return {}
        phi_EOS_data = np.load(f'{Eq} {comp} density.npz')
        phi_exp_data = np.load(f'Ld CoolProp {comp}.npz')

    with open(comp + ' metadata' + '.json') as f:
        metadata = json.load(f)

    #loading vectors
    phi_EOS = np.transpose(phi_EOS_data["Phi_PR"])
    phi_exp = phi_exp_data["Phi_PR"]
    EOS_Error = np.transpose(abs(phi_EOS - phi_exp) * 100 / phi_exp)
    Tgrid = np.array(metadata["grid"]["T"]).tolist()
    Pgrid = np.array(metadata["grid"]["P"]).tolist()

    handle = CoolProp.AbstractState("HEOS", comp)
    pc = handle.p_critical()
    Tc = handle.T_critical()

    Tc_ind = bisect(Tgrid, Tc)
    Pc_ind = bisect(Pgrid, pc)

    Error_compressed_gas = []
    Error_compressed_liquid = []

    #Separating sections that can be separated without the vapor pressure 
    Error_rest = EOS_Error[0:Tc_ind, 0:Pc_ind]
    Error_supercritical_liquid = EOS_Error[0:Tc_ind, Pc_ind:6000].flatten().tolist()
    Error_supercritical_gas = EOS_Error[Tc_ind:6000, 0:Pc_ind].flatten().tolist()
    Error_supercritical = EOS_Error[Tc_ind:6000, Pc_ind:6000].flatten().tolist()
    Error_shape = Error_rest.shape

    for i in range(Error_shape[0]):
        Temp = Tgrid[i]
        try:
            #tries to do it the easy way, by getting the index of vapor pressure, and dividing liquid/gas arrays by that index
            vp = vapor_pressure(Temp, handle)
            Psat_ind = bisect(Pgrid, vp)

            if Psat_ind == 0:
                Error_compressed_liquid.extend(EOS_Error[i, Psat_ind:Error_shape[1]])
            elif Psat_ind > 0:
                Error_compressed_gas.extend(EOS_Error[i, 0:Psat_ind])
                Error_compressed_liquid.extend(EOS_Error[i, Psat_ind:Error_shape[1]])
        except ValueError:
            #vapor_pressure fails constantly, so this function is here to do the same thing, it checks element by element using the last vp index as initial point, and when it sees a 
            #new phase change it divides the arrays
            Error_compressed_gas.extend(EOS_Error[i, 0:Psat_ind])
            for j in range(Psat_ind, Error_shape[1]):
                Phase = CoolProp.CoolProp.PhaseSI("T", Temp, "P", Pgrid[j], comp)
                Clas_error = EOS_Error[i, j]
                match Phase:
                    case "liquid":
                        Error_compressed_liquid.extend(EOS_Error[i, j:Error_shape[1]])
                        break
                    case "gas":
                        Error_compressed_gas.append(Clas_error)

    Med_compressed_liquid = st.median(Error_compressed_liquid)
    Med_compressed_gas = st.median(Error_compressed_gas)
    Med_supercritical_liquid = Median(Error_supercritical_liquid)
    Med_sup_gas = Median(Error_supercritical_gas)
    Med_supercritical = Median(Error_supercritical)

    return {
        f"{Property} Compressed liq": Med_compressed_liquid,
        f"{Property} Compressed gas": Med_compressed_gas,
        f"{Property} Supercritical liquid": Med_supercritical_liquid,
        f"{Property} Supcritical gas": Med_sup_gas,
        f"{Property} Supercritical": Med_supercritical
    }

# Function to calculate saturation properties
def sat_props(Eq, comp, type):
    try:
        if type == "vp":
            if not check_file_exists(f'vp {Eq} {comp}.npz', Eq, comp):
                return {}
            vp_EOS = np.load(f'vp {Eq} {comp}.npz')["vp"]

            try:
                vp_EXP = np.load(f'vp CoolProp {comp}.npz')["vp"]
            except (FileNotFoundError, KeyError, ValueError) as e:
                check_file_exists(f'vp CoolProp {comp}.npz', Eq, comp)
                return {}

            vp_EXP = vp_EXP.astype(np.float64)
            vp_error = abs(vp_EXP - vp_EOS) * 100 / vp_EXP
            vp_errorl = vp_error.tolist()
            return {"vp": st.median(vp_errorl)}

        elif type == "sat_fug":
            if not check_file_exists(f'{Eq}_{comp}_VP.npz', Eq, comp):
                return {}
            VP_error_data = np.load(f'{Eq}_{comp}_VP.npz')["Errors"]
            return {"sat_fug": st.median(VP_error_data.tolist())}

        elif type == "hv":
            if not check_file_exists(f'Hv {Eq} {comp}.npz', Eq, comp) or not check_file_exists(f'Hv CoolProp {comp}.npz', Eq, comp):
                return {}
            hv_EOS = np.load(f'Hv {Eq} {comp}.npz')["Hv"]
            hv_EXP = np.load(f'Hv CoolProp {comp}.npz')["vp"]
            print(hv_EOS,hv_EXP)
            hv_error = abs(hv_EXP - hv_EOS) * 100 / hv_EXP
            hv_error_l = hv_error.tolist()
            return {"hv": st.median(hv_error_l)}

        elif type == "den_liq":
            if not check_file_exists(f'{Eq}_{comp}_density_sat.npz', Eq, comp):
                return {}
            den_error = np.load(f'{Eq}_{comp}_density_sat.npz')["Error_liq_sat"]
            return {"den_liq_sat": st.median(den_error.tolist())}

        elif type == "den_gas":
            if not check_file_exists(f'{Eq}_{comp}_density_sat.npz', Eq, comp):
                return {}
            den_error = np.load(f'{Eq}_{comp}_density_sat.npz')["Error_gas_sat"]
            return {"den_gas_sat": st.median(den_error.tolist())}

    except FileNotFoundError as e:
        check_file_exists(str(e).split("'")[1], Eq, comp)
        return {}

# Initialize dictionary for storing results
Med = {}

# Main loop for processing equations and components
for Eq in CES:
    ES = {}
    for comp in subs:
        try:
            print(f'Processing {comp} in {Eq}')
            ES[comp] = merge_dicts(Medi(Eq, comp, "fugacity"), Medi(Eq, comp, "density"),
                                   sat_props(Eq, comp, "vp"), sat_props(Eq, comp, "sat_fug"),
                                   sat_props(Eq, comp, "hv"), sat_props(Eq, comp, "den_liq"),
                                   sat_props(Eq, comp, "den_gas"))

        except FileNotFoundError:
            pass
    Med[Eq] = ES

# Save the results to JSON files
with open('Med_results.json', 'w') as json_file:
    json.dump(Med, json_file, indent=4)
    print("Med_results.json saved successfully.")

with open('Missing_files.json', 'w') as json_file:
    json.dump(Missing_files, json_file, indent=4)
    print("Missing_files.json saved successfully.")
