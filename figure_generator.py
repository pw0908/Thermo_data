import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import json
import math
import CoolProp
import matplotlib
import thermo
import sys



compounds = ['n-Nonane', 'MethylLinolenate', 'DimethylCarbonate', 'R21', 'DiethylEther', 'trans-2-Butene', 'R245fa', 'ParaDeuterium', 'OrthoDeuterium', 'Isohexane', 'R365MFC', 'n-Dodecane', 'R410A', 'Deuterium', 'D4', 'R13', 'MD2M', 'n-Hexane', 'Methane', 'Ethane', 'CarbonylSulfide', 'EthylBenzene', 'CarbonMonoxide', 'Isopentane', 'Xenon', 'cis-2-Butene', 'R152A', 'Oxygen', 'EthyleneOxide', 'R1234ze(E)', 'n-Octane', 'R404A', 'R236EA', 'CycloHexane', 'n-Heptane', 'R22', 'R113', 'n-Pentane', 'MethylLinoleate', 'R11', 'SulfurDioxide', 'R23', 'Helium', 'R32', 'R227EA', 'R407C', 'HydrogenSulfide', 'Air', 'R245ca', 'Novec649', 'R143a', 'D5', 'R507A', 'R134a', 'Dichloroethane', 'ParaHydrogen', 'R1233zd(E)', 'Acetone', 'n-Decane', 'HeavyWater', 'MethylPalmitate', 'n-Propane', 'R115', 'R1234yf', 'R236FA', 'Ethylene', 'R116', 'MD4M', 'Benzene', 'Methanol', 'SulfurHexafluoride', 'o-Xylene', 'R125', 'Fluorine', 'R1234ze(Z)', 'CarbonDioxide', 'IsoButane', 'n-Butane', 'NitrousOxide', 'DimethylEther', 'RC318', 'Toluene', 'IsoButene', 'MethylStearate', 'Ammonia', 'Argon', 'R218', 'R41', 'Neon', 'Propyne', 'CycloPropane', 'R12', 'Nitrogen', 'Water', 'MethylOleate', 'R161', 'D6', 'SES36', 'HFE143m', 'n-Undecane', 'R123', 'HydrogenChloride', 'm-Xylene', 'R141b', 'R124', '1-Butene', 'Propylene', 'R14', 'p-Xylene', 'Cyclopentane', 'MDM', 'Hydrogen', 'Neopentane', 'Ethanol', 'OrthoHydrogen', 'R114', 'Krypton', 'MD3M', 'R1243zf', 'MM', 'R142b', 'R40', 'R13I1']
#EOS=["ADPCSAFT","BACKSAFT","CKSAFT","CPA","CPPCSAFT","DAPT","GEPCSAFT","HeterogcPCPSAFT","HomogcPCPSAFT","LJSAFT","PCPSAFT","PCSAFT","QPCPSAFT","SAFTVRMie","SAFTVRMie15","SAFTVRQMie","SAFTVRSMie","SAFTVRSW","SAFTgammaMie","gcsPCSAFT","ogSAFT","pharmaPCSAFT","sCKSAFT","sPCSAFT","softSAFT","softSAFT2016","solidsoftSAFT","structSAFTgammaMie","sCPA"]
EOS=["TVTPR","iPCSAFT"]

def replace_none_in_list(data):
    return [math.nan if x is None else x for x in data]

def vapor_pressure(temperature, substance):
    handle = CoolProp.AbstractState("HEOS", substance)
    handle.update(CoolProp.QT_INPUTS, 1, temperature)
    return handle.p()



def graph(subs, k,CES):
    handle = CoolProp.AbstractState("HEOS", subs)
    pc = handle.p_critical()
    Tc = handle.T_critical()
    with open(subs+ ' metadata' + '.json') as f:
        metadata = json.load(f)
    
    phi_PR_data = np.load(CES+" "+subs+" 1.npz")
    phi_exp_data = np.load(f'CoolProp {subs}.npz')
    
    phi_PR=np.transpose(phi_PR_data["Phi_PR"])
    	
    phi_exp=phi_exp_data["Phi_cool"]
    
    
    VP_errpr_data = np.load(CES+"_"+subs+"_VP.npz")
    
    Tvp=VP_errpr_data["Temperature"]
    VP=VP_errpr_data["VaporPressure"]
    errors=VP_errpr_data["Errors"]
    def lowest_T(temperature):
        return vapor_pressure(temperature, subs) - metadata["limits"]["Pmin"]
    
    try:
        Tlow = sp.optimize.fsolve(lowest_T, (Tc + metadata["limits"]["Tmin"]) / 2)
    except ValueError:
        Tlow=metadata["limits"]["Tmin"]
    
    Error=np.transpose(abs(phi_PR-phi_exp)*100/phi_exp)
    P, T = np.meshgrid(np.array(metadata["grid"]["P"]), np.array(metadata["grid"]["T"]))

    levels = np.linspace(0, 30, 11)  # Adjust levels as per your preference
    cmap = plt.get_cmap('RdYlGn_r')
    cmap.set_over('red')
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    
    plt.figure(2**k)
    plt.title(CES+" fugacity coefficient error for " + subs)
    plt.yscale("log")
    contour = plt.contourf(T / Tc, P / pc, Error, levels=levels, cmap=cmap, extend='max')  # Adjust T and P to Tr and Pr
    plt.colorbar(contour, label='Error (%)')
    plt.grid()
    plt.xlabel("Tr ")
    plt.ylabel("Pr ")
    plt.axvline(x=1, linestyle='--', linewidth=3, color="k")
    plt.axhline(y=1, linestyle='--', linewidth=3, color="k")
    
    Tvp2 = np.linspace(Tlow, Tc, 30)
    VP2 = [vapor_pressure(jk, subs) / pc for jk in np.array(Tvp2)]  # Adjust VP to Pr
    plt.gca().set_ylim(bottom=0.01)
    plt.plot(Tvp2 / Tc, VP2, linestyle='-', linewidth=3, color='k')  # Adjust Tvp2 to Tr
    
    scatter_colors = cmap(norm(errors))
    plt.scatter(Tvp / Tc, VP/pc, c=scatter_colors, edgecolor='black', s=50, zorder=2)  # Adjust Tvp to Tr

    plt.savefig(CES+" "+subs+ " big.png")
    plt.close()
    
    plt.figure(3**k)
    plt.title(CES+" Fugacity coefficient error for " + subs)
    contour = plt.contourf(T / Tc, P / pc, Error, levels=levels, cmap=cmap, extend='max')  # Adjust T and P to Tr and Pr
    plt.colorbar(contour, label='Error (%)')
    plt.grid()
    plt.xlabel("Tr ")
    plt.yscale("log")
    plt.ylabel("Pr ")
    plt.gca().set_ylim(bottom=0.01)
    plt.xlim(metadata["limits"]["Tmin"]/Tc,2*(Tc-metadata["limits"]["Tmin"])/Tc)
    plt.axvline(x=1, linestyle='--', linewidth=3, color="k")
    plt.axhline(y=1, linestyle='--', linewidth=3, color="k")
    
    Tvp2 = np.linspace(Tlow, Tc, 30)
    VP2 = [vapor_pressure(jk, subs) / pc for jk in np.array(Tvp2)]  # Adjust VP to Pr

    plt.plot(Tvp2 / Tc, VP2, linestyle='-', linewidth=3, color='k')  # Adjust Tvp2 to Tr
    
    scatter_colors = cmap(norm(errors))
    plt.scatter(Tvp / Tc, VP/pc, c=scatter_colors, edgecolor='black', s=50, zorder=2)  # Adjust Tvp to Tr

    plt.savefig(CES+" "+subs+ " small.png")
    plt.close()

graph(sys.argv[1], 3,sys.argv[2])

