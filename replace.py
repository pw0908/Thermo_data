import os

compounds = ['n-Nonane', 'MethylLinolenate', 'DimethylCarbonate', 'R21', 'DiethylEther', 'trans-2-Butene', 
             'R245fa', 'ParaDeuterium', 'OrthoDeuterium', 'Isohexane', 'R365MFC', 'n-Dodecane', 'R410A', 
             'Deuterium', 'D4', 'R13', 'MD2M', 'n-Hexane', 'Methane', 'Ethane', 'CarbonylSulfide', 
             'EthylBenzene', 'CarbonMonoxide', 'Isopentane', 'Xenon', 'cis-2-Butene', 'R152A', 'Oxygen', 
             'EthyleneOxide', 'R1234ze(E)', 'n-Octane', 'R404A', 'R236EA', 'CycloHexane', 'n-Heptane', 
             'R22', 'R113', 'n-Pentane', 'MethylLinoleate', 'R11', 'SulfurDioxide', 'R23', 'Helium', 
             'R32', 'R227EA', 'R407C', 'HydrogenSulfide', 'Air', 'R245ca', 'Novec649', 'R143a', 'D5', 
             'R507A', 'R134a', 'Dichloroethane', 'ParaHydrogen', 'R1233zd(E)', 'Acetone', 'n-Decane', 
             'HeavyWater', 'MethylPalmitate', 'n-Propane', 'R115', 'R1234yf', 'R236FA', 'Ethylene', 
             'R116', 'MD4M', 'Benzene', 'Methanol', 'SulfurHexafluoride', 'o-Xylene', 'R125', 'Fluorine', 
             'R1234ze(Z)', 'CarbonDioxide', 'IsoButane', 'n-Butane', 'NitrousOxide', 'DimethylEther', 
             'RC318', 'Toluene', 'IsoButene', 'MethylStearate', 'Ammonia', 'Argon', 'R218', 'R41', 
             'Neon', 'Propyne', 'CycloPropane', 'R12', 'Nitrogen', 'Water', 'MethylOleate', 'R161', 
             'D6', 'SES36', 'HFE143m', 'n-Undecane', 'R123', 'HydrogenChloride', 'm-Xylene', 'R141b', 
             'R124', '1-Butene', 'Propylene', 'R14', 'p-Xylene', 'Cyclopentane', 'MDM', 'Hydrogen', 
             'Neopentane', 'Ethanol', 'OrthoHydrogen', 'R114', 'Krypton', 'MD3M', 'R1243zf', 'MM', 
             'R142b', 'R40', 'R13I1','cis-2-Butene']

EOS = ['cPR','ADPCSAFT', 'BACKSAFT', 'Berthelot', 'CKSAFT', 'Clausius', 'CPA', 'CPPCSAFT', 'PR','DAPT', 
       'EPPR78', 'GEPCSAFT' , 'GEPCSAFT' , 'HeterogcPCPSAFT', 'HomogcPCPSAFT', 'iPCSAFT', 'KU', 'LJSAFT',
       'ogSAFT', 'PatelTeja', 'PCPSAFT', 'PCSAFT', 'pharmaPCSAFT', 'PR78','PSRK', 'PTV', 
       'QCPR', 'OPCSAFT', 'RK', 'RKPR','SAFTgammaMie','SAFTVRMie', 'SAFTVRMie15', 'SAFTVRQMie', 
       'SAFTVRSMie', 'SAFTVRSW', 'sCKSAFT','sCPA', 'softSAFT2016','sPCSAFT', 'SRK', 
       'structSAFTgammaMie', 'tcPR', 'tcRK', 'TVTPR', 'gcsPCSAFT','TWUSRK', 'UMRPR', 'vdW', 'VTPR'] 


#compounds = ["Ethanol","HydrogenSulfide","Methanol","Novec649","n-Dodecane"]
#EOS = ["CPPCSAFT","PatelTeja","SAFTVRMie"]


def replace_component_name(template_file, new_component_name, new_eos_name, output_file):
    # Read the content of the template file
    with open(template_file, 'r') as f:
        template_content = f.read()

    # Replace the placeholders with the new component name and EOS
    modified_content = template_content.replace('COMPONENT_NAME', new_component_name).replace('EOS', new_eos_name)

    # Write the modified content to the output file
    with open(output_file, 'w') as f:
        f.write(modified_content)

ind = 1

for i in compounds:
    for j in EOS:
        # Check if the corresponding file exists
        file_name = f"{j} {i} 1.npz"
        # If the file doesn't exist, create and submit the job
        template_file = 'Data.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    for j in EOS:
        # Check if the corresponding file exists
        file_name = f"{j} {i} 1.npz"
        # If the file doesn't exist, create and submit the job
        template_file = 'vp_n.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    for j in EOS:
        # Check if the corresponding file exists
        file_name = f"{j} {i} 1.npz"
        # If the file doesn't exist, create and submit the job
        template_file = 'figures.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    # Check if the corresponding file exists
    file_name = f"CoolProp {i}.npz"
    # If the file doesn't exist, create and submit the job
    template_file = 'New_batch.sh'
    new_component_name = i
    output_file = 'compound_' + str(ind) + '.sh'
    replace_component_name(template_file, new_component_name, "EXP", output_file)
    os.system(f'sbatch {output_file}')
    ind += 1


for i in compounds:
    template_file = 'Properties_EXP.sh'
    new_component_name = i
    output_file = 'compound_' + str(ind) + '.sh'
    replace_component_name(template_file, new_component_name, "EXP", output_file)
    os.system(f'sbatch {output_file}')
    ind += 1


for i in compounds:
    for j in EOS:
        # Define the expected file name for density and vp
        density_file_name = f"{j} {i} density.npz"
        vp_file_name = f"Hv {j} {i}.npz"
        
        # Check if the vp file exists
        
        # If vp file does not exist, proceed with job submission
        template_file = 'Properties_EOS.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
            
            # Replace placeholders and submit the job
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1


for i in compounds:
    for j in EOS:
        # Check if the corresponding file exists
        file_name = f"{j} {i} 1.npz"
        # If the file doesn't exist, create and submit the job
        template_file = 'den_sat.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    for j in EOS:
        # Check if the corresponding file exists
        file_name = f"{j} {i} 1.npz"
        # If the file doesn't exist, create and submit the job
        template_file = 'den_fig.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

