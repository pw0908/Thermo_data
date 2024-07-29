import os

compounds = ['n-Nonane', 'MethylLinolenate', 'DimethylCarbonate', 'R21', 'DiethylEther', 'trans-2-Butene', 'R245fa', 'ParaDeuterium', 'OrthoDeuterium', 'Isohexane', 'R365MFC', 'n-Dodecane', 'R410A', 'Deuterium', 'D4', 'R13', 'MD2M', 'n-Hexane', 'Methane', 'Ethane', 'CarbonylSulfide', 'EthylBenzene', 'CarbonMonoxide', 'Isopentane', 'Xenon', 'cis-2-Butene', 'R152A', 'Oxygen', 'EthyleneOxide', 'R1234ze(E)', 'n-Octane', 'R404A', 'R236EA', 'CycloHexane', 'n-Heptane', 'R22', 'R113', 'n-Pentane', 'MethylLinoleate', 'R11', 'SulfurDioxide', 'R23', 'Helium', 'R32', 'R227EA', 'R407C', 'HydrogenSulfide', 'Air', 'R245ca', 'Novec649', 'R143a', 'D5', 'R507A', 'R134a', 'Dichloroethane', 'ParaHydrogen', 'R1233zd(E)', 'Acetone', 'n-Decane', 'HeavyWater', 'MethylPalmitate', 'n-Propane', 'R115', 'R1234yf', 'R236FA', 'Ethylene', 'R116', 'MD4M', 'Benzene', 'Methanol', 'SulfurHexafluoride', 'o-Xylene', 'R125', 'Fluorine', 'R1234ze(Z)', 'CarbonDioxide', 'IsoButane', 'n-Butane', 'NitrousOxide', 'DimethylEther', 'RC318', 'Toluene', 'IsoButene', 'MethylStearate', 'Ammonia', 'Argon', 'R218', 'R41', 'Neon', 'Propyne', 'CycloPropane', 'R12', 'Nitrogen', 'Water', 'MethylOleate', 'R161', 'D6', 'SES36', 'HFE143m', 'n-Undecane', 'R123', 'HydrogenChloride', 'm-Xylene', 'R141b', 'R124', '1-Butene', 'Propylene', 'R14', 'p-Xylene', 'Cyclopentane', 'MDM', 'Hydrogen', 'Neopentane', 'Ethanol', 'OrthoHydrogen', 'R114', 'Krypton', 'MD3M', 'R1243zf', 'MM', 'R142b', 'R40', 'R13I1']
EOS = ["Berthelot", "Clausius", "KU", "PR", "PTV", "PatelTeja", "RK", "RKPR", "vdW", "EPPR78", "PR78", "PSRK", "QCPR", "SRK", "UMRPR", "cPR", "tcRK"]

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
        template_file = 'Data.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    for j in EOS:
        template_file = 'vp_n.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1

for i in compounds:
    for j in EOS:
        template_file = 'figures.sh'
        new_component_name = i
        new_eos_name = j
        output_file = 'compound_' + str(ind) + '.sh'
        replace_component_name(template_file, new_component_name, new_eos_name, output_file)
        os.system(f'sbatch {output_file}')
        ind += 1
