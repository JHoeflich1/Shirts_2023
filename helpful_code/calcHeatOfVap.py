import subprocess 
import numpy as np
import matplotlib.pyplot as plt
import json

#### Calculate Heat of vaporization. 
# Follow Lipid14 paper for Hvap calc. 60ns sims, cur forst 10 ns, use last 50 ns for sampling
# standard deviation using block ageraves, 5 equal blocks of 10 ns 

avgs_blocks = []
#I have 5, 10ns blocks in both the gas and liquid runs.
# for i in range(1,6):
    
command_gas_Pot = f"echo Potential | gmx energy -f ./Gas/nvt_methylacetate_gas.edr -o ./Gas/Gas_Hvap_Pot.xvg -b 10000 -e 60000"
subprocess.run(command_gas_Pot, shell=True, check=True)
P_gas = np.loadtxt('./Gas/Gas_Hvap_Pot.xvg', comments=['#','@'])

command_liquid_Pot = f"echo Potential| gmx energy -f ./Liquid/npt_methylacetate_liquid.edr -o ./Liquid/Liquid_Hvap_Pot.xvg -b 10000 -e 60000"
subprocess.run(command_liquid_Pot, shell=True, check=True)
P_liquid = np.loadtxt('./Liquid/Liquid_Hvap_Pot.xvg', comments=['#','@'])

block_avg_P_liquid = []
block_avg_P_gas = []
#loop through blocks ot 10 ns and calcualte the average and save values in a list 
#so I have 25_000 total frames in the edr file when cutting out the first 10ns. 
# so each block will be 5_000 ns
for i in range(1, 6):
    Pl_kjm_block = P_liquid[(i-1)*5000:i*5000, 1]
    Pg_kjm_block = P_gas[(i-1)*5000:i*5000, 1]

    #Calcualte the average Pot energy in the gas and liquid sims, and divide by total moelcules in the simulation 
    Avg_P_gas_block = np.average(Pg_kjm_block) / 1
    Avg_P_liq_block = np.average(Pl_kjm_block) / 1000 
    block_avg_P_liquid.append(Avg_P_liq_block)
    block_avg_P_gas.append(Avg_P_gas_block)

print(block_avg_P_liquid, "block_avg_P_liquid")
print(block_avg_P_gas, "block_avg_P_gas")
#Now I have 2 lists, and I can calculate the average P +- std in the gas and liquid phase
Avg_P_liq = np.mean(block_avg_P_liquid)
Avg_P_gas = np.mean(block_avg_P_gas)
print(Avg_P_liq, "Avg_P_liq")
print(Avg_P_gas, "Avg_P_gas")

stdev_P_liq = np.std(block_avg_P_liquid)
stdev_P_gas = np.std(block_avg_P_gas)
print(stdev_P_liq, "stdev_P_liq")
print(stdev_P_gas, "stdev_P_gas")

# Since the volume of the liquid is negligible in comparison with the volume 
# of the gas and using the assumption that the gas is ideal so that the kinetic energies 
# of a molecule in the gas and liquid phases are identical then you can use the formula:
# dHvap= Upot(gas)-Upot(liq) +RT, which is a very reasonable approximation
# dHvap = Ugas - Uliq + RT where U is the potential energy 

R = 8.314/1000 #kJ/mol.K
T = 300 #K 

Hvap = Avg_P_gas - Avg_P_liq  + R*T

prop_error_std = np.sqrt(stdev_P_liq**2 + stdev_P_gas**2)
prop_SE_liq = stdev_P_liq / np.sqrt(len(block_avg_P_liquid))
prop_SE_gas = stdev_P_gas / np.sqrt(len(block_avg_P_gas))
prop_SE = np.sqrt(prop_SE_liq**2 + prop_SE_gas**2)

data = {
    "P_liquid (kJ/mol)": block_avg_P_liquid,
    "P_gas (kJ/mol)": block_avg_P_gas,
    "H_vap (kJ/mol)": Hvap,
    "H_vap_std (kJ/mol)": prop_error_std,
    "H_vap_SE (kJ/mol)": prop_SE
}

json_file = "Hvap_totals.json"
with open(json_file, "w") as f:
    json.dump(data, f, indent=4)

