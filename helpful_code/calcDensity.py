# after a 60 ns npt run, calc density 

import subprocess 
import numpy as np
import matplotlib.pyplot as plt
import json
import argparse 

#### Calculate density. 
# Follow Lipid14 paper for Hvap calc. 60ns sims, cur forst 10 ns, use last 50 ns for sampling
# standard deviation using block ageraves, 5 equal blocks of 10 ns 

def calcDensity(file_path):    
    command_liquid_Pot = f"echo Density | gmx energy -f {file_path} -o ./Liquid_density.xvg -b 10000 -e 60000"
    subprocess.run(command_liquid_Pot, shell=True, check=True)
    P_liquid = np.loadtxt('./Liquid_density.xvg', comments=['#','@'])
    
    block_avg_P_liquid = []
    #loop through blocks ot 10 ns and calcualte the average and save values in a list 
    # dynamically slice into equal-length blocks instead:
    n_blocks = 5
    block_size = len(P_liquid) // n_blocks
    
    for i in range(n_blocks):
        block = P_liquid[i * block_size : (i + 1) * block_size, 1]
        block_avg_P_liquid.append(np.mean(block))

    print(block_avg_P_liquid, "block_avg_P_liquid")
    
    #Now I have 2 lists, and I can calculate the average P +- std in the gas and liquid phase
    Avg_P_liq = np.mean(block_avg_P_liquid)
    print(Avg_P_liq, "Avg_P_liq")
    
    stdev_P_liq = np.std(block_avg_P_liquid)
    print(stdev_P_liq, "stdev_P_liq")
    
    # Since the volume of the liquid is negligible in comparison with the volume 
    # of the gas and using the assumption that the gas is ideal so that the kinetic energies 
    # of a molecule in the gas and liquid phases are identical then you can use the formula:
    # dHvap= Upot(gas)-Upot(liq) +RT, which is a very reasonable approximation
    # dHvap = Ugas - Uliq + RT where U is the potential energy 
    
    
    data = {
        "Density_liquid (kg/m^3)": block_avg_P_liquid,
        "Average density (kg/m^3)": Avg_P_liq,
        "density_std (kg/m^3)": stdev_P_liq,
        # "H_vap_SE (kJ/mol)": prop_SE
    }
    
    json_file = "Density_totals.json"
    with open(json_file, "w") as f:
        json.dump(data, f, indent=4)
    #Plot to check 
    plt.plot(P_liquid[:,0], P_liquid[:,1])
    plt.xlabel("Time (ps)")
    plt.ylabel("Density (kg/m^3)")
    plt.title("Density vs Time")
    plt.savefig("density_plot.png")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fp', '--file_path', type=str, help='Path to production npt.edr file. Note that this assumes simulations are 60 ns')
    args = parser.parse_args()
    calcDensity(args.file_path)
