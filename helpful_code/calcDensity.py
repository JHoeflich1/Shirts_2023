# after a 60 ns npt run, calc density 

import subprocess 
import numpy as np
import matplotlib.pyplot as plt
import json
import argparse 

#### Calculate density. 
# Follow Lipid14 paper for Hvap calc. 60ns sims, cur forst 10 ns, use last 50 ns for sampling
# standard deviation using block ageraves, 5 equal blocks of 10 ns 
molecules = ['pentane','hexane', 'heptane', 'octane', 'decane', 'pentadecane']

def calcDensity():  
    for mol in molecules:
        command_liquid_Density = f"echo Density | gmx energy -f ./{mol}/npt2_{mol}_liquid_f.edr -o ./{mol}/Liquid_density.xvg -b 10000 -e 60000"
        subprocess.run(command_liquid_Density, shell=True, check=True)
        P_liquid = np.loadtxt(f'./{mol}/Liquid_density.xvg', comments=['#','@'])
        
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
        
        
        data = {
            "Density_liquid (kg/m^3)": block_avg_P_liquid,
            "Average density (kg/m^3)": Avg_P_liq,
            "density_std (kg/m^3)": stdev_P_liq,
        }
        
        json_file = f"{mol}_Density_totals.json"
        with open(json_file, "w") as f:
            json.dump(data, f, indent=4)
        #Plot to check 
        plt.plot(P_liquid[:,0], P_liquid[:,1])
        plt.xlabel("Time (ps)")
        plt.ylabel("Density (kg/m^3)")
        plt.title("Density vs Time")
        plt.savefig(f"{mol}_density_plot.png")
        plt.close()

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-fp', '--file_path', type=str, help='Path to Liquid_sims folder. Note that this assumes simulations are 60 ns')
    # args = parser.parse_args()
    calcDensity()
