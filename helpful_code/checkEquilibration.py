

import subprocess 
import numpy as np
import matplotlib.pyplot as plt
import json
import argparse 

molecules = ['hexane', 'heptane', 'octane', 'decane', 'pentadecane']

def plotVolume():  
    for mol in molecules:
        volume = f"echo Volume | gmx energy -f ./{mol}/npt_{mol}_liquid.edr -o ./{mol}/Liquid_volume.xvg -b 10000 -e 60000"
        subprocess.run(volume, shell=True, check=True)
        V_liq = np.loadtxt(f'./{mol}/Liquid_volume.xvg', comments=['#','@'])
        
        #Plot to check 
        plt.plot(V_liq[:,0], V_liq[:,1])
        plt.xlabel("Time (ps)")
        plt.ylabel("Volume (nm^3)")
        plt.title("Volume vs Time")
        plt.savefig(f"{mol}_volume_plot.png")
        plt.close()

def plotDims():  
    for mol in molecules:
        dims_run = f"echo 'Box-X\nBox-Y\nBox-Z' | gmx energy -f ./{mol}/npt_{mol}_liquid.edr -o ./{mol}/dims.xvg -b 10000 -e 60000"
        subprocess.run(dims_run, shell=True, check=True)
        dims = np.loadtxt(f'./{mol}/dims.xvg', comments=['#','@'])

        plt.plot(dims[:,0], dims[:,1], label="X", alpha = 0.2)
        plt.plot(dims[:,0], dims[:,2], label="Y", alpha = 0.2)
        plt.plot(dims[:,0], dims[:,3], label="Z", alpha = 0.2)
        plt.xlabel("Time (ps)")
        plt.ylabel("Dimension (nm)")
        plt.legend()
        plt.title("Dimension vs Time")
        plt.savefig(f"{mol}_Dimension_plot.png")
        plt.close()

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-fp', '--file_path', type=str, help='Path to Liquid_sims folder. Note that this assumes simulations are 60 ns')
    # args = parser.parse_args()
    plotVolume()
    plotDims()
