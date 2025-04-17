import subprocess 
import numpy as np
import json
import argparse 

molecules = ['pentane', 'hexane', 'heptane', 'octane', 'decane', 'pentadecane']

def calcHvap():
    
    for mol in molecules:
        # Extract potential energy from GROMACS .edr files using gmx energy
        command_gas_Pot = f"echo Potential | gmx energy -f ./Gas_sims/{mol}/nvt_{mol}_gas_f.edr -o ./Gas_sims/{mol}/Gas_Hvap_Pot.xvg -b 10000 -e 60000"
        subprocess.run(command_gas_Pot, shell=True, check=True)
        P_gas = np.loadtxt(f'./Gas_sims/{mol}/Gas_Hvap_Pot.xvg', comments=['#','@'])

        command_liquid_Pot = f"echo Potential | gmx energy -f ./Liquid_sims/{mol}/npt2_{mol}_liquid_f.edr -o ./Liquid_sims/{mol}/Liquid_Hvap_Pot.xvg -b 10000 -e 60000"
        subprocess.run(command_liquid_Pot, shell=True, check=True)
        P_liquid = np.loadtxt(f'./Liquid_sims/{mol}/Liquid_Hvap_Pot.xvg', comments=['#','@'])

        block_avg_P_liquid = []
        block_avg_P_gas = []

        n_blocks = 5
        block_size = len(P_gas) // n_blocks
        
        for i in range(n_blocks):
            Pl_kjm_block = P_liquid[i*block_size:(i+1)*block_size, 1]
            Pg_kjm_block = P_gas[i*block_size:(i+1)*block_size, 1]

            Avg_P_gas_block = np.average(Pg_kjm_block) / 1       # adjust for gas molecules
            Avg_P_liq_block = np.average(Pl_kjm_block) / 1000    # adjust for total liquid molecules
            block_avg_P_liquid.append(Avg_P_liq_block)
            block_avg_P_gas.append(Avg_P_gas_block)

        #hard code in mnumber of frames in edr file. above corrects for this hard code value
        # total_samples = 25000  # frames from 10 ns to 60 ns
        # num_blocks = 5
        # block_size = total_samples // num_blocks

        # for i in range(num_blocks):
        #     Pl_kjm_block = P_liquid[i*block_size:(i+1)*block_size, 1]
        #     Pg_kjm_block = P_gas[i*block_size:(i+1)*block_size, 1]

        #     Avg_P_gas_block = np.average(Pg_kjm_block) / 1       # adjust for gas molecules if needed
        #     Avg_P_liq_block = np.average(Pl_kjm_block) / 1000    # adjust for total liquid molecules
        #     block_avg_P_liquid.append(Avg_P_liq_block)
        #     block_avg_P_gas.append(Avg_P_gas_block)

        Avg_P_liq = np.mean(block_avg_P_liquid)
        Avg_P_gas = np.mean(block_avg_P_gas)

        stdev_P_liq = np.std(block_avg_P_liquid)
        stdev_P_gas = np.std(block_avg_P_gas)


        # Since the volume of the liquid is negligible in comparison with the volume 
        # of the gas and using the assumption that the gas is ideal so that the kinetic energies 
        # of a molecule in the gas and liquid phases are identical then you can use the formula:
        # dHvap= Upot(gas)-Upot(liq) +RT, which is a very reasonable approximation
        # dHvap = Ugas - Uliq + RT where U is the potential energy 

        R = 8.314 / 1000  # kJ/mol.K
        T = 300  # K

        Hvap = Avg_P_gas - Avg_P_liq + R * T

        # Propagation of error
        prop_error_std = np.sqrt(stdev_P_liq**2 + stdev_P_gas**2)
        prop_SE_liq = stdev_P_liq / np.sqrt(n_blocks)
        prop_SE_gas = stdev_P_gas / np.sqrt(n_blocks)
        prop_SE = np.sqrt(prop_SE_liq**2 + prop_SE_gas**2)

        data = {
            "P_liquid (kJ/mol)": block_avg_P_liquid,
            "P_gas (kJ/mol)": block_avg_P_gas,
            "H_vap (kJ/mol)": Hvap,
            "H_vap_std (kJ/mol)": prop_error_std,
            "H_vap_SE (kJ/mol)": prop_SE
        }

        with open(f"{mol}_Hvap_totals.json", "w") as f:
            json.dump(data, f, indent=4)

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-fpg', '--file_path_gas', type=str, required=True, help='Path to production GAS npt.edr file (60 ns assumed)')
    # parser.add_argument('-fpl', '--file_path_liquid', type=str, required=True, help='Path to production LIQUID npt.edr file (60 ns assumed)')
    # args = parser.parse_args()
    calcHvap() #args.file_path_gas, args.file_path_liquid)
