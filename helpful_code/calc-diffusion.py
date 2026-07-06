#!/usr/bin/env python3
"""
Diffusion coefficient calculation pipeline.
Combines trajectory processing, MSD calculation, bootstrapping, and fitting
into a single efficient script.

Usage:
    python calc-diffusion.py                    # Full run with all molecules
    python calc-diffusion.py --molecules water  # Single molecule
    python calc-diffusion.py --molecules pentane hexane # Specific molecules only
    python calc-diffusion.py --sizes 512 1024  # Specific sizes only
    python calc-diffusion.py --test             # Quick test with 100 bootstraps
    python calc-diffusion.py --no-plot          # Skip plotting

Requires:
    - Trajectories and TPR files and EDR in current directory (5 ns trajectory, 5 ps time step)
        gmx trjconv -f prod_pentadecane_512.xtc -skip 10 -o prod_cut_pentadecane_512.xtc
Outputs:
    - diffusion_coefficients.csv: Raw diffusion coefficients (columns include `Ds` and `Stderr`)
        * Units: diffusion coefficients (`Ds`) and their standard errors (`Stderr`) are printed in cm^2/s.
    - diffusion_final.pkl: Pickle file with full results DataFrame (same units as above)
    - diffusion_fitting_results.csv: Curve fit parameters
    - Plots in 'plots/' directory for MSDs and bootstrap distributions
"""

import numpy as np
import pandas as pd
import os
import sys
import argparse
import multiprocessing
from pathlib import Path
import scipy
from scipy.optimize import curve_fit, leastsq
from scipy import stats
import matplotlib.pyplot as plt
from datetime import datetime


# CONFIGURATION

MOLECULES = {3: 'water', 17: 'pentane', 20: 'hexane', 23: 'heptane', 
             26: 'octane', 32: 'decane', 47: 'pentadecane'} # dict of alkane atom counts to names
SIZES = [512, 1024, 2048]
MSD_FRAME_COUNT = 1001 # total frames after trajetory is cut
MSD_TIME_PS = 5000
MSD_TIMESTEP_PS = 5
BOOTSTRAP_ITERATIONS = 5000
N_CORES = 8  # GROMACS processes per molecule



# UTILITY FUNCTIONS


def create_directories():
    """Create necessary output directories."""
    for directory in ['msds', 'ndxs', 'plots']:
        Path(directory).mkdir(exist_ok=True)


def log_message(msg, level="INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print("[{}] {}: {}".format(timestamp, level, msg))


def trajectory_already_converted(alkane, alkane_size):
    """Check if trajectory has already been converted."""
    return os.path.exists("prod_cut_test_{}_{}.xtc".format(alkane, alkane_size))

def all_msds_exist(alkane, alkane_size):
    """Check if all MSD files for this molecule/size already exist."""
    for i in range(alkane_size):
        if not os.path.exists(f"msds/msd_{alkane}_{alkane_size}_{i}.xvg"):
            return False
    return True



# STEP 1: TRAJECTORY CONVERSION (Once per molecule/size)

def convert_trajectory(alkane, alkane_size):
    """Convert trajectory with frame skipping. Only runs once per molecule/size."""
    
    trj_short = f"prod_cut_test_{alkane}_{alkane_size}.xtc"
    trj_original = f"prod_{alkane}_{alkane_size}.xtc"

    if os.path.exists(trj_short):
        log_message(f"Using existing converted trajectory: {trj_short}")
        return trj_short
    
    command = f"echo 0 | gmx trjconv -f {trj_original} -skip 10 -o {trj_short} -quiet 2>/dev/null"
    exit_code = os.system(command)
    
    if exit_code != 0:
        log_message(f"Trajectory conversion failed for {alkane}_{alkane_size}", level="ERROR")
        return None
    
    return trj_short

def check_converged(alkane, alkane_size):
    """Plot the density versus time to check for convergece and save plots to /plots."""
    #mkdir density plot directory if it doesn't exist
    Path("densities").mkdir(exist_ok=True)
    command = f"echo Density | gmx energy -f prod_{alkane}_{alkane_size}.edr -o densities/density_{alkane}_{alkane_size}.xvg  -quiet 2>/dev/null"
    density_file = f"densities/density_{alkane}_{alkane_size}.xvg"
    
    if not os.path.exists(density_file):
        log_message(f"Density file not found: {density_file}", level="WARNING")
        return False
    
    time = []
    density = []
    
    with open(density_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            vals = line.split()
            time.append(float(vals[0]))
            density.append(float(vals[1]))
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, density, 'b-', linewidth=2)
    plt.xlabel('Time (ps)')
    plt.ylabel('Density (g/cm³)')
    plt.title(f'Density vs Time for {alkane} {alkane_size}')
    plt.savefig(f'plots/density_{alkane}_{alkane_size}.png', dpi=100, bbox_inches='tight')
    plt.close()
    
    return True
# STEP 2: CREATE INDEX FILES

def create_index_files(alkane, alkane_size, alkane_atoms):
    """Create index files for all particles in parallel."""
    
    def create_single_index(i):
        ndx_file = f"ndxs/ndxs_{alkane}_{alkane_size}_{i}.ndx"
        if os.path.exists(ndx_file):
            return i
        
        with open(ndx_file, "w") as f:
            f.write(f"[ p{i} ]\n")
            n0 = alkane_atoms * i + 1
            f.write(f"{n0}\n")
        
        return i
    
    # Create indices sequentially (fast I/O operation)
    for i in range(alkane_size):
        create_single_index(i)
    
    log_message(f"Created {alkane_size} index files for {alkane}_{alkane_size}")



# STEP 3: CALCULATE MSDs (Parallelized per particle)


def calculate_single_msd(params):
    """Calculate MSD for a single particle."""
    particle_idx, alkane, alkane_size, trj_short, tpr_file = params
    
    msd_file = f"msds/msd_{alkane}_{alkane_size}_{particle_idx}.xvg"
    ndx_file = f"ndxs/ndxs_{alkane}_{alkane_size}_{particle_idx}.ndx"
    
    if os.path.exists(msd_file):
        return particle_idx  # Already computed
    
    if tpr_file is None or not os.path.exists(tpr_file):
        log_message(f"TPR file not found for particle {particle_idx}", level="ERROR")
        return particle_idx
    
    # Use stdin with proper I/O redirection
    command = (
        f"echo 0 | gmx msd -f {trj_short} -s {tpr_file} "
        f"-o {msd_file} -n {ndx_file} -quiet 2>/dev/null"
    )
    
    result = os.system(command)
    
    if result != 0:
        if particle_idx == 0:  # Only log once per batch
            log_message(f"gmx msd had issues (exit code {result}). Checking files...", level="WARNING")
            # Try to see what the actual error was
            os.system(f"echo 0 | gmx msd -f {trj_short} -s {tpr_file} -o {msd_file} -n {ndx_file} -quiet 2>/dev/null")
    
    return particle_idx


def calculate_msds_parallel(alkane, alkane_size, trj_short, alkane_atoms):
    """Calculate all MSDs for a molecule/size combination."""
    
    if all_msds_exist(alkane, alkane_size):
        log_message(f"All MSDs already exist for {alkane}_{alkane_size}")
        return True
    
    tpr_file = f"prod_{alkane}_{alkane_size}.tpr"
    
    # Create index files first
    create_index_files(alkane, alkane_size, alkane_atoms)
    
    # Prepare tasks
    tasks = [(i, alkane, alkane_size, trj_short, tpr_file) for i in range(alkane_size)]
    
    # Parallelize MSD calculations
    log_message(f"Calculating {alkane_size} MSDs for {alkane} (using {N_CORES} cores)")
    pool = multiprocessing.Pool(processes=N_CORES)
    pool.map(calculate_single_msd, tasks)
    pool.close()
    pool.join()
    
    return True



# STEP 4: READ MSD DATA INTO MATRIX


def read_msd_matrix(alkane, alkane_size):
    """Read all MSD files into a numpy matrix [particles x time]."""
    
    matrix = np.zeros((alkane_size, MSD_FRAME_COUNT))
    
    for i in range(alkane_size):
        msd_file = f"msds/msd_{alkane}_{alkane_size}_{i}.xvg"
        
        if not os.path.exists(msd_file):
            log_message(f"MSD file not found: {msd_file}", level="WARNING")
            continue
        
        with open(msd_file, 'r') as f:
            lines = f.readlines()
        
        itv = 0
        for line in lines:
            if line[0] != '#' and line[0] != '@':
                vals = line.split()
                matrix[i, itv] = float(vals[1])
                itv += 1
    
    return matrix



# STEP 5: BOOTSTRAP ANALYSIS


def bootstrap_diffusion_coefficient(matrix, molecule, size, n_bootstraps, plot=True):
    """
    Perform bootstrap analysis on MSD matrix.
    
    Returns:
        (diffusion_coefficient, standard_error)
    """
    
    nparticles, length = np.shape(matrix)
    time_array = np.linspace(0, MSD_TIME_PS, length)
    
    # Calculate average MSD
    avemsd = np.mean(matrix, axis=0)
    
    # Fit average MSD to get slope
    def fit_func(a, mymsd):
        return a * time_array - mymsd
    
    slope_avg = leastsq(fit_func, 0.1, args=avemsd)[0][0]
    D_avg = slope_avg / 6  # 3D diffusion: D = slope / (2*3)
    D_avg_cm2s = (D_avg / (10**7 * 10**7)) * 10**12  # Convert to cm^2/s
    
    # Perform bootstrap
    log_message(f"Running {n_bootstraps} bootstrap iterations for {molecule}_{size}")
    newmsd_data = np.zeros(np.shape(matrix))
    Ds = np.zeros(n_bootstraps)
    
    for n in range(n_bootstraps):
        if (n + 1) % 500 == 0:
            log_message(f"  Bootstrap iteration {n + 1}/{n_bootstraps}")
        
        # Resample particles
        for i in range(nparticles):
            newi = np.random.randint(0, nparticles)
            newmsd_data[i, :] = matrix[newi, :]
        
        # Calculate slope for resampled data
        new_avemsd = np.mean(newmsd_data, axis=0)
        slope = leastsq(fit_func, 0.1, args=new_avemsd)[0][0]
        Ds[n] = (slope / 6) / (10**7 * 10**7) * 10**12  # Convert to cm^2/s
    
    stderr = np.std(Ds)
    
    log_message(f"Result for {molecule}_{size}: D = {D_avg_cm2s:.4g} +/- {stderr:.4g} cm^2/s")
    
    # Optional plotting
    if plot:
        # Plot individual MSDs
        plt.figure(figsize=(10, 6))
        for i in range(min(50, nparticles)):
            plt.plot(matrix[i, :], 'b', alpha=0.02)
        plt.plot(avemsd, 'r', linewidth=2, label='Average MSD')
        plt.xlabel('Time (frames)')
        plt.ylabel('MSD')
        plt.title(f'All MSDs for {molecule} {size}')
        plt.legend()
        plt.savefig(f'plots/allMSD_{molecule}_{size}.png', dpi=100, bbox_inches='tight')
        plt.close()
        
        # Plot bootstrap histogram
        plt.figure(figsize=(10, 6))
        plt.hist(Ds, bins=50, alpha=0.7, edgecolor='black')
        plt.axvline(D_avg_cm2s, color='r', linestyle='--', linewidth=2, label=f'Mean: {D_avg_cm2s:.4g}')
        plt.xlabel('Diffusion Coefficient (cm²/s)')
        plt.ylabel('Frequency')
        plt.title(f'Bootstrap Distribution for {molecule} {size}')
        plt.legend()
        plt.savefig(f'plots/bootstrap_dist_{molecule}_{size}.png', dpi=100, bbox_inches='tight')
        plt.close()
    
    return D_avg_cm2s, stderr, Ds



# STEP 6: CALCULATION


def process_molecule_complete(molecule_name, alkane_atoms, size, plot=True, n_bootstraps=BOOTSTRAP_ITERATIONS):
    """
    Complete processing for a single molecule/size combination:
    1. Convert trajectory
    2. Calculate MSDs
    3. Run bootstrap
    
    Returns:
        (molecule, size, diffusion_D, diffusion_stderr)
    """
    
    try:
        log_message(f"Starting: {molecule_name} (N={size} particles)")
        
        # Step 1: Convert trajectory
        trj_short = convert_trajectory(molecule_name, size)
        if trj_short is None:
            return None
        
        # Step 2 & 3: Calculate MSDs
        log_message(f"DEBUG: Calling calculate_msds_parallel with alkane={molecule_name}, alkane_size={size}, alkane_atoms={alkane_atoms}")
        calculate_msds_parallel(molecule_name, size, trj_short, alkane_atoms)
        
        # Step 4: Read MSD data
        log_message(f"DEBUG: Calling read_msd_matrix with alkane={molecule_name}, alkane_size={size}")
        matrix = read_msd_matrix(molecule_name, size)
        
        # Step 5: Bootstrap
        log_message(f"DEBUG: Calling bootstrap_diffusion_coefficient with molecule={molecule_name}, size={size}")
        D, stderr, Ds = bootstrap_diffusion_coefficient(
            matrix, molecule_name, size, n_bootstraps=n_bootstraps, plot=plot
        )
        
        return (molecule_name, size, D, stderr)
    
    except Exception as e:
        import traceback
        log_message(f"ERROR processing {molecule_name}_{size}: {str(e)}", level="ERROR")
        log_message(f"TRACEBACK:\n{traceback.format_exc()}", level="ERROR")
        return None



# STEP 7: COLLECT RESULTS & CURVE FITTING

def fit_diffusion_models(DS_final, plot=True):
    """Fit diffusion coefficients to model: D = A * (1/L) + D_in"""
    
    # Load box sizes (created during MD setup)
    if not os.path.exists('box_sizes.csv'):
        log_message(f"WARNING: box_sizes.csv not found. Skipping curve fitting.", level="WARNING")
        return None
    
    box = pd.read_csv('box_sizes.csv')
    log_message(f"Loaded box_sizes.csv with {len(box)} rows and columns: {list(box.columns)}")
    if len(box) > 0:
        log_message(f"Sample data: {box.head(1).to_dict('records')[0]}")
    
    fit_results = {}
    diffusion_df = pd.DataFrame(
        columns=['molecule', 'D_in', 'D_inf_err', 'A_coef', 'A_coeff_err']
    )
    
    # Define fitting model
    def model(x, A, D_inf):
        return A * x + D_inf
    
    # Fit for each molecule
    for molecule in box['molecule'].unique():
        molecule_data = box[box['molecule'] == molecule]
        log_message(f"Processing molecule '{molecule}': {len(molecule_data)} entries in box_sizes.csv")
        
        Size_data = []
        DS_data = []
        DS_data_err = []
        
        for _, row in molecule_data.iterrows():
            box_size = row['box_length_avg']
            mol_size = row['size']
            log_message(f"  Checking size {mol_size} (box_size={box_size})")
            
            if (molecule, mol_size) not in DS_final.index:
                log_message(f"    No diffusion data found for ({molecule}, {mol_size})")
                continue
            
            Size_data.append(box_size)
            DS_data.append(DS_final.loc[(molecule, mol_size), 'Ds'])
            DS_data_err.append(DS_final.loc[(molecule, mol_size), 'Stderr'])
            log_message(f"    Found matching data: D={DS_data[-1]:.4g} ± {DS_data_err[-1]:.4g}")
        
        log_message(f"Collected {len(Size_data)} data points for {molecule}")
        
        if len(Size_data) < 2:
            log_message(f"Skipping curve fitting for {molecule}: need at least 2 data points, got {len(Size_data)}", level="WARNING")
            continue
        
        x_data = np.array(Size_data)
        y_data = np.array(DS_data)
        y_err = np.array(DS_data_err)
        
        # Curve fitting
        guess = [1.0, 0.0]
        params, covariance = curve_fit(
            model, 1 / x_data, y_data, sigma=y_err, p0=guess, absolute_sigma=True
        )
        perr = np.sqrt(np.diag(covariance))
        
        fit_results[molecule] = {'params': params, 'errors': perr}
        
        new_row = pd.DataFrame({
            'molecule': [molecule],
            'D_in': [params[1]],
            'D_inf_err': [perr[1]],
            'A_coef': [params[0]],
            'A_coeff_err': [perr[0]]
        })
        diffusion_df = pd.concat([diffusion_df, new_row], ignore_index=True)
        
        if plot:
            plt.figure(figsize=(10, 6))
            plt.errorbar(1 / x_data, y_data, yerr=y_err, fmt='o', markersize=8, label='Data')
            x_fit = np.linspace(min(1 / x_data), max(1 / x_data), 100)
            y_fit = model(x_fit, *params)
            plt.plot(x_fit, y_fit, 'r-', linewidth=2, label='Fit')
            plt.xlabel('1 / L (nm⁻¹)')
            plt.ylabel('D (cm²/s)')
            plt.title(f'Diffusion Fitting for {molecule}')
            plt.legend()
            plt.savefig(f'plots/{molecule}_diffusion_fit.png', dpi=100, bbox_inches='tight')
            plt.close()
            
            log_message(f"{molecule}: D_inf = {params[1]:.4g} ± {perr[1]:.4g} cm²/s")
    
    return diffusion_df



# MAIN


def main():
    """Main execution."""
    
    parser = argparse.ArgumentParser(description='Calculate diffusion coefficients from MD simulations')
    parser.add_argument('--molecules', nargs='+', help='Specific molecules to process')
    parser.add_argument('--sizes', nargs='+', type=int, help='Specific box sizes to process')
    parser.add_argument('--test', action='store_true', help='Quick test with 100 bootstrap iterations')
    parser.add_argument('--no-plot', action='store_true', help='Skip plotting')
    parser.add_argument('--no-fit', action='store_true', help='Skip curve fitting')
    args = parser.parse_args()
    
    # Create directories
    create_directories()
    
    # Configuration
    molecules_to_process = MOLECULES
    sizes_to_process = SIZES
    n_bootstraps = 100 if args.test else BOOTSTRAP_ITERATIONS
    plot = not args.no_plot
    
    if args.molecules:
        # Map molecule names back to atom counts
        molecules_to_process = {k: v for k, v in MOLECULES.items() if v in args.molecules}
    
    if args.sizes:
        sizes_to_process = args.sizes
    
    log_message(f"Starting diffusion calculation pipeline")
    log_message(f"Molecules: {list(molecules_to_process.values())}")
    log_message(f"Sizes: {sizes_to_process}")
    log_message(f"Bootstrap iterations: {n_bootstraps}")
    log_message(f"Plotting: {plot}")
    
    # Process each molecule/size sequentially (better resource usage)
    results = []
    for alkane_atoms, molecule_name in molecules_to_process.items():
        for size in sizes_to_process:
            result = process_molecule_complete(
                molecule_name, alkane_atoms, size, 
                plot=plot, n_bootstraps=n_bootstraps
            )
            if result is not None:
                results.append(result)
    
    # Collect results into DataFrame
    DS_final = pd.DataFrame(results, columns=['Molecule', 'Size', 'Ds', 'Stderr'])
    DS_final = DS_final.set_index(['Molecule', 'Size'])
    
    # Save results
    DS_final.to_csv('diffusion_coefficients.csv')
    DS_final.to_pickle('diffusion_final.pkl')
    log_message(f"Saved diffusion coefficients to diffusion_coefficients.csv and diffusion_final.pkl")
    
    # Curve fitting (only if we have multiple sizes)
    if not args.no_fit and len(sizes_to_process) > 1:
        log_message("Running curve fitting (multiple sizes detected)")
        diffusion_df = fit_diffusion_models(DS_final, plot=plot)
        if diffusion_df is not None:
            diffusion_df.to_csv('diffusion_fitting_results.csv', index=False)
            log_message("Saved fitting results to diffusion_fitting_results.csv")
    else:
        log_message("Skipping curve fitting (single size or disabled)")
    
    log_message("Pipeline complete!")
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(DS_final)


if __name__ == '__main__':
    main()
