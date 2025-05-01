import MDAnalysis as mda
import numpy as np
from collections import Counter
from MDAnalysis.analysis.dihedrals import Dihedral
import pandas as pd
import matplotlib.pyplot as plt

def classify_dihedral(angle_deg):
    """Classify a dihedral angle into conformers based on standard cutoffs that Lipid14 uses
    aka gauche+ 0-120 degrees, trans 120-240 degrees, gaughe- 240-360 degrees.
    And MD anaalysis goes from 0 to 180 degrees - and +, so 
    trans is 120 < theta <= 180 or -180 <= theta < -120 
    gaughe+ is 0 <= theta <= 120
    gauche= is -120 <= theta <= 0 . """
    if 120 < angle_deg <= 180 or -180 <= angle_deg < -120: #notice in MDAnalysis, cis is 0 degrees and trans is  pi radians (180 degrees)
        return "trans"
    elif 0 <= angle_deg <= 120: # gauche is 60 degrees 
        return "gaucheplus"
    elif -120 <= angle_deg < 0: # gauche is 60 degrees 
        return "gaucheminus"
    else: 
        return "NA" 

def analyze_conformers(top, traj, mol):
    '''THis will only work for alkanes, per atom selection language with resname ALK'''
    u = mda.Universe(top,traj)
    alkane_atoms = u.select_atoms("resid 1 and resname ALK and not name H*")


    # Ensure the atoms are ordered properly (assumes topology ordering matches connectivity)
    dihedral_ag_list = []

    # Create overlapping 4-atom groups
    for i in range(len(alkane_atoms) - 3): 
        atom_group = u.atoms[[alkane_atoms[i].index,
                            alkane_atoms[i+1].index,
                            alkane_atoms[i+2].index,
                            alkane_atoms[i+3].index]]
        dihedral_ag_list.append(atom_group)

    # Run Dihedral analysis to calc dihedral angle at all frames
    dih_analysis = Dihedral(dihedral_ag_list).run()

    import pandas as pd


    angles_deg = dih_analysis.results.angles
    df = pd.DataFrame(angles_deg, columns=[f'Dihedral {i+1}' for i in range(angles_deg.shape[1])])

    # Display results
    # print(df)

    conformer_counts = {
        "trans": 0,
        "gauche": 0,
        "NA": 0,
        "eg": 0,   # end gauche
        "gg": 0,   # double gauche
        "gtg": 0   # gtg + gtg'
    }

    for i in range(df.shape[1]): #start index 0, loop through each dihedral angle column 
        dihedrals = df.iloc[:, i] #this is a series of dihedral angles across all frames 
        # print(dihedrals)
        classify = [classify_dihedral(dihedral) for dihedral in dihedrals]
        df[f'dihedral {i+1} classified'] = classify

        # Count trans and gauche conformers
        counter = Counter(classify)
        conformer_counts["trans"] += counter["trans"]
        conformer_counts["gauche"] += counter["gaucheplus"]
        conformer_counts["gauche"] += counter["gaucheminus"]
        conformer_counts["NA"] += counter["NA"]


        # Check for end gauche
        if classify[0] in ("gaucheplus","gaucheminus") or classify[-1] in ("gaucheplus", "gaucheminus"):
            conformer_counts["eg"] += 1
        
        #check for gg
        for j in range(len(classify) -1):
            if classify[j] in ("gaucheplus", "gaucheminus") and classify[j+1] in ("gaucheplus", "gaucheminus"):
                conformer_counts["gg"] += 1

        # CHeck for gtg
        for j in range(len(classify) -2):
            if (classify[j] in ("gaucheplus", "gaucheminus") and 
                classify[j+1] == "trans" and 
                classify[j+2] in ("gaucheplus", "gaucheminus")):
                conformer_counts["gtg"] += 1


    num_dihedrals = angles_deg.shape[1]

    # Store percent gauche for each torsion
    percent_gauche = []

    for i in range(num_dihedrals):
        column_name = f'dihedral {i+1} classified'
        gauche_count = df[column_name].value_counts().get('gaucheplus' or 'gaucheminus', 0)
        trans_count = df[column_name].value_counts().get('trans', 0)
        percent = (gauche_count / len(df)) * 100
        percent_gauche.append(percent)

    # Plotting
    # plt.figure(figsize=(8, 5))
    plt.scatter(range(1, num_dihedrals + 1), percent_gauche)
    plt.xlabel("Torsion Number")
    plt.ylabel("Percent Gauche (%)")
    plt.title(f"Percent Gauche Conformations per Torsion in {mol}")
    plt.xticks(range(1, num_dihedrals + 1))
    plt.ylim(0, 100)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"TG_{mol}.png")

    # Total number of frames and dihedrals
    n_frames = len(df)
    n_dihedrals = angles_deg.shape[1]

    # Total gauche and trans counts
    gauche_total = 0
    trans_total = 0

    for i in range(n_dihedrals):
        col = f'dihedral {i+1} classified'
        value_counts = df[col].value_counts()
        gauche_total += value_counts.get("gaucheplus",0) + value_counts.get("gaucheminus", 0)
        trans_total += value_counts.get("trans", 0)

    # Normalize to number of torsions (so max total = n_dihedrals) Not necessary and I dont think Lipid14 does this in their paper but Im not 100%
    normalized_gauche = gauche_total / n_frames
    normalized_trans = trans_total / n_frames

    # Save or print results
    with open(f"{mol}_TG_analysis.txt", 'w') as f:
        f.write(f'{mol} with files {top} and {traj}')
        for conformer, count in conformer_counts.items():
            f.write(f"{conformer}: {count}\n")
        f.write("/n")
        f.write(f"Normalized gauche: {normalized_gauche:.3f}\n")
        f.write(f"Normalized trans: {normalized_trans:.3f}\n")
        f.write(f"t/g ratio: {normalized_trans/normalized_gauche}\n")
        f.write(f"Check (should sum to {n_dihedrals}): {normalized_gauche + normalized_trans:.3f}\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--gro", help="Topology file (e.g. gro)")
    parser.add_argument("-t","--traj", help="Trajectory file (e.g. xtc)")
    parser.add_argument("-m","--mol", help="molecule name (e.g. pentane)")
    args = parser.parse_args()

    analyze_conformers(args.gro, args.traj, args.mol)
