import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB


def extract_C_alpha(structure: PDB.Structure, model_index=0, chain_id: str = None):
    """Extracts alpha-carbon (Cα) coordinates from a specific model."""
    coords = []
    residues = []

    model = structure[model_index]

    for chain in model:
        if chain_id and chain.id != chain_id:
            continue

        for residue in chain:
            if "CA" in residue:
                atom = residue["CA"]
                coords.append(atom.coord)
                residues.append(residue)

    return np.array(coords), residues


def compute_contact_map(coords: np.ndarray, cutoff: float):
    """Computes the contact map for a given set of atomic coordinates."""
    N = len(coords)
    contact_map = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 2, N):
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < cutoff:
                contact_map[i, j] = 1
                contact_map[j, i] = 1

    return contact_map


def plot_free_energy(free_energy_data, selected_keys=None):
    """Plots the free energy landscapes for selected protein variants with custom colors.

    Parameters:
        free_energy_data (dict): Dictionary containing Q_values and F_values for each protein variant.
        selected_keys (list, optional): List of keys to plot. If None, all are plotted.
    """
    plt.figure(figsize=(15, 10))

    # If no specific keys are provided, plot all
    keys_to_plot = selected_keys if selected_keys else free_energy_data.keys()

    # Color definitions
    color_map = {
        "wt.pdb": "black",
    }

    # Generate shades for groups
    high_variants = [k for k in keys_to_plot if k.startswith("high")]
    low_variants = [k for k in keys_to_plot if k.startswith("low")]
    similar_variants = [k for k in keys_to_plot if k.startswith("similar")]

    high_colors = plt.cm.Blues(np.linspace(0.4, 1, len(high_variants)))
    low_colors = plt.cm.Reds(np.linspace(0.4, 1, len(low_variants)))
    similar_colors = plt.cm.Greens(np.linspace(0.4, 1, len(similar_variants)))

    # Assign colors dynamically
    color_map.update(dict(zip(high_variants, high_colors)))
    color_map.update(dict(zip(low_variants, low_colors)))
    color_map.update(dict(zip(similar_variants, similar_colors)))

    # Plot the free energy landscapes
    for fname in keys_to_plot:
        if fname in free_energy_data:
            data = free_energy_data[fname]
            color = color_map.get(fname, "gray")  # Default to gray if not categorized
            plt.plot(
                data["Q_values"], data["F_values"], label=fname, color=color, alpha=0.7
            )
        else:
            print(f"Warning: {fname} not found in free_energy_data")

    plt.xlabel("Fraction of Native Contacts (Q)")
    plt.ylabel("Free Energy F(Q) (kcal/mol)")
    plt.title("WSME Free Energy Landscape (SSA)")
    plt.legend()
    plt.grid()
    plt.show()


def plot_free_energy_difference(free_energy_data: dict):
    """Plots the free energy difference ΔF(Q) = F_mutant(Q) - F_WT(Q) for each mutant.

    Separates plots for low, similar, and high Tm mutants.

    Parameters:
        free_energy_data (dict): Dictionary containing Q_values and F_values for each protein variant.
    """
    # Extract WT data
    wt_Q = np.array(free_energy_data["wt.pdb"]["Q_values"])
    wt_F = np.array(free_energy_data["wt.pdb"]["F_values"])

    # Categorize mutants
    low_tm_mutants = [m for m in free_energy_data if "low_tm" in m]
    similar_tm_mutants = [m for m in free_energy_data if "similar_tm" in m]
    high_tm_mutants = [m for m in free_energy_data if "high_tm" in m]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    # Plot function to reduce redundancy
    def plot_category(ax, mutants, title):
        ax.set_title(title)
        for mutant in mutants:
            mutant_Q = np.array(free_energy_data[mutant]["Q_values"])
            mutant_F = np.array(free_energy_data[mutant]["F_values"])

            mutant_F_interp = np.interp(wt_Q, mutant_Q, mutant_F)
            delta_F = mutant_F_interp - wt_F
            ax.plot(wt_Q, delta_F, label=mutant)

        ax.axhline(0, linestyle="--", color="gray", linewidth=1)
        ax.set_xlabel("Fraction of Native Contacts (Q)")
        ax.legend()
        ax.grid()

    # Plot for each category
    plot_category(axes[0], low_tm_mutants, "Low tm Mutants")
    plot_category(axes[1], similar_tm_mutants, "Similar tm Mutants")
    plot_category(axes[2], high_tm_mutants, "High tm Mutants")

    axes[0].set_ylabel("ΔF(Q) (kcal/mol)")
    plt.tight_layout()
    plt.show()
