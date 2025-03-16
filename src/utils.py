from typing import Tuple

import numpy as np
from Bio import PDB


def extract_C_alpha(
    structure: PDB.Structure, model_index=0, chain_id: str = None
) -> Tuple[np.ndarray, list]:
    """Extracts alpha-carbon (Cα) coordinates from a specific model of a protein structure.

    Parameters:
    -----------
    structure : PDB.Structure
        A Biopython Structure object containing the protein structure.
    model_index : int, optional (default=0)
        The index of the model to extract Cα atoms from (useful for multi-model structures).
    chain_id : str, optional (default=None)
        If specified, only extracts Cα atoms from the given chain ID. If None, extracts from all chains.

    Returns:
    --------
    np.ndarray
        A NumPy array of shape (N, 3) containing the Cα atom coordinates.
    list
        A list of Residue objects corresponding to the extracted Cα atoms.

    Notes:
    ------
    - This function assumes that the input structure follows the Biopython PDB module format.
    - The function skips residues that do not contain a Cα atom (e.g., non-standard residues).
    """
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


def compute_contact_map(coords: np.ndarray, cutoff: float) -> np.ndarray:
    """Computes the contact map for a given set of atomic coordinates.

    Parameters:
    -----------
    coords : np.ndarray
        A NumPy array of shape (N, 3) representing the atomic coordinates of N residues.
    cutoff : float
        The distance threshold (in Å) below which two residues are considered to be in contact.

    Returns:
    --------
    np.ndarray
        A binary contact map of shape (N, N), where `contact_map[i, j] = 1` if residues `i` and `j`
        are within the cutoff distance, otherwise `0`.

    Notes:
    ------
    - The contact map is symmetric since distances are undirected.
    - The function enforces a minimum sequence separation of 2 residues (`i + 2 < j`)
      to exclude trivial contacts between consecutive residues.
    """
    N = len(coords)
    contact_map = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 2, N):
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < cutoff:
                contact_map[i, j] = 1
                contact_map[j, i] = 1

    return contact_map
