from typing import List, Tuple

import numpy as np
from numba import njit, prange

kB = 1.987e-3  # kcal/mol·K


@njit
def compute_hamiltonian(
    contact_map: np.ndarray, m_state: np.ndarray, contact_energy: float
) -> float:
    """Computes the Hamiltonian (H) for the WSME model with m_{i,j} condition.

    Parameters:
        contact_map (np.ndarray): NxN matrix representing native contacts (1 = contact, 0 = no contact).
        m_state (np.ndarray): Binary array (1 = folded, 0 = unfolded) for each residue.
        contact_energy (float): Energy assigned to each contact (default = -1.5 kcal/mol).

    Returns:
        H (float): Hamiltonian value.
    
    """
    N = contact_map.shape[0]
    H = 0.0

    # Iterate over all pairs (i, j) where i < j
    for i in range(N - 1):
        for j in range(i + 1, N):
            if contact_map[i, j] == 1:  # Check if native contact exists
                m_ij = np.sum(m_state[i : j + 1]) / (j - i + 1)
                H += (
                    contact_energy * m_ij
                )  # Only contribute if all residues in range are folded

    return H


def compute_W(m_state: np.ndarray, entropy_penalty: float) -> float:
    """Computes statistical weight W(m) using SSA.

    Parameters:
        m_state (np.ndarray): Binary array (1 = folded, 0 = unfolded).
        entropy_penalty (float): Entropy penalty per residue.

    Returns:
        W_m (float): Weighting factor.
    """

    S_total = np.sum(entropy_penalty * m_state)
    exponent = S_total / kB
    return np.exp(exponent)


@njit
def compute_partition_function_Q_SSA(
    contact_map: np.ndarray,
    entropy_penalty: float,
    contact_energy: float,
    temperature: float = 310,
) -> np.ndarray:
    """Computes restricted partition function Z(Q) using SSA.

    Parameters:
        contact_map (np.ndarray): NxN matrix for native contacts.
        entropy_penalty (float): Entropy penalty per residue.
        contact_energy (float): Energy per contact.
        temperature (float): Temperature in Kelvin.

    Returns:
        Z_Q (np.ndarray): Restricted partition function per Q.
    """

    N = contact_map.shape[0]
    beta = 1 / (kB * temperature)

    Z_Q = np.zeros(N + 1)  # Initialize partition function storage

    # Iterate over all single contiguous folding segments
    for i in range(N):
        for j in range(i, N):
            m_state = np.zeros(N)
            m_state[i : j + 1] = 1  # Single folding region

            W_m = compute_W(m_state, entropy_penalty)
            H_m = compute_hamiltonian(contact_map, m_state, contact_energy)

            Q = int(np.sum(m_state))  # Fraction of folded residues

            Z_Q[Q] += W_m * np.exp(-beta * H_m)

    # Normalize using vectorized NumPy operation
    Z_total = np.sum(Z_Q)
    if Z_total > 0:
        Z_Q /= Z_total

    return Z_Q


@njit(parallel=True)
def compute_partition_function_Q_DSA(
    contact_map: np.ndarray,
    entropy_penalty: float,
    contact_energy: float,
    temperature=310,
) -> np.ndarray:
    """Computes restricted partition function Z(Q) using Double Sequence Approximation (DSA)
    with parallelization.
    
    Parameters:
        contact_map (np.ndarray): NxN matrix for native contacts.
        entropy_penalty (float): Entropy penalty per residue.
        contact_energy (float): Energy per contact.
        temperature (float): Temperature in Kelvin.

    Returns:
        Z_Q (np.ndarray): Restricted partition function per Q.
    
    """
    N = contact_map.shape[0]
    beta = 1 / (kB * temperature)
    Z_Q = np.zeros(N + 1)

    for i in prange(N):
        if i % 10 == 0:
            print(f"Processing i={i}")

        for j in prange(i, N):
            m_state = np.zeros(N)
            m_state[i : j + 1] = 1

            W_m = np.exp(np.sum(entropy_penalty * m_state) / kB)
            H_m = compute_hamiltonian(contact_map, m_state, contact_energy)
            Q = int(np.sum(m_state))
            Z_Q[Q] += W_m * np.exp(-beta * H_m)

            for k in prange(j + 2, N):
                for m in prange(k, N):
                    m_state_2 = m_state.copy()
                    m_state_2[k : m + 1] = 1

                    W_m2 = np.exp(np.sum(entropy_penalty * m_state_2) / kB)
                    H_m2 = compute_hamiltonian(contact_map, m_state_2, contact_energy)
                    Q_2 = int(np.sum(m_state_2))
                    Z_Q[Q_2] += W_m2 * np.exp(-beta * H_m2)

    # Normalize the partition function
    Z_total = np.sum(Z_Q)
    Z_Q /= Z_total
    return Z_Q


def compute_free_energy(
    Z_Q: np.ndarray, temperature: float = 310
) -> Tuple[List[float], List[float]]:
    """Computes the free-energy landscape F(Q).

    Parameters:
        Z_Q (np.ndarray): Restricted partition function Z(Q) stored as a NumPy array.
        temperature (float): Temperature in Kelvin (default = 310K).

    Returns:
        Q_values (list): Native contact fractions.
        F_values (list): Free energy per Q.
    """
    Z_total = np.sum(Z_Q)
    F_values = []
    Q_values = []

    for q in range(len(Z_Q)):
        if Z_Q[q] > 0:
            F_q = -kB * temperature * np.log(Z_Q[q] / Z_total)
            F_values.append(F_q)
            Q_values.append(q / (len(Z_Q) - 1))

    # Shift the free energy scale so the unfolded state (Q=0) is set to 0
    max_FQ = max(F_values)
    F_values = [F - max_FQ for F in F_values]

    return Q_values, F_values
