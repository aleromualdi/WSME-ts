// partition_function.cpp
#include "partition_function.hpp"
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>
#include <algorithm>
#include <pybind11/iostream.h>  // Include this for pybind11 print support

namespace py = pybind11;

const double kB = 1.987e-3;  // kcal/mol·K

// Compute the Hamiltonian for the WSME model
double compute_hamiltonian(
    const pybind11::array_t<int>& contact_map, 
    const pybind11::array_t<int>& m_state, 
    double contact_energy) {
    
    auto cm = contact_map.unchecked<2>();  // Access the contact map (NxN)
    auto ms = m_state.unchecked<1>();      // Access the m_state array (1D)
    int N = cm.shape(0);  // Assume square contact map

    double H = 0.0;
    // Iterate over all pairs (i, j) where i < j
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (cm(i, j) == 1) {  // Check if contact exists
                // Compute m_ij (fraction of folded residues between i and j)
                double m_ij = 0.0;
                for (int k = i; k <= j; ++k) {
                    m_ij += ms[k];
                }
                m_ij /= (j - i + 1);
                H += contact_energy * m_ij;  // Add to Hamiltonian
            }
        }
    }

    return H;
}

// Compute the statistical weight W(m)
double compute_W(const pybind11::array_t<int>& m_state, double entropy_penalty) {
    auto ms = m_state.unchecked<1>();  // Access m_state array (1D)
    double S_total = 0.0;
    for (ssize_t i = 0; i < ms.shape(0); ++i) {
        S_total += entropy_penalty * ms[i];
    }
    double exponent = S_total / kB;
    return exp(exponent);
}


// Compute the partition function Q using Double Sequence Approximation (DSA)
pybind11::array_t<double> compute_partition_function_Q_DSA(
    const pybind11::array_t<int>& contact_map, 
    double entropy_penalty, 
    double contact_energy, 
    double temperature) {

    auto cm = contact_map.unchecked<2>();  // Access contact map (NxN)
    int N = cm.shape(0);
    double beta = 1.0 / (kB * temperature);
        
    pybind11::array_t<double> Z_Q(N + 1);
    auto Z = Z_Q.mutable_unchecked<1>();  // Access the mutable Z_Q array (1D)

    pybind11::array_t<int> m_state(N);  // Binary array for folded/unfolded state
    auto ms = m_state.mutable_unchecked<1>();  // Access the mutable m_state array (1D)

    // Loop over possible starting points i
    for (int i = 0; i < N; ++i) {
        // Loop over ending points j for the first sequence
        for (int j = i; j < N; ++j) {
            // Set folded residues between i and j (manual fill instead of std::fill)
            for (int k = 0; k < N; ++k) {
                ms[k] = 0;
            }
            for (int k = i; k <= j; ++k) {
                ms[k] = 1;  // Mark residues i to j as folded
            }

            // Calculate the statistical weight and Hamiltonian for the folding state
            double W_m = compute_W(m_state, entropy_penalty);
            double H_m = compute_hamiltonian(contact_map, m_state, contact_energy);
            
            // Calculate Q as the number of folded residues
            int Q = 0;
            for (int k = 0; k < N; ++k) {
                Q += ms[k];
            }

            Z[Q] += W_m * exp(-beta * H_m);

            // Now consider extending the sequence with a second sequence from j+2 to N
            for (int k = j + 2; k < N; ++k) {
                for (int m = k; m < N; ++m) {
                    // Reuse m_state and extend the folded region
                    for (int k = 0; k < N; ++k) {
                        ms[k] = 0;  // Reset to all unfolded
                    }
                    for (int l = k; l <= m; ++l) {
                        ms[l] = 1;  // Mark residues k to m as folded
                    }

                    double W_m2 = compute_W(m_state, entropy_penalty);
                    double H_m2 = compute_hamiltonian(contact_map, m_state, contact_energy);
                    
                    // Calculate Q_2 as the number of folded residues in state 2
                    int Q_2 = 0;
                    for (int l = 0; l < N; ++l) {
                        Q_2 += ms[l];
                    }

                    Z[Q_2] += W_m2 * exp(-beta * H_m2);
                }
            }
        }
    }

    // Normalize the partition function
    double Z_total = 0.0;
    for (int q = 0; q < Z_Q.shape(0); ++q) {
        Z_total += Z[q];
    }

    if (Z_total > 0) {
        for (int q = 0; q < Z_Q.shape(0); ++q) {
            Z[q] /= Z_total;
        }
    }

    return Z_Q;
}
