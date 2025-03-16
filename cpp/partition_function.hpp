//partition_function.hpp
#ifndef PARTITION_FUNCTION_HPP
#define PARTITION_FUNCTION_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

double compute_hamiltonian(const pybind11::array_t<int>& contact_map, const pybind11::array_t<int>& m_state, double contact_energy);
double compute_W(const pybind11::array_t<int>& m_state, double entropy_penalty);
pybind11::array_t<double> compute_partition_function_Q_SSA(const pybind11::array_t<int>& contact_map, double entropy_penalty, double contact_energy, double temperature);
pybind11::array_t<double> compute_partition_function_Q_DSA(const pybind11::array_t<int>& contact_map, double entropy_penalty, double contact_energy, double temperature);

#endif  // PARTITION_FUNCTION_HPP
