//wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "partition_function.hpp"  // Include the header for function declarations

namespace py = pybind11;

PYBIND11_MODULE(partition_function, m) {
    m.def("compute_hamiltonian", &compute_hamiltonian, "Compute the Hamiltonian for the WSME model");
    m.def("compute_W", &compute_W, "Compute the statistical weight W(m)");
    m.def("compute_partition_function_Q_SSA", &compute_partition_function_Q_SSA, "Compute partition function Q using SSA");
    m.def("compute_partition_function_Q_DSA", &compute_partition_function_Q_DSA, "Compute partition function Q using DSA");
}
