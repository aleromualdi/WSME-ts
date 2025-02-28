# Wako-Saitô-Muñoz-Eaton (WSME) model for protein mutant thermostability prediction

This code implements **Wako-Saitô-Muñoz-Eaton-Thermostability** (WSME-ts), a simple physics-based approach to predict the thermostability (as measured by melting point $t_m$) of single-point mutants by comparing their free-energy landscapes with that of the WT protein.

## Background and Motivation

Understanding protein stability and folding is a fundamental challenge in computational biology. While deep learning models have achieved remarkable success in structure prediction, physics-based models remain essential for investigating thermodynamic properties and stability. One such approach is the Wako-Saitô-Muñoz-Eaton (WSME) model [1], a statistical mechanical framework that computes free-energy landscapes from a protein’s native structure. The WSME model provides valuable insights into folding pathways, intermediates, and stability determinants with relatively low computational cost.

The primary goal of this study is to utilize the WSME model to predict the thermostability of single-point mutants. This is achieved by comparing the free-energy landscapes of wild-type (WT) and mutant proteins to evaluate stability changes. The central hypothesis is that mutants with free-energy profiles similar to the WT will maintain comparable thermostability, whereas significant deviations in the landscape may lead to either stabilization or destabilization. This hypothesis is grounded in the fundamental principles of protein folding thermodynamics and energy landscape theory, which state that protein stability is dictated by the free-energy difference between its folded (native) and unfolded states. Mutations can reshape the free-energy landscape, thereby altering the relative stability of these states. If a mutant exhibits a free-energy landscape closely resembling that of the WT, it implies that the mutation has little impact on stability. In contrast, mutations that substantially alter the landscape, such as by modifying the free-energy barrier, can lead to increased or decreased stability. To demonstrate this, I apply the implemented model to a publicly available dataset that includes experimental melting temperature values for single mutants, whose structure I predict using [ESMFold](https://esmatlas.com/resources?action=fold).

---

## Theoretical Foundation: WSME Model

The WSME model is a coarse-grained statistical mechanical approach that describes protein folding as a cooperative process. It assumes that:

- Each residue in the protein exists in either a native-like or unfolded state.

- Native interactions contribute to stability only if all intervening residues between contacting pairs are folded.

- Folding is driven by a balance between enthalpic stabilization (native contacts) and entropic penalties (loss of conformational freedom).

The **Hamiltonian** of the WSME model is given by:

$\hspace{4cm} H(\{m\}) = \sum_{i<j} \epsilon_{ij} \Delta_{ij} m_{i,j}$

where:
- $\epsilon_{ij}$ is the contact energy for residue pair $(i,j)$,
- $\Delta_{ij}$ is 1 if residues $i$ and $j$ are in contact in the native state, otherwise 0,
- $m_{i,j}$ represents cooperative folding, meaning contacts contribute only when all intermediate residues are folded.

The **partition function** is computed as:

$\hspace{4cm}  Z = \sum_{\{m\}} W(\{m\}) e^{-\beta H(\{m\})}$

where $\beta = \frac{1}{k_B T}$ (Boltzmann factor), and $W(\{m\})$ represents the entropy of the system. The partition function  sums over all conformational states of the system. However, in the WSME model, we group conformations are based on their fraction of native contacts $Q$:

$\hspace{4cm}  Z(Q) = \sum_{\{m \mid Q(m) = Q\}} W(m) e^{-\beta H(m)}$,

where $Z(Q)$ is the restricted partition function that considers only states with a specific fraction $Q$ of native contacts.

Under this assumption, the **free-energy landscape** as function of the partition function $Z$ is then derived from the **Helmholtz free energy**:

$\hspace{4cm} F(Q) = - k_B T \ln \frac{Z(Q)}{Z_{\text{total}}}$

where the scaled $\frac{Z(Q)}{Z_{\text{total}}}$ represents the probability of finding the system in a state with a specific native contact fraction $Q$.

---

## Computational Implementation

- In `src/` I implemented the **WSME model** using a **single sequence approximation (SSA)**, and **double sequence approximation (DSA)** (details can be found on the review paper [1]), which allows for efficient computation of the partition function.

- In `single-point-mutants.ipynb` I implement the workflow to demosntrate how to use the WSME-ts model to predict thermostability of single-point mutants. Here, I used the [**Novozymes Kaggle dataset**](https://www.kaggle.com/code/tranminhthuan/novozymes-eda-modelling-protbert-xgboost) for getting the single-point mutant sequences, and their corresponding **$t_m$ values**, and generated mutant structures using [**ESMFold**](https://esmatlas.com/resources?action=fold). The Novozymes dataset includes the sequences and experimentally measured $t_m$ values of single mutants, allowing for model validation.


## References

[1] Ooka K, Liu R, Arai M. The Wako-Saitô-Muñoz-Eaton Model for Predicting Protein Folding and Dynamics. Molecules. (2022) Jul 12;27(14):4460. doi: 10.3390/molecules27144460.
