# WSME-ts

## Background and Motivation

Understanding protein stability and folding is a fundamental challenge in computational biology. While deep learning models have achieved remarkable success in structure prediction, physics-based models remain essential for investigating thermodynamic properties and stability. One such approach is the Wako-Saitô-Muñoz-Eaton (WSME) model [1], a statistical mechanical framework that computes free-energy landscapes from a protein’s native structure. The WSME model provides valuable insights into folding pathways, intermediates, and stability determinants with relatively low computational cost.

The primary goal of this study is to utilize the WSME model to predict the thermostability of single-point mutants. This is achieved by comparing the free-energy landscapes of wild-type (WT) and mutant proteins to evaluate stability changes. The central hypothesis is that mutants with free-energy profiles similar to the WT will maintain comparable thermostability, whereas significant deviations in the landscape may lead to either stabilization or destabilization. This hypothesis is grounded in the fundamental principles of protein folding thermodynamics and energy landscape theory, which state that protein stability is dictated by the free-energy difference between its folded (native) and unfolded states. Mutations can reshape the free-energy landscape, thereby altering the relative stability of these states. If a mutant exhibits a free-energy landscape closely resembling that of the WT, it implies that the mutation has little impact on stability. In contrast, mutations that substantially alter the landscape, such as by modifying the free-energy barrier, can lead to increased or decreased stability. This concept is well-established in statistical mechanical models of protein folding, including the WSME model, where the free-energy profile reflects the probability distribution of different folding states.

I implemented **WSME-ts**, a simple physics-based approach to predict the thermostability (as measured by melting point $t_M$) of single-point mutants by comparing their free-energy landscapes with that of the WT protein. To achieve this, we utilized a publicly available dataset containing experimental melting temperature values of single mutants, coupled with protein structures predicted using [ESMFold](https://esmatlas.com/resources?action=fold).

---

## Theoretical Foundation: WSME Model

The WSME model is a coarse-grained statistical mechanical approach that describes protein folding as a cooperative process. It assumes that:

- Each residue in the protein exists in either a native-like or unfolded state.

- Native interactions contribute to stability only if all intervening residues between contacting pairs are folded.

- Folding is driven by a balance between enthalpic stabilization (native contacts) and entropic penalties (loss of conformational freedom).

The **Hamiltonian** of the WSME model is given by:

$\hspace{4cm} H(\{m\}) = \sum_{i<j} \epsilon_{ij} \Delta_{ij} m_{i,j}$

where:
- $ \epsilon_{ij} $ is the contact energy for residue pair $ (i,j) $,
- $ \Delta_{ij} $ is 1 if residues $ i $ and $ j $ are in contact in the native state, otherwise 0,
- $ m_{i,j} $ represents cooperative folding, meaning contacts contribute only when all intermediate residues are folded.

The **partition function** is computed as:

$\hspace{4cm}  Z = \sum_{\{m\}} W(\{m\}) e^{-\beta H(\{m\})}$

where $ \beta = \frac{1}{k_B T} $ (Boltzmann factor), and $ W(\{m\}) $ represents the entropy of the system. The partition function  sums over all conformational states of the system. However, in the WSME model, we group conformations are based on their fraction of native contacts $Q$:

$\hspace{4cm}  Z(Q) = \sum_{\{m \mid Q(m) = Q\}} W(m) e^{-\beta H(m)}$,

where $Z(Q)$ is the restricted partition function that considers only states with a specific fraction $Q$ of native contacts.

Under this assumption, the **free-energy landscape** as function of the partition function $Z$ is then derived from the **Helmholtz free energy**:

$\hspace{4cm} F(Q) = - k_B T \ln \frac{Z(Q)}{Z_{\text{total}}}$

where the scaled $\frac{Z(Q)}{Z_{\text{total}}}$ represents the probability of finding the system in a state with a specific native contact fraction $Q$.

---

## Computational Implementation

I implemented the **WSME model** using a **single sequence approximation (SSA)**, and **double sequence approximation (DSA)** (details can be found on the review paper [1]), which allows for efficient computation of the partition function while considering native contacts. The workflow is as follows:
1. **Extract alpha-carbon coordinates from the protein structure** (PDB file).
2. **Compute the native contact map** based on a distance cutoff (used 6Å).
3. **Enumerate folding states** under SSA/DSA and compute the partition function.
4. **Calculate the free-energy landscape**  $F(Q)$ for the WT protein.
5. **Repeat steps 1-4 for each single-point mutant**.
6. **Compare free-energy landscapes of WT and mutants** to predict changes in thermostability.

I used the [**Novozymes Kaggle dataset**](https://www.kaggle.com/code/tranminhthuan/novozymes-eda-modelling-protbert-xgboost) for gettinf the single-point mutant sequences, and their corresponding **$t_m$ values**, and generated mutant structures using [**ESMFold**](https://esmatlas.com/resources?action=fold). The Novozymes dataset includes the sequences and experimentally measured $t_m$ values of single mutants, allowing for model validation.


---

## Application to Predicting Mutant Thermostability

To test our hypothesis, I selected mutants with:
- Lowest $t_m$ (high destabilization)
- Highest $t_m$ (high stabilization)
- $t_m$ similar to WT (neutral effect)

By **comparing their free-energy landscapes to that of the wildtype**, I assessed how well the WSME model captures the energetic shifts induced by mutations. My approach provides a simple, physics-based method to **predict thermostability changes** and can be used as a complementary tool for protein engineering and stability optimization.

---

## References

[1] Ooka K, Liu R, Arai M. The Wako-Saitô-Muñoz-Eaton Model for Predicting Protein Folding and Dynamics. Molecules. (2022) Jul 12;27(14):4460. doi: 10.3390/molecules27144460.
