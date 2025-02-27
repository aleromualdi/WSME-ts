import numpy as np
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr, spearmanr


def compute_rmsd(F_wt, F_mutant):
    """Computes RMSD between wild-type and mutant free energy profiles."""
    return np.sqrt(np.mean((np.array(F_wt) - np.array(F_mutant)) ** 2))


def compute_pearson(F_wt, F_mutant):
    """Computes Pearson correlation between wild-type and mutant free energy profiles."""
    return pearsonr(F_wt, F_mutant)[0]  # Only return correlation coefficient


def compute_spearman(F_wt, F_mutant):
    """Computes Spearman correlation between wild-type and mutant free energy profiles."""
    return spearmanr(F_wt, F_mutant)[0]  # Only return correlation coefficient


def compute_dtw(F_wt, F_mutant):
    """Computes Dynamic Time Warping (DTW) distance between two free energy profiles."""
    # Ensure inputs are properly formatted as 1D numpy arrays
    F_wt = np.asarray(F_wt).flatten()
    F_mutant = np.asarray(F_mutant).flatten()

    distance, _ = fastdtw(
        F_wt[:, None], F_mutant[:, None], dist=euclidean
    )  # Convert to Nx1 shape for compatibility
    return distance
