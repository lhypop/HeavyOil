"""
FastPeriodicDBSCAN.py

Implementation of FastPeriodicDBSCAN, a DBSCAN clustering algorithm variant 
that efficiently handles periodic boundary conditions (PBC) commonly encountered 
in molecular simulations and other spatial datasets with periodicity.

The FastPeriodicDBSCAN class inherits from sklearn's DBSCAN, overrides the fitting 
process to compute pairwise distances using a vectorized periodic boundary distance 
function, enabling clustering in periodic spaces.

Author: Haiyang Li
"""


import numpy as np
from sklearn.cluster import DBSCAN

class FastPeriodicDBSCAN(DBSCAN):
    """
    A fast DBSCAN implementation that supports periodic boundary conditions (PBC).
    
    Parameters:
        eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood.
        min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point.
        box_size (array-like or None): Size of the simulation box along each dimension. If None, it will be inferred.
        **kwargs: Additional keyword arguments passed to sklearn's DBSCAN.
    """

    def __init__(self, eps: float = 0.5, min_samples: int = 5, box_size=None, **kwargs):
        super().__init__(eps=eps, min_samples=min_samples, **kwargs)
        self.box_size = box_size

    def fit(self, X: np.ndarray, y=None, sample_weight=None):
        """
        Fit the DBSCAN algorithm using periodic boundary conditions.

        Parameters:
            X (np.ndarray): Input coordinates of shape (n_samples, n_features).
            y: Ignored, for API compatibility.
            sample_weight: Optional sample weights.

        Returns:
            self: Fitted estimator.
        """
        # Check input dimensionality
        if X.ndim != 2:
            raise ValueError("Input X must be a 2D array with shape (n_samples, n_features)")

        # Automatically infer box size if not provided
        if self.box_size is None:
            self.box_size = np.max(X, axis=0) - np.min(X, axis=0)

        # Compute distance matrix with periodic boundary conditions
        dist_matrix = self._periodic_distance_vectorized(X, self.box_size)

        # Use the precomputed distance matrix in DBSCAN
        self.metric = 'precomputed'
        return super().fit(dist_matrix, y, sample_weight)
    

    @staticmethod
    def _periodic_distance_vectorized(X: np.ndarray, box_size) -> np.ndarray:
        """
        Compute the pairwise distance matrix under periodic boundary conditions (PBC).

        Parameters
        ----------
        X : np.ndarray
            Coordinates of shape (n_samples, n_dims).
        box_size : float or np.ndarray
            Box size along each dimension.

        Returns
        -------
        np.ndarray
            Distance matrix of shape (n_samples, n_samples).
        """
        box_size = np.asarray(box_size)
        if box_size.ndim == 0:
            box_size = np.full(X.shape[1], box_size)

        delta = X[:, np.newaxis, :] - X[np.newaxis, :, :]
        delta = np.abs(delta)
        delta = np.where(delta > 0.5 * box_size, box_size - delta, delta)
        return np.sqrt(np.sum(delta ** 2, axis=-1))