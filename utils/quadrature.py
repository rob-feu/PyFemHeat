"""Utility functions for quadrature formula
"""


import numpy as np


def triangleGauss(order=1):
    """Gaussian quadrature's points and associated weights for 2-D 
    isoparametric triangular elements. 
    
    Isoparametric triangle have summits defined as:
    (0, 0) - (1, 0) - (0, 1)
    
    Returns
    -------
    points : ndarray
        quadrature points coordinates in counter clockwise order.
    weights : ndarray
        associated weights.
    """
    # Central point
    if order == 1:
        points = np.array([[1/3, 1/3]])
        weights = np.array([1])
    # Midpoints on segments connecting summits and central points
    elif order == 2:
        points = np.array([[1/3/2, 1/3/2],
                           [2/3, 1/3/2],
                           [1/3/2, 2/3]])
        weights = np.array([[1/3],
                            [1/3],
                            [1/3]])
    return points, weights


def tetraGauss(order=1):
    """Gaussian quadratures's points and weight for 3-D isoparametric
     tetrahedral elements.
     
    Isoparametric tetrahedron have summits defined as:
    (0, 0, 0) - (1, 0, 0) - (0, 1, 0) - (0, 0, 1)
    
    Returns
    -------
    points : ndarray
        quadrature points coordinates in counter clockwise order.
    weights : ndarray
        associated weights.
     """
    if order == 1:
        points = np.array([[1/8, 1/8]])