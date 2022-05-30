import numpy as np


def jacobian(dN, summits):
    """Calculate the inverse Jacobian matrix associated to the affine map F
    going from the global coordinate system to the natural coordinate
    system of the current element.
    
    Parameters
    ----------
    dN : ndarray
        shape functions derivative with respect to the natural 
        elment coordinates.
    summits : ndarray
        coordinates of the element summits in the global coordinates
        system.
        
    Returns
    -------
    Jdet : float
        determinant of the Jacobian matrix.
    Jinv : ndarray
        inverse of the Jacobian matrix.
    """
    # First, calculate the Jacobian :
    J = np.dot(dN, summits)
    
    # Inverse the Jacobian and compute its determinant :
    Jdet = np.linalg.det(J)
    Jinv = np.linalg.inv(J)
    
    return Jdet, Jinv