"""
Robin Feullassier

Compute the linear coefficient of a first order shape function.

"""

import numpy as np
from numpy.linalg import inv

"""Non-isoparametric shape functions"""
class ShapeFunctions:
    def __init__(self, shape):
        self.points = shape.points
        self.delaunay = shape.mesh
        
        self.Nelements = np.shape(shape.mesh.simplices)[0]  # Number of elements
        self.Nnodes = np.shape(shape.mesh.simplices)[1]  # Number of nodes in an element
        dim = np.shape(shape.points)[1]
        self.coeff = np.zeros((self.Nelements, self.Nnodes, dim+1))
        
    def __coeff_N(self, i, element):
        # Basis function coefficients of node i in the current element
        #
        # For a tetrahedra: size(element) = (4, 3)
        # i goes from 0 to 3
    
        nodes = np.shape(element)[0]  # Number of nodes in the element
    
        S = np.append(element, np.ones((nodes, 1)), 1) 
        ibasis_vector = np.zeros((nodes, 1))
        ibasis_vector[i] = 1
    
        N = np.dot(np.linalg.inv(S), ibasis_vector)
        return(N)
        
    def calc_basis(self):
        # Iterate over mesh elements to compute all coeff
        for element in points[self.delaunay.simplices]:
            for node, coord in enumerate(element):
                self.coeff[element, node, :] = self.__coeff_N(node, element)
                
                
    def analytic_tri(self, j, element):
        # Analytical basis functions for 2-D triangle
        
        xj = element[j, 0]
        yj = element[j, 1]
        
        removed = np.delete(element, j, 0)
        xi = removed[0, 0]
        yi = removed[0, 1]
        xk = removed[1, 0]
        yk = removed[1, 1]
        
        alpha_j = (yk - yi) / ( (xj-xk)(yj-ti) - (xj-xi)(yi-yk) )
        beta_j = (xi - xk) / ( (xj-xk)(yj-ti) - (xj-xi)(yi-yk) )
        gamma_j = (xk*yi - xi*yk) / ( (xj-xk)(yj-ti) - (xj-xi)(yi-yk) )
        
        return(alpha_j, beta_j, gamma_j)
    

def triN(x, y):
    """Linear shape functions for isoparametric triangle
    
    Parameters
    ----------
    x : float
        horizontal coordinate.
    y : float
        vertical coordinate.
        
    Returns
    -------
    N : ndarray
        shape functions for each nodes evaluated at (x, y).
    dN : ndarray
        shape functions derivatives at (x, y).
    """
    # N[0] -> shape function for node (0, 0).
    # N[1] -> shape function for node (1, 0).
    # N[2] -> shape function for node (0, 1).
    N = np.array([1 - x - y, x, y])
    
    # dN[0, :] -> shape functions derivative wih respect to x.
    # dN[1, :] -> shape functions derivative wih respect to y.
    dN = np.array([[-1, 1, 0],
                   [-1, 0, 1]])
    
    return N, dN