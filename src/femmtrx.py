import numpy as np

from utils.quadrature import triangleGauss
from utils.shapefunctions import triN
from utils.jacobian import jacobian


"""Functions to calculate matrix and vector involved in the heat
equation for Dirichlet boundary conditions.
"""


def steady_heat_matrix(domain, conductivity, rho, cv, order=1):
    """Functions to calculate M, K matrices and f involved in the 
    time-independent heat equation without sources/sink :
                K T = 0
            ->  Ku T = Kbc Tbar
        with T the temperature at each node to be determined and
        Tbar the temperature fixed by boundary conditions.
    
    Parameters
    ----------
    domain : shapes instance
        mesh : ndarray
            Delaunay triangular mesh.
        points : ndarray
            coordinates of each node (N, 2).
    bound_cond : ndarray
        boundary condition : points coordinates and values.
    conductivity : ndarray
        conductivity matrix
    rho : float
        material density.
    cv : float
        material specific heat.
    order : int
        quadrature order
    
    Returns
    -------
    K : ndarray
        "stiffness" matrix
    f : ndarray
        boundary conditions vector
    """
    dim = np.shape(domain.points)[0]
    K = np.zeros((dim, dim))
    f = np.zeros((dim, 1))
    quad_points, quad_weights = triangleGauss(order)
    
    # Calculate global K #############################################
    
    # Iterate over elements
    for e, element in enumerate(domain.mesh.simplices):
        summits = domain.points[element]  # Summits coordinates
        Ke = np.zeros((len(element), len(element)))
        Me = np.zeros((len(element), len(element)))
        
        # Iterate over quadrature points
        for p, [x, y] in enumerate(quad_points):
            shapef, dN = triN(x, y)  # Shape functions and derivatives of the current element, evaluated on (x, y) being the quadrature point (in isoparametric/natural coordinates).
            
            Jdet, Jinv = jacobian(dN, summits)
            
            glob_dN = np.dot(Jinv, dN)  # Shape functions derivative with respect to global coordinates.
            # Quadrature factor
            quad_factor = quad_weights[p] * Jdet
            
            # Iterate over element's summits
            for i, indi in enumerate(element):
                for j, indj in enumerate(element):
                    dNi = glob_dN[:, i]
                    dNj = glob_dN[:, j]
                    K[indi, indj] += quad_factor * np.dot(np.transpose(dNi), np.dot(conductivity, dNj))
                    Ke += K[indi, indj]
    
    # Assemble Ku and Kbc ############################################
    
    Ku = np.zeros((np.shape(domain.surface_points)[0], 
                   np.shape(domain.surface_points)[0]))
    Kbc = np.zeros((np.shape(domain.surface_points)[0], 
                    np.shape(domain.boundary_points)[0]))
    u_node_nb = np.zeros((np.shape(domain.surface_points)[0], 2))
    bc_node_nb = np.zeros((np.shape(domain.boundary_points)[0], 2))
    ku1, ku2 = 0, 0  # Counter for Ku
    kbc2 = 0  # Counter for Kbc
    points = domain.points[:, 0] + 1j*domain.points[:, 1]
    surface_points = domain.surface_points[:, 0] + 1j*domain.surface_points[:, 1]
    for i in np.unique(domain.mesh.simplices.flatten()):
        if points[i] in surface_points:
            # print("row: ", domain.points[i])
            u_node_nb[ku1] = domain.points[i]
            for j in np.unique(domain.mesh.simplices.flatten()):
                if points[j] in surface_points:
                    # print("Ku column: ", domain.points[j])
                    # print(points[j] in surface_points)
                    Ku[ku1, ku2] = K[i, j]
                    ku2 += 1
                else:
                    # print("Kbc column: ", domain.points[j])
                    Kbc[ku1, kbc2] = K[i, j]
                    bc_node_nb[kbc2] = domain.points[j]
                    kbc2 += 1
            ku1 += 1
            ku2 = 0
            kbc2 = 0
        
    return Ku, Kbc, u_node_nb, bc_node_nb

    

def time_heat_matrix(domain, conductivity, rho, cv, order=1):
    """Functions to calculate M, K matrices and f involved in the 
    time-dependent heat equation without sources/sink :
                M (dT/dt) + K T = f
    
    Parameters
    ----------
    domain : shapes instance
        mesh : ndarray
            
        points : ndarray
            coordinates of each points (N, 2).
    bound_cond : ndarray
        boundary condition : points coordinates and values.
    conductivity : ndarray
        conductivity matrix
    rho : float
        material density.
    cv : float
        material specific heat.
    order : int
        quadrature order
    
    Returns
    -------
    K : ndarray
    M : 
    f : 
    """
    dim = np.shape(domain.points)[0]
    K = np.zeros((dim, dim))
    M = np.zeros((dim, dim))
    f = np.zeros((dim, 1))
    quad_points, quad_weights = triangleGauss(order)
    
    # Iterate over elements
    for element in domain.mesh.simplices:
        summits = domain.points[element]  # Summits coordinates
        Ke = np.zeros((len(element), len(element)))
        Me = np.zeros((len(element), len(element)))
        
        # Iterate over quadrature points
        for p, [x, y] in enumerate(quad_points):
            shapef, dN = triN(x, y)  # Shape functions and derivatives of the current element, evaluated on (x, y) being the quadrature point (in isoparametric/natural coordinates).
            Jdet, Jinv = jacobian(dN, summits)
            
            glob_dN = np.dot(Jinv, dN)  # Shape functions derivative with respect to global coordinates.
            # Quadrature factor
            quad_factor = quad_weights[p] * Jdet
            
            # Iterate over element's summits
            for i, indi in enumerate(element):
                for j, indj in enumerate(element):
                    dNi = glob_dN[:, i]
                    dNj = glob_dN[:, j]
                    K[indi, indj] += quad_factor * np.dot(np.transpose(dNi), np.dot(conductivity, dNj))
                    M[indi, indj] += quad_factor * rho * cv * shapef[i] * shapef[j]
                    Ke += K[indi, indj]
                    Me += M[indi, indj]
        
        #if summits in bound_cond[:, :2]:
         #   for indi in element:
          #      for indj in element:
           #         f[indi] += Ke[0] * 0
    return K, M


    