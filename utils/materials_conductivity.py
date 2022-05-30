from numpy import diag

def conductivity_matrix(ax_conductivity):
    """Thermal conductivity matrix
    
    Output dimension depends on the input dimension.
    Units: W/m/K
    
    Parameters
    ----------
    ax_conductivity : tuple
        diagonal elements
        
            kxx : float
                conductivity along x.
            kyy : float
                conductivity along y.
            kzz : float
                conductivity along z.
    
    Returns
    -------
    K : ndarray
        Thermal conductivity matrix
    """
    K = diag(ax_conductivity)
    return K
        
    