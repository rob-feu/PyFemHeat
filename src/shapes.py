######################################################################
#
#       Define various shapes object
#
#  Robin Feuillassier
#  2022
#
######################################################################

import numpy as np
from numpy.random import uniform
from scipy.spatial import Delaunay


class Shape:
    def __init__(self):
        print('Shape initialisaton')
        self.points = np.array([])
        
    def mesh(self):
        try:
            print('Meshing...')
            self.mesh = Delaunay(self.points)
        except IndexError:
            raise IndexError('Call discretize method before meshing the shape.')
            
            
class Line(Shape):
    """Define a line object.
    
    Parameters
    ----------
    length : float
        total line length 
    """
    def __init__(self, length):
        Shape.__init__(self)
        self.length = length
        
    def discretize(self, N):
        self.points = np.linspace(0, 1, N)
        
        
# 2-D shapes
        

class rectangle(Shape):
    def __init__(self, Lx, Ly):
        Shape.__init__(self)
        self.Lx = Lx
        self.Ly = Ly
        
    def __discretize_boundary(self, delta):
        x = np.linspace(0, self.Lx, round(self.Lx / delta)+1)
        y = np.linspace(0, self.Ly, round(self.Ly / delta)+1)
        
        # Boundaries in counter clockwise order
        bound1 = np.append(x.reshape((x.size, 1)), np.zeros((x.size, 1)), 1)
        bound2 = np.append(np.zeros((y.size, 1)) + self.Lx, y.reshape((y.size, 1)), 1)
        bound3 = np.append(x.reshape((x.size, 1)), np.zeros((x.size, 1)) + self.Ly, 1)
        bound4 = np.append(np.zeros((y.size, 1)), y.reshape((y.size, 1)), 1)
        
        bound = [bound1, bound2, bound3, bound4] 
        points = bound1[:, :]
        for i in bound[1:]:
            points = np.append(points, i, 0)
        return(np.unique(points, axis=0), bound)
    
    def __discretize_domain(self, delta):
        x = np.linspace(delta, self.Lx-delta, round((self.Lx - 2*delta) / delta)+1)
        y = np.linspace(delta, self.Ly-delta, round((self.Ly - 2*delta) / delta)+1)
        x, y = np.meshgrid(x, y)
        points = np.append(x.reshape((x.size, 1)), y.reshape((y.size, 1)), 1)
        return(points)
    
    def discretize(self, delta):
        print('Discratising...')
        self.boundary_points, self.boundary = self.__discretize_boundary(delta)
        self.surface_points = self.__discretize_domain(delta)
        self.points = np.append(self.boundary_points, self.surface_points, 0)

        
class disk(Shape):
    def __init__(self, radius):
        Shape.__init__(self)
        self.radius = radius
        
    def __discretize_boundary(self, N):
        points = np.zeros((N, 2))
        theta = np.linspace(0, 2*np.pi, N)
        points[:, 0] = self.radius * np.cos(theta)  # x
        points[:, 1] = self.radius * np.sin(theta)  # y
        return(points)
    
    def __discretize_domain(self, N):
        points = uniform(-self.radius, self.radius, (N, 2))
        points = points[np.where(points[:, 0]**2+points[:, 1]**2 < self.radius**2)]
        return(points)
    
    def discretize(self, Nbound, Nsurf):
        print('Discratising...')
        self.boundary_points = self.__discretize_boundary(Nbound)
        self.surface_points = self.__discretize_domain(Nsurf)
        self.points = np.append(self.boundary_points, self.surface_points, 0)
        
        
# 3-D shapes
        

class Cylinder(Shape):
    """Define a cylinder object 
    The cylinder is oriented along the z axis
    
    Parameters
    ----------
    radius : float
    length : float
    nbpoints_vol : float
    nbpoints_surf : float
    """
    def __init__(self, radius, length, nbpoints_vol):
        Shape.__init__(self)
        self.radius = radius
        self.length = length
        self.nbpoints_vol = nbpoints_vol
        
        
    def mesh_boundaries(self, nbpoints_surf):
        self.nbpoints_surf = nbpoints_surf
        