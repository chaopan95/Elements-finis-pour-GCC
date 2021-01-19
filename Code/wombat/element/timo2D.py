# -*- coding: utf-8 -*-
"""
Module for class `Beam2D`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *

class Timo2D(Segment):
    """ A 2D Timoshenko beam element
    
        Timo2D is a Segment-type element (2 end nodes) in 2D
        with 3 degrees of freedom/node :\n 
        
        - **Kinematics**: horizontal, vertical displacement \
        and rotation with a linear interpolation inside the element  
        
        .. math:: \{U\}=\\langle u_x^1,u_y^1,\\theta_z^1,u_x^2,u_y^2,\\theta_z^2\\rangle^T 
            
        - **Strains**: 
            + axial strain :math:`\epsilon`  (constant)
            + bending curvature :math:`\chi` (linear)
        - **Stresses**: 
            + normal force :math:`N` (constant)
            + bending moment :math:`M` (constant)
            + shear force :math:`V` (linear)
            
        .. math:: \{\\Sigma\}=\\langle N,M,V^1,V^2\\rangle^T 
    """
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------
        
        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        Segment.__init__(self,node_list,tag)
        # number of degrees of freedom per element
        self.el_dof = 6
        # number of degrees of freedom per node
        self.node_dof = 3
        # number of stress resultants
        self.nb_stresses = 4
    
    def rotation_matrix(self):
        """
        Rotation matrix :math:`[R]` from global to local frame 
        
        shape = (6,6)
        """
        T = self.nodes.coor
        tang = T[1,:]-T[0,:]
        L = self.measure()
        t = tang/L
        
        # rotation matrix [r]{U_x,U_y,Theta_z} = {U_t,U_n,Theta_z}
        r = np.array([[t[0],t[1],0],[-t[1],t[0],0],[0,0,1]])
        # R is such that [R]{U_x1,U_y1,U_x2,U_y2} = {U_t1,U_t2}
        return np.kron(np.eye(2),r) # get a block matrix by repeating [r] twice along the diagonal, 
    
    
    def compute_Be_matrix(self,x):
        """ Strain :math:`[B]` matrix such that
            ?????            
        """
# =============================================================================
#         A COMPLETER
# =============================================================================
        return B
        
        
    def elementary_stiffness(self,mat,sect):
        """ Elementary stiffness matrix :math:`[K_e]` in global frame shape=(12,12)""" 
        L = self.measure()
        R = self.rotation_matrix()
        
        E = mat.Young_modulus
        lamb,mu = mat.compute_lame_coeff()
        S = sect.area
        I = sect.inertia
# =============================================================================
#         A COMPLETER
# =============================================================================
        return Ke_glob
   
    
    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading
        
        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces and couple :math:`(f_x,f_y,c_z)`           
        """
        fx,fy,cz = el_force
            
        L = self.measure()
        R = self.rotation_matrix()
        r = R[:3,:3]   
# =============================================================================
#         A COMPLETER
# =============================================================================
  
    
    def stresses(self,Ue,mat,sect):
        """ Compute generalized stresses
            
            .. math:: \{\Sigma\} = \\langle ?\\rangle^T
        
        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        L = self.measure()
        E = mat.Young_modulus
        lamb,mu = mat.compute_lame_coeff()
        S = sect.area
        I = sect.inertia
# =============================================================================
#         A COMPLETER
# =============================================================================


    def elementary_mass(self,mat,sect,lumped=False): 
# =============================================================================
#         A COMPLETER (EVENTUELLEMENT)
# =============================================================================
        pass
      
        
    def elementary_thermal_vector(self,mat,sect,dilat):        
        pass
    
    def elementary_geometric_stiffness(self,Sig,mat,sect):        
        pass
        
    def deformation(self,Ue):
        """ Interpolation of the deformed element
        
        Parameters
        ----------
        Ue : ndarray
            nodal displacement of the current elements 
        Returns
        -------
        x_def,y_def : ndarray
            returns deformed position of m points along the element
        """
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = x + Ue[::self.node_dof]
        y_def = y + Ue[1::self.node_dof]
        return x_def,y_def
        




        
