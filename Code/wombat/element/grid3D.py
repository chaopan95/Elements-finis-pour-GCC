# -*- coding: utf-8 -*-
"""
Module for class `Grid3D`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *

class Grid3D(Segment):
    """ A 3D grid element combining Euler-Bernoulli out-of-plane bending and St-Venant torsion
    
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
        
# =============================================================================
        r = [[1,0,0], [0,t[0],t[1]],[0,-t[1],t[0]]]
# =============================================================================
        return np.kron(np.eye(2),r) # get a block matrix by repeating [r] twice along the diagonal, 
    
    def shape_functions(self,x):
        """ Cubic shape functions of the Hermite beam element

        Parameters
        ----------
        x : float
            position along the beam axis :math:`x\in[0;L]`
        
        Returns
        -------
        N : ndarray shape=(4,)
            shape functions :math:`[N_i(x)]` evaluated at `x`
        DN : ndarray shape=(4,)
            shape functions first derivatives :math:`[N'_i(x)]` evaluated at `x`
        D2N : ndarray shape=(4,)
            shape functions second derivatives :math:`[N''_i(x)]` evaluated at `x`
            
        """
        L = self.measure()
        N = np.array([1-3*x**2/L**2+2*x**3/L**3,
                      x-2*x**2/L+x**3/L**2,
                      3*x**2/L**2-2*x**3/L**3,
                      -x**2/L+x**3/L**2])
        DN = np.array([6*x/L**2+6*x**2/L**3, 
                       1-4*x/L+3*x**2/L**2, 
                       6*x/L**2-6*x**2/L**3, 
                       -2*x/L+3*x**2/L**2])
        D2N = np.array([6*(-1 + 2*x/L)/L**2,
                        2*(-2 + 3*x/L)/L, 
                        6*(1 - 2*x/L)/L**2,
                        2*(-1 + 3*x/L)/L])
        return N,DN,D2N
    
    def compute_Be_matrix(self,x):
        """ Strain :math:`[B]` matrix such that ???
        
        Parameters
        ----------
        x : float
            position along the beam axis :math:`x\in[0;L]`
            
        """
        L = self.measure()
        N,DN,D2N = self.shape_functions(x)
# =============================================================================
        B = np.zeros((2,6))
        B[1,[2,5]] = np.array([-1,1])/L
        B[0,[0,3]] = D2N[::2]
        B[0,[1,4]] = D2N[1::2]
# =============================================================================
        return B
        
        
    def elementary_stiffness(self,mat,sect):
        """ Elementary stiffness matrix :math:`[K_e]` in global frame shape= ?""" 
        L = self.measure()
        R = self.rotation_matrix()
        E = mat.Young_modulus
        lamb, mu = mat.compute_lame_coeff()
        S = sect.area
        I = sect.inertia
        J = sect.torsion

        
# =============================================================================
        nu=mat.Poisson_coeff
        G = E / (2+2*nu)
        # elementary torsion stiffness matrix in local frame [Ke_torsion]{Omega_1,Omega_2} = {T1,T2}
        Ke_torsion_loc = G*J/L*(np.array([[1,-1],[-1,1]]))
        # elementary bending stiffness matrix in local frame [Ke_bend] = {u_1,Theta_1,U_1,Theta_2}={V1,Mz_1,V2,Mz_2}
        Ke_bend_loc = E*I/L**3*(np.array([[12,6*L,-12,6*L],
                                         [6*L,4*L**2,-6*L,2*L**2],
                                         [-12,-6*L,12,-6*L],
                                         [6*L,2*L**2,-6*L,4*L**2]]))

        Ke_loc = np.zeros((6,6))
        # add torsion stiffness to at rows/columns corresponging to Omega_1,Omega_2
        Ke_loc[np.ix_([2,5],[2,5])] = Ke_torsion_loc
        # add bending moment stiffness to at rows/columns corresponging to u_1,Theta_1,U_1,Theta_2
        Ke_loc[np.ix_([0,1,3,4],[0,1,3,4])] = Ke_bend_loc
        # elementary stiffness matrix in global frame [Ke_glob] = [R]^T*[Ke_loc]*[R]
        Ke_glob = np.dot(np.dot(R.T,Ke_loc),R)
# =============================================================================
        return Ke_glob
        
    
    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading
        
        Parameters
        ----------
        el_force = [fz,cx,cy] : array,list
            contains uniformly distributed forces and couple :math:`(f_z,c_x,c_y)`           
        """    
        fz,cx,cy = el_force
            
        L = self.measure()
        R = self.rotation_matrix()
        r = R[:3,:3]
        
# =============================================================================
        P = np.zeros((3,6))
        P[0,[0,1,3,4]] = np.array([L/2.,L**2/12.,L/2.,-L**2/12.])
        P[1,[0,1,3,4]] = np.array([-1,0,1,0])
        P[2,[2,5]] = np.array([L/2.,L/2.])
        return np.dot(np.dot(np.dot(np.array([fz, cx, cy]), r.T), P), R)   
# =============================================================================
        
    def stresses(self,Ue,mat,sect):
        """ Compute generalized stresses
            
            .. math:: \{\Sigma\} =  ?
                
        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        L = self.measure()
        E = mat.Young_modulus
        lamb, mu = mat.compute_lame_coeff()
        S = sect.area
        I = sect.inertia
        J = sect.torsion
        
        R = self.rotation_matrix()
        
# =============================================================================
        nu = mat.Poisson_coeff
        G = E / (2+2*nu)
        Uloc = np.dot(R,Ue)
        B1 = self.compute_Be_matrix(0)
        B2 = self.compute_Be_matrix(L)
        strain = np.zeros((4,))
        #curvature at x=0 et rotation 
        strain[:2] = np.dot(B1,Uloc)
        # curvature at x=L
        strain[2] = np.dot(B2,Uloc)[0]
        # shear strain
        strain[3] = np.dot(np.array([-12./L**3,-6./L**2,0,12./L**3,-6./L**2,0]),Uloc)
        sig= np.dot(np.diag((E*I,G*J,E*I,E*I)),strain)
        return sig
# =============================================================================
    


    def elementary_mass(self,mat,sect,lumped=False):    
        """ Elementary mass matrix :math:`[M_e]` in global frame shape= ?
        
        Parameters
        ----------
        lumped : bool
            for the `Beam2D` element no lumped mass is defined    
            
        .. note:: Rotation inertia is not taken into account
        """     
        if lumped:
            print ("Warning : no lumped Mass matrix defined for beam elements")
            
# =============================================================================
        L = self.measure()
        S = sect.area
        rho = mat.rho
        
        R = self.rotation_matrix()
        # elementary torsion mass matrix in local frame
        Me_torsion_loc = rho*S*L*(np.array([[0.,0.],[0.,0.]]))
        # elementary bending mass matrix in local frame
        Me_bend_loc = rho*S*L*(np.array([[13/35.,11*L/210.,9/70.,-13*L/420.], \
                                         [11*L/210.,L**2/105.,13*L/420.,-L**2/140.],
                                         [9/70.,13*L/420.,13/35.,-11*L/210.],
                                         [-13*L/420.,-L**2/140.,-11*L/210.,L**2/105.]]))

        Me_loc = np.zeros((6,6))
        # add torsion stiffness to at rows/columns corresponging to U_t1,U_t2
        Me_loc[np.ix_([2,5],[2,5])] = Me_torsion_loc
        # add bending moment stiffness to at rows/columns corresponging to U_n1,Theta_y1,U_n2,Thtea_y2
        Me_loc[np.ix_([0,1,3,4],[0,1,3,4])] = Me_bend_loc

        # elementary stiffness matrix in global frame [Ke_glob] = [R]^T*[Ke_loc]*[R]
        Me_glob = np.dot(np.dot(R.T,Me_loc),R)

        return Me_glob
# =============================================================================
       
        return Me_glob

    def deformation(self,Ue,m=21):
        pass

    def elementary_geometric_stiffness(self,Sig,mat,sect):    
        pass
        
    def elementary_thermal_vector(self,mat,sect,dilat):        
        pass
