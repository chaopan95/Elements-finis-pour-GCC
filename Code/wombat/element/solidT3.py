# -*- coding: utf-8 -*-
"""
Module for class `SolidT3`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *
from .trace_elements import *

class SolidT3(Triangle):
    """ A 2D triangular element for continuum mechanics

        SolidT3 is a Triangle-type element (3 nodes) in 2D
        with 2 degrees of freedom/node :\n

        - **Kinematics**:  horizontal, vertical displacement \
             with a linear interpolation inside the element

        .. math:: \{U\}=\\langle u_x^1,u_y^1,u_x^2,u_y^2,u_x^3,u_y^3\\rangle^T

        - **Strains** :  plane components of :math:`\\underline{\\underline{\\varepsilon}} \
           = (\\underline{\\nabla u} + \\underline{\\nabla u}^T)/2`
           :math:`\\varepsilon_{xx},\\varepsilon_{yy},2\\varepsilon_{xy}` (constant)

        - **Stresses**: :math:`\\{\\sigma_{xx}^e,\\sigma_{yy}^e,\\sigma_{zz}^e,\\sigma_{xy}^e\\}` (constant)

            .. note:: the out-of-plane stress :math:`\\sigma_{zz}` is not necesseary for building the stiffness \
            matrix but is still given as an output (computed from :math:`\\sigma_{zz}= \
            \\lambda(\\varepsilon_{xx}+\\varepsilon_{yy})`)
    """
    trace = TraceSolidT3   # defines the corresponding trace element

    def __init__(self,node_list=[],tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        Triangle.__init__(self,node_list,tag)
        self.elem_type = 'T3'
        self.el_dof = 6
        self.node_dof = 2
        self.nb_stresses = 4  # constant stress in element
        self.kin_field_names = ['U_x','U_y']
        self.strain_field_names = ['eps_xx','eps_yy','2eps_xy']
        self.int_forces_field_names = ['sig_xx','sig_yy','sig_zz','sig_xy']
        self.ext_forces_field_names = ['F_x','F_y']

        self.detJ,self.Jac,self.invJac = self.jacobian()
        self.A = self.detJ/2.
        self.Be = self.compute_Be_matrix(0,0)

    def shape_functions(self,xi,eta):
        """ Returns the shape functions and its derivatives

        Parameters
        -----------
        xi : float
            coordinate of point :math:`\\xi` in the reference space, belongs to :math:`[0;1]`
        eta : float
            coordinate of point :math:`\\eta` in the reference space, belongs to :math:`[0;1]`

        Returns
        --------
        N : ndarray shape = (3,)
            array of shape functions :math:`[N]` evaluated at :math:`(\\xi,\\eta)`
        DN : ndarray shape = (2,3)
            array of shape functions derivatives :math:`[\\nabla N]` evaluated at :math:`(\\xi,\\eta)` (here constant)
        """
        N = np.array([xi,eta,1-xi-eta])
        DN = np.array([[1,0,-1],[0,1,-1]])
        return N,DN

    def jacobian(self):
        """ Computes quantities related to the jacobian of the element

        Returns
        --------
        detJ : float
            determinant of the jacobian matrix (must be strictly positive)
        Jac : ndarray shape (2,2)
            jacobian matrix :math:`[J]`
        invJac : ndarray shape (2,2)
            inverse of the jacobian matrix :math:`[J]^{-1}`
        """
        T = self.node_coor()
        Jac = np.array([T[0,:]-T[2,:],T[1,:]-T[2,:]])
        detJ = np.linalg.det(Jac)
        assert detJ>0, "Jacobian of element %i is negative." % self._id
        invJac = np.linalg.inv(Jac)
        return detJ,Jac,invJac

    def compute_Be_matrix(self,xi,eta):
        """ Local strain matrix :math:`[B_e]` such that

        .. math:: [B_e]\{U_e\} = \\begin{Bmatrix} \\varepsilon_{xx} \\\\ \\varepsilon_{yy} \\\\ 2\\varepsilon_{xy} \end{Bmatrix}

        (here [B_e] is constant, shape = (3,6))
        """
        N,DN = self.shape_functions(xi,eta)
        Be = np.zeros((3,6))
        GN = np.dot(self.invJac,DN)
        Be[0,0::2] = GN[0,:]
        Be[1,1::2] = GN[1,:]
        Be[2,0::2] = GN[1,:]
        Be[2,1::2] = GN[0,:]
        return Be

    def elementary_stiffness(self,mat,sect=None):
        """ Elementary stiffness matrix :math:`[K_e]` shape=(6,6)"""
        Ke = self.A*np.dot(self.Be.T,np.dot(mat.C,self.Be))
        return Ke

    def elementary_mass(self,mat,sect=None,lumped=False):
        """ Elementary mass matrix :math:`[M_e]` shape=(6,6)

        Parameters
        ----------
        lumped : bool
            if False, returns the consistent elementary mass matrix, else the lumped mass matrix
        """
        if lumped:
            Me = mat.density*self.A/3.*np.eye(6)
        else:
            Me = mat.density*self.A/12.*np.kron(np.array([[2,1,1],[1,2,1],[1,1,2]]),np.eye(2))
        return Me

    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading

        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces :math:`(f_x,f_y)`

           .. note:: for the SolidT3 element cz is ignored
        """
        fx,fy,cz = el_force
        return self.A/3.*np.tile([fx,fy],3)

    def deformation(self,Ue):
        """ Interpolation of the deformed element

        Parameters
        ----------
        Ue : ndarray
            nodal displacement of the current elements

        Returns
        -------
        x_def,y_def : ndarray
            returns deformed position of element nodes
        """
        Ux = Ue[0::self.node_dof]
        Uy = Ue[1::self.node_dof]
        s = np.linspace(-1,1,3)
        # repeat twice along axis 0
        s = np.tile(s,(1,1)).T
        S1 = np.concatenate(((1-s)/2.,(1+s)/2.,0.*s),axis=1)
        S2 = np.concatenate((0*s,(1-s)/2.,(1+s)/2.),axis=1)
        S3 = np.concatenate(((1+s)/2.,0*s,(1-s)/2.),axis=1)
        S = np.concatenate((S1,S2,S3),axis=0)
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = np.dot(S,x+Ux)
        y_def = np.dot(S,y+Uy)
        return x_def,y_def

    def stresses(self,Ue,mat,sect=None):
        """ Compute stress state

            .. math:: \{\Sigma\} = \\begin{Bmatrix} \\sigma_{xx} \\\\ \\sigma_{yy} \\\\ \\sigma_{zz} \\\\ \\sigma_{xy} \\end{Bmatrix}

        .. note:: :math:`\\sigma_{zz}` is not used for the computation but computed from the strain

        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        C = mat.C_matrix()

        Eps = np.dot(self.Be,Ue)
        Sig_plane = np.dot(C,Eps)
        Sig = np.array([Sig_plane[0],
                             Sig_plane[1],
                             mat.compute_sigzz(Eps),
                             Sig_plane[2]])
        return Sig

    def internal_forces(self,Sige):
        """ Returns elemental contribution of a stress state Sige to internal forces vector Fint """
        return self.A*np.dot(self.Be.T,Sige[[0,1,3]])


    def internal_forces_nl(self,Ue,Sige,Xe,mat,sect):
        """ Returns elemental contribution to internal forces vector Fint,
            new elemental stresses and internal variable increments
            N : normal force (constant in the element)
        """
        Deps = np.dot(self.Be,Ue)

        sig,dXe,Ct = mat.constitutive_relation(Deps,Sige,Xe)

        fint = self.A*np.dot(self.Be.T,sig[[0,1,3]])

        # elementary stiffness matrix in local frame
        Kte = self.A*np.dot(self.Be.T,np.dot(Ct,self.Be))

        return fint,sig,dXe,Kte
