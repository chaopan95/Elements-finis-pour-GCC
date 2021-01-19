# -*- coding: utf-8 -*-
"""
Module for class `SolidT6`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import *
from .trace_elements import TraceSolidT6

class SolidT6(Triangle6):
    """ A 2D quadratic triangular element for continuum mechanics

        SolidT6 is a Triangle6-type element (6 nodes) in 2D
        with 2 degrees of freedom/node and with **straight sides** (affine interpolation of the geometry):\n

        - **Kinematics**:  horizontal, vertical displacement \
             with a quadratic interpolation inside the element

        .. math:: \{U\}=\\langle u_x^1,u_y^1,\\ldots,u_x^6,u_y^6\\rangle^T

        - **Strains** :  plane components of :math:`\\underline{\\underline{\\varepsilon}} \
           = (\\underline{\\nabla u} + \\underline{\\nabla u}^T)/2`
           :math:`\\varepsilon_{xx},\\varepsilon_{yy},2\\varepsilon_{xy}` (linear)

        - **Stresses**: :math:`\\{\\sigma_{xx}^1,\\sigma_{yy}^1,\\sigma_{zz}^1,\\sigma_{xy}^1,\\ldots,\\sigma_{xy}^3\\}` (linear)

            .. note:: the out-of-plane stress :math:`\\sigma_{zz}` is not necesseary for building the stiffness \
            matrix but is still given as an output (computed from :math:`\\sigma_{zz}= \
            \\lambda(\\varepsilon_{xx}+\\varepsilon_{yy})`)
    """
    trace = TraceSolidT6        # defines the corresponding trace element

    def __init__(self,node_list=[],tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        Triangle6.__init__(self,node_list,tag)
        self.elem_type = 'T6'
        self.el_dof = 12
        self.node_dof = 2
        self.ngauss = 3
        self.nb_stresses = 4*self.ngauss  # 4 stresses at ngauss Gauss points in element

        self.kin_field_names = ['U_x','U_y']
        self.strain_field_names = ['eps_xx','eps_yy','2eps_xy']
        self.int_forces_field_names = ['sig_xx','sig_yy','sig_zz','sig_xy']
        self.ext_forces_field_names = ['F_x','F_y']

        self.detJ,self.Jac,self.invJac = self.jacobian()
        self.A = self.detJ/2.
        self.ag,self.wg = self.gauss_quadrature(self.ngauss)
        self.Be = [self.compute_Be_matrix(a[0],a[1]) for a in self.ag]

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
        N : ndarray shape = (6,)
            array of shape functions :math:`[N]` evaluated at :math:`(\\xi,\\eta)`
        DN : ndarray shape = (2,6)
            array of shape functions derivatives :math:`[\\nabla N]` evaluated at :math:`(\\xi,\\eta)`
        """
        N = np.array([xi*(2*xi-1),
                      eta*(2*eta-1),
                      (1-xi-eta)*(1-2*xi-2*eta),
                      4*xi*eta,
                      4*eta*(1-xi-eta),
                      4*xi*(1-xi-eta)])
        DN = np.array([[4*xi-1,
                        0,
                        -3+4*xi+4*eta,
                        4*eta,
                        -4*eta,
                        4-8*xi-4*eta],[
                        0,
                        4*eta-1,
                        -3+4*xi+4*eta,
                        4*xi,
                        4-4*xi-8*eta,
                        -4*xi]])
        return N,DN

    def jacobian(self):
        """ Computes quantities related to the jacobian of the element (here constant)

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

        evaluated at :math:`(\\xi,\\eta)`, shape = (3,12)
        """
        N,DN = self.shape_functions(xi,eta)
        GN = np.dot(self.invJac,DN)
        Be = np.zeros((3,12))
        Be[0,0::2] = GN[0,:]
        Be[1,1::2] = GN[1,:]
        Be[2,0::2] = GN[1,:]
        Be[2,1::2] = GN[0,:]
        return Be


    def elementary_stiffness(self,mat,sect=None):
        """ Elementary stiffness matrix :math:`[K_e]` shape=(12,12)"""
        C = mat.C_matrix()
        Ke = np.zeros((12, 12))
        detJ = self.detJ
	# =============================================================
	## A COMPLETER
	#
        for i in range(self.ngauss):
            wgi = self.wg[i] # Gauss point weight number i
            Bei = self.Be[i] # [Be] matrix at Gauss point number i
            Ke += wgi*detJ*Bei.T.dot(C.dot(Bei))
	# =============================================================

        return Ke


    def elementary_mass(self,mat,sect=None,lumped=False):
        """ Elementary mass matrix :math:`[M_e]` shape=(12,12)

        Parameters
        ----------
        lumped : bool
            if False, returns the consistent elementary mass matrix, else the lumped mass matrix
        """
        A = self.measure()
        if lumped:
            Me = mat.density*A/3*scipy.linalg.block_diag(np.eye(6),np.zeros((6,6)))
        else:
            Me = np.zeros((12,12))
            for i in range(self.ngauss):
                N,DN = self.shape_functions(self.ag[i,0],self.ag[i,1])
                # change N from 1D array (shape = (6,)) to a 2D array (shape = (1,6))
                N = np.expand_dims(N,axis=0)
                Me += mat.density*self.detJ*self.wg[i]*np.kron(N.T*N,np.eye(2))
        return Me

    def elementary_thermal_vector(self,mat,dilat):
        # TODO
        pass


    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading

        Parameters
        ----------
        el_force = [fx,fy,cz] : array,list
            contains uniformly distributed forces :math:`(f_x,f_y)`

           .. note:: for the SolidT6 element cz is ignored
        """
        fx,fy,cz = el_force
        fe = np.zeros((12,))
        for i in range(self.ngauss):
            N,DN = self.shape_functions(self.ag[i,0],self.ag[i,1])
            fe += self.detJ*self.wg[i]*np.kron(N,np.array([fx,fy]))
        return fe

    def deformation(self,Ue,m=21):
        """ Interpolation of the deformed element

        Parameters
        ----------
        Ue : ndarray
            nodal displacement of the current elements
        m : int
            number of points used to interpolate the deformed configurations

        Returns
        -------
        x_def,y_def : ndarray
            returns deformed position of m points along the element boundary
        """
        Ux = Ue[0::self.node_dof]
        Uy = Ue[1::self.node_dof]
        s = np.linspace(-1,1,m)
        # repeat twice along axis 0
        s = np.tile(s,(1,1)).T
        Nq1 = s*(s-1)/2
        Nq2 = s*(1+s)/2
        Nq3 = 1-s**2
        z = 0*s
        S1 = np.concatenate((Nq1,Nq2,z,Nq3,z,z),axis=1)
        S2 = np.concatenate((z,Nq1,Nq2,z,Nq3,z),axis=1)
        S3 = np.concatenate((Nq2,z,Nq1,z,z,Nq3),axis=1)
        S = np.concatenate((S1,S2,S3),axis=0)
        x = self.node_coor()[:,0]
        y = self.node_coor()[:,1]
        x_def = np.dot(S,x+Ux)
        y_def = np.dot(S,y+Uy)
        return x_def,y_def

    def stresses(self,Ue,mat,sect=None):
        """ Compute stress state evaluated at Gauss points, shape = (4*ngauss,)

            .. math:: \{\Sigma\} = \\begin{Bmatrix} \\Sigma^1 \\\\ \\vdots \\\\ \\Sigma^{ngauss} \end{Bmatrix} \\text{ with } \{\Sigma^g\} = \\begin{Bmatrix} \\sigma_{xx}^g \\\\ \\sigma_{yy}^g \\\\ \\sigma_{zz}^g \\\\ \\sigma_{xy}^g \\end{Bmatrix}

        .. note:: :math:`\\sigma_{zz}` is not used for the computation but computed from the strain

        Parameters
        ----------
        Ue : ndarray
            nodal values of the displacement
        """
        C = mat.C_matrix()
        sig_g = np.zeros((4*self.ngauss,))
        for i in range(self.ngauss):
            Eps = np.dot(self.Be[i],Ue)
            sig_g[[4*i,4*i+1,4*i+3]] =  np.dot(C,Eps)
            sig_g[4*i+2] = mat.compute_sigzz(Eps)
        return sig_g


    def internal_forces_nl(self,Ue,Sige,Xe,mat,sect):
        """ Returns elemental contribution to internal forces vector Fint,
            new elemental stresses and internal variable increments
            N : normal force (constant in the element)
        """
        sig_g = np.zeros((4*self.ngauss,))
        dX_g = np.zeros((5*self.ngauss,))
        fint = np.zeros((12,))
        Kte = np.zeros((12,12))
        for i in range(self.ngauss):
            Deps = np.dot(self.Be[i],Ue)
            sig,dXe,Ct = mat.constitutive_relation(Deps,Sige[4*i:4*(i+1)],
                                                        Xe[5*i:5*(i+1)])
            sig_g[4*i:4*(i+1)] =  sig
            dX_g[5*i:5*(i+1)] =  dXe
            # elementary stiffness matrix in local frame
            Kte += self.detJ*self.wg[i]*np.dot(self.Be[i].T,np.dot(Ct,self.Be[i]))
            fint += self.detJ*self.wg[i]*np.dot(self.Be[i].T,sig[[0,1,3]])

        return fint,sig_g,dX_g,Kte
