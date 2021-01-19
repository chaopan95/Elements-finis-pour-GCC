# -*- coding: utf-8 -*-
"""
Module containing assembling procedures and linear solver
(version with sparse matrix format)

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
from scipy import linalg, compress
import scipy.sparse as sparse
import scipy.sparse.linalg as slin
import time

def assembl_connections(connections,model):
    """ Assembly procedure for relations between dofs

        Parameters
        ----------
        connections
            :class:`Connections <connections.Connections>` object

        model
            :class:`Model <model.Model>` object

        Returns
        -------
        L : sparse lil matrix
            :math:`[L]` connection matrix of shape (Nl,Nd) where
            Nl number of relations and Nd number of dofs
        Ud : ndarray
            :math:`\{U_d\}` imposed displacement vector of shape (Nl,)
    """
    mesh = model.mesh
    ndof = mesh.node_dof
    Nd = ndof*mesh.Nno
    Nl = connections.nb_relations
    L =  np.zeros((Nl,Nd))
    Ud = np.zeros((Nl,))
    buff = 0
    for i,master in enumerate(connections.master_list):
        for j,slave in enumerate(connections.slave_list[i]):
            slave_dof = slave.get_dof(ndof)
            for k,comp in enumerate(connections.components_list[i]):
                if master is not None:
                    master_dof = master.get_dof(ndof)
                    L[buff,master_dof[comp]] = connections.lin_rela_list[i][0]
                if len(connections.lin_rela_list[i])==2:
                    jrela = 1
                else:
                    jrela = j+1
                L[buff,slave_dof[comp]] = connections.lin_rela_list[i][jrela]
                Ud[buff] = connections.imposed_value_list[i][k]
                buff += 1

    sorted_idx = np.lexsort(L.T)
    sorted_data = L[sorted_idx,:]

    # Get unique row mask
    row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))

    # Get unique rows
    return sparse.coo_matrix(sorted_data[row_mask,:]),Ud[sorted_idx[row_mask]]

def assembl_stiffness_matrix(model):
    """ Assembly procedure of the global stiffness matrix

        Parameters
        ----------
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        K : sparse coo matrix
            :math:`[K]` global stiffness matrix of shape=(Nd,Nd) where Nd number of dofs
    """
    mesh = model.mesh
    Nd = mesh.node_dof*mesh.Nno
    Nb_non_zero = mesh.Nel*mesh.el_dof**2
    K = sparse.coo_matrix((Nd,Nd))
    K.data = np.zeros((Nb_non_zero,),dtype='float64')
    K.row = np.zeros((Nb_non_zero,),dtype='int32')
    K.col = np.zeros((Nb_non_zero,),dtype='int32')
    buff=0

    for e in mesh.elem_list:
        dofe = e.get_dof()

        # call elementary stiffness method
        Ke = sparse.coo_matrix(e.elementary_stiffness(e.mat,e.sect))

        nnz = Ke.data.shape[0]
        K.data[buff:buff+nnz]=Ke.data
        K.row[buff:buff+nnz]=dofe[Ke.row]
        K.col[buff:buff+nnz]=dofe[Ke.col]
        buff+=nnz

    return K


def assembl_external_forces(forces,model):
    """ Assembly procedure of the external forces vector

        Parameters
        ----------
        forces
            :class:`ExtForce <forces.ExtForce>` object
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        F : ndarray
            :math:`\{F\}` vector of equivalent external forces shape=(Nd,)
            where Nd number of dofs
    """
    mesh = model.mesh
    ndof = mesh.node_dof
    Nd = ndof*mesh.Nno
    F = np.zeros((Nd,))
    for index,node in enumerate(forces.node_list):
        Fx = forces.Fx_n[index]
        Fy = forces.Fy_n[index]
        Cz = forces.Cz_n[index]
        node_dof = node.get_dof(ndof)
        if Fx is not None:
            F[node_dof[0]] += Fx
        if Fy is not None:
            F[node_dof[1]] += Fy
        if ndof==3 and Cz is not None:
            F[node_dof[2]] += Cz

    for i,e in enumerate(forces.el_list):
        fx_el = forces.fx_e[i]
        fy_el = forces.fy_e[i]
        cz_el = forces.cz_e[i]

        dofe = e.get_dof()

        F[dofe] += e.elementary_distributed_forces([fx_el,fy_el,cz_el])

    return F


def assembl_thermal_strains(dilat,model):
    """ Assembly procedure of the thermal strain vector

        Parameters
        ----------
        dilat : ndarray
            array of shape (Nel,) containing thermal dilatation of each element
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        Fther : ndarray
            :math:`\{F_{ther}\}` vector of equivalent thermal forces shape=(Nd,)
            where Nd number of dofs
    """
    mesh = model.mesh
    ndof = mesh.node_dof
    Nd = ndof*mesh.Nno
    Fther = np.zeros((Nd,))
    for i,e in enumerate(mesh.elem_list):
        dofe = e.get_dof()

        Fther[dofe] += e.elementary_thermal_vector(e.mat,e.sect,dilat[i])

    return Fther


def assembl_mass_matrix(model,lumped=False):
    """ Assembly procedure of the global mass matrix

        Parameters
        ----------
        model
            :class:`Model <model.Model>` object
        lumped : bool
            if True compute the lumped mass matrix (when available), otherwise compute the consistent mass matrix (default)

        Returns
        -------
        M : sparse coo matrix
            :math:`[M]` global mass matrix of shape=(Nd,Nd) where Nd number of dofs
    """
    mesh = model.mesh
    Nd = mesh.node_dof*mesh.Nno
    Nb_non_zero = mesh.Nel*mesh.el_dof**2
    M = sparse.coo_matrix((Nd,Nd))
    M.data = np.zeros((Nb_non_zero,),dtype='float64')
    M.row = np.zeros((Nb_non_zero,),dtype='int32')
    M.col = np.zeros((Nb_non_zero,),dtype='int32')
    buff=0

    for e in mesh.elem_list:
        dofe = e.get_dof()

        # call elementary mass method
        Me = sparse.coo_matrix(e.elementary_mass(e.mat,e.sect,lumped))

        nnz = Me.data.shape[0]
        M.data[buff:buff+nnz]=Me.data
        M.row[buff:buff+nnz]=dofe[Me.row]
        M.col[buff:buff+nnz]=dofe[Me.col]
        buff+=nnz

    return M

def assembl_geometric_stiffness_matrix(U,model):
    """ Assembly procedure of the global geometric stiffness matrix

        Parameters
        ----------
        U : ndarray
            global displacement vector on which geometrical stiffness effects are evaluated
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        KG : sparse coo matrix
            :math:`[K_G]` global geometric stiffness matrix of shape=(Nd,Nd) where Nd number of dofs
    """
    mesh = model.mesh
    Nd = mesh.node_dof*mesh.Nno
    Nb_non_zero = mesh.Nel*mesh.el_dof**2
    KG = sparse.coo_matrix((Nd,Nd))
    KG.data = np.zeros((Nb_non_zero,),dtype='float64')
    KG.row = np.zeros((Nb_non_zero,),dtype='int32')
    KG.col = np.zeros((Nb_non_zero,),dtype='int32')
    buff=0

    for e in mesh.elem_list:
        dofe = e.get_dof()

        # compute internal forces
        Fint_e = e.stresses(U[dofe],e.mat,e.sect)
        # call elementary geometric stiffness method
        KGe = sparse.coo_matrix(e.elementary_geometric_stiffness(Fint_e,e.mat,e.sect))

        nnz = KGe.data.shape[0]
        KG.data[buff:buff+nnz]=KGe.data
        KG.row[buff:buff+nnz]=dofe[KGe.row]
        KG.col[buff:buff+nnz]=dofe[KGe.col]
        buff+=nnz

    return KG


def solve(K,F,L=[],Ud=[],print_info=False):
    """ Resolution of the global finite element linear system using Lagrange multipliers

        :math:`\\begin{bmatrix} K & L^T \\\\ L & 0 \end{bmatrix}\\begin{Bmatrix} U \\\\ \lambda \end{Bmatrix}=
        \\begin{Bmatrix} F \\\\ U_d \end{Bmatrix}`

        Parameters
        ----------
        K : sparse matrix
            global stiffness matrix :math:`[K]` shape (Nd,Nd)
        F : ndarray
            global forces vector :math:`\{F\}` shape (Nd,)
        L : sparse matrix
            connection matrix :math:`[L]` shape (Nl,Nd)
        Ud : ndarray
            imposed displacement vector :math:`\{U_d\}` shape (Nl,)

        Returns
        -------
        U : ndarray
            :math:`\{U\}` displacement vector shape (Nd,)
        lamb : ndarray
            :math:`{\lambda}\}` Lagrange multiplier vector shape (Nl,)
    """
    n = L.shape[0]
    Z = np.zeros((n,n))
    # Building complete system [A] = [[K, L^T],[L,0]]
    A = sparse.vstack((sparse.hstack((K,L.T)),sparse.hstack((L,Z)))).tocsc()
    # Building complete right-hand-side {b} = {F,Ud}
    b = np.hstack((F,Ud)).T
    # solving linear system [A]{x} = {b}
    tic = time.clock()
    fact = slin.splu(A)
    x = fact.solve(b)
#    x = slin.spsolve(A,b)
    toc = time.clock()
    if print_info:
        print("Linear solver time : %f s" % (toc-tic)    )
    return x[0:K.shape[0]],x[K.shape[0]:]

def stresses(U,model):
    """ Compute generalized stresses (elastic behaviour)

    Parameters
    ----------
    U : ndarray
        displacement vector solution :math:`\{U\}`
    model
        :class:`Model <model.Model>` object

    Returns
    -------
    Sig : ndarray
        vector of generalized stresses :math:`\{\Sigma\}` (depends on the element)
    """
    mesh = model.mesh
    nS = mesh.nb_stresses
    Sig = np.zeros((nS*mesh.Nel,))
    for (i,e) in enumerate(mesh.elem_list):
        Sig[nS*i:nS*(i+1)] = e.stresses(U[e.get_dof()],e.mat,e.sect)
    return Sig

def assembl_initial_state(Sig,model):
    """ Assembly procedure of the initial state internal force vector :math:`\{F^{int,0}\}`

        Parameters
        ----------
        Sig : ndarray
            array containing initial stress state
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        Fint : ndarray
            :math:`\{F^{int,0}\}` vector of equivalent internal forces shape=(Nd,)
            where Nd number of dofs
    """
    mesh = model.mesh
    nS = mesh.nb_stresses
    Nd = mesh.node_dof*mesh.Nno
    Fint = np.zeros((Nd,))

    for (i,e) in enumerate(mesh.elem_list):
        dofe = e.get_dof()
        Fint[dofe] += e.internal_forces(Sig[nS*i:nS*(i+1)])

    return Fint

def eigenmodes(K,B,L=[],nmodes=1,mode="normal",shift=None,compute_error=False):
    """ Extraction of eigenmodes for the generalized eigenvalue problem

        :math:`[K]\{\\xi\\} = \\pm\\lambda[B]\{\\xi\}`

            sign :math:`+` for mode="normal", :math:`-` for mode="buckling"

        Parameters
        ----------
        K : sparse matrix
            global stiffness matrix
        B : sparse matrix
            :math:`[B]=[M]` mass matrix (with mode="normal") or :math:`[B]=[K_G]`
            geometric stiffness matrix (with mode="buckling")
        L : sparse matrix
            connection matrix :math:`[L]` shape (Nl,Nd)
        nmodes : int
            number of modes with smallest eigenvalues to extract, must be less than matrix rank
        mode : {"normal","buckling"}
            + "normal" : eigenvalue problem for dynamic modal analysis :math:`[K]\{\\xi\\} = \\lambda[M]\{\\xi\}`
            + "buckling" : eigenvalue problem for buckling analysis :math:`[K]\{\\xi\\} = -\\lambda[K_G]\{\\xi\}`
        shift : float
            value around which eigenvalues are looked for (improves accuracy, especially for "buckling" mode)
        compute_error : bool
            if True, prints the residual norm :math:`\|[K]\{\\xi\}\mp \\lambda[B]\{\\xi\}\|_2`

        Returns
        -------
        lamb : ndarray shape = (nmodes,)
            array of eigenvalues :math:`\\lambda`
        xi : ndarray shape = (Nd,nmodes)
            array of eigenmodes :math:`\{\\xi\}`
    """
    tic = time.clock()
    Z = sparse.csc_matrix(null(L.toarray()))
    r = Z.shape[1]
    Kred = Z.T.dot(K.dot(Z))
    Bred = Z.T.dot(B.dot(Z))
    if mode == "normal":
        sign = 1.
        if shift == None:
            shift = 0
    elif mode == "buckling":
        sign = -1.
        if shift == None:
            # Needs an estimate of the minimum eigenvalue
            shift = slin.norm(B)/slin.norm(K,ord=-np.inf)
    if nmodes < r:
        sig,v = slin.eigsh(Kred,nmodes,sign*Bred,sigma=shift,which='LM',mode=mode)
    elif nmodes == r:
        sig,v = linalg.eigh(Kred.toarray(),sign*Bred.toarray())
    else:
        raise ValueError("Requested number of modes "+str(nmodes)+" is larger than number of dofs "+str(r))
    reorder = np.argsort(np.abs(sig))
    sig = sig[reorder]
    v = v[:,reorder]

    if compute_error:
            for i in range(nmodes):
                print("Error on mode %i: %e"%(i,np.linalg.norm(Kred.dot(v[:,i])-sign*sig[i]*Bred.dot(v[:,i]))/np.linalg.norm(sig[i]*Bred.dot(v[:,i]))))
    toc = time.clock()
    print("Eigensolver time : %f s\n" % (toc-tic))
    return sig,Z.dot(v)

def null(A, eps=1e-12):
    """ Compute null space of matrix A up to a tolerance eps """
    u, s, vh = linalg.svd(A)
    padding = max(0,np.shape(A)[1]-np.shape(s)[0])
    null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
    null_space = compress(null_mask, vh, axis=0)
    return null_space.T

def compute_Rayleigh_damping(K,M,alpha_K=0,alpha_M=0):
    """ Returns Rayleigh damping matrix :math:`[C] = \\alpha_K[K]+\\alpha_M[M]` """
    return sparse.coo_matrix(alpha_K*K+alpha_M*M)

def internal_forces(DU,Sig,X,model,tangent_matrix=False):
    mesh = model.mesh
    nS = mesh.nb_stresses
    nV = mesh.ngauss*model.mat.nvar
    new_Sig = np.zeros((nS*mesh.Nel,))
    dX = np.zeros((nV*mesh.Nel,))
    Nd = mesh.node_dof*mesh.Nno
    Fint = np.zeros((Nd,))
    if tangent_matrix:
        Nb_non_zero = mesh.Nel*mesh.el_dof**2
        Kt = sparse.coo_matrix((Nd,Nd))
        Kt.data = np.zeros((Nb_non_zero,),dtype='float64')
        Kt.row = np.zeros((Nb_non_zero,),dtype='int32')
        Kt.col = np.zeros((Nb_non_zero,),dtype='int32')
        buff=0

    for (i,e) in enumerate(mesh.elem_list):
        dofe = e.get_dof()
        fint_el,sig_el,dX_el,Kte = e.internal_forces_nl(DU[e.get_dof()],Sig[nS*i:nS*(i+1)],X[nV*i:nV*(i+1)],e.mat,e.sect)
        if tangent_matrix:
            Kte = sparse.coo_matrix(Kte)
            nnz = Kte.data.shape[0]
            Kt.data[buff:buff+nnz]=Kte.data
            Kt.row[buff:buff+nnz]=dofe[Kte.row]
            Kt.col[buff:buff+nnz]=dofe[Kte.col]
            buff+=nnz
        new_Sig[nS*i:nS*(i+1)] = sig_el
        dX[nV*i:nV*(i+1)] = dX_el
        Fint[dofe] += fint_el
    if tangent_matrix:
        return Fint,new_Sig,dX,Kt
    else:
        return Fint,new_Sig,dX,0
