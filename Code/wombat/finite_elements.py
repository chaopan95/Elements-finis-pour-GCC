# -*- coding: utf-8 -*-
"""
Module containing assembling procedures and linear solver
(version with dense matrix format)

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
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
        L : ndarray
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

    return sorted_data[row_mask,:],Ud[sorted_idx[row_mask]]


def assembl_stiffness_matrix(model):
    """ Assembly procedure of the global stiffness matrix

        Parameters
        ----------
        model
            :class:`Model <model.Model>` object

        Returns
        -------
        K : ndarray
            :math:`[K]` global stiffness matrix of shape=(Nd,Nd) where Nd number of dofs
    """
    mesh = model.mesh
    Nd = mesh.node_dof*mesh.Nno
    K = np.zeros((Nd,Nd))

    for e in mesh.elem_list:
        dofe = e.get_dof()

        # call elementary stiffness method
        Ke = e.elementary_stiffness(e.mat,e.sect)

        # assembling step
        K[np.ix_(dofe,dofe)] += Ke

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


def solve(K,F,L=[],Ud=[]):
    """ Resolution of the global finite element linear system using Lagrange multipliers

        :math:`\\begin{bmatrix} K & L^T \\\\ L & 0 \end{bmatrix}\\begin{Bmatrix} U \\\\ \lambda \end{Bmatrix}=
        \\begin{Bmatrix} F \\\\ U_d \end{Bmatrix}`

        Parameters
        ----------
        K : ndarray
            global stiffness matrix :math:`[K]` shape (Nd,Nd)
        F : ndarray
            global forces vector :math:`\{F\}` shape (Nd,)
        L : ndarray
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
    A = np.vstack((np.hstack((K,L.T)),np.hstack((L,Z))))
    # Building complete right-hand-side {b} = {F,Ud}
    b = np.hstack((F,Ud))
    # solving linear system [A]{x} = {b}
    tic = time.clock()
    x = np.linalg.solve(A,b)
    toc = time.clock()
    print "Linear solver time : %f s" % (toc-tic)

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