# -*- coding: utf-8 -*-
"""
Module for nonlinear computations

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
from .utils import value_at_next_step
from .finite_elements_sparse import solve as ssolve
from .finite_elements_sparse import assembl_stiffness_matrix, internal_forces

def nonlinear_solver(evol,F_list,L,Ud_list,model,Sig=0,X=0,method="newton-raphson",tol=1e-6,nitermax=100,print_info=False):
    """ Nonlinear solver using a Newton method

    Parameters
    ----------
    evol : list, ndarray
        value of the loading amplitude at the corresponding time
    F_list : ndarray or list
        external reference force vector, amplitude is multiplied by evol
        :math:`F(t) = F_{list}[0]\\cdot evol(t)+F_{list}[1]` if two elements in F_list
    L : sparse matrix
        connection matrix
    Ud_list : ndarray or list
        imposed displacements reference vector, amplitude is multiplied by evol
        :math:`Ud(t) = Ud_{list}[0]\\cdot evol(t)+Ud_{list}[1]` if two elements in Ud_list
    model
        :class:`Model <model.Model>` object
    Sig : ndarray
        initial stress state (default is natural state)
    X : ndarray
        initial internal variable vector (default is zero)
    method : {"newton-raphson","modified-newton"}
        method type (use of tangent matrix or not)
    tol : float
        stopping criterion tolerance (default = 1e-6)
    nitermax : int
        maximum number of iterations before stopping an increment (default = 100)
    print_info : bool
        print iteration informations if True

    Returns
    -------
    U_list : ndarray
        displacement vector at all Nsteps time steps shape=(Nd,Nsteps)
    lamb_list : ndarray
        vector of Lagrange multipliers at all time steps
    Sig_list : ndarray
        vector of stress state at all time steps
    lX_list : ndarray
        vector of internal variables at all time steps
    info_list : list
        list of algorithm information at all time steps
    """
    solver_parameters = {"method":method, "print_info":print_info, "tol":tol, "max_iterations":nitermax}
    mesh = model.mesh
    Nd = mesh.node_dof*mesh.Nno

    nincr = len(evol)
    U_list = np.zeros((Nd,nincr))
    lamb = np.zeros((L.shape[0],))
    if not(isinstance(Sig,np.ndarray)):
        Sig = np.zeros((mesh.nb_stresses*mesh.Nel,))
    if not(isinstance(X,np.ndarray)):
        X = np.zeros((mesh.ngauss*mesh.Nel*model.mat.nvar,))

    K = assembl_stiffness_matrix(model)
    Fint = np.zeros((Nd,))
    Fint,new_Sig,dX,_ = internal_forces(np.zeros((Nd,)),Sig,X,model,tangent_matrix=True)

    last_Ud = 0*lamb
    Sig_list = np.repeat(Sig[:,np.newaxis],nincr,axis=1)
    X_list = np.repeat(X[:,np.newaxis],nincr,axis=1)
    lamb_list = np.repeat(lamb[:,np.newaxis],nincr,axis=1)
    info_list = []
    for (i,t) in enumerate(evol[1:]):
        new_F = value_at_next_step(F_list,t)
        new_Ud = value_at_next_step(Ud_list,t)
        DU,lamb,Sig,X,Fint,info = solve_increment(K,new_F,L,new_Ud-last_Ud,lamb,Sig,X,Fint,model,solver_parameters)
        U_list[:,i+1] = U_list[:,i] + DU
        Sig_list[:,i+1] = Sig
        X_list[:,i+1] = X
        lamb_list[:,i+1] = lamb
        info_list.append(info)
        last_Ud = new_Ud
        if print_info:
            print("Increment %i : Converged in %i iterations" % (i+1,info["niter"]))
    return U_list,lamb_list,Sig_list,X_list,info_list


def solve_increment(K,Fext,L,Ud,lamb,Sig,X,Fint,model,solver_parameters={"max_iterations":100,"method":"newton-raphson","print_info":False,"tol":1e-6}):
    """
    Solves a nonlinear increment using the Newton method

    Parameters
    ----------
    K : sparse matrix
        stiffness matrix
    Fext : ndarray
        external force vector at current time step
    L : sparse matrix
        connection matrix
    Ud : ndarray
        vector of imposed displacement at given time step
    lamb : ndarray
        Lagrange multiplier at previous time step
    Sig : ndarray
        stress state at previous time step
    X : ndarray
        internal variables at previous time step
    Fint : ndarray
        vector of internal forces at previous time step
    model
        :class:`Model <model.Model>` object
    solver_parameters
        dictionnary containing the solver parameters
    """
    sp = solver_parameters
    nitermax = sp["max_iterations"]
    if sp["method"] == "newton-raphson":
        tangent_matrix = True
    else:
        tangent_matrix = False

    DU = np.zeros((K.shape[0],))
    DX = 0*X
    Res = Fint+L.T.dot(lamb)-Fext
    nRes0 = max(np.linalg.norm(Res),np.linalg.norm(Ud))
    nRes = nRes0
    new_Sig = Sig
    flag_Ud = 1
    niter = 0
    info = {"niter":0,"residual":nRes0}
    K_newton = K
    while ((nRes/nRes0>sp["tol"]) and (niter<nitermax)):
        niter += 1
        dU,dlamb = ssolve(K_newton,-Res,L,flag_Ud*Ud)
        DU += dU
        lamb += dlamb
        Fint,new_Sig,dX,Kt = internal_forces(DU,Sig,X,model,tangent_matrix=tangent_matrix)
        if tangent_matrix:
            K_newton = Kt
        else:
            K_newton = K
        DX = dX
        Res = Fint+L.T.dot(lamb)-Fext
        nRes = np.linalg.norm(Res)
        info["niter"] += 1
        info["residual"] = nRes/nRes0
        if sp["print_info"]:
            print("     Iteration %i Normalized residual norm %e" % (niter,nRes/nRes0))
        flag_Ud = 0
    if niter == nitermax:
        raise ValueError("Newton-Raphson method did not converge in less than %i iterations\n Residual norm is %e" % (nitermax,nRes/nRes0))
    return DU,lamb,new_Sig,X+DX,Fint,info