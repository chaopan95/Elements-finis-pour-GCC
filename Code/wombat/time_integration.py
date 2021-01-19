# -*- coding: utf-8 -*-
"""
Module for time integration

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
from .finite_elements_sparse import solve as ssolve
from .utils import value_at_next_step


class TimeIntegrator:
    def __init__(self,time,evol,beta=1/4.,gamma=1/2.):
        """ Time integrator using the Newmark method

        Parameters
        ----------
        time : ndarray
            values of discrete times
        evol : ndarray
            value of the loading amplitude at the corresponding time
        beta : float
            :math:`\\beta` parameter of the Newmark method (default = 0.25)
        gamma : float
            :math:`\\gamma` parameter of the Newmark method (default = 0.5)
        """
        assert gamma>=1/2., "gamma < 1/2 => unstable scheme"
        self.time = time
        self.time_step = time[-1]/len(time)
        self.evolution = evol
        self.newmark_beta = beta
        self.newmark_gamma = gamma
        self.curr_time = evol[0]
        self.nb_time_steps = evol.shape[0]-1

    def next_step(self,Uold,Vold,Aold,K,C,M,F,L,Ud):
        """ Performs one step of the Newmark method using values at previous step

        Parameters
        -----------
        Uold : ndarray
            displacement at previous time step
        Vold : ndarray
            velocity at previous time step
        Aold : ndarray
            acceleration at previous time step
        K : sparse matrix
            stiffness matrix
        C : sparse matrix
            damping matrix
        M : sparse matrix
            mass matrix
        F : ndarray
            external force vector
        L : sparse matrix
            connection matrix
        Ud : ndarray
            imposed displacements vector

        Returns
        -------
        Unew : ndarray
            displacement at current time step
        Vold : ndarray
            velocity at current time step
        Aold : ndarray
            acceleration at current time step

        """
        beta = self.newmark_beta
        gamma = self.newmark_gamma
        dt = self.time_step
        Upred = Uold+dt*Vold+dt**2/2.*(1-2*beta)*Aold
        Vpred = Vold+dt*(1-gamma)*Aold
        R = F-K.dot(Upred)-C.dot(Vpred)
        S = M+gamma*dt*C+beta*dt**2*K
        Anew,lambdanew = ssolve(S,R,L,beta*dt**2*Ud)
        Unew = Upred + dt**2*beta*Anew
        Vnew = Vpred + dt*gamma*Anew
        return Unew, Vnew, Anew

    def solve_evolution(self,U0,V0,K,C,M,F_list,L,Ud_list):
        """ Solve evolution problem for all time steps

        Parameters
        ----------
        K : sparse matrix
            stiffness matrix
        C : sparse matrix
            damping matrix
        M : sparse matrix
            mass matrix
        F_list : ndarray or list
            external reference force vector, amplitude is multiplied by evol
            :math:`F(t) = F_{list}[0]\\cdot evol(t)+F_{list}[1]` if two elements in F_list
        L : sparse matrix
            connection matrix
        Ud_list : ndarray or list
            imposed displacements reference vector, amplitude is multiplied by evol
            :math:`Ud(t) = Ud_{list}[0]\\cdot evol(t)+Ud_{list}[1]` if two elements in Ud_list

        Returns
        -------
        U : ndarray
            displacement vector at all Nsteps time steps shape=(Nd,Nsteps)
        """
        Uold = U0
        Vold = V0
        F = value_at_next_step(F_list,self.evolution[0])
        Ud = value_at_next_step(Ud_list,self.evolution[0])
        Aold = ssolve(M,F-K.dot(U0)-C.dot(V0),L,Ud)[0]
        U = np.zeros((Uold.shape[0],self.nb_time_steps+1))
        U[:, 0] = U0
        V = 0*U
        V[:, 0] = V0
        for i in range(self.nb_time_steps):
            self.curr_time = self.time[i+1]
            # loading at next time step (sum of a variable loading and a fixed loading)
            F = value_at_next_step(F_list,self.evolution[i+1])
            # same for displacements
            Ud = value_at_next_step(Ud_list,self.evolution[i+1])
            print('Time step %i over %i : t = %f\n' % (i+1,self.nb_time_steps+1,self.curr_time))
            Unew, Vnew, Anew = self.next_step(Uold, Vold, Aold, K, C, M, F, L, Ud)
            U[:, i+1] = Unew
            V[:, i+1] = Vnew
            Uold, Vold, Aold = Unew, Vnew, Anew
        return U, V
