# -*- coding: utf-8 -*-
"""
Module for class `Material`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
from math import sqrt

class Material:
    """ Abstract class for material properties"""
    pass

#class LinearElastic(Material):
#    """ Linear elastic material for bars
#
#    Attributes
#    ----------
#    Young_modulus : float
#        Material Young modulus :math:`E`
#    """
#    def __init__(self,E=1e6):
#        self.Young_modulus = E
#

class LinearElastic(Material):
    """ Linear elastic material

    Attributes
    ----------
    Young_modulus : float
        material Young modulus :math:`E`
    Poisson_coeff : float
        material Poisson coefficient :math:`\\nu` (with :math:`-1<\\nu<1/2`),
        ignored for :class:`Bar2D <bar2D.Bar2D>`  and :class:`Beam2D <beam2D.Beam2D>`  elements
    rho : float
        material volumetric mass density :math:`\\rho`
    model : {'plane_strain','plane_stress'}
        type of 2D model
    C : ndarray
        elasticity matrix :math:`[C]` shape=(3,3)
    """
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain"):
        assert (nu<=0.5) and (nu>-1), "Wrong Poisson coefficient"
        self.Young_modulus = E
        self.Poisson_coeff = nu
        self.rho = rho
        self.model = model
        self.C = self.C_matrix()

    def compute_lame_coeff(self):
        """Returns Lamé coefficients :math:`\lambda,\mu`"""
        E = self.Young_modulus
        nu = self.Poisson_coeff
        lamb = E*nu/(1+nu)/(1-2*nu)
        mu = E/2./(1+nu)
        return lamb,mu

    def from_lame(self,lamb,mu):
        self.Young_modulus = mu*(3*lamb+2*mu)/(lamb+mu)
        self.Poisson_coeff = lamb/2./(lamb+mu)

    def C_matrix(self):
        """Compute elasticity matrix :math:`[C]`"""
        if self.model == "plane_strain":
            lamb,mu = self.compute_lame_coeff()
            return np.array([[lamb+2*mu,lamb,0],[lamb,lamb+2*mu,0],[0,0,mu]])
        elif self.model == "plane_stress":
            E = self.Young_modulus
            nu = self.Poisson_coeff
            return E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2.]])

    def compute_sigzz(self,Eps):
        if self.model == "plane_strain":
            lamb,mu = self.compute_lame_coeff()
            return lamb*(Eps[0]+Eps[1])
        elif self.model == "plane_stress":
            return 0.

class ReinforcedConcrete(LinearElastic):
    def __init__(self,Ec,nuc,Es,Phi=[],es=[]):
        self.Ec = Ec
        self.nuc = nuc
        self.Es = Es
        self.ns = len(Phi)
        self.Phi = Phi
        self.es = np.array(es)
        self.C = self.C_matrix()

    def C_matrix(self):
        E = self.Ec
        nu = self.nuc
        C = E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2.]])
        for i in range(self.ns):
            c = self.es[i,0]
            s = self.es[i,1]
            R = np.array([[c**2,s**2,0],[s**2,c**2,0],[0,0,1]])
            C_steel_loc = self.Phi[i]*self.Es*np.diag([1,0,0])
            C_steel_glob = R.T.dot(C_steel_loc.dot(R))
            C += C_steel_glob
        return C

    def compute_sigzz(self,Eps):
        return 0.

class ElastoPlastic(LinearElastic):
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain",Et=0,fy=1e3):
        LinearElastic.__init__(self,E,nu,rho,model)
        self.tangent_modulus = Et
        self.hardening_modulus = E*Et/(E-Et)
        self.yield_stress = fy
        self.nvar = 2

    def constitutive_relation(self,dEps,Sig,X):
        E = self.Young_modulus
        Et = self.tangent_modulus
        H = self.hardening_modulus
        fy = self.yield_stress
        Sig_elas = Sig+E*dEps
        p = X[1]
        fyield = abs(Sig_elas) - (fy + H*p)
        if fyield < 0:
            new_Sig = Sig_elas
            dX = 0
            Ct = E
        else:
            sign = dEps/abs(dEps)
#            dp = fyield/(H+E)
            dp = fyield*(E-Et)/E**2
            dEps_p = sign*dp
            new_Sig = Sig_elas - E*dEps_p
            dX = [dEps_p,dp]
            Ct = max(Et,1e-6*E)
        return new_Sig,dX,Ct


class vonMises_2D(ElastoPlastic):
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain",Et=0,fy=1e3):
        ElastoPlastic.__init__(self,E,nu,rho,model,Et,fy)
        self.nvar = 5

#    def constitutive_relation(self,dEps,Sig,X):
#        C = self.C_matrix()
#        lamb,mu = self.compute_lame_coeff()
#        H = self.hardening_modulus
#        fy = self.yield_stress
#        # projector to deviatoric space in 2D (sigxx,sigyy,sigzz,sigxy)
#        K = np.array([[2/3.,-1/3.,-1/3.,0],
#                      [-1/3.,2/3.,-1/3.,0],
#                        [-1/3.,-1/3.,2/3.,0],
#                        [0,0,0,1/2.]])
#        norm = lambda x: sqrt(x[0]**2+x[1]**2+x[2]**2+2*x[3]**2)
#        sig_eq = lambda x: sqrt(1./2.)*norm(x)
#
#        dSig_plane = np.dot(C,dEps)
#        Sig_elas = Sig+np.array([dSig_plane[0],dSig_plane[1],lamb*sum(dEps[0:2]),dSig_plane[2]])
#        s_elas = Sig_elas - 1/3.*np.sum(Sig_elas[:3])*np.array([1,1,1,0])
#
#        alp = X[-1]
#        fyield = sig_eq(s_elas) - (fy/sqrt(3) + H*alp)
#        if fyield < 0:
#            new_Sig = Sig_elas
#            dX = np.zeros((self.nvar,))
#            Ct = C
#        else:
#            dalp = fyield/(H+mu)
#            n_el = s_elas/sqrt(2)/sig_eq(s_elas)
#
#            dEps_p = dalp*n_el/sqrt(2.)
#            new_Sig = Sig_elas - 2*mu*dEps_p
#            dEps_p[3] *= 2
#
#            dX = np.append(dEps_p,dalp)
#            D = 2*mu**2/(mu+H)*np.outer(n_el,n_el) + \
#                2*mu**2*dalp/sig_eq(s_elas)*(K-np.outer(n_el,n_el))
#            Ct = C - D[np.ix_([0,1,3],[0,1,3])]
#        return new_Sig,dX,Ct

    def constitutive_relation(self,dEps,Sig,X):
        C = self.C_matrix()
        lamb,mu = self.compute_lame_coeff()
        # projector to deviatoric space in 2D (sigxx,sigyy,sigzz,sigxy)
        P = np.array([[2/3.,-1/3.,-1/3.,0],
                      [-1/3.,2/3.,-1/3.,0],
                        [-1/3.,-1/3.,2/3.,0],[0,0,0,1/2.]])
        norm = lambda x: (x[0]**2+x[1]**2+x[2]**2+2*x[3]**2)**0.5
        sig_eq = lambda x: (3./2.)**0.5*norm(x)

        H = self.hardening_modulus
        fy = self.yield_stress
        dSig_plane = np.dot(C,dEps)
        Sig_elas = Sig+np.array([dSig_plane[0],dSig_plane[1],lamb*sum(dEps[0:2]),dSig_plane[2]])
        s_elas = Sig_elas - 1/3.*np.sum(Sig_elas[:3])*np.array([1,1,1,0])

        p = X[-1]
        fyield = sig_eq(s_elas) - (fy + H*p)
        if fyield < 0:
            new_Sig = Sig_elas
            dX = np.zeros((self.nvar,))
            Ct = C
        else:
            dp = fyield/(H+3*mu)
            beta = 3*mu*dp/sig_eq(s_elas)
            gamma = 3*mu/(H+3*mu)
            dEps_p = beta/2./mu*s_elas
            new_Sig = Sig_elas - 2*mu*dEps_p
            dEps_p[3] *= 2

            dX = np.append(dEps_p,dp)
            n_elas = s_elas/sig_eq(s_elas)
            D = 3*mu*(gamma-beta)*np.outer(n_elas,n_elas)+2*mu*beta*P
            Ct = C - D[np.ix_([0,1,3],[0,1,3])]
        return new_Sig,dX,Ct

class Cable(LinearElastic):
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain"):
        LinearElastic.__init__(self,E,nu,rho,model)
        self.nvar = 1

    def constitutive_relation(self,dEps,Sig,X=0):
        E = self.Young_modulus
        Sig_elas = Sig+E*dEps
        if Sig_elas <= 0:
            new_Sig = 0
            Ct = 0*E
        else:
            new_Sig = Sig_elas
            Ct = E
        return new_Sig,dEps,Ct


class ElastoPlastic_bar(LinearElastic):
    def __init__(self,E=1e6,nu=0.,rho=0.,model="plane_strain",fy=1e3):
        LinearElastic.__init__(self,E,nu,rho,model)
        self.yield_stress = fy
        self.nvar = 1

    def constitutive_relation(self,dEps,Sig,X=0):
        E = self.Young_modulus
        fy = self.yield_stress
        Sig_elas = Sig+E*dEps

        if abs(Sig_elas) < fy:
            new_Sig = Sig_elas
            Ct = E
        else:
            sign = Sig_elas/abs(Sig_elas)
            new_Sig = sign*fy
            Ct = 0*E
        return new_Sig,0*X,Ct