# -*- coding: utf-8 -*-
"""
Module for trace elements

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .generic_element import Segment,Segment3
import numpy as np

class TraceSolidT3(Segment):
    """ Trace elements of the :class:`SolidT3 <solidT3.SolidT3>` element

        corresponds to a segment with 2 nodes
    """
    def __init__(self,node_list=[],tag=1):
        Segment.__init__(self,node_list,tag)
        self.el_dof = 4
        self.node_dof = 2
        self.nb_stresses = 2  # 2 stress vector at 2 Gauss points in element

        self.ngauss = 2
        self.ag,self.wg = self.gauss_quadrature(self.ngauss)


    def shape_functions(self,xi):
        """ Return shape function along the reference edge coordinate :math:`\\xi\\in[-1;1]` """
        N = np.array([(xi+1)/2.,
                      (1-xi)/2.])
        return N

    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading on the edge"""
        L = self.measure()
        fx,fy,cz = el_force
        fe = np.zeros((4,))
        for i in range(self.ngauss):
            N = self.shape_functions(self.ag[i])
            fe += L/2.*self.wg[i]*np.kron(N,np.array([fx,fy]))
        return fe

class TraceSolidT6(Segment3):
    """ Trace elements of the :class:`SolidT6 <solidT6.SolidT6>` element

        corresponds to a segment with 3 nodes
    """
    def __init__(self,node_list=[],tag=1):
        Segment3.__init__(self,node_list,tag)
        self.el_dof = 6
        self.node_dof = 2
        self.nb_stresses = 4  # 2 stress vector at 2 Gauss points in element

        self.ngauss = 2
        self.ag,self.wg = self.gauss_quadrature(self.ngauss)


    def shape_functions(self,xi):
        """ Return shape function along the reference edge coordinate :math:`\\xi\\in[-1;1]` """
        N = np.array([xi*(xi+1)/2.,
                      xi*(xi-1)/2.,
                      1-xi**2])
        return N

    def elementary_distributed_forces(self,el_force):
        """ Elementary force vector for uniform distributed loading on the edge"""
        L = self.measure()
        fx,fy,cz = el_force
        fe = np.zeros((6,))
        for i in range(self.ngauss):
            N = self.shape_functions(self.ag[i])
            fe += L/2.*self.wg[i]*np.kron(N,np.array([fx,fy]))
        return fe
