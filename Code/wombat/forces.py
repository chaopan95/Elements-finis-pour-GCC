# -*- coding: utf-8 -*-
"""
Module for class `ExtForce`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .element import ElementGroup
from .utils import append_component

class ExtForce:
    """ External force object for imposing loading conditions

    Attributes
    ----------
    node_list : list
        list of node groups on which concentrated forces are applied
    Fx_n : list
        list of corresponding horizontal concentrated forces components
    Fy_n : list
        list of corresponding vertical concentrated forces components
    Cz_n : list
        list of corresponding concentrated couples (for :class:`Beam2D <beam2D.Beam2D>` only)
    el_list : list
        list of element groups on which distributed forces are applied
    fx_e : list
        list of corresponding horizontal distributed forces components
    fy_e : list
        list of corresponding vertical distributed forces components
    cz_e : list
        list of corresponding distributed couples  (for :class:`Beam2D <beam2D.Beam2D>` only)

    """
    def __init__(self):
        self.node_list = []
        self.Fx_n = []
        self.Fy_n = []
        self.Cz_n = []
        self.el_list = []
        self.fx_e = []
        self.fy_e = []
        self.cz_e = []

    def add_distributed_forces(self,el_list,fx=[],fy=[],cz=[]):
        """ Adds a uniformly distributed force to a list of elements

        Parameters
        ----------
        el_list : :class:`Element <generic_element.Element>`, list of :class:`Elements <generic_element.Element>`, :class:`ElementGroup <generic_element.ElementGroup>`
            element(s) on which distributed forces are applied
        fx : float, list
            imposed value of horizontal distributed force :math:`f_x`
        fy : float, list
            imposed value of vertical distributed force :math:`f_y`
        cz : float, list
            imposed value of distributed couple :math:`c_z` (only for :class:`Beam2D <beam2D.Beam2D>` elements)

            .. note :: if one value only is provided, the same applies to all elements of the list

        """
        if isinstance(el_list,ElementGroup):
            el_list = el_list.elem_list
        if not isinstance(el_list, list): el_list = [el_list]
        self.el_list += el_list
        self.fx_e += append_component(fx,len(el_list))
        self.fy_e += append_component(fy,len(el_list))
        self.cz_e += append_component(cz,len(el_list))

    def add_concentrated_forces(self,node_list,Fx=[],Fy=[],Cz=[]):
        """ Adds a concentrated force to a list of nodes

        Parameters
        ----------
        node_list : :class:`Node <node.Node>`, list of  :class:`Nodes <node.Node>`, :class:`NodeGroup <node.NodeGroup>`
            node(s) on which concentrated forces are applied
        Fx : float, list
            imposed value of horizontal concentrated force :math:`F_x`
        Fy : float, list
            imposed value of vertical concentrated force :math:`F_y`
        Cz : float, list
            imposed value of concentrated couple :math:`C_z` (only for :class:`Beam2D <beam2D.Beam2D>` elements)

            .. note :: if one value only is provided, the same applies to all nodes of the list

        """
        if not isinstance(node_list, list): node_list = [node_list]
        self.node_list += node_list
        self.Fx_n += append_component(Fx,len(node_list))
        self.Fy_n += append_component(Fy,len(node_list))
        self.Cz_n += append_component(Cz,len(node_list))
