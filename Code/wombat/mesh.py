# -*- coding: utf-8 -*-
"""
Module for class `Mesh`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .node import NodeGroup
from .element import ElementGroup
import numpy as np

class Mesh(ElementGroup):
    """ Contains a group of elements

    Attributes
    ----------
    el_list : list
        list of elements composing the mesh
    nodes : :class:`NodeGroup <node.NodeGroup>`
        list of all nodes in the mesh
    Nno
        total number of nodes
    Nel
        total number of elements
    connec : ndarray
        connectivity matrix shape=(Nel,node/element)
    coor : ndarray
        coordinate matrix shape=(Nno,dim)
    """
    def __init__(self,el_list=[],elem_type=None):
        """
        Parameters
        ----------
        el_list : list
            list of elements composing the mesh
        elem_type : `Element`, optional
            if not None, converts all elements of `el_list` into elements of
            type `elem_type`
        """
        if isinstance(el_list,ElementGroup):
            el_list = el_list.elem_list
        if elem_type != None:
            el_list = [elem_type(x.nodes.node_list,x.physical_group) for x in el_list]
        elif len(el_list) > 0:
            elem_type = el_list[0].__class__
        else:
            elem_type = None
        ElementGroup.__init__(self,el_list)

        self.nodes = NodeGroup(self.get_nodes())
        # renumber all nodes inside the mesh
        for i,n in enumerate(self.nodes.node_list):
            n._id = i
        self.Nno = self.nodes.nb_nodes
        self.Nel = self.nb_elem
        self.connec = np.array([x.nodes.get_id() for x in self.elem_list])
        self.coor = np.array([x.coor for x in self.nodes.node_list])

    def volume(self):
        """Gives the total area/volume of a mesh"""
        return np.sum(self.measure())

    def print_info(self):
        print("Number of nodes : %i\nNumber of elements : %i\n" % (self.Nno,self.Nel))




