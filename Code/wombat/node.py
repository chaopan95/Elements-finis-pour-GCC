# -*- coding: utf-8 -*-
"""
.. module:: node

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

.. moduleauthor:: Jeremy Bleyer, Ecole des Ponts ParisTech, \
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205) <jeremy.bleyer@enpc.fr>
"""
import numpy as np

# id counter used to label nodes from 0
last_id = 0
class Node:
    """Abstract `Node` object (cannot be used as such)

    Attributes
    ----------
    _id : int
        node id number
    coor : ndarray   
        node coordinates
    """
    def __init__(self,coor=np.array([0,0])):
        if type(self) is Node:
            raise Exception("<Node> must be subclassed by Node2D or Node3D.")
        global last_id
        coor = np.array(coor)
        self._id = last_id
        last_id +=1
        self.coor = coor

    def get_dof(self,ndof):
        """ Return the dof list of the corresponding nodes
        
        Parameters
        ----------
        ndof : int
            number of dof/node
        """
        node_id = self._id
        return ndof*node_id + np.array(range(0,ndof))
    
class Node2D(Node):
    """ Two-dimensional node 

        Attributes
        ----------
        coor : ndarray
            coor = [x,y] with :math:`x` and :math:`y` coordinates 
        
           can init with a list for `coor`
    """
    def __init__(self,coor=np.array([0,0])):
        dim = len(coor)
        if dim != 2:
            raise ValueError('Coordinates of Node2D must be of length 2')
        Node.__init__(self,coor)
        self.dim = 2    
        
    def __str__(self):
        return "Noeud %i : coordonnes [%f,%f]" % (self._id,self.coor[0],self.coor[1])

    def copy(self,n=0):
        """ Creates a new node with the same coordinates but different id
        
              - if n=0      returns a single instance
              - otherwise   returns a list of length n
        """
        if n == 0:
            return Node2D(self.coor)
        else:
            copies = []
            for i in range(n):
                copies.append(Node2D(self.coor))
            return copies
        
class NodeGroup:
    """ Group of nodes 
    
    Attributes
    ----------
    nb_nodes : int
        number of nodes in the group
    node_list : list
        list of nodes
    coor : ndarray
        array of coordinates for the whole group, shape=(nb_nodes,2)
    """
    def __init__(self,node_list=[]):
        self.nb_nodes = len(node_list)
        self.node_list = node_list
        self.coor = np.zeros((self.nb_nodes,2))
        self.dim = 2
        for i in range(self.nb_nodes):
            assert self.dim==self.node_list[i].dim
            self.coor[i,:] = np.array(self.node_list[i].coor)        
    
    def add_node(self,node):
        """ Add `node` to the group """
        self.node_list.append(node)
        self.nb_nodes += 1
    
    def get_dof(self,ndof):
        """ Get an array of all node dofs, ndof = number of dof/node"""
        return np.array([n.get_dof(ndof) for n in self.node_list]).flatten()
        
    def get_id(self):
        """ Get an array of all nodes ids """
        id_array = np.zeros((self.nb_nodes,),dtype=int)
        for i in range(self.nb_nodes):
            id_array[i] = self.node_list[i]._id
        return id_array