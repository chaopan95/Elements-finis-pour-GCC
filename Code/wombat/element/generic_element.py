# -*- coding: utf-8 -*-
"""
Module for generic element classes

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from wombat.node import NodeGroup
from wombat.utils import uniquify
import numpy as np
from math import sqrt

last_id = 0
class Element:
    """Abstract `Element` object (cannot be used as such)

    Attributes
    ----------
    _id : int
        element id number
    nodes : :class:`NodeGroup <node.NodeGroup>`
        list of nodes attached to the element
    nb_nodes : int
        number of nodes
    physical_group : int,str
        physical group to which the element belongs
    el_dof : int
        number of dof/element
    node_dof : int
        number of dof/node
    nb_stresses : int
        number of stress variables/element
    """
    def __init__(self,node_list=[],tag=1):
        global last_id
        if type(self) is Element:
            raise Exception("<Element> must be subclassed (Bar2D,Triangle,...).")
        self._id = last_id
        last_id +=1

        self.nodes = NodeGroup(node_list)
        self.nb_nodes = self.nodes.nb_nodes
        self.physical_group = "".join(str(tag).splitlines()) # remove unwanted characters when reading with meshio
        self.el_dof = 0
        self.node_dof = 0
        self.nb_stresses = 0


    def get_dof(self):
        """ Returns list of dofs for nodes in the element """
        node = self.nodes.node_list
        ndof = self.node_dof
        nb_nodes = len(node)
        assert self.el_dof == ndof*nb_nodes, "Element has internal dofs"
        return np.array([ndof*node[i]._id+j for i in range(nb_nodes) for j in range(ndof) ])

    def barycenter(self):
        """ Compute coordinates of the element barycenter """
        return np.mean(self.node_coor(),axis=0)

    def node_coor(self):
        """ Returns array of nodal coordinates """
        return self.nodes.coor

    def __str__(self):
        return "Element ("+self.__class__.__name__ +") numero %i Tag %s\n   Nombre de noeuds %i\n   Barycentre [%f,%f]\n" % (self._id,str(self.physical_group),self.nb_nodes,self.barycenter()[0],self.barycenter()[1])


class ElementGroup:
    """ Group of elements (must be of the same kind)

    Attributes
    ----------
    nb_elem : int
        number of elements in the group
    eleme_list : list
        list of elements
    physical_group : int,str
        list of physical groups for each element
    el_dof : int
        number of dof/element
    node_dof : int
        number of dof/node
    nb_stresses : int
        number of stress variables/element
    """
    def __init__(self,el_list=[]):
        self.nb_elem = len(el_list)
        self.elem_list = el_list
        self.physical_group = []
        if len(el_list) != 0:
            assert all([isinstance(e,el_list[0].__class__) for e in el_list]), "Element group contains different kinds of elements"
            for i,e in enumerate(self.elem_list):
                self.physical_group.append(e.physical_group)
            self.el_dof = self.elem_list[0].el_dof
            self.node_dof = self.elem_list[0].node_dof
            self.nb_stresses = self.elem_list[0].nb_stresses
            if hasattr(self.elem_list[0],"ngauss"):
                self.ngauss = self.elem_list[0].ngauss
            else:
                self.ngauss = 1

    def get_elem_from_tag(self,tag):
        """ Returns elements with physical_group given by `tag`"""
        return ElementGroup([x for x in self.elem_list if x.physical_group==str(tag)])

    def get_nodes(self):
        """ Returns list of nodes inside the group """
        return uniquify([n for e in self.elem_list for n in e.nodes.node_list])
#        for e in self.elem_list:
#            for n in e.nodes.node_list:
#                if n._id not in [x._id for x in nodes]:
#                    nodes.append(n)
#        return nodes



    def __add__(self,other):
        return ElementGroup(self.elem_list+other.elem_list)

    def __str__(self):
        s = ""
        for i in range(self.nb_elem):
            s += self.elem_list[i].__str__()+"\n"
        return s

    def measure(self):
        """ Returns list of element measures """
        measure=[]
        for e in self.elem_list:
            measure.append(e.measure())
        return measure

    def add_element(self,el):
        """ Adds an element to the group """
        self.elem_list.append(el)
        self.nb_elem += 1

    def edge_list(self):
        node2edges = {}
        el2edges = []
        edge2nodes = []
        edg_map = [[0,1],[1,2],[2,0]]
        num_edge = 0
        for e in self.elem_list:
            el2edges.append([])
            if isinstance(e, Triangle):
                for i in range(3):  # loop over edges
                    na = e.nodes.node_list[edg_map[i][0]]
                    nb = e.nodes.node_list[edg_map[i][1]]
                    if na._id > nb._id: # cognvention for ordering element in edges is first one with lower _id
                        nb, na = na, nb
                    label = str(na._id)+"-"+str(nb._id)
                    if not node2edges.has_key(label):
                        node2edges[label] = num_edge
                        num_edge += 1
                        edge2nodes.append([na,nb])
                    el2edges[-1].append(node2edges[label])
            else:
                print("Ignoring edges for non-Triangle elements")
        return el2edges,node2edges,edge2nodes

class Segment(Element):
    """ A uniaxial element with two end-nodes """
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing two nodes
        tag : int,str
            tag of physical group
        """
        assert len(node_list)==2, "Segment must have 2 nodes."
        Element.__init__(self,node_list,tag)

    def measure(self):
        """Segment length"""
        return np.linalg.norm(self.node_coor()[1,:]-self.node_coor()[0,:])

    def hsize(self):
        """Segment length"""
        return self.measure()

#    def trace(self):
#        """No trace element defined for Segments"""
#        pass
#
    def gauss_quadrature(self,ngauss):
        """ Gauss quadrature on a segment

        Parameters
        ----------
        ngauss : {1,2}
            number of Gauss points

        Returns
        -------
        ag : ndarray
            position of Gauss points in the reference segment shape=(`ngauss`,)
        wg : ndarray
            weights of corresponding Gauss points shape=(`ngauss`,)
        """
        if (ngauss==1):
            ag = np.array([0])
            wg = np.array([2.])
        elif (ngauss==2):
            ag = np.array([-1/sqrt(3),1/sqrt(3)])
            wg = np.array([1.,1.])
        else:
            raise(ValueError,"Wrong number of Gauss points inside a Triangle")
        return ag,wg

class Segment3(Segment):
    """ A uniaxial element with two end-nodes + one mid-node """
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing three nodes, mid-node is the last one
        tag : int,str
            tag of physical group
        """
        assert len(node_list)==3, "Segment3 must have 3 nodes."
        Element.__init__(self,node_list,tag)

class Triangle(Element):
    """ A triangular element with three nodes """
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing three nodes, mid-node is the last one
        tag : int,str
            tag of physical group
        """
        assert len(node_list)==3, "Triangle must have 3 nodes."
        Element.__init__(self,node_list,tag)

    def measure(self):
        """Triangle area"""
        T = self.node_coor()
        J = np.array([T[1,:]-T[0,:],T[2,:]-T[0,:]])
        return abs(np.linalg.det(J)/2)

    def hsize(self,rad="circ"):
        """Typical size of the element computed as

        - if rad ="circ" : circumscribed circle radius
        - if rad ="insc" : inscribed circle radius
        """
        T = self.node_coor()
        a = np.linalg.norm(T[1,:]-T[0,:])
        b = np.linalg.norm(T[2,:]-T[1,:])
        c = np.linalg.norm(T[0,:]-T[2,:])
        semi_perimeter = (a+b+c)/2.
        if rad == "insc":
            return self.measure()/semi_perimeter
        elif rad == "circ":
            return a*b*c/((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))**0.5

#    def trace(self):
#        """ Trace element for :class:`Triangle` is :class:`Segment`"""
#        return Segment
#
    def gauss_quadrature(self,ngauss):
        """ Gauss quadrature on a triangle

        Parameters
        ----------
        ngauss : {1,3,6}
            number of Gauss points

        Returns
        -------
        ag : ndarray
            position of Gauss points in the reference triangle shape=(`ngauss`,2)
        wg : ndarray
            weights of corresponding Gauss points shape=(`ngauss`,)
        """
        if (ngauss==1):
            ag = np.array([[1/3.,1/3.]])
            wg = np.array([1/2.])
        elif (ngauss==3):
            ag = np.array([[1/6.,1/6.],[2/3.,1/6.],[1/6.,2/3.]])
            wg = np.array([1/6.,1/6.,1/6.])
        elif (ngauss==4):
            ag = np.array([[1/3.,1/3.],[3/5.,1/5.],[1/5.,3/5.],[1/5.,1/5.]])
            wg = np.array([-27/96.,25/96.,25/96.,25/96.])
        elif (ngauss==6):
            a1 = 0.445948490915965
            a2 = 0.091576213509771
            ag = np.array([[a1,a1],
                           [1-2*a1,a1],
                           [a1,1-2*a1],
                           [a2,a2],
                           [1-2*a2,a2],
                           [a2,1-2*a2]])
            w1 = 0.111690794839005
            w2 = 0.054975871827661
            wg = np.concatenate((w1*np.ones((3,)),w2*np.ones((3,))))
        else:
            raise(ValueError,"Wrong number of Gauss points inside a Triangle")
        return ag,wg

class Triangle6(Triangle):
    def __init__(self,node_list,tag=1):
        assert len(node_list)==6, "Triangle6 must have 6 nodes."
        Element.__init__(self,node_list,tag)
        self.trace = Segment3
    def measure(self):
        return Triangle.measure(self)

    def hsize(self,rad="circ"):
        return Triangle.hsize(self,rad)

class Quadrangle(Element):
    """ A quadrangular element with four nodes"""
    def __init__(self,node_list,tag=1):
        """
        Parameters
        ----------

        node_list : list
            list containing four nodes
        tag : int,str
            tag of physical group
        """
        assert len(node_list)==4, "Quadrangle must have 4 nodes."
        Element.__init__(self,node_list,tag)

    def measure(self):
        """Quadrangle area"""
        T = self.node_coor()
        J1 = np.array([T[1,:]-T[0,:],T[3,:]-T[0,:]])
        J2 = np.array([T[3,:]-T[2,:],T[1,:]-T[2,:]])
        return abs((np.linalg.det(J1)+np.linalg.det(J2))/2)

    def hsize(self,rad="circ"):
        """Typical size of the element computed as

        - if rad ="circ" : circumscribed circle radius
        - if rad ="insc" : inscribed circle radius
        """
        T = self.node_coor()
        a = np.linalg.norm(T[1,:]-T[0,:])
        b = np.linalg.norm(T[2,:]-T[1,:])
        c = np.linalg.norm(T[0,:]-T[2,:])
        semi_perimeter = (a+b+c)/2.
        if rad == "insc":
            return self.measure()/semi_perimeter
        elif rad == "circ":
            return a*b*c/((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))**0.5

#    def trace(self):
#        """ Trace element for :class:`Triangle` is :class:`Segment`"""
#        return Segment
#
    def gauss_quadrature(self,ngauss):
        """ Gauss quadrature on a triangle

        Parameters
        ----------
        ngauss : {1,3,6}
            number of Gauss points

        Returns
        -------
        ag : ndarray
            position of Gauss points in the reference triangle shape=(`ngauss`,2)
        wg : ndarray
            weights of corresponding Gauss points shape=(`ngauss`,)
        """
        if (ngauss==1):
            ag = np.array([[0.,0.]])
            wg = np.array([4.])
        elif (ngauss==4):
            ag = np.array([[-1/sqrt(3),-1/sqrt(3)],
                            [1/sqrt(3),-1/sqrt(3)],
                            [1/sqrt(3),1/sqrt(3)],
                            [-1/sqrt(3),1/sqrt(3)]])
            wg = np.ones((4,))
        elif (ngauss==9):
            ag = np.array([[-sqrt(3/5.),-sqrt(3/5.)],
                            [0,-sqrt(3/5.)],
                            [sqrt(3/5.),-sqrt(3/5.)],
                            [-sqrt(3/5.),0],
                            [0,0],
                            [sqrt(3/5.),0],
                            [-sqrt(3/5.),sqrt(3/5.)],
                            [0,sqrt(3/5.)],
                            [sqrt(3/5.),sqrt(3/5.)]])
            wg1d = np.array([5./9.,8./9.,5./9.])
            wg = np.kron(wg1d,wg1d)
        else:
            raise(ValueError,"Wrong number of Gauss points inside a Quadrangle")
        return ag,wg

class Quadrangle9(Quadrangle):
    def __init__(self,node_list,tag=1):
        assert len(node_list)==9, "Quadrangle9 must have 9 nodes."
        Element.__init__(self,node_list,tag)
        self.trace = Segment3
    def measure(self):
        return Quadrangle.measure(self)

    def hsize(self,rad="circ"):
        return Quadrangle.hsize(self,rad)