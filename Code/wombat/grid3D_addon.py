# -*- coding: utf-8 -*-
"""
Additional functions for using `Grid3D` elements

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .utils import append_component
from .element.generic_element import ElementGroup
from .node import NodeGroup
from .connections import Connections
from .forces import ExtForce
from .geometric_caract import BeamSection

class GridSection(BeamSection):
    """ Represents geometric properties of a bar/beam cross-section

    Attributes
    ----------
    area : float
        area :math:`S` of the cross-section
    inertia : float
        bending inertia :math:`I` for a torsional beam
    torsion : float
        torsional inertia :math:`J` for a torsional beam
    """
    def __init__(self,S=1.,I=0,J=0):
        self.area = S
        self.inertia = I
        self.torsion = J

class GridConnections(Connections):
    """ `Connections` objects are used to apply displacement boundary conditions
    and relations between degrees of freedom

    Attributes
    ----------
    nb_relations : int
        total number of relations between degrees of freedom
    master_list : list
        list of master elements for each relation (max one per relation)
    slave_list : list
        list of all slave elements for each relation (one or more per element)
    components_list : list
        list of degrees of freedom component affected by each relation
    lin_rela_list : list
        list of linear relation coefficients for each relation
    imposed_value_list : list
        list of imposed values for each relation

    See Also
    --------
    add_relation : for more information on the structure of one relation
    """

    def add_imposed_displ(self,location,w=[],thetax=[],thetay=[]):
        """ Imposes a given displacement to a list of nodes

        Parameters
        ----------
        location : :class:`Node`, list of Nodes, :class:`NodeGroup`, :class:`ElementGroup`
            node(s) on which displacement conditions are applied
        ux : float, list
            imposed value of horizontal displacement :math:`u_x`
        uy : float, list
            imposed value of vertical displacement :math:`u_y`
        thetaz : float, list
            imposed value of rotation :math:`\\theta_z` (only for `Beam2D` elements)

        .. note:: if one value only is provided, the same applies to all elements of the list

            use None to let the corresponding dof free
        """
        if isinstance(location,ElementGroup):
            node_list = location.get_nodes()
        elif isinstance(location,NodeGroup):
            node_list = location.node_list
        elif not isinstance(location, list):
            node_list = [location]
        else:
            node_list = location
        self.Ux_n = append_component(w,len(node_list),None)
        self.Uy_n = append_component(thetax,len(node_list),None)
        self.Thetaz_n = append_component(thetay,len(node_list),None)
        for i,node in enumerate(node_list):
            if self.Ux_n[i] is not None:
                self.add_relation(None,node,comp=0,lin_relation=[0,1.],imposed_value=self.Ux_n[i])
            if self.Uy_n[i] is not None:
                self.add_relation(None,node,comp=1,lin_relation=[0,1.],imposed_value=self.Uy_n[i])
            if self.Thetaz_n[i] is not None:
                self.add_relation(None,node,comp=2,lin_relation=[0,1.],imposed_value=self.Thetaz_n[i])


class GridExtForce(ExtForce):
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
    
    def add_distributed_forces(self,el_list,fz=[],cx=[],cy=[]):
        """ Adds a uniformly distributed force to a list of elements

        Parameters
        ----------
        el_list : :class:`Element <generic_element.Element>`, list of :class:`Elements <generic_element.Element>`, :class:`ElementGroup <generic_element.ElementGroup>`
            element(s) on which distributed forces are applied
        fz : float, list
            imposed value of transversal distributed force :math:`f_z`
        cx : float, list
            imposed value of x-distributed couple :math:`c_x`
        cy : float, list
            imposed value of y-distributed couple :math:`c_y` 

            .. note :: if one value only is provided, the same applies to all elements of the list
        
        """
        if isinstance(el_list,ElementGroup):
            el_list = el_list.elem_list
        if not isinstance(el_list, list): el_list = [el_list]
        self.el_list += el_list
        self.fx_e += append_component(fz,len(el_list))
        self.fy_e += append_component(cx,len(el_list))
        self.cz_e += append_component(cy,len(el_list))
        
    def add_concentrated_forces(self,node_list,Fz=[],Cx=[],Cy=[]):
        """ Adds a concentrated force to a list of nodes

        Parameters
        ----------
        node_list : :class:`Node <node.Node>`, list of  :class:`Nodes <node.Node>`, :class:`NodeGroup <node.NodeGroup>`
            node(s) on which concentrated forces are applied
        Fz : float, list
            imposed value of transversal concentrated force :math:`F_z`
        Cx : float, list
            imposed value of concentrated couple along x :math:`C_x`
        Cy : float, list
            imposed value of concentrated couple along y :math:`C_y` 

            .. note :: if one value only is provided, the same applies to all nodes of the list
        
        """
        if not isinstance(node_list, list): node_list = [node_list]
        self.node_list += node_list
        self.Fx_n += append_component(Fz,len(node_list))
        self.Fy_n += append_component(Cx,len(node_list))
        self.Cz_n += append_component(Cy,len(node_list))

