# -*- coding: utf-8 -*-
"""
Module for class `Connections`

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

class Connections:
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
    def __init__(self):
        self.nb_relations = 0
        self.master_list = []
        self.slave_list = []
        self.components_list = []
        self.lin_rela_list = []
        self.imposed_value_list = []

    def add_relation(self,master,slave,comp=[0,1],lin_relation=[1.,-1.],imposed_value=0.):
        """ Add a relation between degrees of freedom to the `Connections` instance

            general form :math:`a_0\\underline{u}^{master}_j+
            \sum_{i=1}a_iu^{slave_i}_j = u^{imposed}_j` for component :math:`j`

            default is :math:`\\underline{u}_{master}=\\underline{u}_{slave}`

        Parameters
        ----------
        master : Node
            master node (maximum one)
        slave : Node, list of Nodes
            slave nodes
        comp : list
            components of master and slaves affected by the relation
            (default is [0,1]),

            0 : :math:`u_x`, 1 : :math:`u_y`, 2 : :math:`\\theta_z` in 2D
        lin_relation : list
            coefficients :math:`[a_0,a_1,\ldots,a_k]` of the linear relation
            between degrees of freedom,

            the first one corresponds to the master node, the others to the
            slave elements,

            if lin_rela contains only two values, the second is applied to all
            slave elements

            Example : is slave=[s1,s2]

                - if lin_relation = [a0,a1,a2] :math:`a_0 u^{master}+a_1 u^{s1} + a_2 u^{s2} = u_{imposed}`
                - if lin_relation = [b0,b1] :math:`b_0 u^{master}+b_1 u^{s1} + b_1 u^{s2} = u^{imposed}`

        imposed_value : float or list
            value of the relation right hand side :math:`u^{imposed}`

            if float, the same value is imposed for all components,
            otherwise the list must matches the list of components
        """
        if not isinstance(comp, list):
            comp = [comp]
        if not isinstance(slave, list):
            slave = [slave]
        if not isinstance(imposed_value, list):
            imposed_value = [imposed_value]*len(comp)
        self.master_list.append(master)
        self.slave_list.append(slave)
        self.components_list.append(comp)
        self.lin_rela_list.append(lin_relation)
        self.imposed_value_list.append(imposed_value)
        self.nb_relations += len(comp)*len(slave)

    def add_imposed_displ(self,location,ux=[],uy=[],thetaz=[]):
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
        self.Ux_n = append_component(ux,len(node_list),None)
        self.Uy_n = append_component(uy,len(node_list),None)
        self.Thetaz_n = append_component(thetaz,len(node_list),None)
        for i,node in enumerate(node_list):
            if self.Ux_n[i] is not None:
                self.add_relation(None,node,comp=0,lin_relation=[0,1.],imposed_value=self.Ux_n[i])
            if self.Uy_n[i] is not None:
                self.add_relation(None,node,comp=1,lin_relation=[0,1.],imposed_value=self.Uy_n[i])
            if self.Thetaz_n[i] is not None:
                self.add_relation(None,node,comp=2,lin_relation=[0,1.],imposed_value=self.Thetaz_n[i])
