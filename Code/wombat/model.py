# -*- coding: utf-8 -*-
"""
Module for class `Model`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .material import Material
from .geometric_caract import Section
from .element.generic_element import ElementGroup

class Model:
    """ Contains a mesh, material and (possibly) section properties

    Attributes
    ----------
    mesh
        :class:`Mesh <mesh.Mesh>` object
    mat
        :class:`Material <material.Material>` object
    sect
        :class:`Section <geometric_caract.Section>` object (optional)
    """
    def __init__(self,mesh,mat,sect=Section()):
        self.mesh = mesh
        self.mat = mat
        self.sect = sect
        self.affect_property(mat)
        self.affect_property(sect)

    def affect_property(self,prop,el_list=[],tag=None):
        """ Affect property (material or section) to a specific list of elements

        Parameters
        ----------
        prop : {:class:`Material <material.Material>`, :class:`Section <geometric_caract.Section>`}
            property to affect, either a material or a section
        el_list : list, :class:`ElementGroup`
            list of elements on which `prop` is affected

            if list is empty, all elements of the mesh are selected
        tag : None, int, string
            if tag is not None, only the elements of `el_list`
            with `tag` as a physical region marker are selected
        """
        if isinstance(el_list,ElementGroup):
            el_list = el_list.elem_list

        assert ((len(el_list) == 0) or (tag is None)), "Element list and tag cannot be affected at the same time"

        if len(el_list)==0 and tag is None:
            to_do_list = self.mesh.elem_list
        elif len(el_list)!=0:
            to_do_list = el_list
        elif tag is not None:
            to_do_list = self.mesh.get_elem_from_tag(tag).elem_list

        for e in to_do_list:
            if isinstance(prop,Material):
                e.mat = prop
            elif isinstance(prop,Section):
                e.sect = prop
            else:
                raise ValueError('Unknown property type')


