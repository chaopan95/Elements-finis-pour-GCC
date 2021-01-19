# -*- coding: utf-8 -*-
"""
Module for class `Section`

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
class Section:
    """ Abstract class for section properties"""
    pass


class BeamSection(Section):
    """ Represents geometric properties of a bar/beam cross-section

    Attributes
    ----------
    area : float
        area :math:`S` of the cross-section
    inertia : float
        bending inertia :math:`I` for a planar beam (unused for :class:`Bar2D <bar2D.Bar2D>` elements)
    """
    def __init__(self,S=1.,I=0):
        self.area = S
        self.inertia = I

    def rect(self,b,h=None):
        """ Builds a rectangular cross-section of dimensions :math:`b\\times h`

        :math:`S=bh` and :math:`I=bh^3/12`

        assumes a square cross-section if only one parameter is given
        """
        if h is None:  # square section
            h = b
        self.area = b*h
        self.inertia = b*h**3./12
        return self

class ShellSection(Section):
    def __init__(self,thick=0):
        self.thickness = thick