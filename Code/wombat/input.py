# -*- coding: utf-8 -*-
"""
Module for generating and reading Gmsh mesh files

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
from .node import Node2D, NodeGroup
from .element import ElementGroup, Segment,Segment3,Triangle, Triangle6, Quadrangle, Quadrangle9
from .mesh import Mesh
import numpy as np
import os
import platform
try:
    import meshio
except:
    pass

def call_gmsh(geofile, elem_type):
    """ Calls ``Gmsh`` to generate a mesh

    Parameters
    ----------
    geofile : str
        path to GEO file, must end in ".geo"
    elem_type : {:class:`SolidT3`,:class:`SolidT6`}
        type of generated elements

    Returns
    -------
    regions
        a :class:`Mesh <mesh.Mesh>` object
    """
    if issubclass(elem_type, Triangle6) or issubclass(elem_type, Quadrangle9):
        order = 2
    else:
        order = 1
    if platform.system()=="Windows":
        commands = [r'"gmsh.exe" ', r'"C:/Program Files (x86)/gmsh/gmsh.exe" ']
    else:
        commands = ['gmsh ', '"/opt/gmsh-4.2.2-Linux64/bin/gmsh" ']
    for gmsh_command in commands:
        command = gmsh_command + '-format msh2 -order '+str(order)+' '+geofile+' -2'
        s = os.system(command)
        if s == 0:
            break
    print("Finished meshing!\n")
    mshfile = geofile[:-3]+'msh'
    return read_msh(mshfile,elem_type)

def read_msh(mshfile,elem_type):
    """Reads ``Gmsh`` 2.0 mesh files using ``meshio`` package

    Parameters
    ----------
    mshfile : str
        path to MSH file, must end in ".msh"

    Returns
    -------
    regions
        a :class:`Mesh <mesh.Mesh>` like object with abstract geometrical entities (Segment, Triangles, etc.) containing
        mesh cells, facets and nodes
    """
#    mesh = Mesh()

    points, cells, point_data, cell_data, field_data = meshio.read(mshfile)
    v = field_data.values()
    region_names = [None]*max(v)
    for keys, values in field_data.items():
        region_names[values-1] = keys
    def get_name(v):
        name = region_names[v-1]
        if name is not None:
            return name
        else:
            return v

    Nno = points.shape[0]
    nodes = NodeGroup([Node2D(points[i,:2]) for i in range(Nno)])
    if "line" in cells:
        cell_l = cells["line"]
        cell_data_l = cell_data["line"]["physical"]
    else:
        cell_l = np.empty(shape=(0,0))
    if "line3" in cells:
        cell_l3 = cells["line3"]
        cell_data_l3 = cell_data["line3"]["physical"]
    else:
        cell_l3 = np.empty(shape=(0,0))
    if "triangle" in cells:
        cell_T3 = cells["triangle"]
        cell_data_T3 = cell_data["triangle"]["physical"]
    else:
        cell_T3 = np.empty(shape=(0,0))
    if"triangle6" in cells:
        cell_T6 = cells["triangle6"]
        cell_data_T6 = cell_data["triangle6"]["physical"]
    else:
        cell_T6 = np.empty(shape=(0,0))
    if "quad" in cells:
        cell_Q4 = cells["quad"]
        cell_data_Q4 = cell_data["quad"]["physical"]
    else:
        cell_Q4 = np.empty(shape=(0,0))
    if "quad9" in cells:
        cell_Q9 = cells["quad9"]
        cell_data_Q9 = cell_data["quad9"]["physical"]
    else:
        cell_Q9 = np.empty(shape=(0,0))

    if issubclass(elem_type,Triangle) and not issubclass(elem_type,Triangle6):
        if (len(cell_T6)>0 or len(cell_l3)>0):
            raise(ValueError,"Trying to use 3-noded elements with 6-noded triangles or 3-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_T3[i,j]] for j in range(3)], tag=get_name(cell_data_T3[i])) for i in range(cell_T3.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l[i,j]] for j in range(2)], tag=get_name(cell_data_l[i])) for i in range(cell_l.shape[0])]
    elif issubclass(elem_type,Triangle6):
        if len(cell_T3)>0 or len(cell_l)>0:
            raise(ValueError,"Trying to use 6-nodes elements with 3-noded triangles or 2-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_T6[i,j]] for j in range(6)], tag=get_name(cell_data_T6[i])) for i in range(cell_T6.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l3[i,j]] for j in range(3)], tag=get_name(cell_data_l3[i])) for i in range(cell_l3.shape[0])]
    if issubclass(elem_type,Quadrangle) and not issubclass(elem_type,Quadrangle9):
        if (len(cell_Q9)>0 or len(cell_l3)>0):
            raise(ValueError,"Trying to use 4-noded elements with 9-noded quadranlges or 3-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_Q4[i,j]] for j in range(4)], tag=get_name(cell_data_Q4[i])) for i in range(cell_Q4.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l[i,j]] for j in range(2)], tag=get_name(cell_data_l[i])) for i in range(cell_l.shape[0])]
    elif issubclass(elem_type,Quadrangle9):
        if len(cell_Q4)>0 or len(cell_l)>0:
            raise(ValueError,"Trying to use 9-noded elements with 4-noded quadrangles or 2-noded lines in .msh file")
        else:
            elem_list = [elem_type([nodes.node_list[cell_Q9[i,j]] for j in range(9)], tag=get_name(cell_data_Q9[i])) for i in range(cell_Q9.shape[0])]
            bd_list = [elem_type.trace([nodes.node_list[cell_l3[i,j]] for j in range(3)], tag=get_name(cell_data_l3[i])) for i in range(cell_l3.shape[0])]
    if issubclass(elem_type,Segment):
            elem_list = [elem_type([nodes.node_list[cell_l[i,j]] for j in range(2)], tag=get_name(cell_data_l[i])) for i in range(cell_l.shape[0])]
            bd_list = []

    mesh = Mesh(elem_list)
    boundary = ElementGroup(bd_list)

    return mesh,boundary

