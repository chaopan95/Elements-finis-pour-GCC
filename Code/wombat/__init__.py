"""
WomBat package

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr

Submodules
-----------
node
    manages `Node` objects
element
    package containing all finite element classes (Bar2D, Beam2D, ...)
mesh
    manages `Mesh` objects
material
    material constitutive relations
geometric_caract
    manages `Section` objects for `Bar` and `Beam` elements
model
    regroups mesh, material and section properties
connections
    manages imposed displacement and relations between nodes
forces
    manages loading conditions
finite_elements
    manages assembling procedures and linear solver
finite_elements_sparse
    same with sparse matrix format and eigenvalue solver
post_process
    manages plotting functions
input
    manages Gmsh file formats and converts to `Mesh` object
time_integration
    transient analysis with time integration schemes
nonlinear
    iterative procedures for nonlinear constitutive relations
"""
from .node import *
from .element import *
from .mesh import *
from .material import *
from .geometric_caract import *
from .model import *
from .connections import *
from .forces import *
from .finite_elements_sparse import *
from .post_process import *

from .input import *

from .time_integration import *

from .nonlinear import *