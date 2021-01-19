# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:22:23 2020

@author: ENPC\ZIJIE LI et Chao PAN
"""

from wombat import *
from wombat.element.grid3D import *
from wombat.grid3D_addon import *
import re

L_plaque = 1
H_plaque = 1
E= 210e9
nu = 0.3
rho = 2700 
q = 10e6

b=0.1
h=0.1
I = b*h**3/12
J = b*h**3/12+b**3*h/12
S = b*h

Nx=Ny=2
with open(r"geometry/rect_grid.geo","r+") as file:
    text = file.read()
    text = re.sub("\nNy = \d*;","\nNy = %s;"%Ny,text)
    file.seek(0)
    file.write(text)
    file.close()

sect = GridSection(S, I, J)   
mat = LinearElastic(E, nu, model="plane_stress", rho=rho)   

mesh, boundary = call_gmsh('geometry/rect_grid.geo',Grid3D)
horz = mesh.get_elem_from_tag("horz")
vert = mesh.get_elem_from_tag("vert")
bottom = mesh.get_elem_from_tag("bottom")
right = mesh.get_elem_from_tag("right")
top = mesh.get_elem_from_tag("top")
left = mesh.get_elem_from_tag("left")
nF = right.get_nodes()[-1]
model = Model(mesh,mat,sect)

appuis = GridConnections()
appuis.add_imposed_displ(bottom, w=0)
appuis.add_imposed_displ(right, w=0)
appuis.add_imposed_displ(top, w=0)
appuis.add_imposed_displ(left, w=0)

forces = GridExtForce()
forces.add_distributed_forces(horz, fz=-q)
forces.add_distributed_forces(vert, fz=-q)
K = assembl_stiffness_matrix(model)
L,Ud = assembl_connections(appuis,model)
F = assembl_external_forces(forces,model)

U,lamb = solve(K,F,L,Ud)
w = U[::3]
thetax = U[1::3]
fleche_EF = max(abs(w))
fig = Figure(1,"Displacement w")
fig.plot_bc(mesh, appuis)
fig.plot_field_lines(mesh, w, scale=1.)
fig.show()
print("\nFlèche aux points conjoints de grillage: %f\n" %fleche_EF)

Sigma = stresses(U,model)
M1 = Sigma[0::4]
T = Sigma[1::4]
M2 = Sigma[2::4]
V = Sigma[3::4]
fig = Figure(2,r"$\sigma{xx}$ stress")
fig.plot_field_lines(mesh,M1)
fig.show()
