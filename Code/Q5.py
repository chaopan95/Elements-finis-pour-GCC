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
q = 1e6

b=0.1
h=0.1
I = b*h**3/12
J = b*h**3/12+b**3*h/12
S = b*h

Nx=Ny=10
with open(r"geometry/rect_grid.geo","r+") as file:
    text = file.read()
    text = re.sub("\nNy = \d*;","\nNy = %s;"%Nx,text)
    file.seek(0)
    file.write(text)
    file.close()

sect = GridSection(S, I, J)   
mesh, boundary = call_gmsh('geometry/rect_grid.geo',Grid3D) #nombre de poutre égale à 9
horz = mesh.get_elem_from_tag("horz")
vert = mesh.get_elem_from_tag("vert")
bottom = mesh.get_elem_from_tag("bottom")
right = mesh.get_elem_from_tag("right")
top = mesh.get_elem_from_tag("top")
left = mesh.get_elem_from_tag("left")
nF = right.get_nodes()[-1]
mat = LinearElastic(E, nu, model="plane_stress", rho=rho)   
model = Model(mesh,mat,sect)

'''
Valider la solution sur un cas-test simple pour lequel la flexion se
produit uniquement dans une direction. 
On fixe le déplacement vertical du côté left et right, et on n'applique qu'une 
force distribuée verticale sur les éléments "horz", "top" et "bottom".
Dans ce cas-là, cette plaque est équivalente à une poutre simplement appuyée. 
De ce fait, on peut vérifer notre solution en comparant le déplacement au 
mi-travée et l'angle d'un extrémité.
'''

appuis = GridConnections()
appuis.add_imposed_displ(right, w=0)
appuis.add_imposed_displ(left, w=0)

forces = GridExtForce()
forces.add_distributed_forces(horz, fz=-q)
forces.add_distributed_forces(top, fz=-q)
forces.add_distributed_forces(bottom, fz=-q)

K = assembl_stiffness_matrix(model)
L,Ud = assembl_connections(appuis,model)
F = assembl_external_forces(forces,model)

U,lamb = solve(K,F,L,Ud)
w = U[::3]

thetaX = U[1::3]
thetaY = U[2::3]
fig = Figure(1,"Displacement w - cas simple")
fig.plot_bc(mesh, appuis)
fig.plot_field_lines(mesh, w, scale=1.)
fig.show()

fleche_EF = max(abs(w))
fleche_theorique = 5*q*L_plaque**4/(384*E*I)
error = (fleche_EF-fleche_theorique)/fleche_theorique
print("\nFlèche à mi-travée de grillage: %f" % (fleche_EF))
print("Flèche à mi-travée de poutre simplement appuyée: %f\n" % (fleche_theorique))
print("Erreur de la flèche à mi-travée: %f%%\n" % (error))
# Les erreur sont quasiment 0, donc, on peut dire que notre code est bon