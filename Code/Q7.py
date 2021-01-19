# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:22:23 2020

@author: ENPC\ZIJIE LI et Chao PAN
"""

from wombat import *
from wombat.element.grid3D import *
from wombat.grid3D_addon import *
import re 
from math import *
import numpy as np
# Paramètres de plaque
h_plaque = 0.1 # Epaisseur de plaque
q_plaque = 1e6
E= 210e9
nu = 0.3
rho = 2700.

def grid_fleche(Nx, L_plaque=1, H_plaque=1):
    """
    L_plaque: Largeur de plaque
    H_plaque: Longueur de plaque
    """
    Ny = H_plaque/L_plaque*Nx
    with open(r"geometry/rect_grid_Q7.geo","r+") as file:
        text = file.read()
        text = re.sub("\nNy = \d*;","\nNy = %s;"%int(Ny),text)
        text = re.sub("\nNx = \d*;","\nNx = %s;"%int(Nx),text)
        text = re.sub("\nL = \d*.;","\nL = %s;"%int(L_plaque),text)
        text = re.sub("\nH = \d*.;","\nH = %s;"%int(H_plaque),text)
        file.seek(0)
        file.write(text)
        file.close()
    # Determination des paramètres de grillage selon les paramètres de plaque
    b=L_plaque/Nx # Largeur de poutre
    h=h_plaque # Hauteur de poutre 
    S = b*h
    I = b*h**3/12
    J = b*h**3/12+b**3*h/12
    sect = GridSection(S, I, J)
    
    q = q_plaque/(Nx/L_plaque+Ny/H_plaque)
    
    mesh, boundary = call_gmsh('geometry/rect_grid_Q7.geo',Grid3D)
    horz = mesh.get_elem_from_tag("horz")
    vert = mesh.get_elem_from_tag("vert")
    bottom = mesh.get_elem_from_tag("bottom")
    right = mesh.get_elem_from_tag("right")
    top = mesh.get_elem_from_tag("top")
    left = mesh.get_elem_from_tag("left")
    nF = right.get_nodes()[-1]
    mat = LinearElastic(E, nu, model="plane_stress", rho=rho)   
    model = Model(mesh,mat,sect)
    
    appuis = GridConnections()
    appuis.add_imposed_displ(bottom, w=0)
    appuis.add_imposed_displ(top, w=0)
    appuis.add_imposed_displ(right, w=0)
    appuis.add_imposed_displ(left, w=0)
    
    forces = GridExtForce()
    forces.add_distributed_forces(horz, fz=-q)
    forces.add_distributed_forces(vert, fz=-q)
    K = assembl_stiffness_matrix(model)
    L,Ud = assembl_connections(appuis,model)
    F = assembl_external_forces(forces,model)
    U,lamb = solve(K,F,L,Ud)    
    return U, mesh, appuis

def plaque_fleche(a, b, h=h_plaque, nu=0.3):
    '''Calculer la fleche de la plaque'''
    D = E*h**3/(12*(1-nu**2))
    coef1 = 4*q_plaque*a**4/D/pi**5
    coef2 = 0
    for i in range(1, 10):
        m = 2*i-1
        am = m*pi*b/(2*a)
        coef2+=(1-(2+am*tanh(am))/2/cosh(am)*cosh(0))*sin(m*pi/2)/m**5
    f_plaque = coef1*coef2
    return f_plaque

U, mesh, appuis = grid_fleche(10)
w = U[::3]
fig = Figure(1,"Displacement w")
fig.plot_bc(mesh, appuis)
fig.plot_field_lines(mesh, w, scale=1.)
fig.show()
print("Flèche de grillage carré: %f\n" % (max(abs(w))))
print(plaque_fleche(1, 1, h=h_plaque, nu=0))

