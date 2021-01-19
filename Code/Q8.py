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
from operator import itemgetter
# Paramètres de plaque
h_plaque = 0.1 # Epaisseur de plaque
E= 210e9
nu = 0.3
rho = 2700.
Nmodes = 4

def grid_frequence(Nx, Nmodes, L_plaque=1, H_plaque=1):
    """
    L_plaque: Longueur de plaque
    H_plaque: Largeur de plaque
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
    S = 0.5*b*h
    I = b*h**3/12
    J = b*h**3/12+b**3*h/12
    
    #c = 1/3*(1-0.63*h/b+0.052*(h/b)**5)
    #J = c*b*h**3
    sect = GridSection(S, I, J)

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
    
    K = assembl_stiffness_matrix(model)
    M = assembl_mass_matrix(model)
    L,Ud = assembl_connections(appuis,model)
    omega2, xi = eigenmodes(K,M,L,Nmodes)
    return omega2, xi, mesh

def plaque_frequence(L_plaque=1, H_plaque=1):
    '''Calculer la fréquence de la plaque, dont l\'équation est montrée dans le
    rapport. Enfin, on offre les premières modes 4 par défaut'''
    D = E*h_plaque**3/12/(1-nu**2)
    M = rho*h
    coef = pi**2*(D/M)**0.5
    modes = {}
    for m in range(1, 5):
        for n in range(1, 5):
            f = coef*((m/L_plaque)**2+(n/H_plaque)**2)/(2*pi)
            modes[(m, n)] = f
    return sorted(modes.items(), key=itemgetter(1), reverse=False)[:Nmodes]

if __name__ == "__main__":
    '''Comparer la fréquence de la grille et celle de la plaque avec deux cas,
    soit carré, soit rectangulaire'''
    for L_plaque, H_plaque in [(1, 1), (1, 5)]:
        print("########## Longueur={} et largeur={} ##########".format(L_plaque,
                                                                       H_plaque))
        omega2, xi, mesh = grid_frequence(8, Nmodes, L_plaque, H_plaque)
        for i in range(Nmodes):
            fig = Figure(4+i, "Eigenmodes: %d" % (i+1))
            fig.plot_field_lines(mesh, xi[::3,i], 1.)
            fig.show()
        print(omega2**0.5/2/pi)
        pf = plaque_frequence(L_plaque, H_plaque)
        print(np.array(pf)[:, 1])
        print("Erreur:")
        print(abs(np.array(pf)[:, 1]-omega2**0.5/2/pi)/np.array(pf)[:, 1])
