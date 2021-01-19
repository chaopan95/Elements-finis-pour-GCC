# -*- coding: utf-8 -*-
"""
Module for `Figure` class

This file is part of the **WomBat** finite element code
used for the *Civil Engineering Finite Element Course*
of *Ecole des Ponts ParisTech* 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""
import numpy as np
from .element import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Circle
from matplotlib import cm


lin_shape_functions = lambda xi,eta: np.array([xi,eta,1-xi-eta])
quad_shape_functions = lambda xi,eta: np.array([(1-xi)*(1-eta)/4.,
                                                (1+xi)*(1-eta)/4.,
                                                (1+xi)*(1+eta)/4.,
                                                (1-xi)*(1+eta)/4.])


class Figure:
    """ Figure object for plotting mesh, results, etc. using
    matplotlib.pyplot """

    def __init__(self,n=1, title="",ax=[1,1]):
        """
        Parameters
        ----------
        n : int
            figure number
        title : string
            figure title
        ax : list
            generate p x q subplots with ax=[p,q]
        """
        self.fig = plt.figure(n)
        self.ax_size = ax
        self.ax0 = self.fig.add_subplot(ax[0],ax[1],1, aspect='equal')
        plt.title(title,fontsize=16)

    def show(self,block=True):
        """ Draw the figure

        Parameters
        ----------
        block : bool
            if block = True, execution is suspended until figure is closed,
            otherwise execution continues and figure remains open
        """
        plt.show(block=block)
    def clear(self):
        """ Clear the figure
        """
        self.fig.clf()

    def add_subplot(self,i,aspect='equal'):
        """ Add subplot and set it as current

        Parameters
        ----------
        i : int
            number of the current subplot
        """
        if i>0:
            self.fig.add_subplot(self.ax_size[0],self.ax_size[1],i+1, aspect=aspect,sharex=self.ax0)
        else:
            self.ax0.set_aspect(aspect)

    def plot(self,mesh,nodelabels=False,elemlabels=False,nodesymbols=True):
        """ Plot mesh and node/element informations

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        nodelabels : bool
            if True, plot node numbers on the mesh
        elemlabels : bool
            if True, plot element numbers on the mesh
        nodesymbols : True
            if True, plot nodes as small circles on the mesh
        """
        hsize = [e.hsize() for e in mesh.elem_list]
        mesh_bbox = np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0)
        node_symb_rad = min(0.2*min(hsize),0.015*max(mesh_bbox))

        ax = self.fig.gca()
        if all([isinstance(e,Triangle) for e in mesh.elem_list]):
            plt.triplot(mesh.coor[:,0], mesh.coor[:,1], mesh.connec[:,:3],"-k",linewidth=0.5)
        else:
            for el in mesh.elem_list:
                c = el.node_coor()
                if isinstance(el,Triangle6): # reorder nodes to be in consecutive orders
                    c = c[[0,3,1,4,2,5],:]
                P = Polygon(c,True,fill=True,ec="black",fc='lightgray',linewidth=1)
                ax.add_patch(P)

                cG = c.mean(axis=0)
                if elemlabels:
                    plt.text(cG[0],cG[1],str(el._id),
                    ha='center',va='center',color='r',fontsize=16,fontweight='bold',
                    backgroundcolor='w')

                for node in el.nodes.node_list:
                    if nodesymbols:
                        C = Circle((node.coor[0],node.coor[1]),node_symb_rad,ec="black",facecolor='w')
                        ax.add_patch(C)
                    if nodelabels:
                        plt.text(node.coor[0]+2*node_symb_rad,node.coor[1]+2*node_symb_rad,str(node._id),
                        color='b',fontsize=16,fontweight='bold')

        ax.axis('equal')
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)


    def plot_bc(self,mesh,imposed_displ):
        """ Plot boundary conditions, relations between master/slave dofs
        and imposed displacement to a non-zero value are not represented.

        - vertical green lines : fixed horizontal displacement
        - horizontal green lines : fixed vertical displacement
        - magental crosses : fixed rotation

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        imposed_displ
            imposed displacement :class:`Connections <connections.Connections>` object
        """
#        lgth_max = np.linalg.norm(np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0))
#        bc_symbol_ratio = 0.04
#        bc_symbol_size = bc_symbol_ratio*lgth_max
        hsize = [e.hsize() for e in mesh.elem_list]
        mesh_bbox = np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0)
        bc_symbol_size = max(min(0.5*min(hsize),0.04*max(mesh_bbox)),0.02*max(mesh_bbox))

        ax = self.fig.gca()
        lw = 2.
        for (j,master) in enumerate(imposed_displ.master_list):
            if master is None:
                for i,node in enumerate(imposed_displ.slave_list[j]):
                    comp = imposed_displ.components_list[j]
                    val = imposed_displ.imposed_value_list[j]
                    if 0 in comp and (val[comp.index(0)] == 0):
                        astyle = '-'
                        acolor = 'g'
                        beg_arrow1 = node.coor-bc_symbol_size*np.array([0.5,1])
                        end_arrow1 = node.coor+bc_symbol_size*np.array([-0.5,1])
                        beg_arrow2 = node.coor-bc_symbol_size*np.array([-0.5,1])
                        end_arrow2 = node.coor+bc_symbol_size*np.array([0.5,1])

                        ax.plot([beg_arrow1[0],end_arrow1[0]],[beg_arrow1[1],end_arrow1[1]],astyle+acolor,linewidth=lw)
                        ax.plot([beg_arrow2[0],end_arrow2[0]],[beg_arrow2[1],end_arrow2[1]],astyle+acolor,linewidth=lw)

                    if 1 in comp and (val[comp.index(1)] == 0):
                        astyle = '-'
                        acolor = 'g'
                        beg_arrow1 = node.coor-bc_symbol_size*np.array([1,0.5])
                        end_arrow1 = node.coor+bc_symbol_size*np.array([1,-0.5])
                        beg_arrow2 = node.coor-bc_symbol_size*np.array([1,-0.5])
                        end_arrow2 = node.coor+bc_symbol_size*np.array([1,0.5])

                        ax.plot([beg_arrow1[0],end_arrow1[0]],[beg_arrow1[1],end_arrow1[1]],astyle+acolor,linewidth=lw)
                        ax.plot([beg_arrow2[0],end_arrow2[0]],[beg_arrow2[1],end_arrow2[1]],astyle+acolor,linewidth=lw)


                    if 2 in comp and (val[comp.index(2)] == 0):
                        astyle = '-'
                        acolor = 'm'
                        beg_arrow1 = node.coor-bc_symbol_size*np.array([0.5,0.5])
                        end_arrow1 = node.coor+bc_symbol_size*np.array([0.5,0.5])
                        beg_arrow2 = node.coor-bc_symbol_size*np.array([0.5,-0.5])
                        end_arrow2 = node.coor+bc_symbol_size*np.array([0.5,-0.5])

                        ax.plot([beg_arrow1[0],end_arrow1[0]],[beg_arrow1[1],end_arrow1[1]],astyle+acolor,linewidth=lw)
                        ax.plot([beg_arrow2[0],end_arrow2[0]],[beg_arrow2[1],end_arrow2[1]],astyle+acolor,linewidth=lw)

#        ax.axis('equal')
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)


    def plot_def(self,mesh,U,ampl=1.,undef=True, color = "b"):
        """ Plot deformed mesh shape

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        U : ndarray
            displacement solution used to deform the mesh shape
        ampl : float
            amplification factor
        undef : bool
            if True, plots also undeformed mesh
        """
        ax = self.fig.gca()
        if all([isinstance(e,Triangle) and not(isinstance(e,Triangle6)) for e in mesh.elem_list]):
            if undef:
                plt.triplot(mesh.coor[:,0], mesh.coor[:,1], mesh.connec[:,:3],"-k",linewidth=0.5,alpha=0.5)
#            plt.tripcolor(mesh.coor[:,0]+ampl*U[::2], mesh.coor[:,1]+ampl*U[1::2], mesh.connec[:,:3],facecolors=np.ones((mesh.Nel,)),alpha=0.8,cmap=cm.Blues)
            plt.triplot(mesh.coor[:,0]+ampl*U[::2], mesh.coor[:,1]+ampl*U[1::2], mesh.connec[:,:3],"-" + color,linewidth=0.5, alpha=0.5)
        else:
            el0 = mesh.elem_list[0]
            dofe = el0.get_dof()
            Ue = U[dofe]
            x_def,y_def = el0.deformation(ampl*Ue)
            m = x_def.shape[0]
            X_def = np.zeros((m,mesh.Nel))
            Y_def = np.zeros((m,mesh.Nel))
            X = np.zeros((m,mesh.Nel))
            Y = np.zeros((m,mesh.Nel))
            for (i,el) in enumerate(mesh.elem_list):
                dofe = el.get_dof()
                Ue = U[dofe]
                x_def,y_def= el.deformation(ampl*Ue)
                X_def[:,i] = x_def
                Y_def[:,i] = y_def
                if undef:
                    x,y= el.deformation(0*Ue)
                    X[:,i] = x
                    Y[:,i] = y
            if isinstance(el,Segment):
                thick = 2
            else:
                thick = 0.5
            ax.plot(X_def,Y_def,'-'+color,linewidth=thick,alpha=0.5)
            if undef:
                ax.plot(X,Y,'-k',linewidth=0.5,alpha=0.5)
#        ax.axis('equal')
        ax.relim()
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)

    def plot_field(self,mesh,V,extrapolation=True,U=[]):
        """ Plot field defined on a triangular mesh (for :class:`Bar2D <element.bar2D.Bar2D>`  or :class:`Beam2D <element.beam2D.Beam2D>`
        use :func:`plot_field_diagrams` instead)

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        V : ndarray
            field to plot, can be of shape :

            - (Nno,) : nodal field

                + for T3 triangles linear interpolation between nodes
                + for T6 triangles, each quadratic triangle is splitted \
                    into 4 subtriangles overwhich linear interpolation is performed
            - (Nel,) : constant field over elements (discontinuous plot)
            - (3*Nel,) : element field with linear variation, values are expressed \
                at Gauss points and extrapolated at nodes without smoothing (discontinuous plot)
        extrapolation : bool
            if True (default) extrapolation is performed, otherwise values are plotted as is at vertices
        """
        ax = self.fig.gca()

        e0 = mesh.elem_list[0]
        try:
            ngauss = e0.ngauss
        except:
            ngauss = 1
        if V.shape[0] == mesh.Nel:      # case of constant field/elem
            if len(U) > 0:
                Ux = U[::2]
                Uy = U[1::2]
            else:
                Ux = 0
                Uy = 0
            plt.tripcolor(mesh.coor[:,0] + Ux, mesh.coor[:,1] + Uy, triangles=mesh.connec[:,:3], facecolors=V, shading='flat')
        elif mesh.connec.shape[1]==6 and V.shape[0] == mesh.Nno:    # case of T6 element
            print("Warning: Quadratic triangles are approximated by 4 piecewise-linear triangles")
            tri = np.zeros((4*mesh.Nel,3))
            for (j,e) in enumerate(mesh.elem_list):
                nodes = e.nodes.get_id()
                tri[4*j:4*(j+1),:] = np.array([nodes[[0,3,5]],
                                               nodes[[3,1,4]],
                                               nodes[[4,2,5]],
                                               nodes[[3,4,5]]])
            if len(U) > 0:
                Ux = U[::2]
                Uy = U[1::2]
            else:
                Ux = 0
                Uy = 0
            plt.tripcolor(mesh.coor[:,0]+Ux, mesh.coor[:,1]+Uy, tri, V, shading='gouraud')
        elif mesh.connec.shape[1]==4 and V.shape[0] == ngauss*mesh.Nel:  # case of quadrangles
            print("Warning: Values are extrapolated at nodes from Gauss points")
            V_node = np.zeros((6*mesh.Nel,))    # value at vertices of the element
            tri = np.reshape(np.arange(0,6*mesh.Nel),(2*mesh.Nel,3))
            coor = np.zeros((6*mesh.Nel,2))
            for (j,e) in enumerate(mesh.elem_list):
                nodes = e.nodes.get_id()
                coor[6*j:6*j+3,:] = e.node_coor()[[0,1,3],:]
                coor[6*j+3:6*j+6,:] = e.node_coor()[[1,2,3],:]
                if extrapolation:
                    N = np.array([quad_shape_functions(e.ag[i,0],e.ag[i,1]) for i in range(ngauss)])
                    V_quad = np.linalg.lstsq(N,V[ngauss*j:ngauss*(j+1)], rcond=None)[0]
                else:
                    V_quad = V[ngauss*j:ngauss*(j+1)]
                V_node[6*j:6*j+3] = V_quad[[0,1,3]]
                V_node[6*j+3:6*j+6] = V_quad[[1,2,3]]
            plt.tripcolor(coor[:,0], coor[:,1], tri, V_node, shading='gouraud')
        elif V.shape[0] == ngauss*mesh.Nel and ngauss > 1:  # case of non uniform field/elem
            print("Warning: Values are extrapolated at nodes from Gauss points")
            V_node = np.zeros((3*mesh.Nel,))    # value at vertices of the element
            tri = np.reshape(np.arange(0,3*mesh.Nel),(mesh.Nel,3))
            coor = np.zeros((3*mesh.Nel,2))
            Udef = np.zeros((3*mesh.Nel,2))
            for (j,e) in enumerate(mesh.elem_list):
                coor[3*j:3*(j+1),:] = e.node_coor()[:3,:]
                if extrapolation:
                    N = np.array([lin_shape_functions(e.ag[i,0],e.ag[i,1]) for i in range(ngauss)])
                    V_node[3*j:3*(j+1)] = np.linalg.lstsq(N,V[ngauss*j:ngauss*(j+1)], rcond=None)[0]
                else:
                   V_node[3*j:3*(j+1)] = V[ngauss*j:ngauss*(j+1)]
                if len(U) > 0:
                    nodes = e.nodes.get_id()
                    Udef[3*j:3*(j+1),0] = U[2*nodes[:3]]
                    Udef[3*j:3*(j+1),1] = U[2*nodes[:3]+1]
            plt.tripcolor(coor[:,0] + Udef[:,0], coor[:,1] + Udef[:,1], tri, V_node, shading='gouraud')
        else: # case of nodal field for T3 or constant field/element
            if len(U) > 0:
                Ux = U[::2]
                Uy = U[1::2]
            else:
                Ux = 0
                Uy = 0
            plt.tripcolor(mesh.coor[:,0] + Ux, mesh.coor[:,1] + Uy, mesh.connec, V, shading='gouraud')

        cbar = plt.colorbar(format="%.2e")
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)
        return cbar

    def plot_field_diagrams(self,mesh,V,scale=1.,unsigned=False):
        """ Plot field defined on a structural elements (only for :class:`Bar2D <element.bar2D.Bar2D>` or :class:`Beam2D <element.beam2D.Beam2D>`)
        using piecewise linear diagrams.

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        V : ndarray
            field to plot, can be of shape :

            - (Nno,) : nodal field, linear interpolation is used between nodes
            - (Nel,) : constant field over elements (discontinuous plot)
            - (2*Nel,) : element field with linear variation (discontinuous plot)
        scale : float
            use to scale the size of the diagrams
        """
        Vm = max(abs(V))
        hsize = [e.hsize() for e in mesh.elem_list]
        mesh_bbox = np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0)
        ratio = 0.2*max(min(hsize),max(mesh_bbox))

        ax = self.fig.gca()
        ax.set_aspect('equal')
        X = np.zeros((5*mesh.Nel,))
        Y = np.zeros((5*mesh.Nel,))
        connec = np.zeros((2*mesh.Nel,3))
        C = np.zeros((5*mesh.Nel,))

        for index,el in enumerate(mesh.elem_list):
            node = el.nodes.node_list
            if len(V) == mesh.Nel:
                Ve = np.repeat(V[index],2)
            elif len(V) == mesh.Nno:
                Ve = V[[node[0]._id,node[1]._id]]
            elif len(V) == 2*mesh.Nel:
                Ve = V[2*index:2*index+2]
            unit_scale = ratio*Ve/Vm
            c = el.node_coor()
            if unit_scale[0]*unit_scale[1] < 0:
                s = -np.mean(unit_scale)/np.diff(unit_scale)
            else:
                s = 0
            x = np.array([-0.5,0.5,-0.5,0.5,s])
            if unsigned:
                y = np.array([-unit_scale[0],-unit_scale[1],unit_scale[0],unit_scale[1],0])*scale/2.
            else:
                y = np.array([0,0,scale*unit_scale[0],scale*unit_scale[1],0])
            tang = c[1,:]-c[0,:]
            t = tang
            X[5*index+np.arange(5)] = x*t[0]-y*t[1]/np.linalg.norm(tang)+np.mean(c[:,0])
            Y[5*index+np.arange(5)] = x*t[1]+y*t[0]/np.linalg.norm(tang)+np.mean(c[:,1])
            C[5*index+np.arange(4)] = np.tile(Ve,2)
            if unit_scale[0]*unit_scale[1] >= 0:    # case of same sign on the element
                connec[2*index,:] = np.array([0,1,2])+5*index
                connec[2*index+1,:] = np.array([2,3,1])+5*index
            else:
                connec[2*index,:] = np.array([0,4,2])+5*index
                connec[2*index+1,:] = np.array([4,3,1])+5*index
#            plt.plot(c[:,0],c[:,1],'-k',linewidth=1)

        plt.tripcolor(X, Y, connec, C, shading='gouraud')
        ax.axis('equal')
        plt.colorbar()
        #ax.relim()
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)


    def plot_field_lines(self,mesh,V,scale=1.):
        """ Plot field defined on a structural elements (only for :class:`Bar2D <element.bar2D.Bar2D>` or :class:`Beam2D <element.beam2D.Beam2D>`)
        using piecewise linear diagrams.

        Parameters
        ----------
        mesh
            :class:`Mesh <mesh.Mesh>` object
        V : ndarray
            field to plot, can be of shape :

            - (Nno,) : nodal field, linear interpolation is used between nodes
            - (Nel,) : constant field over elements (discontinuous plot)
            - (2*Nel,) : element field with linear variation (discontinuous plot)
        scale : float
            use to scale the size of the diagrams
        """
        Vm = max(abs(V))
        hsize = [e.hsize() for e in mesh.elem_list]
        mesh_bbox = np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0)
        ratio = 0.2*max(min(hsize),max(mesh_bbox))

        ax = self.fig.gca()
        ax.set_aspect('equal')
        X = np.zeros((5*mesh.Nel,))
        Y = np.zeros((5*mesh.Nel,))
        connec = np.zeros((2*mesh.Nel,3))
        C = np.zeros((5*mesh.Nel,))

        for index,el in enumerate(mesh.elem_list):
            node = el.nodes.node_list
            if len(V) == mesh.Nno:
                Ve = V[[node[0]._id,node[1]._id]]
            elif len(V) == mesh.Nel:
                Ve = np.repeat(V[index],2)
            elif len(V) == 2*mesh.Nel:
                Ve = V[2*index:2*index+2]
            unit_scale = ratio/2.
            c = el.node_coor()
            x = np.array([-0.5,0.5,-0.5,0.5,0])
            y = np.array([-1,-1,1,1,0])*scale*unit_scale
            tang = c[1,:]-c[0,:]
            t = tang
            X[5*index+np.arange(5)] = x*t[0]-y*t[1]/np.linalg.norm(tang)+np.mean(c[:,0])
            Y[5*index+np.arange(5)] = x*t[1]+y*t[0]/np.linalg.norm(tang)+np.mean(c[:,1])
            C[5*index+np.arange(4)] = np.tile(Ve,2)
            connec[2*index,:] = np.array([0,1,2])+5*index
            connec[2*index+1,:] = np.array([2,3,1])+5*index
            plt.plot(c[:,0],c[:,1],'-k',linewidth=0.5)

        plt.tripcolor(X, Y, connec, C, shading='gouraud')
        ax.axis('equal')
        plt.colorbar()
        #ax.relim()
        ax.autoscale_view(tight=True)
        ax = axes_rescale(ax)

def axes_rescale(ax,padx=0.1,pady=0.1):
    """ Rescale axis with pad values"""
    Xmin,Xmax = ax.get_xlim()
    Ymin,Ymax = ax.get_ylim()
    DX = Xmax-Xmin
    DY = Ymax-Ymin
    Xc = (Xmax+Xmin)/2
    Yc = (Ymax+Ymin)/2

    ax.set_xlim((Xc-(1+padx)*DX/2,Xc+(1+padx)*DX/2))
    ax.set_ylim((Yc-(1+pady)*DY/2,Yc+(1+pady)*DY/2))

    return ax


#def plot_loading(self,mesh,forces,ampl,value=False):
#    lgth_max = np.linalg.norm(np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0))
#    ax = self.fig.gca()
#    for i,node in enumerate(forces.node_list):
#        Fx = forces.Fx_n[i]
#        Fy = forces.Fy_n[i]
#        if Fx is None:
#            Fx = 0
#        if Fy is None:
#            Fy = 0
#        norm = (Fx**2+Fy**2)**0.5
#        end_arrow = node.coor - np.array([Fx,Fy])*ampl*0.1*lgth_max
#        label = ''
#        if norm>0 and value:
#            label = str(norm)+' kN'
#        ax.annotate(label,xy=(node.coor[0],node.coor[1]), xycoords='data',
#        xytext=(end_arrow[0],end_arrow[1]), textcoords='data', ha = 'center',va = 'center',
#        fontsize = 14, color='r',
#        arrowprops=dict(arrowstyle="->",
#                        connectionstyle="arc3",linewidth=2,color='k'),
#        )
#
#    ax.axis('equal')
#    ax.autoscale_view(tight=True)
#    ax = axes_rescale(ax)
#
#def plot_reactions(self,mesh,imposed_displ,R,ampl):
#    lgth_max = np.linalg.norm(np.max(mesh.coor,axis=0)-np.min(mesh.coor,axis=0))
#    ax = self.fig.gca()
#    for index,node in enumerate(imposed_displ.node_list):
#        Udx = imposed_displ.get_ux_value(node)
#        Udy = imposed_displ.get_uy_value(node)
#        Rx = 0
#        Ry = 0
#        if Udx is not None:
#            Rx = R[2*index]
#        if Udy is not None:
#            Ry = R[2*index+1]
#        norm = (Rx**2+Ry**2)**0.5
#        end_arrow = node.coor - np.array([Rx,Ry])*ampl*0.1*lgth_max
#        label = ''
#        if norm>1e-10:
#            label = '{:.1} kN'.format(norm)
#        print Rx,Ry
#        ax.annotate(label,xy=(node.coor[0],node.coor[1]), xycoords='data',
#        xytext=(end_arrow[0],end_arrow[1]), textcoords='data', ha = 'center',va = 'center',
#        fontsize = 14, color='g',
#        arrowprops=dict(arrowstyle="->",
#                        connectionstyle="arc3",linewidth=2,color='g'),
#        )
#
#    ax.axis('equal')
#    ax.relim()
#    ax.autoscale_view(tight=True)
#    ax = axes_rescale(ax)
#