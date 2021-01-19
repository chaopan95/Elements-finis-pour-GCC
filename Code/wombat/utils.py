# -*- coding: utf-8 -*-
"""
Contains utility functions common to different submodules

This file is part of the WomBat finite element code
used for the Civil Engineering Finite Element Course 
of Ecole des Ponts ParisTech 2017-2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech, 
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.fr
"""

def append_component(comp,node_list_length,default_value=0):
    """ Appends the given component/list of components to the list of imposed displacements
    
    Parameters
    ----------
    comp : float,list
        list of components or single value
    node_list_length : int
        length of the list of corresponding nodes
    default_value 
        default value of comp when empty list is given
    """
    if not isinstance(comp, list): 
        comp = [comp]
    if len(comp)>0:
        if len(comp)==1:
            comp = comp*node_list_length
        if len(comp)!=node_list_length:
            raise ValueError("Given components list must be the same length as the given node list")
    else:
        comp = [default_value]*node_list_length
    return comp

def uniquify(seq):
    """ Remove duplicate elements from a list while preserving order"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]

def value_at_next_step(X_list,ampl):
    """ Compute value of a given step
    
    Parameters
    ----------
    X_list : ndarray or list of ndarray
        if type is ndarray, the full vector is time-dependent and its value is given by ampl*X_list
        if type is list, only the first element is time-dependent, the second one is fixed in time,
        the value is then ampl*X_list[0]+X_list[1]
    ampl : float
        amplitude of the loading at the given step
    """
        
    if not isinstance(X_list, list):
        X = ampl*X_list
    else:
        X = ampl*X_list[0]
        if len(X_list)>1:
            X += X_list[1]
    return X