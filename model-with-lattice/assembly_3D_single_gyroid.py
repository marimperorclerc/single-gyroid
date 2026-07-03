r"""
single gyroid model for finite size crystal with (na,nb,nc) unit cells

Definition
----------

TO BE WRITTEN

Orientation average: The 1D form factor corresponds to the orientation average with all the possible orientations having the same probability.
Instead of rotating the shape through all the possible orientations,
it is equivalent to integrate the 3D scattering vector over a sphere of radius q with the shape in its reference orientation.

The sphere is sampled using Fibonacci quadrature to provide a quasi-uniform distribution of points on the unit sphere.
The distribution of the N points is computed using the golden ratio (see fibonacci.py). 
Each point of the quadrature on the unit sphere corresponds to a vector :math:`\mathbf{u}_{j}`.
In the sum, all weights :math:`w_j` are taken identical and equal to :math:`\frac{1}{N}`.

.. math::

    P(q) =  \sum_{j=1}^{N} w_j I(q\mathbf{u}_{j}, n, R, L)

.. figure:: img/fibonacci_sphere.png

    Fibonacci sphere using N=5810 points.

Validation
----------

TO DO

References
----------

TO DO

1. Marcone, J., Trazo, J. G., Nag, R., Goldmann, C., Ratel-Ramond, N., Hamon, C., & Impéror-Clerc, M. (2025).
   Form factor of prismatic particles for small-angle scattering analysis.
   *Journal of Applied Crystallography*, 58(2), 543‑552. https://doi.org/10.1107/S1600576725000676


Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr) 3d July 2026

* **Last Modified by:** MIC **Date:** 11 December 2025

* **Last Reviewed by:** SM **Date:** 17 April 2026

"""
import numpy as np
from numpy import inf
import scipy.special as sp

from sasmodels.special import sas_sinx_x
from sasmodels.special.fibonacci import fibonacci_sphere

name = "assembly_3D_single_gyroid"
title = "single gyroid model for finite size crystal with (na,nb,nc) unit cells"
description = """
        single gyroid model for finite size crystal with (na,nb,nc) unit cells"""

# assembly category: to be created
category = "shape-independent"

#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 50., [-inf, inf], "sld",
               "matrix scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 0., [-inf, inf], "sld",
               "solvent scattering length density"],
              ["t", "0.", 0.1, [0., 0.5], "nodal surface parameter",
               ""],
              ["a_lattice", "Ang", 3540., [0., inf], "volume",
               "a_lattice"],
              ["n_a", "", 10, [1, 50], "volume",
               "assembly number along a_lattice"],
              ["n_b", "", 10, [1, 50], "volume",
               "assembly number along b_lattice"],
              ["n_c", "", 10, [1, 50], "volume",
               "assembly number along c_lattice"]
             ]

def f(x,y,z,a):
    """ 
    nodal surface equation
    a: cubic cell parameter
    """ 
    term1=np.cos(2*np.pi*x/a)*np.sin(2*np.pi*y/a)
    term2=np.cos(2*np.pi*y/a)*np.sin(2*np.pi*z/a)
    term3=np.cos(2*np.pi*z/a)*np.sin(2*np.pi*x/a)
    return term1+term2+term3

def Funitcell(qa,qb,qc,a,t,nmax=50):
    """ 
    returns the complex form factor amplitude of the cubic unit cell
    a: cubic cell parameter
    t: nodal surface value
    nmax: number of steps for integration
    """ 

    delta = a/nmax
    x0 = y0 = z0 = -a/2
    sumcos = 0.
    sumsin = 0.
    for nx in range(nmax):
        for ny in range(nmax):
            for nz in range(nmax):
                x = x0 + nx*delta
                y = y0 + ny*delta
                z = z0 + nz*delta
                if f(x,y,z,a)>t:
                    sumcos += np.cos(qa*x+qb*y+qc*z)*delta*delta*delta
                    sumsin += np.sin(qa*x+qb*y+qc*z)*delta*delta*delta
    return (sumcos+1J*sumsin)/a**3

def Iunitcell(qa,qb,qc,a,t,nmax=50):
    """ 
    returns Iunitcell the intensity for one unit cell
    a: cubic cell parameter
    t: nodal surface value
    nmax: number of steps for integration
    """ 
    A=Funitcell(qa,qb,qc,a,t,nmax)
    return (np.abs(A))**2 


def vol_fraction(a,t,nmax=50):
    """ 
    a: cubic cell parameter
    t: nodal surface value
    nmax: number of steps for integration
    returns the matrix volume fraction
    """ 
    
    vf=abs(Funitcell(0,0,0,a,t,nmax))

    return vf

def unit_cell(a_lattice, b_lattice, c_lattice): #volume of assembly for orthomrhombic symmetry only
    return a_lattice*b_lattice*c_lattice


def form_volume(a_lattice,b_lattice,c_lattice,n_a,n_b,n_c):
    """
    returns the volume of the 3D assembly
    """
    assembly=n_a*n_b*n_c*unit_cell(a_lattice,b_lattice,c_lattice)

    return assembly


#Structure factor for 1D assembly
   
def laue_cheb(x,nsize):  #Laue function usign Chebichev polynomia
    val=sp.eval_chebyu(nsize-1, np.cos(x/2))**2    
    return(val/nsize) 
#Normalisation choice is to divide by the number of sites nsize. 
#Peak intensity at the maximum is proportionnal to nsize and not to nsize**2 or 1.


def ILaueqabc(qa,qb,qc,a_lattice,b_lattice,c_lattice,n_a,n_b,n_c): # orthorhombic only
    n_a = int(n_a)
    n_b = int(n_b)
    n_c = int(n_c)
    a_vect=[a_lattice,0,0]
    b_vect=[0,b_lattice,0]
    c_vect=[0,0,c_lattice]
    laue_a=laue_cheb(qa*a_vect[0]+qb*a_vect[1]+qc*a_vect[2],n_a)**2
    laue_b=laue_cheb(qa*b_vect[0]+qb*b_vect[1]+qc*b_vect[2],n_b)**2
    laue_c=laue_cheb(qa*c_vect[0]+qb*c_vect[1]+qc*c_vect[2],n_c)**2
    return laue_a*laue_b*laue_c

def Iqabc(qa,qb,qc,a_lattice,t,n_a,n_b,n_c,nmax:int=50): # proportionnal to the volume**2
    """
    Parameters
    ----------
    qa, qb, qc :  np.ndarray shape (N,)
        components of the scattering vector q
    a_lattice : float
        cubic unit cell parameter
    n_a : int
        Number of unit cells along a_vect (from cubic unit cell)
    n_b : int
        Number of unit cells along b_vect (from cubic unit cell)
    n_c : int
        Number of unit cells along c_vect (from cubic unit cell)
    Returns
    -------
    Iqabc : np.ndarray of float, shape (N,)
        Scattered intensity for single gyroid unit cell at the specific three dimensional q components.
        multiplied by Laue function for 3D assembly
    """

    n_a = int(n_a)
    n_b = int(n_b)
    n_c = int(n_c)

    coeff=n_a*n_b*n_c

    laue=ILaueqabc(qa,qb,qc,a_lattice,a_lattice,a_lattice,n_a,n_b,n_c) #because of cubic symmetry

    unit = Iunitcell(qa,qb,qc,a_lattice,t,nmax) 

    intensity = coeff*unit*laue  # intensity is proportional to volume
    return intensity



def Iq(q, sld, sld_solvent, a_lattice:float, t:float, n_a:int, n_b:int, n_c:int, nmax:int=50, npoints_fibonacci:int= 1000):
    """
    Computes the scattering intensity I(q) of nanoprisms averaged over all orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci. Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float ou array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    n_sides, radius_average, length :
        Geometrical parameters of the prism
    npoints_fibonacci : int
        Number of Fibonacci points on the sphere, set to 500 by default
        (higher number increases accuracy but also computation time, 500 is usually a good compromise)
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    """
    n_a = int(n_a)
    n_b = int(n_b)
    n_c = int(n_c)

    q = np.atleast_1d(q)
    q_unit,w = fibonacci_sphere(npoints_fibonacci)   # shape (npoints, 3)
    # Projections
    qa = q[:, np.newaxis] * q_unit[:, 0][np.newaxis, :]
    qb = q[:, np.newaxis] * q_unit[:, 1][np.newaxis, :]
    qc = q[:, np.newaxis] * q_unit[:, 2][np.newaxis, :]
    # # Compute intensity

    intensity = Iqabc(qa.ravel(), qb.ravel(), qc.ravel(),a_lattice, t, n_a, n_b, n_c,nmax).reshape(qa.shape)
    # Uniform average over the sphere
    integral = np.mean(intensity, axis=1)
    return (integral) * (sld - sld_solvent)**2 * 10**-4

Iq.vectorized = True

# TO ADAPT LATER:
tests = [
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     0.01, 5.62789],
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [5.62789, 0.73696]],
]
