r"""
This model adds a q**-4 Porod regime and Bragg peaks for single gyroid structure.
It correspond to the cubic space group I 4_1 3_2.
This structure is encountered in many butterfly wings made of chitin networks with structural colors. 


Definition
----------

Bragg peaks are modeled using pseudo-Voigt peak function.

The peak positions are computed for a cubic lattice with the I 4_1 3_2 extinctions rules.
q_hkl =l (2\pi/a_lattice) * sqrt(h**2 + k**2 + l**2)
with
h**2+k**2+l**2 = 2, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24 ...
i.e. the (110), (211), (220), (310), (222), (321), (400), (411)+(330), (420), (332), (422) ...

(411) and (330) have the same q_hkl modulus and are treated together as a single peak for 1D data.

Validation
----------


References
----------

1. G Porod. *Kolloid Zeit* 124 (1951) 83

2. L A Feigin, D I Svergun, G W Taylor 
   Structure Analysis by Small-Angle X-ray and Neutron Scattering 
   Springer (1987)

3. Aaron L. Stancik, Eric B. Brauns
   A simple asymmetric lineshape for fitting infrared absorption spectra
   Vibrational Spectroscopy 47 (2008) 66-69

   
Authorship and Verification
----------------------------

* **Authors:** Marianne Imperor-Clerc (marianne.imperor@cnrs.fr) 30 June 2026

* **Last Modified by:** MIC **Date:** 30 June 2026

* **Last Reviewed by:** **Date:**

"""

import numpy as np
from numpy import inf, errstate

name = "porod_peak_voigt_single_gyroid"
title = "Porod function with added single gyroid Bragg peaks"
description = """\
          I(q) = scale(scale_Porod/q^4+sum_of_single_gyroid_peaks) + background
"""

category = "shape-independent"

parameters = [["scale_Porod", "", 0.05, [0, inf], "", "Scale factor for Porod"],
              ["a_cell", "Ang", 900, [20, 10000], "", "cubic cell parameter"],
              ["w_f", "", 0.8, [0, 1], "", "lorentzian/gaussian weighting factor"],
              ["hwhm_q110", "1/Ang", 0.001, [0, 1], "", "HWHM of 110 peak"],
              ["hwhm_q211", "1/Ang", 0.001, [0, 1], "", "HWHM of 211 peak"],
              ["hwhm_q220", "1/Ang", 0.001, [0, 1], "", "HWHM of 220 peak"],
              ["hwhm_q310", "1/Ang", 0.001, [0, 1], "", "HWHM of 310 peak"],
              ["hwhm_q222", "1/Ang", 0.001, [0, 1], "", "HWHM of 222 peak"],
              ["hwhm_q321", "1/Ang", 0.001, [0, 1], "", "HWHM of 321 peak"],
              ["hwhm_q400", "1/Ang", 0.001, [0, 1], "", "HWHM of 400 peak"],
              ["hwhm_q411", "1/Ang", 0.001, [0, 1], "", "HWHM of 411+330 peak"],
              ["hwhm_q420", "1/Ang", 0.001, [0, 1], "", "HWHM of 420 peak"],
              ["hwhm_q332", "1/Ang", 0.001, [0, 1], "", "HWHM of 332 peak"],
              ["hwhm_q422", "1/Ang", 0.001, [0, 1], "", "HWHM of 422 peak"],
              ["scale_q110", "", 1, [0,inf], "", "Scale for 110 peak"],
              ["scale_q211", "", 1, [0,inf], "", "Scale for 211 peak"],
              ["scale_q220", "", 1, [0,inf], "", "Scale for 220 peak"],
              ["scale_q310", "", 1, [0,inf], "", "Scale for 310 peak"],
              ["scale_q222", "", 1, [0,inf], "", "Scale for 222 peak"],
              ["scale_q321", "", 1, [0,inf], "", "Scale for 321 peak"],
              ["scale_q400", "", 1, [0,inf], "", "Scale for 400 peak"],
              ["scale_q411", "", 1, [0,inf], "", "Scale for 411+330 peak"],
              ["scale_q420", "", 1, [0,inf], "", "Scale for 420 peak"],
              ["scale_q332", "", 1, [0,inf], "", "Scale for 332 peak"],
              ["scale_q422", "", 1, [0,inf], "", "Scale for 422 peak"]]



def Ipeak(q, wf, q0, hwhm):
    """
    When $w_f$ = 1 a Lorentzian peak is returned, and when $w_f$ = 0 a
    Gaussian peak is returned.

    The peak is taken to be centered at $q_0$ with a HWHM (half-width
    half-maximum) for the Lorentzian and $sigma$=HWHM/1.17741, where $sigma$ is the standard deviation
    of the Gaussian. In other words, the widths of the Lorentzian and the
    Gaussian have been coupled for convenience of parameterisation.
    
    1.17741=np.sqrt(2*np.ln(2))
    """
    #cste=np.sqrt(2*np.ln(2))
    cste=1.17741
    sigma=hwhm/cste
    intensity = (wf*(1/(1+((q-q0)**2.0/hwhm**2.0))))+((1.0-wf)*np.exp((-0.5*(q-q0)**2.0)/(sigma**2.0)))
    return intensity


def Iq(q,scale_Porod,a_cell, w_f, hwhm_q110, hwhm_q211, hwhm_q220, hwhm_q310, hwhm_q222, hwhm_q321, hwhm_q400, hwhm_q411, hwhm_q420, hwhm_q332, hwhm_q422, scale_q110, scale_q211, scale_q220, scale_q310, scale_q222, scale_q321, scale_q400, scale_q411, scale_q420, scale_q332, scale_q422):
    """
    scale_Porod: scale coefficent for Porod term
    a_cell: cubic unit cell parameter
    w_f: weighting coefficient in the pseudo Voigt peak function 
    $w_f=1$ for Lorentzian and $w_f=0$ for Gaussian peak.
    hwhm_q110: HWHM of 110 peak
    scale_q110: Scale factor for 110 peak
    hwhm_q211: HWHM of 211 peak
    scale_q211: Scale factor for 211 peak
    hwhm_q220: HWHM of 220 peak
    scale_q220: Scale factor for 220 peak
    and so on ... until 422 peak
    411 and 330 peaks have the same scattering modulus and are treated as a single peak (1D data)
    """
    
    with errstate(divide='ignore'):
        q110=np.sqrt(2)*2*np.pi/a_cell
        q211=np.sqrt(6)*2*np.pi/a_cell
        q220=np.sqrt(8)*2*np.pi/a_cell
        q310=np.sqrt(10)*2*np.pi/a_cell
        q222=np.sqrt(12)*2*np.pi/a_cell
        q321=np.sqrt(14)*2*np.pi/a_cell
        q400=np.sqrt(16)*2*np.pi/a_cell
        q411=np.sqrt(18)*2*np.pi/a_cell #includes 330 peak
        q420=np.sqrt(20)*2*np.pi/a_cell
        q332=np.sqrt(22)*2*np.pi/a_cell
        q422=np.sqrt(24)*2*np.pi/a_cell
        porod = scale_Porod/q**4
        L2 = Ipeak(q,w_f,q110,hwhm_q110)
        L6 = Ipeak(q,w_f,q211,hwhm_q211)
        L8 = Ipeak(q,w_f,q220,hwhm_q220)
        L10 = Ipeak(q,w_f,q310,hwhm_q310)
        L12 = Ipeak(q,w_f,q222,hwhm_q222)
        L14 = Ipeak(q,w_f,q321,hwhm_q321)
        L16 = Ipeak(q,w_f,q400,hwhm_q400) 
        L18 = Ipeak(q,w_f,q411,hwhm_q411) #includes 330 peak
        L20 = Ipeak(q,w_f,q420,hwhm_q420)
        L22 = Ipeak(q,w_f,q332,hwhm_q332)
        L24 = Ipeak(q,w_f,q422,hwhm_q422)

        return porod+scale_q110*L2+scale_q211*L6+scale_q220*L8+scale_q310*L10+scale_q222*L12+scale_q321*L14+scale_q400*L16+scale_q411*L18+scale_q420*L20+scale_q332*L22+scale_q422*L24

Iq.vectorized = True  # Iq accepts an array of q values

tests = [
    [{"scale": 0.00001, 
      "background":0.01,
      "scale_Porod": 1.,
      "scale_q110": 0.,
      "scale_q211": 0.,
      "scale_q220": 0.,
      "scale_q310": 0.,
      "scale_q222": 0.,
      "scale_q321": 0.,
      "scale_q400": 0.,
      "scale_q411": 0.,
      "scale_q420": 0.,
      "scale_q332": 0.,
      "scale_q422": 0.,
      }, 
      0.04, 3.916250],
]
