import numpy as np

def f(x,y,z,a):
    """ 
    nodal surface equation
    a: cubic cell parameter
    """ 
    term1=np.cos(2*np.pi*x/a)*np.sin(2*np.pi*y/a)
    term2=np.cos(2*np.pi*y/a)*np.sin(2*np.pi*z/a)
    term3=np.cos(2*np.pi*z/a)*np.sin(2*np.pi*x/a)
    return term1+term2+term3

def volume_fraction(a,t,nmax=50):
    """ 
    a: cubic cell parameter
    t: nodal surface value
    nmax: number of steps for integration
    """ 
    delta = a/nmax
    x0 = y0 = z0 = -a/2
    sum=0
    for nx in range(nmax):
        for ny in range(nmax):
            for nz in range(nmax):
                x = x0 + nx*delta
                y = y0 + ny*delta
                z = z0 + nz*delta
                if f(x,y,z,a)>t:
                    sum += delta*delta*delta
    return sum/(a*a*a)

def Fhkl(a,t,h,k,l,nmax=50):
    """ 
    returns Fhklcos,Fhklsin the real and imaginary parts of Fhkl
    a: cubic cell parameter
    t: nodal surface value
    hkl: Miller indices
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
                    sumcos += np.cos(2*np.pi*(h*x+k*y+l*z)/a)*delta*delta*delta
                    sumsin += np.sin(2*np.pi*(h*x+k*y+l*z)/a)*delta*delta*delta
    return sumcos/(a*a*a), sumsin/(a*a*a)

def Ihkl(a,t,h,k,l,nmax=50):
    """ 
    returns Ihkl for a single crystal
    a: cubic cell parameter
    t: nodal surface value
    hkl: Miller indices
    nmax: number of steps for integration
    """ 
    return (Fhkl(a,t,h,k,l,nmax)[0]**2+Fhkl(a,t,h,k,l,nmax)[1]**2)


def vol_fraction(a,t,nmax=50):
    """ 
    alternative volume fraction calculation for h=k=l=0
    a: cubic cell parameter
    t: nodal surface value
    nmax: number of steps for integration
    """ 
    return Fhkl(a,t,0,0,0,nmax)[0]