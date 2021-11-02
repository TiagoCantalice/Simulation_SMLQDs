import numpy as np
import math

##VARIABLES MEANING
# r0 : radius of the spherical island
# mu : surface chemical potential for InAs
# mu_0 : surface chemical potential for unstressed InAs surface
# sigma_t : tangential stress contribution of surface curvature 
# r: distance from the center of the island
# B_InAs : bulk modulus of InAs
# E_GaAs : Young's Modulus of GaAs
# E_InAs : Young's Modulus of InAs
# nu_GaAs : Poisson's ratio of GaAs 
# d : size of the island's base
# e0 : natural lattice mismatch between GaAs and InAs
# e_xxGaAs : strain component on the surface (z = z_s), with additional correction from an image island due to the existence of a free surface 
# l: distance between center of two InAs islands
# D: surface diffusion coefficient of In adatoms
# T: growth temperature
# F: In flux
# W: strain strength
# Omega_InAs: atomic volume of InAs 
# Omega_GaAs: atomic volume of GaAs
# B_InAs : Bulk modulus of InAs
# C_11InAs : Elastic Constant of InAs
# Ld: diffusion length
# kb: Boltzmann constant
# ls: distance between center of InAs islands and inflation of potential
# tal: average lifetime for incorporation of In adatoms into lattice in the region outside the islands
# P: island pairing probability
# M1 : Material arriving in the region -d/2 < x < d/2
# K : the net probability  of In adatoms being sucked into the boundaries -d/2 < x < d/2
# T : growth temperature
# M_{t}*F*l : total material delivered in unit time within length l
# D: is the surface diffusion coefficient of the In adatoms
# D/(kb*T): is the mobility of the In adatmos
# tal : time for incorporation outside the islands
# n(+-d/2) : adatom concentration at the boundaries x = +- d/2
# z0 : characteristic spacer layer thickness below which a vertically self-organized growth occurs

def Q(l, d, Ld):
    return (l-d)/(2*Ld)

def Vd(tal):
    return Ld(D, tal)/tal

def n(F):
    numerator = F*Ld(D, tal)*(math.sinh(Q(l, d, Ld)))
    denominator = math.cosh(Q(l, d, Ld)) + math.cosh(Q(l, d, Ld))
    return numerator/denominator    

def e_xxGaAs():
    numerator = 2*A()*e0*r0**3
    denominator = (x**2+zs**2)**(3/2)
    return numerator/denominator

def A(B_InAs,E_GaAs, nu_GaAs):
    numerator = 3*B_InAs
    denominator = 3*B_InAs + 2*E_GaAs/(1+nu_GaAs)
    return numerator/denominator

def W(Omega_InAs, C_11InAs, e0, E_InAs):
    numerator   = Omega_InAs*(C_11InAs*e0)**2
    denominator = 2*E_InAs
    return numerator/denominator  

def Ld(D, tal):
    return np.sqrt(D*tal)

def adatom_concentration(F, Ld, Q, K, Vd):
    return 0

def z0(r0, Ld, W, l, kb, T, A):
    numerator = r0*(   8*Ld*W * A  )**(1/3)
    denominator = (l*kb*T)**(1/3)
    return numerator/denominator

def Omega_InAs(lat_par_InAs):
    #volume = (4/3)*math.pi*(Angstrom_to_cm(lat_par_InAs)/4)**3
    volume = (Angstrom_to_cm(lat_par_InAs)/2)**3
    return volume

def e0(lat_par_InAs, lat_par_GaAs):
    return (lat_par_InAs-lat_par_GaAs)/lat_par_GaAs

# UNIT CONVERTERS

def GPa_to_dyn_over_cm2(x):
    return x*10**10

def dyn_times_cm_to_eV(x):
    return x*624150636309.4

def Angstrom_to_ML_GaAs(x):
    return x*2/5.6533

def Angstrom_to_ML_InAs(x):
    return x*2/6.0583

def Angstrom_to_cm(x):
    return x*10**-8

if __name__ == '__main__':
    #parameters for A
    B_InAs  = 5.8*10**11 #dyn cm^-2 || reference: http://www.ioffe.ru/SVA/NSM/Semicond/InAs/mechanic.html
    nu_GaAs = 0.31  #no unit        || reference: https://www.memsnet.org/material/galliumarsenidegaasbulk/
    E_GaAs  = 85.5  #GPa            || reference: https://www.memsnet.org/material/galliumarsenidegaasbulk/

    print('A (no unit) :', round(A(B_InAs,GPa_to_dyn_over_cm2(E_GaAs), nu_GaAs), 2))

    #parameters for W
    lat_par_InAs = 6.0583 #A
    lat_par_GaAs = 5.6533 #A
    Omega_InAs   = Omega_InAs(lat_par_InAs) #cm^3
    C_11InAs     = 8.03*10**11  # dyn*cm^2 || reference: http://www.ioffe.ru/SVA/NSM/Semicond/InAs/mechanic.html
    e0           = e0(lat_par_InAs, lat_par_GaAs) #no unit
    E_InAs       = 5.14*10**11 # dyn*cm^2 || reference: http://www.ioffe.ru/SVA/NSM/Semicond/InAs/mechanic.html

    print('W (strain strength (eV)) :', round(dyn_times_cm_to_eV(W(Omega_InAs, C_11InAs, e0, E_InAs)), 5))

    #parameters for z0
    W    = round(dyn_times_cm_to_eV(W(Omega_InAs, C_11InAs, e0, E_InAs)), 5) #eV
    Ld   = 0.28          #mu meters
    r0   = 37            #A
    kb   = 8.62*10**-5   #eV/K 
    l    = 0.0535        #mu meters
    T    = 773           #K (500 C)
    A    = round(A(B_InAs,GPa_to_dyn_over_cm2(E_GaAs), nu_GaAs), 2)          #

    print('z0 (characteristic length (A)) :', round(z0(r0, Ld, W, l, kb, T, A), 2))
    print('z0 (characteristic length (ML):', round(Angstrom_to_ML_GaAs(z0(r0, Ld, W, l, kb, T, A)), 2))

