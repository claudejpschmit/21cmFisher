from hmf import MassFunction
from hmf import cosmo
from astropy.cosmology import LambdaCDM
import numpy as np
import sys
import argparse

############## Parsing input ##############
descr = 'This program runs ARES for the Fisher code.'

parser = argparse.ArgumentParser(description=descr)
parser.add_argument('--version', action='version', 
        version='%(prog)s v0.1')
parser.add_argument('--omega_m_0', metavar = 'omega_m_0', 
        type = float, default = 0.3089)
parser.add_argument('--omega_b_0', metavar = 'omega_b_0', 
        type = float, default = 0.0486)
parser.add_argument('--omega_l_0', metavar = 'omega_l_0',
        type = float, default = 0.6911)
parser.add_argument('--hubble_0', metavar = 'hubble_0', 
        type = float, default = 67.74)
parser.add_argument('--cmb_temp_0', metavar = 'cmb_temp_0', 
        type = float, default = 2.7255)
parser.add_argument('--redshift', metavar = 'redshift', 
        type = float, default = 0.0)
parser.add_argument('--Mmin', metavar = 'Mmin', 
        type = int, default = 5)
parser.add_argument('--Mmax', metavar = 'Mmax', 
        type = int, default = 15)


args = parser.parse_args()

new_model = LambdaCDM(H0=args.hubble_0, Om0 = args.omega_m_0, Tcmb0 = args.cmb_temp_0,\
        Ob0 = args.omega_b_0, Ode0 = args.omega_l_0)
my_cosmo = cosmo.Cosmology(cosmo_model = new_model)
hmf = MassFunction()
print ("redshift called: ", args.redshift)
hmf.update(Mmin = args.Mmin, Mmax = args.Mmax, z=args.redshift, hmf_model="ST")

#print "Matter density: ", my_cosmo.cosmo.Om0
#print "Hubble constant: ", my_cosmo.cosmo.H0
#print "Dark Energy density: ", my_cosmo.cosmo.Ode0
#print "Baryon density: ",  my_cosmo.cosmo.Ob0
#print "Curvature density: ", my_cosmo.cosmo.Ok0

f = open("hmf.dat", "w")

for n in range(0,len(hmf.M)):
    a = hmf.M[n]
    b = hmf.dndm[n]
    f.write(str(a) + " " + str(b) + "\n")

f.close()
