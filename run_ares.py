import ares 
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
        type = float, default = 0.6774)
parser.add_argument('--helium_by_mass', metavar = 'helium_by_mass', 
        type = float, default = 0.2453)
parser.add_argument('--cmb_temp_0', metavar = 'cmb_temp_0', 
        type = float, default = 2.7255)
parser.add_argument('--sigma_8', metavar = 'sigma_8', 
        type = float, default = 0.8159)
parser.add_argument('--primordial_index', metavar = 'primordial_index', 
        type = float, default = 0.9667)
parser.add_argument('--fstar', metavar = 'fstar', 
        type = float, default = 0.1)
parser.add_argument('--Tmin', metavar = 'Tmin', 
        type = float, default = 10000)
parser.add_argument('--Nion', metavar = 'Nion', 
        type = float, default = 4000)
parser.add_argument('--fesc', metavar = 'fesc', 
        type = float, default = 0.1)
parser.add_argument('--Nlw', metavar = 'Nlw', 
        type = float, default = 9690)
parser.add_argument('--cX', metavar = 'cX', 
        type = float, default = 3.4E40)
parser.add_argument('--fX', metavar = 'fX', 
        type = float, default = 0.2)

args = parser.parse_args()

#This file should be able to read parameters from either command line or a file.
#Then at the end it should output the dTb(z) to some file which can be read in by
# the c++ code.

sim = ares.simulations.Global21cm(verbose = False, omega_m_0 = args.omega_m_0,\
        omega_b_0 = args.omega_b_0, omega_l_0 = args.omega_l_0,\
        hubble_0 = args.hubble_0, helium_by_mass = args.helium_by_mass,\
        cmb_temp_0 = args.cmb_temp_0, sigma_8 = args.sigma_8,\
        primordial_index = args.primordial_index, fstar = args.fstar,\
        Tmin = args.Tmin, Nion = args.Nion, fesc = args.fesc, Nlw = args.Nlw,\
        cX = args.cX, fX = args.fX)
sim.run()

with open("dTb_ares.dat", 'w') as file:
    for i in range(0,len(sim.history['z'])):
        z = sim.history['z'][i]
        dTb_z = sim.history['dTb'][i]
        file.write(str(z) + " " + str(dTb_z) + "\n")
