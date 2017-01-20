# Mina Jafari
# 19-01-2017
# This script reads optimized files, finds unique ones, and sorts them based on
# energy. TODO
# Step 4 of surface ZStruct: reading in optimized structures and generating driving coordinates
# and initial####.xyz files containing binding sites for SE-GSM calculations.

'''
Python modules
'''
import sys 
sys.path.append('/export/zimmerman/paulzim/ase')
# this is for parallel gsm, so VASP does not read numOfThreads from submit script
import os
os.environ['OMP_NUM_THREADS'] = '1' 
import shutil
from os import listdir
from os.path import isfile, join
from string import Template
from subprocess import call
import numpy as np
#import pandas as pd
import csv

'''
ASE modules
'''
from ase import Atoms,Atom
from ase.calculators.emt import EMT 
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write, read
from ase.lattice.surface import surface, add_adsorbate, fcc100, fcc110, fcc111,\
        hcp0001, bcc100, bcc110, bcc111

'''
functions and main()
'''
print ("")

'''
set up slab and adsorbates by reading it from input file

slab = fcc100('Cu', size=(4,3,2), vacuum=14)
adsorbate= Atoms('NH3OH2')
add_adsorbate(slab, adsorbate, 1.7, 'ontop')
#current position read in
#slabatoms = read("output-46-55-0-0.xyz")
slabatoms = read("output-25-28-0-0-in.xyz")
slab.set_positions(slabatoms.get_positions())

mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
calc = Vasp(xc='PBE', lreal='Auto', kpts=[1,1,1], ismear=1, sigma=0.2, algo='fast', istart=0, npar=2, encut=300)
#calc = EMT()
slab.set_calculator(calc)


dyn = QuasiNewton(slab, trajectory='trajectory.traj')
dyn.run(fmax=0.05)

wriSlabAtomste("${outputName}", slab)
'''

# returns a list of all the files in a directory
def listFiles(folderPath):
    alignedStructures = [file for file in listdir(folderPath) if isfile(join(folderPath, file))]
    return (alignedStructures)

def readInputFile():
    #inputFile = np.genfromtxt("INPUT", dtype='str', delimiter="\t", usecols = range(0,2))
    inputFile = np.genfromtxt("INPUT", dtype='str', delimiter=" +")#, usecols=(0,1))
    
    findSites = int(inputFile[0].split()[1])
    slabFile = inputFile[1].split()[1]
    numOfAdsorbates = int(inputFile[2].split()[1])
    slabIndex1 = int(inputFile[3].split()[1])
    radius1 = float(inputFile[4].split()[1])
    adsorbFile1 = inputFile[5].split()[1]
    adsorbIndex1 = int(inputFile[6].split()[1])
    reactive_1Index1 = int(inputFile[7].split()[1])
    reactive_1Index2 = int(inputFile[7].split()[2])

    slabIndex2 = int(inputFile[9].split()[1])
    radius2 = float(inputFile[10].split()[1])
    adsorbFile2 = inputFile[11].split()[1]
    adsorbIndex2 = int(inputFile[12].split()[1])
    reactive_2Index1 = int(inputFile[13].split()[1])
    reactive_2Index2 = int(inputFile[13].split()[2])

    addMoves = int(inputFile[15].split()[1])
    breakMoves = int(inputFile[16].split()[1])

    return (findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
            adsorbIndex1, reactive_1Index1, reactive_1Index2, slabIndex2, radius2,\
            adsorbFile2, adsorbIndex2, reactive_2Index1, reactive_2Index2, addMoves,\
            breakMoves)



def main():
    # read INPUT file
    findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
    adsorbIndex1, reactive_1Index1, reactive_1Index2, slabIndex2, radius2,\
    adsorbFile2, adsorbIndex2, reactive_2Index1, reactive_2Index2, addMoves,\
    breakMoves = readInputFile()

    # read binding sites
    binding_sites = read("bindingSites.xyz")
    # set tangs = -3 for binding sites
    binding_sites.set_tags(-3)
    # list all unique files in unique-structures directory
    folderName = "unique-structures/"
    files = listFiles(folderName)

    # define slab
    slab = fcc100('Cu', size=(4,3,2), vacuum=14)
    adsorbate= Atoms('NH3OH2')
    add_adsorbate(slab, adsorbate, 1.7, 'ontop')

    # update slab coordinates from optimized files
    for file in files:
        #current position read in
        slab_atoms = read(folderName + file)
        slab.set_positions(slab_atoms.get_positions())
        # check if all tags are zero (when slab is not defined and only read form
        # xyz file)
        if (all(v == 0 for v in slab.get_tags())):
            print "111"
            # call a function to tag atoms
        # uni-molecular reaction
        print numOfAdsorbates
        if (numOfAdsorbates == 1):
            print "stuff"
        # bi-molecular reaction
        elif (numOfAdsorbates == 2):
            print "stuff"


if __name__ == "__main__":
    main()


'''


#print (slab[0].get_symbol())

print (slab[0].symbol)

print (slab.get_tags())
binding_sites.set_tags(-1)
print (binding_sites.get_tags())
#slab_2 = fcc111('Al', size=(2,2,3), vacuum=10.0)
#print (slab_2.get_tags())
vasp_slab = vasp.read_vasp('POSCAR')
print (vasp_slab.get_tags())


list files in the unique geometry directory
open files and create slab object
if all tags == 0
set tags based on INPUT file
open binding site file
if name == X, set tag = -3
if (tag1 != tag2)
if num of adsorbates == 2
    generate add and break moves (bimolecular)
elif num of adsorbates == 1
    find nearby sites, how????????
    generate combinations



find number of slab atoms in x, y, z
'''










