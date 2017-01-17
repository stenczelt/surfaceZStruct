# Mina Jafari
# 17-01-2017
# This script reads all the output files generated after aligning structures
# and runs geometry optimization on them

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

'''
ASE modules
'''
from ase.lattice.surface import surface
from ase import Atoms,Atom
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write
from ase.io import read
from ase.lattice.surface import add_adsorbate, fcc100, fcc110, fcc111, hcp0001,\
        bcc100, bcc110, bcc111

'''
functions and main()
'''
# returns a list of all the files in a directory
def listFiles(folderPath):
    alignedStructures = [file for file in listdir(folderPath) if isfile(join(folderPath, file))]
    return (alignedStructures)

#def runOptimization()

#def findUniques
#def filterUniques

def main():

    # list of files in aligned-structures folder
    folderPath = "aligned-structures"
    alignedStructures = listFiles(folderPath)
    cwd = os.getcwd()

    #count number of binding sites
    numOfBSites = 0
    fh = open("aligned-structures/" + alignedStructures[0], "r")
    lines = fh.readlines()
    fh.close()
    for line in lines:
        if line[0] == "X":
            numOfBSites += 1
    numOfBSites = int(lines[0]) - numOfBSites


    for element in alignedStructures:
        #open file
        fh = open("aligned-structures/" + element, "r")
        lines = fh.readlines()
        fh.close()
        #create new directory
        newFolder = "optimization/" + element.split(".")[0] + "/"
        if not os.path.exists(newFolder):
            os.makedirs(newFolder)
        #remove X and write to new directory
        fileName = newFolder + element.split(".")[0] + "-in.xyz"
        fh = open(fileName,"w")
        fh.write( str(numOfBSites) )
        fh.write("\n")
        for i in range(1, len(lines)):
            if lines[i][0] != "X":
                fh.write(lines[i])
        fh.close()

    '''
    generate all the necessary input files for ASE/VASP optimization. Enable
    EMT calculator too.
    '''
        
        #read file using ASE
        #fileName = "aligned-structures/noBS/" + element
        #slab = read(fileName, format='xyz')
        #set up ASE objects
        #submit optimizations in new folders in parallel

if __name__ == "__main__":
    main()
