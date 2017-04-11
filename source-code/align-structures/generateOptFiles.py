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
from string import Template
from subprocess import call

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
# returns a list of all the files in a directory
def listFiles(folderPath):
    alignedStructures = [file for file in listdir(folderPath) if isfile(join(folderPath, file))]
    return (alignedStructures)

#def runOptimization()

#def findUniques
#def filterUniques

def main():

    # list of files in aligned-structures folder
    folderPath = "aligned_structures"
    alignedStructures = listFiles(folderPath)
    cwd = os.getcwd()

    #count number of binding sites
    numOfBSites = 0
    fh = open("aligned_structures/" + alignedStructures[0], "r")
    lines = fh.readlines()
    fh.close()
    for line in lines:
        if line[0] == "X":
            numOfBSites += 1
    numOfAtoms = int(lines[0]) - numOfBSites


    for element in alignedStructures:
        #open file
        fh = open("aligned_structures/" + element, "r")
        lines = fh.readlines()
        fh.close()
        #create new directory
        newFolder = "optimization/" + element.split(".")[0] + "/"
        if not os.path.exists(newFolder):
            os.makedirs(newFolder)
        #remove X and write to new directory
        fileName = newFolder + element.split(".")[0] + "-in.xyz"
        fh = open(fileName,"w")
        fh.write( str(numOfAtoms) )
        fh.write("\n")
        for i in range(1, len(lines)):
            if lines[i][0] != "X":
                fh.write(lines[i])
        fh.close()

        '''
        generate all the necessary input files for ASE/VASP optimization (.xyz file,
        run.py, and submit.qsh). These files are generated using templates.
        http://stackoverflow.com/questions/6385686/python-technique-or-simple-templating-system-for-plain-text-output
        '''
        # setting up run.py files
        fh = open("run.py")
        templateFile = Template(fh.read())
        fh.close()
        input = element.split(".")[0] + "-in.xyz"
        outputName = element
        myDictionary={'input':input, 'outputName':outputName}
        #do the substitution
        result = templateFile.substitute(myDictionary)
        pyFile = newFolder + "/run.py"
        fh = open(pyFile,"w")
        fh.write(result)
        fh.close()

        # setting up submit.qsh files
        fh = open("submit.qsh")
        templateFile2 = Template(fh.read())
        fh.close()
        # user can replace this with a meaningful name
        jobName = "test"
        fileNumber = (element.split(".")[0]).replace("output-", "")
        myDictionary2={'jobName':jobName, 'fileNumber':fileNumber}
        result2 = templateFile2.substitute(myDictionary2)
        PBSfile = newFolder + "/submit.qsh"
        fh = open(PBSfile,"w")
        fh.write(result2)
        fh.close()

        #submit optimizations
        os.chdir(newFolder)
        call(["qsub", "submit.qsh"])
        os.chdir(cwd)

if __name__ == "__main__":
    main()
