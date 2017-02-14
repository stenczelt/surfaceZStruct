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
import math

'''
ASE modules
'''
from ase import Atoms,Atom
from ase.calculators.emt import EMT 
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write, read
from ase.build import surface, add_adsorbate, fcc100, fcc110, fcc111,\
        hcp0001, bcc100, bcc110, bcc111

'''
Open Babel modules
'''
import openbabel as OB

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
    reactiveIndices_1 = []
    reactiveIndices_2 = []

    # reading max of 10 indices
    print ("RRR ", len(inputFile[7].split()))
    if ( len(inputFile[7].split()) > 11 ):
        print ("ERROR: Maximum of 10 reactive atoms on each adsorbate are allowed!")
        sys.exit(-1)
    else:
        for i in range(1, len(inputFile[7].split()) ):
            reactiveIndices_1.append( int(inputFile[7].split()[i]) )
    #reactive_1Index1 = int(inputFile[7].split()[1])
    #reactive_1Index2 = int(inputFile[7].split()[2])

    #if (numOfAdsorbates == 2):
    slabIndex2 = int(inputFile[9].split()[1])
    radius2 = float(inputFile[10].split()[1])
    adsorbFile2 = inputFile[11].split()[1]
    adsorbIndex2 = int(inputFile[12].split()[1])
    if ( len(inputFile[13].split()) > 11 ):
        print ("ERROR: Maximum of 10 reactive atoms on each adsorbate are allowed!")
        sys.exit(-1)
    else:
        for i in range(1, len(inputFile[13].split()) ):
            reactiveIndices_2.append( int(inputFile[13].split()[i]) )
    #reactive_2Index1 = int(inputFile[13].split()[1])
    #reactive_2Index2 = int(inputFile[13].split()[2])

    addMoves = int(inputFile[15].split()[1])
    breakMoves = int(inputFile[16].split()[1])

    return (findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
            adsorbIndex1, reactiveIndices_1, slabIndex2, radius2,\
            adsorbFile2, adsorbIndex2, reactiveIndices_2, addMoves, breakMoves)

'''
recursive 
def generateIsomers(numOfBreakMoves, numOfAddMoves, reactiveIndices_1, reactiveIndices_2,\
        currentSet):
    for i in range(0, numOfBreakMoves):
        generateIsomers(numOfBreakMoves, numOfAddMoves, reactiveIndices_1[1:],\
                (reactiveIndices_2[:i] + reactiveIndices_2[(i+1):]),\
                (currentSet.append(reactiveIndices_1[0]).append(reactiveIndices_2[i])))
#        currentSet.append(reactiveIndices_1[0], reactiveIndices_2[i])
        if i == 2:
            return currentSet

'''

# creates a list of lists. Each list has indices of two atoms that wil be added
def generateIsomerPair(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
        numOfAds1Atoms, numOfAds2Atoms):
    # check reactive indices are <= than number of atoms
    assert( reactiveIndices_1[0] <= numOfAds1Atoms and reactiveIndices_1[1] <= numOfAds1Atoms  )
    assert( reactiveIndices_2[0] <= numOfAds2Atoms and reactiveIndices_2[1] <= numOfAds2Atoms  )
    listOfAdds = []
    #for i in range(0, len(reactiveIndices_1)):
    for index_1 in reactiveIndices_1:
        #for j in range(0, len(reactiveIndices_2)):
        for index_2 in reactiveIndices_2:
            list = []
            if (index_1 != 0 and index_2 != 0):
                list.append(index_1 + numOfSlabAtoms - 1)
                list.append(index_2 + numOfSlabAtoms + numOfAds1Atoms - 1)
                listOfAdds.append(list)
                # TODO
                '''
                if (index_1 . coordination() == max_coordination ):
                    BREAK index_1[0] index_1[1]
                    if index_1[0] not connected to index_1[1]:
                        raise Error
                '''
    return listOfAdds

# create lists of two pairs of add moves
def generate2IsomerPairs(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
        numOfAds1Atoms, numOfAds2Atoms):
    pairs = generateIsomerPair(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
            numOfAds1Atoms, numOfAds2Atoms)
    lists = []
    # add unique pairs to the list
    for i in range(0, len(pairs)):
        for j in range(i, len(pairs)):
            if pairs[i][0] != pairs[j][0] and\
               pairs[i][1] != pairs[j][1]:
                list = []
                list.append(pairs[i])
                list.append(pairs[j])
                lists.append(list)
    return lists

# find sites that are within MAX_DISTANCE of a given adsorbate index
def findNearbySites(slab, atomIndex, bindingSiteFile):
    fh = open(bindingSiteFile, 'r')
    lines = fh.readlines()
    fh.close()

    listOfSites = []

    MAX_DISTANCE = 2.0 # Angstroms
    THRESH = 0.1

    atomX = slab[atomIndex - 1].x
    atomY = slab[atomIndex - 1].y
    atomZ = slab[atomIndex - 1].z
    lines.pop(0)
    lines.pop(0)

    wrongIndices = []
    for i in range(0, len(lines)):
        inX, inY, inZ = lines[i].split()[1:4]
        inX = float(inX)
        inY = float(inY)
        inZ = float(inZ)

        distance = (atomX - inX)**2 + (atomY - inY)**2 + (atomZ - inZ)**2
        if ( math.sqrt(distance) <= MAX_DISTANCE and math.sqrt(distance) > 0.01 ):
            for atom in slab:
                if ( (atom.x+THRESH >= inX and atom.x-THRESH <= inX)\
                        and (atom.y+THRESH >= inY and atom.y-THRESH <= inY)\
                        and (atom.z+THRESH >= inZ and atom.z-THRESH <= inZ) ):
                    wrongIndices.append(i)

                else:
                    listOfSites.append(i)

    listOfSites = set(listOfSites)
    for element in wrongIndices:
        if element in listOfSites:
            listOfSites.remove(element)
    return listOfSites


# find two indices from two input lists that are farthest from eachother
def findFarthestIndices(list1, list2, bindingSiteFile):
    # read as ASE Atoms object
    binding_sites = read(bindingSiteFile)
    outIndex = []
    distance = 0.0
    for index_1 in list1:
        for index_2 in list2:
            x_1 = float(binding_sites[index_1].x)
            y_1 = float(binding_sites[index_1].y)
            z_1 = float(binding_sites[index_1].z)

            x_2 = float(binding_sites[index_2].x)
            y_2 = float(binding_sites[index_2].y)
            z_2 = float(binding_sites[index_2].z)

            new_distance = math.sqrt( (x_1 - x_2)**2 + (y_1 - y_2)**2 + (z_1 - z_2)**2 )
            if (new_distance > distance):
                distance = new_distance
                outIndex = [ index_1, index_2 ]

    return outIndex

def areConnected(inputFile, twoIndices):
    obConversion = OB.OBConversion()
    obConversion.SetInFormat("xyz")
    adsorbate = OB.OBMol()
    obConversion.ReadFile(adsorbate, inputFile)

    index_1 = twoIndices[0]
    index_2 = twoIndices[1]
    #print (adsorbate.GetAtom( twoIndices[0] ).GetType())
    #print (adsorbate.GetAtom( twoIndices[1] ).GetType())
    #if (adsorbate.GetBond(index_1, index_2).GetBondOrder() != None):
    if (adsorbate.GetBond(index_1, index_2) != None):
        return True;
    else:
        return False;
    #print (adsorbate.NumBonds())




# TODO write function to find number of atoms in each layer, x, y


def main():
    # read INPUT file
    findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
    adsorbIndex1, reactiveIndices_1, slabIndex2, radius2,\
    adsorbFile2, adsorbIndex2, reactiveIndices_2, addMoves, breakMoves = readInputFile()

    # combine reactive indices into two lists for each adsorbate
    #reactiveIndices_1 = [ reactive_1Index1, reactive_1Index2 ]
    #reactiveIndices_2 = [ reactive_2Index1, reactive_2Index2 ]

    # read number of slab atoms
    with open(slabFile, 'r') as fh:
        numOfSlabAtoms = int(fh.readline().strip())

    # read number of adsorbate 1 atoms
    with open(adsorbFile1, 'r') as fh:
        numOfAds1Atoms = int(fh.readline().strip())

    # read number of adsorbate 1 atoms
    if (numOfAdsorbates == 2):
        with open(adsorbFile2, 'r') as fh:
            numOfAds2Atoms = int(fh.readline().strip())

    # read number of binding site atoms
    with open("bindingSites.xyz", 'r') as fh:
        numOfBSAtoms = int(fh.readline().strip())

    # read binding sites
    binding_sites = read("bindingSites.xyz")
    # set tangs = -3 for binding sites
    binding_sites.set_tags(-3)
    # list all unique files in unique-structures directory
    folderName = "unique-structures/"
    files = listFiles(folderName)

    # define slab
    slab = fcc100('Cu', size=(4,3,2), vacuum=14)
    #adsorbate= Atoms('NH3OH2')
    adsorbate = Atoms('NH3')
    add_adsorbate(slab, adsorbate, 1.7, 'ontop')

    # current working directory
    cwd = os.getcwd()

    # update slab coordinates from optimized files
    for file in files:
        #current position read in
        slab_atoms = read(folderName + file)
        slab.set_positions(slab_atoms.get_positions())
        # check if all tags are zero (when slab is not defined and only read form
        # xyz file)
        if (all(v == 0 for v in slab.get_tags())):
            print ("111")
            # TODO call a function to tag atoms
            # TODO ads1 = -1, ads2 = -2, BS = -3. Find ads1 and ads2 atoms from INPUT or numberOfAds1/2Atoms

        # uni-molecular reaction
        # TODO check number of adds and breaks match number of reactive atoms. Raise error if not.
        if (numOfAdsorbates == 1):
            # TODO make sure the site is empty
            '''
            if one adsorbate, we have two reactive atoms at a time. Break the bond
            between them and add each fragment to a new site.
            - break bond
            - find nearby sites
            - make sure the site is empty: at this stage of development we're skipping this step
            - find smaller fragment ( needed ? )
            - attach two fragments to sites that are > 2.0 A away
            - push two fragments in opposite directions, satisfied by above step
            - if two fragments have similar size, always move one of them
            - generate driving coordinates
            '''
            # create slab to be written to initial000# file
            slab_out = slab[0:numOfSlabAtoms]
            #binding_sites = read("bindingSites.xyz", format='xyz')
            for i in range(0, len(binding_sites)):
                slab_out.append(binding_sites[i])
            for i in range(0+numOfSlabAtoms, numOfAds1Atoms+numOfSlabAtoms):
                slab_out.append(slab[i])

            # find all combinations of reactive indices
            breakCombos = []
            for i in range(0, len(reactiveIndices_1)):
                for j in range(i+1, len(reactiveIndices_1)):
                    temp = []
                    temp.append(reactiveIndices_1[i])
                    temp.append(reactiveIndices_1[j])
                    breakCombos.append(temp)

            # find sites nearby a given adsorbate index
            i = 1
            for element in breakCombos:
                connected = areConnected(adsorbFile1, element)
                if (connected == True):
                    indexIn = element[0] + numOfSlabAtoms
                    listOfSitesFrag_1 = findNearbySites(slab, indexIn, "bindingSites.xyz")
                    indexIn = element[1] + numOfSlabAtoms
                    listOfSitesFrag_2 = findNearbySites(slab, indexIn, "bindingSites.xyz")
                    # find two sites that are farthest from two lists of binding sites
                    farthestSites = findFarthestIndices(listOfSitesFrag_1,\
                        listOfSitesFrag_2, "bindingSites.xyz")
                    # TODO add to different binding site types
                    # write initial### and ISOMERS file
                    #i = 1
                    folder = "se_gsm_cals/" + file.split(".")[0] + "/" +\
                        str(i).zfill(4) + "/scratch/"
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                    fh.write("NEW\n")
                    firstNum = element[0] + numOfSlabAtoms + numOfBSAtoms
                    secondNum = element[1] + numOfSlabAtoms + numOfBSAtoms
                    fh.write("BREAK " + str(firstNum) + "  " + str(secondNum) + "\n")
                    fh.write("ADD   " + str(farthestSites[0] + 1 + numOfSlabAtoms) + "  " +\
                        str(firstNum) + "\n")
                    fh.write("ADD   " + str(farthestSites[1] + 1 + numOfSlabAtoms) + "  " +\
                        str(secondNum) + "\n")
                    fh.close()

                    slab_out.write(folder + "initial" + str(i).zfill(4) + ".xyz")

                    # submit SE-GSM
                    # files needed for SE-GSM calculation: inpfileq, grad.py,
                    # status, gfstringq.exe, scratch/submit_gsm.qsh, scratch/initial000.xyz,
                    # and scratch/ISOMERS000
                    folder = "se_gsm_cals/" + file.split(".")[0] + "/" +\
                        str(i).zfill(4)
                    src = "inputs_se_gsm/"
                    src_files = listFiles(src)
                    for file_name in src_files:
                        full_file_name = os.path.join(src, file_name)
                        #if (os.path.isfile(full_file_name)):
                        shutil.copy(full_file_name, folder)

                    # setting up runGSM.qsh files
                    fh = open("inputs_se_gsm/scratch/runGSM.qsh")
                    templateFile = Template(fh.read())
                    # user can replace this with a meaningful name
                    jobName = "test" + str(i).zfill(4)
                    print (jobName)
                    #fileNumber = (element.split(".")[0]).replace("output-", "") 
                    #myDictionary = {'jobName':jobName, 'fileNumber':fileNumber}
                    myDictionary = {'jobName':jobName, 'jobID':i}
                    result = templateFile.substitute(myDictionary)
                    PBSfile = folder + "/scratch/runGSM.qsh"
                    fh = open(PBSfile,"w")
                    fh.write(result)
                    fh.close()

                    os.chdir(folder)
                    call(["qsub", "scratch/runGSM.qsh"])
                    os.chdir(cwd)
                    i += 1


        # bi-molecular reaction
        # TODO check number of adds and breaks match number of reactive atoms. Raise error if an add
        # move requires breaking a bond but number of breaks is zero.
        elif (numOfAdsorbates == 2):
            # generate breaks
            # generate adds
            pairsSet = generateIsomerPair(reactiveIndices_1, reactiveIndices_2,\
                    numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms)

            if addMoves == 2:
                twoPairsSet = generate2IsomerPairs(reactiveIndices_1, reactiveIndices_2,\
                        numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms)

            '''
            write combinations to file:
            create new directories for each input geometry file and isomer
            combination.
            '''
            # TODO convert this to function
            i = 1
            for element in pairsSet:
                folder = "se_gsm_cals/" + file.split(".")[0] + "/" +\
                        str(i).zfill(4) + "/scratch/"
                if not os.path.exists(folder):
                    os.makedirs(folder)
                fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                fh.write("NEW\n")
                fh.write("ADD  " + str(element[0]) + "  "+ str(element[1]) )
                fh.close()

                slab.write(folder + "initial" + str(i).zfill(4) + ".xyz")

                i += 1

            if addMoves == 2:
                for element in twoPairsSet:
                    print ("element", element)
                    print (element[0][0], "      ", element[0][1])
                    print (element[1][0], "      ", element[1][1])
                    folder = "se_gsm_cals/" + file.split(".")[0] + "/" +\
                            str(i).zfill(4) + "/scratch/"
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                    fh.write("NEW\n")
                    fh.write("ADD  " + str(element[0][0]) + "  "+ str(element[0][1]) + "\n" )
                    fh.write("ADD  " + str(element[1][0]) + "  "+ str(element[1][1]) + "\n" )
                    fh.close()
                    slab.write(folder + "initial" + str(i).zfill(4) + ".xyz")
                    i += 1
            '''
            for i in range(0, breakMoves):
                print "BREAK    ", 
            for i in range(0, addMoves):
                for index_1 in reactiveIndices_1:
                    if (index_1 != 0):
                        index_1 += numOfSlabAtoms
                    else:
                        continue
                    for index_2 in reactiveIndices_2:
                        if (index_2 != 0):
                            index_2 += (numOfSlabAtoms + numOfAds1Atoms)
                            print "ADD  ", index_1, "      ", index_2
                        else:
                            continue
                            '''




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










