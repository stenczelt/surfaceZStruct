# Mina Jafari
# 19-01-2017
# This script reads the unique structures and creates driving coordinates and runs SE-GSM.
# Step 5 of surface ZStruct: reading in optimized structures and generating driving coordinates
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
import csv
import math

from CoordinationNumbers import CoordinationNumbers

'''
ASE modules
'''
from ase import Atoms,Atom
from ase.calculators.emt import EMT 
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write, read
import ase.io.vasp
#from ase.build import surface, add_adsorbate, fcc100, fcc110, fcc111,\
from ase.lattice.surface import surface, add_adsorbate, fcc100, fcc110, fcc111,\
        hcp0001, bcc100, bcc110, bcc111

'''
Open Babel modules
'''
import openbabel as OB

'''
functions and main()
'''

# returns a list of all the files in a directory
def listFiles(folderPath):
    alignedStructures = [file for file in listdir(folderPath) if isfile(join(folderPath, file))]
    return (alignedStructures)

def readInputFile():
    inputFile = np.genfromtxt("INPUT", dtype='str', delimiter=" +")
    
    findSites = int(inputFile[0].split()[1])
    slabFile = inputFile[1].split()[1]
    numOfAdsorbates = int(inputFile[2].split()[1])
    slabIndex1 = int(inputFile[3].split()[1])
    radius1 = float(inputFile[4].split()[1])
    adsorbFile1 = inputFile[5].split()[1]
    adsorbIndex1 = int(inputFile[6].split()[1])
    reactiveIndices_1 = []
    reactiveIndices_2 = []
    surfaceReactiveIndices = []

    # reading max of 10 indices 11 = 10 indices + first column of the file
    if ( len(inputFile[7].split()) > 11 ):
        print ("ERROR: Maximum of 10 reactive atoms on each adsorbate are allowed!")
        sys.exit(-1)
    else:
        for i in range(1, len(inputFile[7].split()) ):
            reactiveIndices_1.append( int(inputFile[7].split()[i]) )

    #if (numOfAdsorbates == 2):
    slabIndex2 = int(inputFile[10].split()[1])
    radius2 = float(inputFile[11].split()[1])
    adsorbFile2 = inputFile[12].split()[1]
    adsorbIndex2 = int(inputFile[13].split()[1])
    # reading max of 10 indices 11 = 10 indices + first column of the file
    if ( len(inputFile[14].split()) > 11 ):
        print ("ERROR: Maximum of 10 reactive atoms on each adsorbate is allowed!")
        sys.exit(-1)
    else:
        for i in range(1, len(inputFile[14].split()) ):
            reactiveIndices_2.append( int(inputFile[14].split()[i]) )

    addMoves = int(inputFile[17].split()[1])
    breakMoves = int(inputFile[18].split()[1])
    # reactive slab indices
    if ( len(inputFile[19].split()) > 5 ):
        print ("ERROR: Maximum of 4 reactive atoms on surface is allowed!")
        sys.exit(-1)
    else:
        for i in range(1, len(inputFile[19].split()) ):
            surfaceReactiveIndices.append( int(inputFile[19].split()[i]) )
    reactionType = int(inputFile[20].split()[1])

    return (findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
            adsorbIndex1, reactiveIndices_1, slabIndex2, radius2,\
            adsorbFile2, adsorbIndex2, reactiveIndices_2, addMoves, breakMoves,\
            surfaceReactiveIndices, reactionType)

# creates a list of lists. Each list has indices of two atoms that wil be added
# this function is for the case with reactive indices from surface
def generateIsomerPair_1(reactiveIndices_1, surfaceReactiveIndices, numOfSlabAtoms,\
        numOfAds1Atoms):
    # check reactive indices are <= than number of atoms
    assert( reactiveIndices_1[0] <= numOfAds1Atoms and reactiveIndices_1[1] <= numOfAds1Atoms  )
    assert( surfaceReactiveIndices[0] <= numOfSlabAtoms and surfaceReactiveIndices[1] <= numOfSlabAtoms )
    listOfAdds = []
    for index_1 in reactiveIndices_1:
        for index_2 in surfaceReactiveIndices:
            list = []
            if (index_1 != 0 and index_2 != 0):
                list.append(index_1 + numOfSlabAtoms)
                list.append(index_2)
                listOfAdds.append(list)

    return listOfAdds

# creates a list of lists. Each list has indices of two atoms that wil be added
def generateIsomerPair_2(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
        numOfAds1Atoms, numOfAds2Atoms):
    # check reactive indices are <= than number of atoms
    assert( reactiveIndices_1[0] <= numOfAds1Atoms and reactiveIndices_1[1] <= numOfAds1Atoms  )
    assert( reactiveIndices_2[0] <= numOfAds2Atoms and reactiveIndices_2[1] <= numOfAds2Atoms  )
    listOfAdds = []
    for index_1 in reactiveIndices_1:
        for index_2 in reactiveIndices_2:
            list = []
            if (index_1 != 0 and index_2 != 0):
                list.append(index_1 + numOfSlabAtoms)
                list.append(index_2 + numOfSlabAtoms + numOfAds1Atoms)
                listOfAdds.append(list)

    return listOfAdds

# create lists of two pairs of add moves
def generate2IsomerPairs(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
        numOfAds1Atoms, numOfAds2Atoms, reactionType):
    if (reactionType == 1):
        pairs = generateIsomerPair_1(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
                numOfAds1Atoms)
    else:
        pairs = generateIsomerPair_2(reactiveIndices_1, reactiveIndices_2, numOfSlabAtoms,\
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
    return list(listOfSites)


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
    if (adsorbate.GetBond(index_1, index_2) != None):
        return True;
    else:
        return False;

def getCoordinationNum(inputFile, indexIn):
    obConversion = OB.OBConversion()
    obConversion.SetInFormat("xyz")
    adsorbate = OB.OBMol()
    obConversion.ReadFile(adsorbate, inputFile)

    anAtom = adsorbate.GetAtom(indexIn)
    return anAtom.GetValence()

def getAtomicNumber(inputFile, indexIn):
    obConversion = OB.OBConversion()
    obConversion.SetInFormat("xyz")
    adsorbate = OB.OBMol()
    obConversion.ReadFile(adsorbate, inputFile)

    anAtom = adsorbate.GetAtom(indexIn)
    return anAtom.GetAtomicNum()

def createTemplateFiles(file, folder, index, slab, slabType):
    slab.write(folder + "initial" + str(index).zfill(4) + ".xyz")

    # Read slab type form slab file input (surface.xyz)
    fh = open(folder + "initial" + str(index).zfill(4) + ".xyz", 'r')
    lines = fh.readlines()
    lines[1] = slabType + '\n'
    fh.close()

    # read slab size and number of frozen atoms from run.py
    fh = open("run.py", 'r')
    run_lines = fh.readlines()
    fh.close()
    atomsInX = 0
    atomsInY = 0
    atomsInZ = 0
    frozenLayers = 0
    for line in run_lines:
        if "size" in line:
            splitted_line = line.split(',', -1)
            atomsInX = int( splitted_line[1].split('(')[1].strip() )
            atomsInY = int( splitted_line[2].strip() )
            atomsInZ = int( splitted_line[3].split(')')[0].strip() )
        if "atom.tag" in line:
            frozenLayers = int( line.split('>')[1].split('f')[0].strip() )
    numOfFrozenAtoms = atomsInX * atomsInY * (atomsInZ - frozenLayers)
    
    for k in range(2, numOfFrozenAtoms+2):
        lines[k] = lines[k].rstrip()
        lines[k] = lines[k] + " *\n"

    fh = open(folder + "initial" + str(index).zfill(4) + ".xyz", 'w')
    fh.writelines(lines)
    fh.close()

    # setting up runGSM.qsh files
    fh = open("inputs_se_gsm/scratch/runGSM.qsh")
    templateFile = Template(fh.read())
    fh.close()
    jobName = file.split(".")[0].split('-', 1)[1] + str(index).zfill(4)
    myDictionary = {'jobName':jobName, 'jobID':index}
    result = templateFile.substitute(myDictionary)
    PBSfile = folder + "/runGSM.qsh"
    fh = open(PBSfile,"w")
    fh.write(result)
    fh.close()

def submitSE_GSM(file, extension, index, cwd):
    # submit SE-GSM
    folder = "se_gsm_cals_" + str(extension) + "/" + file.split(".")[0] + "/" +\
            str(index).zfill(4)
    src = "inputs_se_gsm/"
    src_files = listFiles(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        shutil.copy(full_file_name, folder)

    os.chdir(folder)
    call(["qsub", "scratch/runGSM.qsh"])
    #call(["sbatch", "scratch/runGSM.qsh"])
    os.chdir(cwd)

# read input files and set slab size--number of atoms in X, Y, Z directions
# supported file types are xyz and VASP POSCAR
def findSlabSize(inFile, numOfSlabAtoms):
    numOfAtomsInX = 0
    numOfAtomsInY = 1
    TOLERANCE = 0.15
    fileType = "xyz"

    fh = open(inFile)
    lines = fh.readlines()
    fh.close()

    numOfAtoms = 0
    try:
        numOfAtoms = int(lines[0])
    except:
        TypeError
        fileType = "POSCAR"
        print("ERROR: File type should be xyz")
        exit(-1)

    #if (fileType == "xyz"):
    xCoordinate = float(lines[2].split()[1])
    yCoordinate = float(lines[2].split()[2])
    zCoordinate = float(lines[2].split()[3])
    currentYcoordinate = float(lines[2].split()[2])
    for i in range(2, len(lines)):
        if (float(lines[i].split()[3]) <= zCoordinate + TOLERANCE and
                float(lines[i].split()[3]) >= zCoordinate - TOLERANCE):
            # atoms in the same z layer
            if (float(lines[i].split()[2]) <= yCoordinate + TOLERANCE and
                    float(lines[i].split()[2]) >= yCoordinate - TOLERANCE):
                numOfAtomsInX += 1

            newYcoordinate = float(lines[i].split()[2])
            if (newYcoordinate >= currentYcoordinate + TOLERANCE or
                    newYcoordinate <= currentYcoordinate - TOLERANCE):
                numOfAtomsInY += 1
                currentYcoordinate = float(lines[i].split()[2])

    numOfAtomsInZ = numOfSlabAtoms / (numOfAtomsInX * numOfAtomsInY)
    return (numOfAtomsInX, numOfAtomsInY, numOfAtomsInZ)



# TODO
#def tagAtoms():
    # find atoms in each layer


def main():
    # read INPUT file
    findSites, slabFile, numOfAdsorbates, slabIndex1, radius1, adsorbFile1,\
    adsorbIndex1, reactiveIndices_1, slabIndex2, radius2,\
    adsorbFile2, adsorbIndex2, reactiveIndices_2, addMoves, breakMoves,\
    surfaceReactiveIndices, reactionType = readInputFile()

    assert(addMoves < 3)
    assert(addMoves > 0)
    assert(breakMoves < 3)
    assert(breakMoves > 0)

    if (reactionType == 1 and numOfAdsorbates == 1):
        if (len(surfaceReactiveIndices) < 1 or surfaceReactiveIndices[0] == 0):
            print ("ERROR: For inter-molecular reaction with one added adsorbate\
                    'reactive-atoms-surface' should be set in INPUT file.")
            exit (-1)

    slabType = ""
    numOfSlabAtoms = 0
    # read number of slab atoms
    with open(slabFile, 'r') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                numOfSlabAtoms = int(line.strip())
            if i == 1:
                slabType = line.strip()

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
    folderName = "unique_structures/"
    files = listFiles(folderName)

    # read slab set up from slab_setup.py
    fh = open("slab_setup.py")
    slab_lines = fh.readlines()
    fh.close()
    POSCAR = False

    slab = fcc111('Cu', size=(1,1,1), vacuum=15)
    adsorbate = Atoms('X')
    for line in slab_lines:
        if "slab_in" in line:
            if '#' not in line:
                slab_type = line.split('=')[1].split('(')[0].strip()
                slab_atom = line.split('=')[1].split("'")[1].strip()
                slab_size_x = int(line.split('=')[2].split('(')[1].split(',')[0])
                slab_size_y = int(line.split('=')[2].split('(')[1].split(',')[1])
                slab_size_z = int(line.split('=')[2].split('(')[1].split(',')[2].split(')')[0])
                slab_vac = float(line.split('=')[3].split(')')[0])
                slab = eval(slab_type)(slab_atom, size=(slab_size_x,slab_size_y,slab_size_z), vacuum=slab_vac)
        if "adsorbate_in" in line:
            if '#' not in line:
                temp = line.split('=')[1].split("'")[1]
                adsorbate = Atoms(temp)
        
        if "slab_POSCAR_file_name" in line:
            if '#' not in line:
                POSCAR = line.split('=')[1].strip().replace("'", "")
                slab = ase.io.vasp.read_vasp(POSCAR)
                POSCAR = True

    try:
        slab
    except NameError:
        print "slab is not defined! Check slab_setup.py file."
        exit (-1)
    if (POSCAR == False):
        try:
            adsorbate
        except NameError:
            print "adsorbate is not defined! Check slab_setup.py file."
            exit (-1)
        else:
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
            findSlabSize(folderName + file, numOfSlabAtoms)
            # TODO call a function to tag atoms
            # TODO ads1 = -1, ads2 = -2, BS = -3. Find ads1 and ads2 atoms from INPUT or numberOfAds1/2Atoms


        ###############################################
        # uni-molecular reaction
        ###############################################
        if (numOfAdsorbates == 1 and reactionType == 0):
            '''###
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
            '''###
            if (addMoves != 2):
                print ("ERROR: You need at least 2 add moves for dissociation of one adsorbate.")
                exit(-1)
            # find all combinations of reactive indices
            breakCombos = []
            for i in range(0, len(reactiveIndices_1)):
                for j in range(i+1, len(reactiveIndices_1)):
                    if ( (reactiveIndices_1[i] != 0) and (reactiveIndices_1[j] != 0) ):
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
                    # TODO make sure the site is empty
                    farthestSites = findFarthestIndices(listOfSitesFrag_1,\
                        listOfSitesFrag_2, "bindingSites.xyz")
                    # TODO add to different binding site types input by the user
                    # write initial### and ISOMERS file
                    folder = "se_gsm_cals_1/" + file.split(".")[0] + "/" +\
                        str(i).zfill(4) + "/scratch/"
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                    fh.write("NEW\n")
                    firstNum = element[0] + numOfSlabAtoms #+ numOfBSAtoms
                    secondNum = element[1] + numOfSlabAtoms #+ numOfBSAtoms
                    fh.write("BREAK " + str(firstNum) + "  " + str(secondNum) + "\n")
                    fh.write("ADD   " + str(farthestSites[0] + numOfSlabAtoms + numOfAds1Atoms + 1)\
                            + "  " + str(firstNum) + "\n")
                    fh.write("ADD   " + str(farthestSites[1] + numOfSlabAtoms + numOfAds1Atoms + 1)\
                            + "  " + str(secondNum) + "\n")
                    fh.close()

                    index = i
                    createTemplateFiles(file, folder, index, slab, slabType)

                    # submit SE-GSM
                    # files needed for SE-GSM calculation: inpfileq, grad.py,
                    # status, gfstringq.exe, scratch/submit_gsm.qsh, scratch/initial000.xyz,
                    # and scratch/ISOMERS000
                    extension = 1
                    submitSE_GSM(file, extension, index, cwd)

                    i += 1



        ###############################################
        ## ligand transfer or ligand exchange reactions
        ###############################################
        #TODO ligand transfer
        elif (numOfAdsorbates == 1 and reactionType == 1):
            if (addMoves != 2 or breakMoves != 2):
                print ("ERROR: You need at least 2 add and 2 break moves for this type of reaction.")
                exit(-1)
            if (len(reactiveIndices_1) < 2 or len(surfaceReactiveIndices) < 2):
                    print ("ERROR: You need at least 2 reactive atoms on each adsorbate for ligand exchange.")
                    exit(-1)

            # find all combinations of reactive indices 1
            breakCombos_1 = []
            for i in range(0, len(reactiveIndices_1)):
                for j in range(i+1, len(reactiveIndices_1)):
                    if ( (reactiveIndices_1[i] != 0) and (reactiveIndices_1[j] != 0) ):
                        temp = []
                        temp.append(reactiveIndices_1[i])
                        temp.append(reactiveIndices_1[j])
                        breakCombos_1.append(temp)
            # find all combinations of slab reactive indices
            breakCombos_surface = []
            for i in range(0, len(surfaceReactiveIndices)):
                for j in range(i+1, len(surfaceReactiveIndices)):
                    if ( (surfaceReactiveIndices[i] != 0) and (surfaceReactiveIndices[j] != 0) ):
                        temp = []
                        temp.append(surfaceReactiveIndices[i])
                        temp.append(surfaceReactiveIndices[j])
                        breakCombos_surface.append(temp)

            for element_1 in breakCombos_1:
                for element_2 in breakCombos_surface:
                    connected_1 = areConnected(adsorbFile1, element_1)
                    connected_2 = areConnected(slabFile, element_2)
                    if (connected_1 == True and connected_2 == True):
                        numOfAds2Atoms = 0
                        twoPairsSet = generate2IsomerPairs(element_1, element_2,\
                                numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms, reactionType)

                        for kk in range(0, len(twoPairsSet)):
                            # write initial### and ISOMERS file
                            folder = "se_gsm_cals_2/" + file.split(".")[0] + "/" +\
                                    str(i).zfill(4) + "/scratch/"
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                            fh.write("NEW\n")
                            firstNum_1 = element_1[0] + numOfSlabAtoms
                            secondNum_1 = element_1[1] + numOfSlabAtoms
                            firstNum_2 = element_2[0]
                            secondNum_2 = element_2[1]
                            fh.write("BREAK " + str(firstNum_1) + "  " + str(secondNum_1) + "\n")
                            fh.write("BREAK " + str(firstNum_2) + "  " + str(secondNum_2) + "\n")
                            fh.write("ADD   " + str(twoPairsSet[kk][0][0]) + "  " + str(twoPairsSet[kk][0][1]) + "\n")
                            fh.write("ADD   " + str(twoPairsSet[kk][1][0]) + "  " + str(twoPairsSet[kk][1][1]) + "\n")
                            fh.close()

                            index = i
                            createTemplateFiles(file, folder, index, slab, slabType)
                            extension = 2
                            #submitSE_GSM(file, extension, index, cwd)

                            i += 1



        ###############################################
        # bi-molecular reaction
        ###############################################
        elif (numOfAdsorbates == 2):
            # find all combinations of reactive indices 1
            breakCombos_1 = []
            for i in range(0, len(reactiveIndices_1)):
                for j in range(i+1, len(reactiveIndices_1)):
                    if ( (reactiveIndices_1[i] != 0) and (reactiveIndices_1[j] != 0) ):
                        temp = []
                        temp.append(reactiveIndices_1[i])
                        temp.append(reactiveIndices_1[j])
                        breakCombos_1.append(temp)
            # find all combinations of reactive indices 2
            breakCombos_2 = []
            for i in range(0, len(reactiveIndices_2)):
                for j in range(i+1, len(reactiveIndices_2)):
                    if ( (reactiveIndices_2[i] != 0) and (reactiveIndices_2[j] != 0) ):
                        temp = []
                        temp.append(reactiveIndices_2[i])
                        temp.append(reactiveIndices_2[j])
                        breakCombos_2.append(temp)

            # for addMove == 1 and breakMove == 1
            i = 1
            for element_1 in breakCombos_1:
                connected_1 = areConnected(adsorbFile1, element_1)
                if (connected_1 == True):
                    pairsSet = generateIsomerPair(element_1, reactiveIndices_2,\
                            numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms)

                    for kk in range(0, len(pairsSet)):
                        # check coordination number
                        coordNum = getCoordinationNum(adsorbFile2, (pairsSet[kk][1]-\
                                numOfSlabAtoms - numOfAds1Atoms))
                        atomicNumber = getAtomicNumber(adsorbFile2, (pairsSet[kk][1]-\
                                numOfSlabAtoms - numOfAds1Atoms))
                        #CoordinationNumbers()
                        MAX_COORDINATION_NUM = CoordinationNumbers().getMaxCoordNum(atomicNumber)
                        if (coordNum < MAX_COORDINATION_NUM):
                            # write initial### and ISOMERS file
                            folder = "se_gsm_cals_2/" + file.split(".")[0] + "/" +\
                                    str(i).zfill(4) + "/scratch/"
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                            fh.write("NEW\n")
                            firstNum_1 = element_1[0] + numOfSlabAtoms
                            secondNum_1 = element_1[1] + numOfSlabAtoms
                            fh.write("BREAK " + str(firstNum_1) + "  " + str(secondNum_1) + "\n")
                            fh.write("ADD   " + str(pairsSet[kk][0]) + "  " + str(pairsSet[kk][1]) + "\n")
                            fh.close()

                            index = i
                            createTemplateFiles(file, folder, index, slab, slabType)
                            extension = 2
                            submitSE_GSM(file, extension, index, cwd)
                            
                            i += 1

            for element_2 in breakCombos_2:
                connected_2 = areConnected(adsorbFile2, element_2)
                if (connected_2 == True):
                    pairsSet = generateIsomerPair(element_2, reactiveIndices_1,\
                            numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms)
                    for kk in range(0, len(pairsSet)):
                        # check coordination number
                        coordNum = getCoordinationNum(adsorbFile1, (pairsSet[kk][0] -\
                                numOfSlabAtoms))
                        atomicNumber = getAtomicNumber(adsorbFile1, (pairsSet[kk][0]-\
                                numOfSlabAtoms))
                        #CoordinationNumbers()
                        MAX_COORDINATION_NUM = CoordinationNumbers().getMaxCoordNum(atomicNumber)
                        if (coordNum < MAX_COORDINATION_NUM):
                            # write initial### and ISOMERS file
                            folder = "se_gsm_cals_2/" + file.split(".")[0] + "/" +\
                                    str(i).zfill(4) + "/scratch/"
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                            fh.write("NEW\n")
                            firstNum_2 = element_2[0] + numOfSlabAtoms + numOfAds1Atoms
                            secondNum_2 = element_2[1] + numOfSlabAtoms + numOfAds1Atoms
                            fh.write("BREAK " + str(firstNum_2) + "  " + str(secondNum_2) + "\n")
                            fh.write("ADD   " + str(pairsSet[kk][1]) + "  " + str(pairsSet[kk][0]) + "\n")
                            fh.close()

                            index = i
                            createTemplateFiles(file, folder, index, slab, slabType)
                            extension = 2
                            submitSE_GSM(file, extension, index, cwd)

                            i += 1


            # for addMove == 2 and breakMove == 2
            if (addMoves != 2):
                if (breakMoves != 2):
                    print ("ERROR: You need at least 2 add and break moves for ligand exchange.")
                    exit(-1)
            if (len(reactiveIndices_1) < 2):
                if (len(reactiveIndices_2) < 2):
                    print ("ERROR: You need at least 2 reactive atoms on each adsorbate for ligand exchange.")
                    exit(-1)

            # if (addMoves == 2 and breakMoves == 2):
            for element_1 in breakCombos_1:
                for element_2 in breakCombos_2:
                    connected_1 = areConnected(adsorbFile1, element_1)
                    connected_2 = areConnected(adsorbFile2, element_2)
                    if (connected_1 == True and connected_2 == True):
                        #pairsSet = generateIsomerPair(element_1, reactiveIndices_2,\
                        #        numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms)
                        twoPairsSet = generate2IsomerPairs(element_1, element_2,\
                                numOfSlabAtoms, numOfAds1Atoms, numOfAds2Atoms, reactionType)

                        for kk in range(0, len(twoPairsSet)):
                            # check coordination number
                            # Do we need this here? don't think so.
                            '''###
                            coordNum = getCoordinationNum(adsorbFile2, (pairsSet[kk][1]-\
                                    numOfSlabAtoms - numOfAds1Atoms))
                            atomicNumber = getAtomicNumber(adsorbFile2, (pairsSet[kk][1]-\
                                    numOfSlabAtoms - numOfAds1Atoms))
                            #CoordinationNumbers()
                            MAX_COORDINATION_NUM = CoordinationNumbers().getMaxCoordNum(atomicNumber)
                            if (coordNum < MAX_COORDINATION_NUM):
                            '''###
                            #if (True):
                            # write initial### and ISOMERS file
                            folder = "se_gsm_cals_2/" + file.split(".")[0] + "/" +\
                                    str(i).zfill(4) + "/scratch/"
                            if not os.path.exists(folder):
                                os.makedirs(folder)
                            fh = open(folder + "ISOMERS"+str(i).zfill(4), 'w')
                            fh.write("NEW\n")
                            firstNum_1 = element_1[0] + numOfSlabAtoms
                            secondNum_1 = element_1[1] + numOfSlabAtoms
                            firstNum_2 = element_2[0] + numOfSlabAtoms + numOfAds1Atoms
                            secondNum_2 = element_2[1] + numOfSlabAtoms + numOfAds1Atoms
                            fh.write("BREAK " + str(firstNum_1) + "  " + str(secondNum_1) + "\n")
                            fh.write("BREAK " + str(firstNum_2) + "  " + str(secondNum_2) + "\n")
                            fh.write("ADD   " + str(twoPairsSet[kk][0][0]) + "  " + str(twoPairsSet[kk][0][1]) + "\n")
                            fh.write("ADD   " + str(twoPairsSet[kk][1][0]) + "  " + str(twoPairsSet[kk][1][1]) + "\n")
                            fh.close()

                            index = i
                            createTemplateFiles(file, folder, index, slab, slabType)
                            extension = 2
                            submitSE_GSM(file, extension, index, cwd)

                            i += 1



if __name__ == "__main__":
    main()

