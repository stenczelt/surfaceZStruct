# Mina Jafari
# 03-10-2017
# This script compares optimized structures and copies unique ones into
# "unique-structures" folder. It also sorts them based on energy.

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
import math
import operator


# function to calculate root mean square deviation of two chemical structures
def calculateRMSD(structure_1, structure_2):
    try:
        fh_1 = open(structure_1)
        lines_1 = fh_1.readlines()
        fh_1.close()
        fh_2 = open(structure_2)
        lines_2 = fh_2.readlines()
        fh_2.close()

        # assert number of atoms are the same
        assert(int(lines_1[0]) == int(lines_2[0]))
        # assert atom names are the same
        for i in range(2, int(lines_1[0])+2):
            assert(lines_1[i].split()[0] == lines_2[i].split()[0])


        sum = 0.0
        for i in range(2, len(lines_1)):
            x_1, y_1, z_1 = lines_1[i].split()[1:]
            x_2, y_2, z_2 = lines_2[i].split()[1:]
            x_1 = float(x_1)
            x_2 = float(x_2)
            y_1 = float(y_1)
            y_2 = float(y_2)
            z_1 = float(z_1)
            z_2 = float(z_2)

            sum += (x_1 - x_2)**2 + (y_1 - y_2)**2 + (z_1 - z_2)**2

        rmsd = 0.0
        rmsd = math.sqrt( sum / (len(lines_1) - 2) )
        return rmsd

    except IOError:
        print ("ERROR: File ", structure_1, " or ", structure_2, " does not exist")
        return -1

def sortEnergy(listOfStructuresPath):
    dictionary = {}

    for elem in listOfStructuresPath:
        try:
            temp = elem.split('/')
            folderPath = temp[0] + '/' + temp[1]
            fileName = folderPath + "/OUTCAR"
            # read OUTCAR file in each folder
            fh = open(fileName, 'r')
            lines = fh.readlines()
            fh.close()
            finished = False
            for line in reversed(lines):
                if "Total CPU time used" in line:
                    finished = True
                if (finished == True):
                    if ("free energy    TOTEN" in line):
                        energy = float( line.split('=')[1].split('e')[0].strip() )
                        dictionary.update( {elem : energy} )
                        break
        except IOError:
            print("ERROR: No OUTCAR file in folder ", elem)
            continue

    # sort dictionary based on value/Energy and return as a sorted list
    sorted_E = sorted(dictionary.items(), key=operator.itemgetter(1))
    return sorted_E



def main():
    # list files in optimization directory
    folderPath = "optimization"
    filePaths = []
    for elem in os.listdir(folderPath):
        filePaths.append("optimization/" + elem + "/" + elem + ".xyz")

    # threshold for comparing similarity of structures
    THRESHOLD = 0.01
    # list of unique structures
    uniqueFileNames = []
    loopRange_1 = range(0, len(filePaths))

    # identify unique files based on RMSD calculation
    for i in loopRange_1:
        for j in range(i+1, len(filePaths)):
            rmsd = calculateRMSD(filePaths[i], filePaths[j])
            if (rmsd >= THRESHOLD):
                # two strucures are not identical
                uniqueFileNames.append(filePaths[i])
                uniqueFileNames.append(filePaths[j])
            elif (rmsd <= THRESHOLD and rmsd > 0.0):
                # two strucures are identical
                if j in loopRange_1:
                    loopRange_1.remove(j)
            elif (rmsd == -1):
                print ("ERROR: File does not exist!")
                exit(-1)
    uniqueFileNames = list( set(uniqueFileNames) )

    # sort uniqueFileNames based on energy and pick top 5 lowest energy structures
    sorted_files = sortEnergy(uniqueFileNames)

    # overwrite folder if exists and create it if does not exist
    dst = "unique_structures"
    if os.path.exists(dst):
        shutil.rmtree(dst)
        os.makedirs(dst)
    if not os.path.exists("unique_structures"):
        os.makedirs(dst)

    numOfFiles = len(sorted_files)
    print sorted_files[0][0]
    if (numOfFiles > 5):
        # only use 5 structures with lowest energy
        for i in range(0, 5):
            shutil.copy(sorted_files[i][0], dst)
    else:
        # use all structures
        for i in range(0, numOfFiles):
            shutil.copy(sorted_files[i][0], dst)




if __name__ == "__main__":
    main()
