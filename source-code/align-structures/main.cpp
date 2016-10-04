// Mina Jafari, August 2016
#include "align.h"
#include "icoord.h"
#include "utils.h"
#include "BindingSiteClass.h"
#include "SurfaceClass.h"
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>


// ic1, atom1, .. : Surface
// ic2, atom2, .. : Adsorbate

bool populateArrayFromVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms, 
                 std::string* inSymbols);
bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray,
                  int &secondStructureIndex, std::vector<double> &angles);

int main(int argc, char* argv[])
{
    std::vector<std::string> xyzFileSlab; // vector to store the input file
    if (argc != 3)
    {   
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " <slab file> " << "adsorbate file" << std::endl;
        return (1);
    }   
    else // reading the input file
    std::cout << "Reading input files" << std::endl;
    {   
        std::ifstream inFile (argv[1]);
        std::string newLine;
        while (std::getline(inFile, newLine))
        {   
            xyzFileSlab.push_back(newLine);
        }   
    }   
    int numOfSlabAtoms = std::stoi(*(xyzFileSlab.begin()));
    //int numOfAdsorbateAtoms = std::stoi(*(xyzFileAdsorbate.begin()));
    std::string surfaceType = xyzFileSlab[1];
    if (xyzFileSlab[1].empty())
    {   
        std::cout << "ERROR: Set surface type in the input file" << std::endl;
    }   

    int slabSize = 3 * numOfSlabAtoms;
    //int sizeAdsorbate = 3*numOfAdsorbateAtoms;
    double* slabCartesianCoords = new double[slabSize];
    //double* adsorbateCartesianCoords = new double[sizeAdsorbate];
    std::string* slabAtomicSymbols = new std::string[numOfSlabAtoms];
    //std::string* adsorbateAtomicSymbols = new std::string[numOfAdsorbateAtoms];

    populateArrayFromVector(xyzFileSlab, slabCartesianCoords, numOfSlabAtoms, slabAtomicSymbols);
    SurfaceClass aSurface;
    std::string outFName = "bindingSites.xyz";
    if (aSurface.setSurfaceType(surfaceType))
    {   
        if (!aSurface.setAtoms(numOfSlabAtoms, slabCartesianCoords, slabAtomicSymbols))
        {   
            std::cout << "ERROR: " << std::endl;
            return (2);
        }   
        aSurface.findAllSites();
        aSurface.writeToFile(outFName);
    } 

    ICoord slab, adsorbate1;
    slab.init(outFName);
    //ic1.init(argv[1]);
    adsorbate1.init(argv[2]);

    std::vector<ICoord> adsorbates;
    adsorbates.push_back(adsorbate1);
    Align totalSystem(slab, adsorbates);
//    alignObj.init(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,ic2.natoms,ic2.anames,ic2.anumbers,ic2.coords);

    int numOfAdd = 0;
    int addArray[4] = {};
    int secondStructureIndex = 0;
    std::vector<double> angleSet;
    std::string orientationIn = "horiz";
    readFromFile("INPUT", numOfAdd, addArray, secondStructureIndex, angleSet);
    totalSystem.add_align(numOfAdd, addArray, secondStructureIndex, angleSet, orientationIn/*default is horiz*/);
//    std::string outFileName = "alignedStr.xyz";
//    alignObj.writeToFile(outFileName);
    std::cout << "\n***********************************\n";
    std::cout << "\nOutput is written to aligned.xyz file\n";
    std::cout << "\n***********************************\n";
    
    delete [] slabCartesianCoords;
    delete [] slabAtomicSymbols;
    return (0);
}

bool populateArrayFromVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms,
                 std::string* inSymbols)
{
    bool success = true;
    int k = 0;
    for (auto i=inVec.begin()+2; i != inVec.end(); ++i)
    {
        std::string unwanted = "";
        double temp[3];
        std::istringstream iss(*i);
        iss >> unwanted >> temp[0] >> temp[1] >> temp[2] >> unwanted;
        for (int j=0; j<3; ++j)
        {
            inArr[3*k+j] = temp[j];
        }
        ++k;
    }
    int m = 0;
    for (auto l=inVec.begin()+2; l != inVec.end(); ++l)
    {
        std::string unwanted = "";
        std::string symbol;
        std::istringstream iss(*l);
        iss >> symbol >> unwanted >> unwanted >> unwanted >> unwanted;
        inSymbols[m] = symbol;
        ++m;
    }
    return (success);
}

bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray,
                  int &secondStructureIndex, std::vector<double> &angles)
{
    bool success = false;
    //int numOfAdd = numOfAddIn;
    //int* addArray = addArrayIn;
    //int secondStructureIndex = secondStructureIndexIn;
    std::vector<std::string> inputFile;
    std::string newLine;
    ifstream inFile;
    inFile.open(inFileName);
    //std::ofstream inFile = inFileName;
    while (std::getline(inFile, newLine))
    {
        inputFile.push_back(newLine);
    }

    std::stringstream ss(*inputFile.begin());
    std::string unwanted = "";
    std::string temp[6] = {};
    //if (ss.peek() != "#")
    ss >> unwanted >> numOfAdd;
    if (numOfAdd == 1)
    {
        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+1));
        ss >> unwanted >> unwanted >> temp[0] >> temp[1];
        addArray[0] = std::stoi(temp[0]);
        addArray[1] = std::stoi(temp[1]);
        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+2));
        ss >> unwanted >> temp[3];
        secondStructureIndex = std::stoi(temp[3]);

        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+3));
        ss >> unwanted >> temp[4];
        for(int i=0; i<std::stoi(temp[4]); i++)
        {
            ss >> temp[5];
            angles.push_back(std::stod(temp[5]));
        }
        success = true;
    }
    else if (numOfAdd == 2)
    {
        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+1));
        ss >> unwanted >> unwanted >> temp[0] >> temp[1];
        addArray[0] = std::stoi(temp[0]);
        addArray[1] = std::stoi(temp[1]);
        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+2));
        ss >> unwanted >> unwanted >> temp[2] >> temp[3];
        addArray[2] = std::stoi(temp[2]);
        addArray[3] = std::stoi(temp[3]);
        ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+3));
        ss >> unwanted >> temp[5];
        secondStructureIndex = std::stoi(temp[5]);
        success = true;
    }

    return success;
}

