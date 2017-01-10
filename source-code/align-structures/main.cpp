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
#include <iostream>
#include <fstream>
#include <stdio.h>

// slab, atom1, .. : Surface
// adsorbate1, atom2, .. : Adsorbate-1

bool populateArrayFromVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms, 
                 std::string* inSymbols);
//bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray);
//bool readFromFile(std::string inFileName, vector<std::string>& xyzFileSlab);
//bool readFromFile(std::string inFileName);
int readFromFile(std::string inFileName, int &numOfAdsorbates, int* slabIndices, double* radius,
                  std::string* adsorbateFiles, int* adsorbateIndices, int* reactiveIndex1,
                  int* reactiveIndex2, int &numOfAdd, int &numOfBreak);

//bool readSlabFile(std::string slabFileName, vector<std::string>& xyzFileSlab)
bool readSlabFile(std::string &slabFileName);

int main(int argc, char* argv[])
{
    //if (argc < 3 || argc > 4)


    if (argc != 1)
    {   
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " " << std::endl;
        return (1);
    }
    // reading input files
    else
    {
        std::cout << "Reading input files" << std::endl;
        int numOfAdsorbates = 0;
        int slabIndices[2] = {};
        double radius[2] = {};
        std::string adsorbateFiles[2] = {};
        int adsorbateIndices[2] = {};
        int reactiveIndex1[4] = {0, 0, 0, 0}; //maximum of 4 reactive atoms on each atom
        int reactiveIndex2[4] = {0, 0, 0, 0};
        int numOfAdd = 0;
        int numOfBreak = 0;

        int returnVal = 0;
        returnVal = readFromFile("INPUT", numOfAdsorbates, slabIndices, radius, adsorbateFiles, 
                     adsorbateIndices, reactiveIndex1, reactiveIndex2,
                     numOfAdd, numOfBreak);
        int addArray[4] = {};
        if (returnVal == 0)
        {
            // create objects from input structures
            ICoord slab, adsorbate1, adsorbate2;
            std::vector<ICoord> adsorbates;

            slab.init("bindingSites.xyz");
            adsorbate1.init(adsorbateFiles[0]);
            adsorbates.push_back(adsorbate1);
            // list of adsorbate and slab indices
            addArray[0] = slabIndices[0] -1 ;
            addArray[1] = adsorbateIndices[0] + slab.natoms - 1;
            if (numOfAdsorbates == 2)
            {
                adsorbate2.init(adsorbateFiles[1]);
                adsorbates.push_back(adsorbate2);
                //addArray[0] = slabIndices[0];
                //addArray[1] = adsorbateIndices[0];
                addArray[2] = slabIndices[1] - 1;
                addArray[3] = adsorbateIndices[1] + slab.natoms + adsorbate1.natoms - 1;
                std::cout << "*********** " << addArray[0] << " " <<addArray[1] << "    " << addArray[2] << "    "
                    << addArray[3] << std::endl;
            }

            Align totalSystem(slab, adsorbates);

            // read parameters from INPUT file
            //std::string orientationIn = "horiz";
            // numOfAdsorbates = numOfAdd passed to add_align
            totalSystem.add_align(numOfAdsorbates, addArray, radius); /*orientationIn//default is horiz*/
            std::cout << "\n***************************************\n";
            std::cout << "\nOutput is written to output-*.xyz file\n";
            std::cout << "\n***************************************\n";
        }
    }
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

int readFromFile(std::string inFileName, int &numOfAdsorbates, int* slabIndices, double* radius,
                  std::string* adsorbateFiles, int* adsorbateIndices, int* reactiveIndex1,
                  int* reactiveIndex2, int &numOfAdd, int &numOfBreak)
{
    // populate vector by input file lines
    std::vector<std::string> inputFile;
    std::string newLine;
    ifstream inFile;
    inFile.open(inFileName);
    while (std::getline(inFile, newLine))
    {
        inputFile.push_back(newLine);
        //std::cout << newLine << std::endl;
    }

    // start parsing input file
    std::stringstream ss(*inputFile.begin());
    std::string unwanted = "";
    std::string temp[2] = {};
    
    std::string findSites = "";
    ss >> unwanted >> findSites >> unwanted;
    std::string slabFileName = "";

    if (std::stoi(findSites) == 1)
    {
        // stage 1: find binding sites
        std::cout << "** Finding binding sites" << std::endl;
        ss.str(*(inputFile.begin()+1));
        ss >> unwanted >> slabFileName >> unwanted;
        readSlabFile(slabFileName);
        return 1;
    }
    else if (std::stoi(findSites) == 0)
    {
        // stage 2: add adsorbates to binding sites
        std::cout << "** reading input file information" << std::endl;
        std::cout << "** WARNING: slab is read from bindingSites.xyz" << std::endl;
        // line 3 of INPUT file
        ss.clear();
        ss.str(*(inputFile.begin()+2));
        ss >> unwanted >> temp[0] >> unwanted;
        numOfAdsorbates = std::stoi(temp[0]);
        // line 4
        ss.str(*(inputFile.begin()+3));
        ss >> unwanted >> temp[0] >> unwanted;
        slabIndices[0] = std::stoi(temp[0]);
        // line 5
        ss.str(*(inputFile.begin()+4));
        ss >> unwanted >> temp[0] >> unwanted;
        radius[0] = std::stod(temp[0]);
        // line 6
        ss.str(*(inputFile.begin()+5));
        ss >> unwanted >> adsorbateFiles[0] >> unwanted;
        // line 7
        ss.str(*(inputFile.begin()+6));
        ss >> unwanted >> temp[0] >> unwanted;
        adsorbateIndices[0] = std::stoi(temp[0]); // TODO + num of slab atoms
        // line 8
        ss.str(*(inputFile.begin()+7));
        ss >> unwanted >> temp[0] >> unwanted;
        reactiveIndex1[0] = std::stoi(temp[0]);

        if (numOfAdsorbates == 2)
        {
            // lines 9-14
            ss.str(*(inputFile.begin()+9));
            ss >> unwanted >> temp[0] >> unwanted;
            slabIndices[1] = std::stoi(temp[0]);
            ss.str(*(inputFile.begin()+10));
            ss >> unwanted >> temp[0] >> unwanted;
            radius[1] = std::stod(temp[0]);
            ss.str(*(inputFile.begin()+11));
            ss >> unwanted >> adsorbateFiles[1] >> unwanted;
            ss.str(*(inputFile.begin()+12));
            ss >> unwanted >> temp[0] >> unwanted;
            adsorbateIndices[1] = std::stoi(temp[0]); // TODO + num of slab+adsorbate1 atoms
            ss.str(*(inputFile.begin()+13));
            ss >> unwanted >> temp[0] >> temp[1] >> unwanted; //TODO
            reactiveIndex2[0] = std::stoi(temp[0]);
            reactiveIndex2[1] = std::stoi(temp[1]);
        }

        ss.str(*(inputFile.begin()+15));
        ss >> unwanted >> temp[0] >> unwanted;
        numOfAdd = std::stoi(temp[0]);
        ss.str(*(inputFile.begin()+16));
        ss >> unwanted >> temp[0] >> unwanted;
        numOfBreak = std::stoi(temp[0]);

        std::cout << numOfAdsorbates << "   " << slabIndices[0] << "    " << slabIndices[1] 
            << "    " << radius[0] << "     " << radius[1] << "     " << adsorbateFiles[0]
            << "    " << adsorbateFiles[1] << " " << adsorbateIndices[0] << "   " <<
            adsorbateIndices[1] << "    " << numOfAdd
            << "    " << numOfBreak << "   " << reactiveIndex1[0] << "     " << reactiveIndex2[0] << std::endl;
    return 0;
    }
}

bool readSlabFile(std::string &slabFileName)
{
    bool success = false;
    std::vector<std::string> xyzFileSlab; // vector to store the input file
    std::ifstream inFile;
    inFile.open(slabFileName);
    if (!inFile.is_open())
    {
        std::cout << "ERROR: File " << slabFileName << " does not exist!" << std::endl;
    }
    std::string newLine;
    while (std::getline(inFile, newLine))
    {   
        xyzFileSlab.push_back(newLine);
    }   

    int numOfSlabAtoms = std::stoi(*(xyzFileSlab.begin()));
    std::string surfaceType = xyzFileSlab[1];
    if (xyzFileSlab[1].empty())
    {   
        std::cout << "ERROR: Set surface type in the input file" << std::endl;
    }   

    int slabSize = 3 * numOfSlabAtoms;
    double* slabCartesianCoords = new double[slabSize];
    std::string* slabAtomicSymbols = new std::string[numOfSlabAtoms];

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
        success = true;
    }

    delete [] slabCartesianCoords;
    delete [] slabAtomicSymbols;

    return success;
}





