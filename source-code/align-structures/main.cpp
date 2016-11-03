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
bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray);
//                  std::vector<double> &angles);

int main(int argc, char* argv[])
{
    std::vector<std::string> xyzFileSlab; // vector to store the input file
    if (argc < 3 || argc > 4)
    {   
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " <slab file> " 
                  << "<adsorbate-1 file> " << "<adsorbate-2 file>" << std::endl;
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
    } 
    // create objects from input structures
    ICoord slab, adsorbate1, adsorbate2;
    std::vector<ICoord> adsorbates;
    std::cout << "iIn main: HEREEEEEEEEE  11\n";
    if (argc == 3)
    {
        slab.init(outFName);
        adsorbate1.init(argv[2]);
        adsorbates.push_back(adsorbate1);
    }
    else if (argc == 4)
    {
        slab.init(outFName);
        adsorbate1.init(argv[2]);
        adsorbate2.init(argv[3]);
        adsorbates.push_back(adsorbate1);
        adsorbates.push_back(adsorbate2);
    }
    std::cout << "iIn main: HEREEEEEEEEE\n";
    Align totalSystem(slab, adsorbates);

    // read parameters from INPUT file
    int numOfAdd = 0;
    int addArray[4] = {};
    //std::vector<double> angleSet;
    std::string orientationIn = "horiz";

    //readFromFile("INPUT", numOfAdd, addArray, angleSet);
    readFromFile("INPUT", numOfAdd, addArray);
    //totalSystem.add_align(numOfAdd, addArray, angleSet, orientationIn/*default is horiz*/);
    totalSystem.add_align(numOfAdd, addArray, orientationIn/*default is horiz*/);
    std::cout << "\n***************************************\n";
    std::cout << "\nOutput is written to aligned-*.xyz file\n";
    std::cout << "\n***************************************\n";
    
    delete [] slabCartesianCoords;
    delete [] slabAtomicSymbols;
    //slab.freemem();
    //adsorbate1.freemem();
    //adsorbate2.freemem();
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

//bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray,
//                  std::vector<double> &angles)
bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray)
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
    std::string temp[2] = {};
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

        /*ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+2));
        ss >> unwanted >> temp[0];
        for(int i=0; i<std::stoi(temp[0]); i++)
        {
            ss >> temp[1];
            angles.push_back(std::stod(temp[1]));
        }
        success = true;*/
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
        ss >> unwanted >> unwanted >> temp[0] >> temp[1];
        addArray[2] = std::stoi(temp[0]);
        addArray[3] = std::stoi(temp[1]);

        /*ss.str("");
        ss.clear();
        ss.str(*(inputFile.begin()+3));
        ss >> unwanted >> temp[0];
        for(int i=0; i<std::stoi(temp[0]); i++)
        {
            ss >> temp[1];
            angles.push_back(std::stod(temp[1]));
        }*/
        success = true;
    }

    return success;
}

