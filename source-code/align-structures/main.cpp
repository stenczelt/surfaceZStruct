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

bool populateVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms, 
                 std::string* inSymbols);

int main(int argc, char* argv[])
{
    //std::ifstream inputFile1(argv[1]);
    //std::ifstream inputFile2(argv[2]);
    //std::string drivingCoord = argv[3];
    //std::string inputFile1 = argv[1];
    //std::string inputFile2 = argv[2];
    
    
    // Input cartesians of slab and adsorbate: Done, in command line
    // find surface sites
    // pass new structures to align code
    // Find binding sites on the slab and return coordinates for one of them
    // Align and add adsorbate to that binding site
    

    std::vector<std::string> xyzFileSlab; // vector to store the input file
    //std::vector<std::string> xyzFileAdsorbate;
    std::vector<std::string> xyzFileSites;
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

    int sizeSlab = 3*numOfSlabAtoms;
    //int sizeAdsorbate = 3*numOfAdsorbateAtoms;
    double* slabCartesianCoords = new double[sizeSlab];
    //double* adsorbateCartesianCoords = new double[sizeAdsorbate];
    std::string* slabAtomicSymbols = new std::string[numOfSlabAtoms];
    //std::string* adsorbateAtomicSymbols = new std::string[numOfAdsorbateAtoms];

    populateVector(xyzFileSlab, slabCartesianCoords, numOfSlabAtoms, slabAtomicSymbols);
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

        /*std::ifstream inFile2 (outFName);
        std::string newLine2;
        while (std::getline(inFile2, newLine2))
        {   
            xyzFileSites.push_back(newLine2);
        } */  

        /*int sizeSlabSites = 3*numOfSlabAtoms;
        double* slabCartesianCoords = new double[sizeSlab];
        std::string* slabAtomicSymbols = new std::string[numOfSlabAtoms];
        populateVector(xyzFileSlab, slabCartesianCoords, numOfSlabAtoms, slabAtomicSymbols);*/
    } 

    ICoord ic1, ic2;
    ic1.init(outFName);
    ic2.init(argv[2]);

    Align alignObj;
    alignObj.init(ic1.natoms,ic1.anames,ic1.anumbers,ic1.coords,ic2.natoms,ic2.anames,ic2.anumbers,ic2.coords);
    std::cout << "Before alignment: \n";
//    alignObj.print_xyz();

    /*
    int siteNumber = 0;
    std::cout << "Enter the binding site number:" << std::endl;
    std::cin >> siteNumber;

    double xSite = ic1.coords[3*siteNumber];
    double ySite = ic1.coords[3*siteNumber+1];
    double zSite = ic1.coords[3*siteNumber+2];
    std::cout << "z site: " << zSite << "\n";
    */
    int nadd = 1;
    int add[2] = {29, 60};
    alignObj.add_align(nadd, add);
    std::cout << "After alignment: \n\n"; 
    alignObj.print_xyz();
    

    delete [] slabCartesianCoords;
    delete [] slabAtomicSymbols;
    return (0);
}

bool populateVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms,
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


