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
#include <algorithm>
#include <fenv.h> // for NaN debugging

// slab, atom1, .. : Surface
// adsorbate1, atom2, .. : Adsorbate-1

bool populateArrayFromVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms, 
                 std::string* inSymbols);
//bool readFromFile(std::string inFileName, int &numOfAdd, int* addArray);
//bool readFromFile(std::string inFileName, vector<std::string>& xyzFileSlab);
//bool readFromFile(std::string inFileName);
int readFromFile(std::string inFileName, int &numOfAdsorbates, int* slabIndices, double* radius,
                  std::string* adsorbateFiles, int* adsorbateIndices, int* reactiveIndex1,
                  int* reactiveIndex2, int &numOfAdd, int &numOfBreak, std::string &slabFileName);

//bool readSlabFile(std::string slabFileName, vector<std::string>& xyzFileSlab)
bool readSlabFileAndWrite(std::string &slabFileName);
bool readSlabFile(std::string &slabFileName, SurfaceClass &aSurface);

int main(int argc, char* argv[])
{
    //feenableexcept(FE_INVALID | FE_OVERFLOW); // for nan debugging
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
        std::string slabFileName = "";
        int numOfAdsorbates = 0;
        int slabIndices[2] = {};
        double radius[2] = {};
        std::string adsorbateFiles[2] = {};
        int adsorbateIndices[2] = {};
        int reactiveIndex1[4] = {0, 0, 0, 0}; //maximum of 4 reactive atoms on each atom
        int reactiveIndex2[4] = {0, 0, 0, 0};
        int numOfAdd = 0;
        int numOfBreak = 0;

        int returnVal = 1;
        returnVal = readFromFile("INPUT", numOfAdsorbates, slabIndices, radius, adsorbateFiles, 
                     adsorbateIndices, reactiveIndex1, reactiveIndex2,
                     numOfAdd, numOfBreak, slabFileName);
        int addArray[4] = {};
        if (returnVal == 0)
        {
            SurfaceClass aSurface;
            readSlabFile(slabFileName, aSurface);
            aSurface.findAllSites();

          // returns all sites including the input index
            std::vector<int> allSites1 = 
                aSurface.findNearbySites(slabIndices[0], radius[0], "all");

          std::cout << " slab index 0 " << slabIndices[0] << std::endl;
          std::cout << " slab index 1 " << slabIndices[1] << std::endl;


          // removes the element that equals the input index
            //allSites1.erase(std::remove(allSites1.begin(), allSites1.end(), 
            //            slabIndices[0]-aSurface.getNumOfAtoms()-1), allSites1.end());
                //aSurface.findNearbySites(slabIndices[0] - aSurface.getNumOfAtoms(), radius[0], "all");

            // create objects from input structures
            ICoord slab, adsorbate1, adsorbate2;
            std::vector<ICoord> adsorbates;

            slab.init("bindingSites.xyz");
            adsorbate1.init(adsorbateFiles[0]);
            adsorbates.push_back(adsorbate1);
            // list of adsorbate and slab indices
            addArray[0] = slabIndices[0] - 1;
            addArray[1] = adsorbateIndices[0] + slab.natoms - 1;
            //std::vector<int> allSitesCombined = allSites1;
            std::vector<int> allSites2;

            if (numOfAdsorbates == 2)
            {
//                allSites2 = 
//                    aSurface.findNearbySites(slabIndices[1], radius[1], "all");
                // removes the element that equals the input index
//                allSites2.erase(std::remove(allSites2.begin(), allSites2.end(), 
//                            slabIndices[1]-aSurface.getNumOfAtoms()-1), allSites2.end());
                    //aSurface.findNearbySites(slabIndices[1] - aSurface.getNumOfAtoms(), radius[1], "all");
//                std::cout << "input index  " << slabIndices[1] - aSurface.getNumOfAtoms() << std::endl;
                adsorbate2.init(adsorbateFiles[1]);
                adsorbates.push_back(adsorbate2);
                //addArray[0] = slabIndices[0];
                //addArray[1] = adsorbateIndices[0];
                addArray[2] = slabIndices[1] - 1;
                addArray[3] = adsorbateIndices[1] + slab.natoms + adsorbate1.natoms - 1;
                // http://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
                //allSitesCombined.reserve( allSites1.size() + allSites2.size() );
                //allSitesCombined.insert( allSitesCombined.end(), allSites1.begin(), allSites1.end() );
                //allSitesCombined.insert( allSitesCombined.end(), allSites2.begin(), allSites2.end() );
                std::cout << "*********** " << addArray[0] << " " <<addArray[1] << "    " << addArray[2] << "    "
                    << addArray[3] << std::endl;
            }
//            std::cout << " before ~~~~~~~~~~@@@@@@@ vector " << allSites1.size() << std::endl;
//            std::cout << " before ~~~~~~~~~~@@@@@@@ vector " << allSites2.size() << std::endl;
//            for (int i=0; i<allSites1.size(); i++)
//            {
//                std::cout << allSites1[i]+aSurface.getNumOfAtoms() << "   \n";
//                std::cout << allSites2[i]+aSurface.getNumOfAtoms() << "   \n";
//            }

            Align totalSystem(slab, adsorbates);

            // read parameters from INPUT file
            //std::string orientationIn = "horiz";
            // numOfAdsorbates = numOfAdd passed to add_align
            int numOfAdds = numOfAdsorbates;
//            std::cout << "after ~~~~~~~~~~@@@@@@@ vector " << allSites1.size() << std::endl;
            //std::cout << "~~~~~~~~~~@@@@@@@ size " << allSitesCombined.size() << std::endl;
//            for (int i=0; i<allSites1.size(); i++)
//            {
                //std::cout << allSites1[i]+slab.natoms << "   ";
//                std::cout << "allSites1 " <<  allSites1[i]+aSurface.getNumOfAtoms() << "   ";
//                std::cout << "allSites2 " <<  allSites2[i]+aSurface.getNumOfAtoms() << "   ";
                //std::cout << "slab.natoms " << slab.natoms << std::endl;
//            }
//            std::cout << std::endl;
                    std::cout << "```````````*********** add array " << addArray[0] << " " <<addArray[1] 
                        << "    " << addArray[2] << "    " << addArray[3] << " " << numOfAdds << std::endl;

                    /*
            totalSystem.add_align(numOfAdds, addArray); //, allSitesCombined); //orientationIn//default is horiz
                    std::cout << "TTTSSS " << addArray[0] << std::endl;
                    std::cout << "TTTSSS " << addArray[0] << std::endl;
                    */

//            int outerLoopMax = std::min(allSites1.size(), allSites2.size());
//            int innerLoopMax = std::max(allSites1.size(), allSites2.size());
//            for (int k=0; k<outerLoopMax; k++)
//            {
//                for (int l=0; l<innerLoopMax; l++)
//                {
                    //Align totalSystem(slab, adsorbates);
//                    std::cout << "#####  k: " << k << "  l: " << l << std::endl;
                    // TODO match iterator with sites1 or 2
                    // TODO match iterator with index input by the user
//                    addArray[0] = allSites1[k] + aSurface.getNumOfAtoms(); //TODO do we need -1?
//                    addArray[2] = allSites2[l] + aSurface.getNumOfAtoms(); //TODO if numOfAdds > 1
//                    std::cout << "```````````*********** add array " << addArray[0] << " " <<addArray[1] 
//                        << "    " << addArray[2] << "    " << addArray[3] << " " << numOfAdds << std::endl;
//                    totalSystem.add_align(numOfAdds, addArray);
//                }
//            }

                for (int l = 0; l < allSites1.size() ; l++)
                {
                  //Align totalSystem(slab, adsorbates);
                    //std::cout << "#####  k: " << k << "  l: " << l << std::endl;
                  // TODO match iterator with sites1 or 2
                  // TODO match iterator with index input by the user
                    addArray[0] = allSites1[l] + aSurface.getNumOfAtoms(); 
                    std::cout << "TTTSSS " << addArray[0] << std::endl;
                    std::cout << "TTTSSS " << addArray[1] << std::endl;
                    totalSystem.add_align(numOfAdds, addArray);
                }

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
        int* reactiveIndex2, int &numOfAdd, int &numOfBreak, std::string &slabFileName)
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

    if (std::stoi(findSites) == 1)
    {
        // stage 1: find binding sites
        std::cout << "** Finding binding sites" << std::endl;
        ss.str(*(inputFile.begin()+1));
        std::string name = "";
        ss >> unwanted >> slabFileName >> unwanted;
        readSlabFileAndWrite(slabFileName);
        return 1;
    }
    else if (std::stoi(findSites) == 0)
    {
        // stage 2: add adsorbates to binding sites
        std::cout << "** reading input file information" << std::endl;
        std::cout << "** WARNING: slab is read from bindingSites.xyz" << std::endl;
        ss.clear();
        ss.str(*(inputFile.begin()+1));
        std::string name = "";
        ss >> unwanted >> slabFileName >> unwanted;
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

        std::cout << slabFileName << "--" << numOfAdsorbates << "   " << slabIndices[0] << "    " << slabIndices[1] 
            << "    " << radius[0] << "     " << radius[1] << "     " << adsorbateFiles[0]
            << "    " << adsorbateFiles[1] << " " << adsorbateIndices[0] << "   " <<
            adsorbateIndices[1] << "    " << numOfAdd
            << "    " << numOfBreak << "   " << reactiveIndex1[0] << "     " << reactiveIndex2[0] << std::endl;
        return 0;
    }
}

bool readSlabFileAndWrite(std::string &slabFileName)
{
    bool success = false;
    std::vector<std::string> xyzFileSlab; // vector to store the input file
    std::ifstream inFile;
    inFile.open(slabFileName);
    if (!inFile.is_open())
    {
        std::cout << "ERROR: File 1111 " << slabFileName << " does not exist!" << std::endl;
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

bool readSlabFile(std::string &slabFileName, SurfaceClass &aSurface)
{
    bool success = false;
    std::vector<std::string> xyzFileSlab; // vector to store the input file
    std::ifstream inFile;
    inFile.open(slabFileName);
    if (!inFile.is_open())
    {
        std::cout << "ERROR: File 222 " << slabFileName << " does not exist!" << std::endl;
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
    if (aSurface.setSurfaceType(surfaceType))
    {   
        if (!aSurface.setAtoms(numOfSlabAtoms, slabCartesianCoords, slabAtomicSymbols))
        {   
            std::cout << "ERROR: " << std::endl;
            return (2);
        }   
        //aSurface.findAllSites();
        success = true;
    }

    delete [] slabCartesianCoords;
    delete [] slabAtomicSymbols;

    return success;
}




