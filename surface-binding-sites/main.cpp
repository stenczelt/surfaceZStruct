#include "BindingSiteClass.h"
#include "SurfaceClass.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <string>

bool parseVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms, 
                 std::string* inSymbols);

int main(int argc , char* argv[])
{
    std::vector<std::string> xyzFile; // vector to store the input file
    if (argc != 3)
    {   
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " <surface type> <input file name> " << std::endl;
        return (1);
    }
    else // reading the input file
    {
        std::ifstream inFile (argv[2]);
        std::string newLine;
        while (std::getline(inFile, newLine))
        {
            xyzFile.push_back(newLine);
        }
    }
    int numOfAtoms = std::stoi(*(xyzFile.begin()));

    int size = 3*numOfAtoms;
    double* xyz = new double[size];
    std::string* atomicSymbols = new std::string[numOfAtoms];

    parseVector(xyzFile, xyz, numOfAtoms, atomicSymbols);
    SurfaceClass aSurface;
    std::string outFName = "bindingSites.xyz";
    if (aSurface.setSurfaceType(argv[1]))
    {
        aSurface.setAtoms(numOfAtoms, xyz, atomicSymbols);

        int numOfAtomsInLayer = aSurface.getSurfaceWidth() * aSurface.getSurfaceLength(); 
        int index = 0;
        while (index > aSurface.getNumOfAtoms() || 
               index < aSurface.getNumOfAtoms() - (numOfAtomsInLayer-1))
        {
            std::cout << "Enter the atom number between " << 
                aSurface.getNumOfAtoms() - (numOfAtomsInLayer-1) << " and "
                << aSurface.getNumOfAtoms() << std::endl;
            std::cin >> index;
        }
        double radius = 0.0;
        std::cout << "Enter the radius: " << std::endl;
        std::cin >> radius;
        std::string type = " ";
        std::cout << "Enter the binding site type: " << std::endl;
        std::cin >> type;

        aSurface.findNearbySites(index, radius, type);
        aSurface.writeToFile(outFName);
    }
    delete [] xyz;
    delete [] atomicSymbols;
    return (0);
}

bool parseVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms,
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
    int m =0;
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
