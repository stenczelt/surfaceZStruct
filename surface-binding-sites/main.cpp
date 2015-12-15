#include "BindingSiteClass.h"
#include "SurfaceClass.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cstring>

bool parseVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms);

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

    parseVector(xyzFile, xyz, numOfAtoms);
    SurfaceClass aSurface;
    if (aSurface.setSurfaceType(argv[1]))
    {
        aSurface.setAtoms(numOfAtoms, xyz);
    }
    delete [] xyz;
    return (0);
}

bool parseVector(std::vector<std::string> inVec, double* inArr, int numOfAtoms)
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
    return (success);
}
