// Mina Jafari
// 12-02-2015

#include "fcc100.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

int main(int argc, char* argv[])
{
    std::vector<std::string> xyzFile; // vector to store the input file
    if (argc != 2)
    {   
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " <input file name> " << std::endl;
        return (1);
    }
    else // reading the input file
    {
        std::ifstream inFile (argv[1]);
        std::string newLine;
        while (std::getline(inFile, newLine))
        {
            xyzFile.push_back(newLine);
        }
    }

    fcc100 aSlab;
    aSlab.setAtoms(xyzFile);
    aSlab.findHollow();
    aSlab.findAtop();
    aSlab.findBridge();

    return (0);
} // main()
