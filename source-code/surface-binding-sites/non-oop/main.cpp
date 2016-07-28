// Mina Jafari
// 12-07-2015

#include "fcc100.h"
#include "fcc110.h"
#include "fcc111.h"
#include "bcc100.h"
#include "bcc110.h"
#include "bcc111.h"
#include "hcp0001.h"
#include "diamond100.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>

static struct option long_options[] =
{
    {"fcc100", required_argument, nullptr, 'a'},
    {"fcc110", required_argument, nullptr, 'b'},
    {"fcc111", required_argument, nullptr, 'c'},
    {"bcc100", required_argument, nullptr, 'd'},
    {"bcc110", required_argument, nullptr, 'e'},
    {"bcc111", required_argument, nullptr, 'f'},
    {"hcp0001", required_argument, nullptr, 'g'},
    {"diamond100", required_argument, nullptr, 'i'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0} // what is this??
};


int main(int argc, char* argv[])
{
    int index = 0;
    char c;
    std::vector<std::string> xyzFile; // vector to store the input file

    if (argc != 3)
    {
        std::cout << "ERROR - wrong number of parameters" << std::endl;
        std::cout << "Usage: " << argv[0] << " <--surface type>  <input file> " << std::endl;
        return (1);
    }
    while ((c = getopt_long(argc, argv, "a:b:c:h", long_options, &index)) != -1)
    {
        switch (c)
        {
            case 'a': // fcc100
                {
                    fcc100 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findHollow();
                    aSlab.findAtop();
                    aSlab.findBridge();
                    break;
                }

            case 'b': // fcc110
                {
                    fcc110 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findHollow();
                    aSlab.findAtop();
                    aSlab.findLongBridge();
                    aSlab.findShortBridge();
                    break;
                }

            case 'c': // fcc111
                {
                    fcc111 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findAtop();
                    aSlab.findBridge();
                    aSlab.findFcc();
                    aSlab.findHcp();
                    break;
                }

            case 'd': // bcc100
                {
                    bcc100 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findHollow();
                    aSlab.findAtop();
                    aSlab.findBridge();
                    break;
                }

            case 'e': // bcc110
                {
                    bcc110 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findHollow();
                    aSlab.findAtop();
                    aSlab.findLongBridge();
                    aSlab.findShortBridge();
                    break;
                }

            case 'f': // bcc111
                {
                    bcc111 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findAtop();
                    aSlab.findFcc();
                    aSlab.findHcp();
                    aSlab.findShallowBridge();
                    aSlab.findDeepBridge();
                    break;
                }

            case 'g': // hcp0001
                {
                    hcp0001 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findAtop();
                    aSlab.findBridge();
                    aSlab.findFcc();
                    aSlab.findHcp();
                    break;
                }

            case 'i': // diamond100
                {
                    diamond100 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findAtop();
                    aSlab.findFirstBridge();
                    aSlab.findSecondBridge();
                    aSlab.findThirdBridge();
                    break;
                }

            case 'h':
                {
                    std::cout << "Usage: " << argv[0] << " <--surface type>  <input file> " << std::endl;
                    std::cout << "Surface types: fcc100, fcc110, fcc111\n";
                    break;
                }

            case ':':
            default:
                {
                    std::cout << "ERROR - " << std::endl;
                    return (2);
                }

        } // switch
    } // while
    return (0);
} // main()
