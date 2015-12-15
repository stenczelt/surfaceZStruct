// Mina Jafari
// 12-07-2015

#include "Termination100.h"
#include "Termination111.h"
#include "fcc110.h"
#include "bcc110.h"
#include "bcc111.h"
//#include "diamond100.h"
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
//    {"diamond100", required_argument, nullptr, 'i'},
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
    while ((c = getopt_long(argc, argv, "a:b:c:d:e:f:g:h", long_options, &index)) != -1)
    {
        switch (c)
        {
            case 'a': // fcc100
            case 'd': // bcc100
                {
                    Termination100 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findHollow(1, 1);
                    aSlab.findAtop(1, 1);
                    aSlab.findBridge(1, 1);
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
                    aSlab.findHollow(1, 1);
                    aSlab.findAtop(1, 1);
                    aSlab.findLongBridge(1, 1);
                    aSlab.findShortBridge(1, 1);
                    break;
                }

            case 'c': // fcc111
            case 'g': // hcp0001
                {
                    Termination111 aSlab;
                    std::ifstream inFile (optarg);
                    std::string newLine;
                    while (std::getline(inFile, newLine))
                    {
                        xyzFile.push_back(newLine);
                    }
                    aSlab.setAtoms(xyzFile);
                    aSlab.findAtop(1, 1);
                    aSlab.findBridge(1, 1);
                    aSlab.findFcc(1, 1);
                    aSlab.findHcp(1, 1);
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
                    aSlab.findHollow(1, 1);
                    aSlab.findAtop(1, 1);
                    aSlab.findLongBridge(1, 1);
                    aSlab.findShortBridge(1, 1);
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
                    aSlab.findAtop(1, 1);
                    aSlab.findFcc(1, 1);
                    aSlab.findHcp(1, 1);
                    aSlab.findShallowBridge(1, 1);
                    aSlab.findDeepBridge(1, 1);
                    break;
                }

/*            case 'i': // diamond100
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
*/
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
