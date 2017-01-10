// Mina Jafari, August 2016
#ifndef ALIGN_H
#define ALIGN_H

#include "icoord.h"
//#include "stringtools.h"
#include "utils.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>

void print_xyz_gen(int numOfAtoms, string* atomicNames, double* coords);

class Align
{
    private:
        // slab defined as a ICoord object
        ICoord mSlab;
        // a vector containing adsorbate molecuels
        std::vector<ICoord> mAdsorbates;
        // keep track of adsorbate 1 or 2
        int mAdsorbateNum = 0;
        
        int mNumOfAtomsCombined = 0;
        std::string* mAtomicNamesCombined;
        int* mAtomicNumbersCombined;
        double* mCoordinatesCombined;

        // checks to see if the atoms to be added are in the same fragment
        int check_frag(int atom1, int atom2);
        //int check_frag_3(int atom1, int atom2);
        // get number of atoms bonded to atom 1 (atom1)
        int get_bonds(int atom1, ICoord ic1, int* bonded);
        void align_to_z(int numOfAtoms, int t1, int t2, double* xyz, string* atomicNames);
        //void align_to_z(int numOfAtoms, int t1, int t2, double* xyz, string* atomicNames,
        //                int sign, double offset);
        // rotates fragment around X axis
        //void rotate_around_z(int numOfAtoms1, int numOfAtoms2, double torv, double* xyz);
        void rotate_around_z(int numOfAtoms2, double torv, double* xyz, int numAdsorbate);
        void linear_right(double* v1, int atom1, int* bonded, double* xyz, std::string orientationIn);
        void planar_cross(double* v1, int atom1, int* bonded, double* xyz);
        void align_v1(int nvf, double* v1);
        // not implemented!
        //void point_out(double* v1, int numOfAtoms, double* xyz);
        void vdw_vector_opt(double* v1, ICoord icp);
        //void vdw_vector_opt(int numOfAtoms1, int numOfAtoms2, double* v1, ICoord icp, int atom1, int atom2);
        void unifyStructures();

    public:
        int inited;
        //Constructor
        Align(ICoord slab, std::vector<ICoord> adsorbates);
        //void add_third(int numOfAtoms3i, string* atomicNames3i, int* anumbers3i, double* xyz3i);
        //void align_zero();
        void add_align(int nadd1, int* add1, double* radius); /*std::string orientationIn*/
        //int add_align_v(int nadd1, int* add1, int wtm, double* aprv);
        //void shuttle_align(int nadd1, int* add1);
        void print_xyz();
        //void print_xyz_3();
        double norm(double* x, int size);
        void print_xyz_gen(int natoms, string* anames, double* coords);
        bool writeToFile(std::string &outFile);
        void moveToOrigin(int atom2);
        void moveToBindingSite(int atom1, int atom2, int numAdsorbate);
        void applyRotationMatrix(int numOfAtoms, double* xyz, double** rotationMatrix);
};

#endif

