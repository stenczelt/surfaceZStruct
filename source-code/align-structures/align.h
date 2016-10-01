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
        // Atomic symbols of the fragments
        string* atomicNames1;
        string* atomicNames2;
        // Atomic numbers of the fragments
        int* anumbers1;
        int* anumbers2;
        // Cartesian coordinates of the fragments
        double* xyz1s;
        double* xyz2s;
        // ???
        double* xyz1;
        double* xyz2;
        // ???
        int numOfAtoms3;
        string* atomicNames3;
        int* anumbers3;
        double* xyz3;
        // set of angles being sampled during geometry alignment
        // TODO: for now max. is 4 angles
        int angleSet[4] = {};

        // checks to see if the atoms to be added are in the same fragment
        int check_frag(int atom1, int atom2);
        int check_frag_3(int atom1, int atom2);
        // get number of atoms bonded to atom 1 (atom1)
        int get_bonds(int atom1, ICoord ic1, int* bonded);
        void align_to_z(int numOfAtoms, int t1, int t2, double* xyz, string* atomicNames,
                        int sign);
        //void align_to_z(int numOfAtoms, int t1, int t2, double* xyz, string* atomicNames,
        //                int sign, double offset);
        // rotates fragment around X axis
        //void rotate_around_z(int numOfAtoms1, int numOfAtoms2, double torv, double* xyz);
        void rotate_around_z(int numOfAtoms2, double torv, double* xyz);
        void linear_right(double* v1, int atom1, int* bonded, double* xyz, std::string orientationIn);
        void planar_cross(double* v1, int atom1, int* bonded, double* xyz);
        void align_v1(int nvf, double* v1);
        // not implemented!
        void point_out(double* v1, int numOfAtoms, double* xyz);
        void vdw_vector_opt(int numOfAtoms1, int numOfAtoms2, double* v1, ICoord icp);
        //void vdw_vector_opt(int numOfAtoms1, int numOfAtoms2, double* v1, ICoord icp, int atom1, int atom2);

    public:
        int inited;
        double* xyzAugmented;
        double* xyzAugmented3;
        int numOfAtoms1;
        int numOfAtoms2;
        int avec1; //final rotation alignment, not implemented
        int avec2;

        void init(int numOfAtoms1i, string* atomicNames1i, int* anumbers1i, double* xyz1i, 
                  int numOfAtoms2i, string* atomicNames2i, int* anumbers2i, double* xyz2i);
        void add_third(int numOfAtoms3i, string* atomicNames3i, int* anumbers3i, double* xyz3i);
        void align_zero();
        //void add_align(int nadd1, int* add1, double zBindingSite);
        void add_align(int nadd1, int* add1, int atomIndex2, vector<double> angleSet, std::string orientationIn);
        int add_align_v(int nadd1, int* add1, int wtm, double* aprv);
        void shuttle_align(int nadd1, int* add1);

        void print_xyz();
        void print_xyz_3();

        void freemem();
        double norm(double* x, int size);
        //void get_rotation_matrix(double** rotMat, double* thetas);
        void print_xyz_gen(int natoms, string* anames, double* coords);
        bool writeToFile(std::string &outFile);
        //void moveToOrigin(int nadd1, int* add1, int atomIndex2);
        void moveToOrigin(int atom2, int atomIndex2, int i, int nadd1);
        void moveToBindingSite(int atom1, int atom2, int atomIndex2, int i, int nadd1);
        void applyRotationMatrix(int numOfAtoms, double* xyz, double** rotationMatrix);
};

#endif

