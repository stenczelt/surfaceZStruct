// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //
// Mina Jafari, August 2016


#include <math.h>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include "align.h"
using namespace std;
using namespace Eigen;

//TODO: implement point_out function for planes

int Align::check_frag(int atom1, int atom2)
{
    int f1, f2;
    if (atom1 < numOfAtoms1) f1 = 1;
    else f1 = 2;
    if (atom2 < numOfAtoms1) f2 = 1;
    else f2 = 2;

    int same_frag = 1;
    if (f1!=f2) same_frag = 0;

    return same_frag;
}

/*int Align::check_frag_3(int atom1, int atom2)
{
    int f1,f2;
    if (atom1<numOfAtoms1+numOfAtoms2) f1 = 1;
    else f1 = 2;
    if (atom2<numOfAtoms1+numOfAtoms2) f2 = 1;
    else f2 = 2;

    int same_frag = 1;
    if (f1!=f2) same_frag = 0;

    return same_frag;
}*/

int Align::get_bonds(int atom1, ICoord ic1, int* bonded)
{
    int nfound = 0;

    for (int i=0;i<ic1.natoms;i++)
        if (i!=atom1 && ic1.bond_exists(i,atom1))
            bonded[nfound++] = i;
#if 0
    for (int i=0;i<nfound;i++)
        printf(" %2i bonded to %2i \n",atom1+1,bonded[i]+1);
#endif
    return nfound;
}

void Align::applyRotationMatrix(int numOfAtoms, double* xyz, double** rotationMatrix)
{
    double* xyzTemp = new double[3*numOfAtoms];
    for (int i=0;i<3*numOfAtoms;i++) xyzTemp[i] = 0.;
    for (int i=0;i<numOfAtoms;i++)
        for (int j=0;j<3;j++)
            for (int k=0;k<3;k++)
                xyzTemp[3*i+j] += rotationMatrix[j][k]*xyz[3*i+k];
    for (int i=0;i<3*numOfAtoms;i++)
        xyz[i] = xyzTemp[i];

    delete [] xyzTemp;
}

void Align::align_to_z(int numOfAtoms, int t1, int t2, double* xyz, string* atomicNames, int sign)
{
    // xyz1 temporary Cartesian coordinates of vectors pointing from 
    // the arbitrary origin to each atom
    double* xyz1 = new double[3*numOfAtoms];
    for (int i=0;i<numOfAtoms;i++)
    {
        xyz1[3*i+0] = xyz[3*i+0] - xyz[3*t1+0];
        xyz1[3*i+1] = xyz[3*i+1] - xyz[3*t1+1];
        xyz1[3*i+2] = xyz[3*i+2] - xyz[3*t1+2];
    }

    // Rotation matrix, check Wikipedia!
    double** rotm = new double*[3];
    rotm[0] = new double[3];
    rotm[1] = new double[3];
    rotm[2] = new double[3];

    //unit vectors representing angles wrt the xyz planes
    double* u1 = new double[3];
    double* x1 = new double[3];

    // A vector pointing from the arbitrary origin to the free side 
    // of the molecule (t2 coordinates are defined based on averagedBondVector2)
    x1[0] = xyz1[3*t2+0] - xyz1[3*t1+0];
    x1[1] = xyz1[3*t2+1] - xyz1[3*t1+1];
    x1[2] = xyz1[3*t2+2] - xyz1[3*t1+2];

    double n1 = sqrt(x1[0]*x1[0]+x1[2]*x1[2]);

    if (n1>0.)
        for (int i=0;i<3;i++) u1[i] = x1[i]/n1;

    double* angles = new double[3];
    angles[0] = angles[1] = angles[2] = 0.;

    angles[1] = acos(u1[2]);
    if (u1[2]<0.) angles[1] = angles[1] * -1;

    Utils::get_rotation_matrix(rotm,angles);

    applyRotationMatrix(numOfAtoms, xyz1, rotm);

    /*
    double* xyzn = new double[3*numOfAtoms];
    for (int i=0;i<3*numOfAtoms;i++) xyzn[i] = 0.;
    for (int i=0;i<numOfAtoms;i++)
        for (int j=0;j<3;j++)
            for (int k=0;k<3;k++)
                xyzn[3*i+j] += rotm[j][k]*xyz1[3*i+k];
    for (int i=0;i<3*numOfAtoms;i++)
        xyz1[i] = xyzn[i];
        */

    //start second rotation

    n1 = sqrt(x1[2]*x1[2]+x1[1]*x1[1]); //Mina
    if (n1>0.0)
        for (int i=0;i<3;i++) u1[i] = x1[i]/n1;


    angles[0] = angles[1] = angles[2] = 0.;
    angles[0] = acos(u1[2]); //Mina

    Utils::get_rotation_matrix(rotm,angles);
    applyRotationMatrix(numOfAtoms, xyz1, rotm);


    if (sign==-1)
    {
        angles[0] = angles[1] = angles[2] = 0.;
        angles[0] = PI;

        Utils::get_rotation_matrix(rotm,angles);
        applyRotationMatrix(numOfAtoms, xyz1, rotm);

    }


#if 0
    //checking final geom
    double THRESH = 0.1;
    if (abs(xyzn[3*t1+0]) > THRESH 
            || abs(xyzn[3*t1+1]) > THRESH
            || abs(xyzn[3*t1+2]) > THRESH
            || abs(xyzn[3*t2+1]) > THRESH
            || abs(xyzn[3*t2+2]) > THRESH)
    {
        printf("\n WARNING: align did not set geom to 0,0,0/a,0,0 \n");
        //align_to_z(numOfAtoms,t1,t2,xyzn,atomicNames);
    }
#endif

    for (int i=0;i<3*numOfAtoms;i++)
        xyz[i] = xyz1[i];

    delete [] angles;
    for (int i=0;i<3;i++)
        delete [] rotm[i];
    delete [] rotm;

    //delete [] xyzn;
    delete [] xyz1;

    return;
}

//void Align::rotate_around_z(int numOfAtoms1, int numOfAtoms2, double torv, double* xyz)
void Align::rotate_around_z(int numOfAtoms2, double torv, double* xyz)
{
    //printf("  rotating fragment around x axis: %3.2f (degrees) \n",torv);

    //int numOfAtoms = numOfAtoms1 + numOfAtoms2;
    int numOfAtoms = numOfAtoms2;

    double* xyzTemp = new double[3*numOfAtoms];
    for (int i=0;i<3*numOfAtoms;i++) xyzTemp[i] = 0.;

    double* angles = new double[3];
    double** rotm = new double*[3];
    rotm[0] = new double[3];
    rotm[1] = new double[3];
    rotm[2] = new double[3];

    //angles[0] = torv * PI/180.;
    //angles[1] = angles[2] = 0.;
    //angles[0] = torv * 180./PI;
    //angles[1] = angles[2] = 0.;
    angles[0] = angles[1] = 0.;
    if (torv >= 0)
    {
        // clockwise rotation for positive torv
        //angles[2] = torv * -180./PI;
        angles[2] = -1 * torv;
    }
    else if (torv < 0)
    {
        // counter-clockwise rotation for negative torv
        //angles[2] = torv * 180./PI;
        angles[2] = torv;
    }
    //angles[2] = torv * PI/180.;

    Utils::get_rotation_matrix(rotm, angles);
    applyRotationMatrix(numOfAtoms, xyz, rotm);

    delete [] xyzTemp;
    delete [] angles;
    delete [] rotm[0];
    delete [] rotm[1];
    delete [] rotm[2];
    delete [] rotm;

    return;
}

void Align::linear_right(double* v1, int atom1, int* bonded, double* xyz, std::string orientationIn="horiz")
{
    double* x1 = new double[3];
    for (int i=0;i<3;i++)
        x1[i] = xyz[3*atom1+i] - xyz[3*bonded[0]+i];

    std::string orientation = orientationIn;
    double* x2 = new double[3];
    if (orientation == "horiz")
    {
        x2[0] = x2[1] = 0.;
        x2[2] = 1.;
    }
    else if (orientation == "vertical")
    {
        //TODO
    }

    cross(v1,x1,x2);

    double n1 = norm(v1,3);
    if (n1<0.001)
    {
        x2[0] = x2[2] = 0.;
        x2[1] = 1.0;
        cross(v1,x1,x2);
    }

    delete [] x1;
    delete [] x2;

    return;
}

void Align::planar_cross(double* v1, int atom1, int* bonded, double* xyz)
{
    double* x1 = new double[3];
    double* x2 = new double[3];

    int b2 = bonded[0];
    int b3 = bonded[1];
    for (int i=0;i<3;i++)
        x1[i] = xyz[3*b2+i] - xyz[3*atom1+i];
    for (int i=0;i<3;i++)
        x2[i] = xyz[3*b3+i] - xyz[3*atom1+i];

    double n1 = norm(x1,3);
    for (int i=0;i<3;i++)
        x1[i] = x1[i] / n1;
    n1 = norm(x2,3);
    for (int i=0;i<3;i++)
        x2[i] = x2[i] / n1;

    cross(v1,x1,x2);

    delete [] x1;
    delete [] x2;

    return;
}

void Align::align_v1(int nvf, double* v1)
{
    double dot;
    for (int i=1;i<nvf;i++)
    {
        dot = 0.;
        for (int j=0;j<3;j++)
            dot += v1[j] * v1[3*i+j];
        //printf("  dot (%i to 0): %4.3f \n",i,dot);
        if (dot<0.0)
            for (int j=0;j<3;j++)
                v1[3*i+j] = - v1[3*i+j];
    }

    return;
}

void Align::point_out(double* v1, int numOfAtoms, double* xyz)
{
    return;
}

//void Align::vdw_vector_opt(int numOfAtoms1, int numOfAtoms2, double* v1, ICoord icp, int atom1, int atom2)
void Align::vdw_vector_opt(int numOfAtoms1, int numOfAtoms2, double* v1, ICoord icp)
{
    //printf(" in vdw_vector_opt \n");
    double Energy = 0.;
    double potentialE = 0.;
    double THRESH = 0.001;
    double energyMin = 1000000.;
    int numOfSteps = 100;
    for (int i=0;i<numOfSteps;i++)
    {
        potentialE = Energy;
        Energy = icp.mm_energy();
        //Energy = icp.mm_energy(atom1, atom2);
        //printf(" Energy: %6.5f \n",Energy);

        if (Energy<energyMin) energyMin = Energy;
        if (Energy>potentialE+THRESH && Energy>energyMin) break;

        for (int j=0;j<numOfAtoms2;j++)
            for (int k=0;k<3;k++)
                icp.coords[3*numOfAtoms1+3*j+k] -= 0.1*v1[k];
                //printf(" v1[3]: %6.5f \n",v1[3]);

        //std::cout << "after each step:\n";
        //print_xyz_gen(icp.natoms,icp.anames,icp.coords);
    }
    return;
}

void Align::init(int numOfAtoms1i, string* atomicNames1i, int* anumbers1i, double* xyz1i,
                 int numOfAtoms2i, string* atomicNames2i, int* anumbers2i, double* xyz2i)
{
    //printf(" Align init: %i %i \n",numOfAtoms1i, numOfAtoms2i);
    inited = 1;

    numOfAtoms1 = numOfAtoms1i;
    numOfAtoms2 = numOfAtoms2i;
    atomicNames1 = new string[numOfAtoms1];
    atomicNames2 = new string[numOfAtoms2];
    anumbers1 = new int[numOfAtoms1];
    anumbers2 = new int[numOfAtoms2];
    xyz1 = new double[3*numOfAtoms1];
    xyz2 = new double[3*numOfAtoms2];
    xyz1s = new double[3*numOfAtoms1];
    xyz2s = new double[3*numOfAtoms2];
    xyzAugmented = new double[3*(numOfAtoms1+numOfAtoms2)];

    for (int i=0;i<numOfAtoms1;i++)
    {
        atomicNames1[i] = atomicNames1i[i];
        anumbers1[i] = anumbers1i[i];
    }
    for (int i=0;i<numOfAtoms2;i++)
    {
        atomicNames2[i] = atomicNames2i[i];
        anumbers2[i] = anumbers2i[i];
    }
    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1s[i] = xyz1[i] = xyz1i[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2s[i] = xyz2[i] = xyz2i[i];
    for (int i=0;i<3*numOfAtoms1;i++)
        xyzAugmented[i] = xyz1i[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyzAugmented[3*numOfAtoms1+i] = xyz2i[i];

    atomicNames3 = NULL;
    anumbers3 = NULL;
    xyz3 = NULL;
    xyzAugmented3 = NULL;

    return;
}

void Align::add_third(int numOfAtoms3i, string* atomicNames3i, int* anumbers3i, double* xyz3i)
{
    if (!inited)
    {
        printf(" WARNING: cannot add_third without initing first \n");
        return;
    }
    if (inited==3)
    {
        if (atomicNames3!=NULL) delete [] atomicNames3;
        if (anumbers3!=NULL) delete [] anumbers3;
        if (xyz3!=NULL) delete [] xyz3;
        if (xyzAugmented3!=NULL) delete [] xyzAugmented3;
    }
    inited = 3;

    numOfAtoms3 = numOfAtoms3i;
    atomicNames3 = new string[numOfAtoms3];
    anumbers3 = new int[numOfAtoms3];
    xyz3 = new double[3*numOfAtoms3];
    xyzAugmented3 = new double[3*(numOfAtoms1+numOfAtoms2+numOfAtoms3)];

    for (int i=0;i<numOfAtoms3;i++)
    {
        atomicNames3[i] = atomicNames3i[i];
        anumbers3[i] = anumbers3i[i];
    }
    for (int i=0;i<3*numOfAtoms3;i++)
        xyz3[i] = xyz3i[i];

    return;
}

void Align::align_zero()
{
  if (inited==0) return;

  //printf(" in align_zero \n");

 //get vector pointed along COM's
 //move close


 //full pair geometry
  for (int i=0;i<3*numOfAtoms1;i++)
    xyzAugmented[i] = xyz1[i];
  for (int i=0;i<3*numOfAtoms2;i++)
    xyzAugmented[3*numOfAtoms1+i] = xyz2[i];

  return;
}

//void Align::add_align(int nadd1, int* add1, double zBindingSite)
void Align::add_align(int nadd1, int* add1, int atomIndex2, vector<double> angleSet, std::string orientaionIn="horiz")
{
    if (inited==0)
    {
        printf(" WARNING: in add_align, not init'd \n");
        return;
    }

    //get v1's for all add atoms
    // TODO: magic numbers??
    double* averagedBondVector1 = new double[3*nadd1];
    for (int i=0;i<3*nadd1;i++) averagedBondVector1[i] = 0.;
    double* averagedBondVector2 = new double[3*nadd1]; //x, y, z, coordinate
    for (int i=0;i<3*nadd1;i++) averagedBondVector2[i] = 0.;
    //double* c1 = new double[6]; //center points

//    int* bonded1 = new int[8];
    int* bonded2 = new int[8]; //size of 8 bc of the octohedral geometry
//    int* bonded3 = new int[8];

    // TODO: There is a Align::init which is called before calling "add_align"
    // function, so why repeating it here?
    ICoord ic1, ic2;
    ic1.init(numOfAtoms1,atomicNames1,anumbers1,xyz1);
    ic1.make_frags();
    ic1.bond_frags();

    ic2.init(numOfAtoms2,atomicNames2,anumbers2,xyz2);
    ic2.make_frags();
    ic2.bond_frags();

    std::string orientaion = orientaionIn;

    //int npf1 = 0;
    int npf2 = 0;
    //int nvf1 = 0;
    int nvf2 = 0;
    int found = 0;
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        if (atom1>atom2) {std::swap(atom1, atom2);}

        int same_frag = check_frag(atom1,atom2);
        //printf(" atom1: %i atom2: %i same_frag: %i \n",atom1+1,atom2+1,same_frag);
        atom2 -= numOfAtoms1;
        if (!same_frag)
        {
            found = 1;
            //int nbondsatom1, numOfBondedToAtom2;
            int numOfBondedToAtom2 = 0;
            // TODO: we can use this to see if the binding site is available
            //nbondsatom1 = get_bonds(atom1,ic1,bonded1);
            numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);

            //first atom's vector (surface)
            averagedBondVector1[0] = averagedBondVector1[1] = 0.0;
            averagedBondVector1[2] = 1.0;

            //second atom's vector
            //Move atom2 to origin so the calculated vector will start from origin
            //Other way is to calculate the vector first and then add its coordinates
            //to the coordinates of atom3
            moveToOrigin(atom2, atomIndex2, i, nadd1);
            for (int i=0;i<3*nadd1;i++) averagedBondVector2[i] = 0.; //TODO, not sure if this line should be here
            for (int j=0;j<numOfBondedToAtom2;j++)
            {
                averagedBondVector2[0] += xyz2[3*atom2+0] - xyz2[3*bonded2[j]+0];
                averagedBondVector2[1] += xyz2[3*atom2+1] - xyz2[3*bonded2[j]+1];
                averagedBondVector2[2] += xyz2[3*atom2+2] - xyz2[3*bonded2[j]+2];
            }

            if (numOfBondedToAtom2==2)
            {
                double anglev = ic2.angle_val(bonded2[0],atom2,bonded2[1]);
                if (anglev>175.)
                {
                    //printf("  using ool v1 (linear angle found) \n");
                    //linear_right(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz2, orientaion);
                    linear_right(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz2);
                }
            }

            if (numOfBondedToAtom2==3) //TODO: Does this work for planar molecules like benzene?
            {
                double improperTorsionAngle = ic2.torsion_val(bonded2[0],atom2,bonded2[1],bonded2[2]);
                //printf("  improperTorsionAngle: %4.2f \n",improperTorsionAngle);
                if (fabs(improperTorsionAngle)>175.)
                {
                    //printf("  using oop v1 (planar atom found) \n");
                    planar_cross(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz2);
                    npf2++;
                }
            }
            if (npf2>1)
            {
                align_v1(nvf2+1,averagedBondVector2); //TODO: what does this do?
            }
            //TODO double n1 = norm(&averagedBondVector2[3*nvf2],3);
            double n1 = norm(&averagedBondVector2[0],3);
            for (int j=0;j<3;j++)
                averagedBondVector2[3*nvf2+j] = averagedBondVector2[3*nvf2+j] / n1; // TODO: not sure if this does anything

            nvf2++;

        } //if !same_frag
    } //loop i over nadd

    if (!found)
    {
        //printf("  couldn't find anything to align \n");
        return;
    }

    //averaging over averagedBondVector1 and averagedBondVector2 to get v2
    double* v2 = new double[6];
    for (int i=0;i<6;i++) v2[i] = 0.;
    for (int i=0;i<nadd1;i++) //TODO: why averaging over all adds ??
    {
        for (int j=0;j<3;j++)
        {
            //v2[j]   += averagedBondVector1[3*i+j];
            v2[j]   = averagedBondVector1[3*i+j];
            //v2[3+j] += averagedBondVector2[3*i+j];
            v2[3+j] = averagedBondVector2[3*i+j]; //TODO: changed += to =, nothing changed for one add
        }
    }

//    double n1 = norm(v2,3);
    double n2 = norm(&v2[3],3);
/*    if (n1<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        if (atom1>atom2) {std::swap(atom1, atom2);}
        int nbondsatom1 = get_bonds(atom1,ic1,bonded1);
        linear_right(v2,atom1,bonded1,xyz1);
        n1 = norm(v2,3);
    }
*/
    if (n2<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        if (atom1>atom2) {std::swap(atom1, atom2);}
        atom2 -= numOfAtoms1;
        int numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);
        linear_right(&v2[3],atom2,bonded2,xyz2);
        n2 = norm(&v2[3],3);
        //    printf(" new norm: %8.6f \n",n2);
    }
    //for (int i=0;i<3;i++)
    //    v2[i] = v2[i] / n1;
    // normalize v2
    for (int i=3;i<6;i++)
        v2[i] = v2[i] / n2;

    // get ready for rotation
    double* xyz2a = new double[3*(numOfAtoms2+2)]; // TODO: +2 ??
    string* atomicNames2a = new string[numOfAtoms2+2];
    for (int i=0;i<numOfAtoms2;i++) atomicNames2a[i] = atomicNames2[i];
    atomicNames2a[numOfAtoms2] = "X";
    atomicNames2a[numOfAtoms2+1] = "X";

    for (int i=0;i<numOfAtoms2;i++)
    {
        for (int j=0;j<3;j++)
        {
            xyz2a[3*i+j] = xyz2[3*i+j];
        }
    }
    xyz2a[3*(numOfAtoms2)+0] = xyz2a[3*(numOfAtoms2)+1] = xyz2a[3*(numOfAtoms2)+2] = 0.;
    xyz2a[3*(numOfAtoms2+1)+0] = v2[3];
    xyz2a[3*(numOfAtoms2+1)+1] = v2[4];
    xyz2a[3*(numOfAtoms2+1)+2] = v2[5];

    // align to Z
    int t1 = numOfAtoms2; // arbitrary origin
    int t2 = numOfAtoms2+1; // head of v2 vector
    align_to_z(numOfAtoms2+2,t1,t2,xyz2a,atomicNames2a,-1);

    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyz2a[i];

    // move atom2 to the origin of Cartesians, and obviously move its attached atoms too!
    for (int i=0; i<nadd1; i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        atom2 -= numOfAtoms1;
        moveToOrigin(atom2, atomIndex2, i, nadd1);

        // assuming atom2 is the central atom. Find the vector from atom2 to one of
        // the attached atoms. Find the cross product of that vector with averagedBondVector2, the 
        // result should equal the surface normal
        double surfNormal[3] = {};
        double bondVector[3] = {};
        bondVector[0] = xyz2[3*bonded2[0]+0] - xyz2[3*atom2+0];
        bondVector[1] = xyz2[3*bonded2[0]+1] - xyz2[3*atom2+1];
        bondVector[2] = xyz2[3*bonded2[0]+2] - xyz2[3*atom2+2];

        double zDir[3] = {0., 0., -1.};
        cross(surfNormal, bondVector, zDir);
        // find the angle between surface normal and x axis
        double xAxis[3] = {1., 0., 0.};
        double dotProduct = surfNormal[0] * xAxis[0] + surfNormal[1] * xAxis[1] + surfNormal[2] * xAxis[2];
        double lengthOf1st = sqrt(pow(surfNormal[0], 2) + pow(surfNormal[1], 2) + pow(surfNormal[2], 2));
        double lengthOf2nd = 1.0; // bc it's a unit vector along X axis
        double cosTheta = dotProduct / (lengthOf1st * lengthOf2nd);
        double angleToX = acos(cosTheta)*180./PI;

        // align surface normal to x
        rotate_around_z(numOfAtoms2, angleToX, xyz2);

        double xyz2AtZeroDegree[3*numOfAtoms2] = {};
        for (int i=0;i<3*numOfAtoms2;i++) xyz2AtZeroDegree[i] = xyz2[i];
        moveToBindingSite(atom1, atom2, atomIndex2, i, nadd1);
        std::string outFileName = "alignedStr-0.xyz";
        this->writeToFile(outFileName);


        for (int j=0; j<angleSet.size(); j++)
        {
            for (int i=0;i<3*numOfAtoms2;i++) xyz2[i] = xyz2AtZeroDegree[i];
            //xyzTemp = xyz2
            //moveToOrigin(atom2, atomIndex2, i, nadd1);
            rotate_around_z(numOfAtoms2, angleSet[j]*PI/180., xyz2);
            moveToBindingSite(atom1, atom2, atomIndex2, i, nadd1);
            std::string outFileName = "alignedStr-" + std::to_string((int)angleSet[j]) + ".xyz";
            this->writeToFile(outFileName);
        }
    }

        //face each other at X Angstroms
        //double X = 8.;
        //for (int i=0;i<3*numOfAtoms1;i++)
        //xyz1[i] = xyz1a[i];
        //for (int i=0;i<3*numOfAtoms2;i++)
        //    xyz2[i] = xyz2Displaced[i];
        //for (int i=0;i<numOfAtoms2;i++)
        //xyz2[3*i+0] += X;

    //prepare for vdw opt
    string atomicNamesCombined[numOfAtoms1+numOfAtoms2];
    int atomicNumbersCombined[numOfAtoms1+numOfAtoms2];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNamesCombined[i] = atomicNames1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNamesCombined[numOfAtoms1+i] = atomicNames2[i];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNumbersCombined[i] = anumbers1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNumbersCombined[numOfAtoms1+i] = anumbers2[i];
    //full pair geometry, before vdw opt
    for (int i=0;i<3*numOfAtoms1;i++)
        xyzAugmented[i] = xyz1[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyzAugmented[3*numOfAtoms1+i] = xyz2[i];

    ICoord icp;
    icp.init(numOfAtoms1+numOfAtoms2,atomicNamesCombined,atomicNumbersCombined,xyzAugmented);

    //print_xyz();

    //rotate about central axis
    /*if (nadd1>1)
    {
        int atom1 = add1[0];
        int atom2 = add1[1];
        //int tmp1 = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp1; }
        if (atom1>atom2) { std::swap(atom1, atom2); }
        int atom3 = add1[2];
        int atom4 = add1[3];
        //int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }
        if (atom3>atom4) { std::swap(atom3, atom4); }

        int same_frag = check_frag(atom1,atom2) + check_frag(atom3,atom4);
        if (!same_frag && atom1!=atom3 && atom2!=atom4)
        {
            double torv = icp.torsion_val(atom1,atom3,atom4,atom2);
            //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,atom1+1,a3+1,a4+1,atom2+1);
            rotate_around_z(numOfAtoms1,numOfAtoms2,-torv,icp.coords); //TODO
        }
    }*/

    double* vx = new double[3];
    vx[0] = vx[1] = 0.;
    vx[2] = 1.0;
    // Suppose atom1tmp is always a surface atom
    vdw_vector_opt(numOfAtoms1, numOfAtoms2, vx, icp);
    delete [] vx;

    for (int i=0;i<3*(numOfAtoms1+numOfAtoms2);i++)
        xyzAugmented[i] = icp.coords[i];

    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyzAugmented[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyzAugmented[3*numOfAtoms1+i];

    //print_xyz_gen(numOfAtoms1+numOfAtoms2,icp.atomicNames,xyzAugmented);

    int badgeom = 0;
    for (int i=0;i<3*numOfAtoms1;i++)
        if (xyz1[i]!=xyz1[i])
        {
            badgeom = 1;
            break;
        }
    for (int i=0;i<3*numOfAtoms2;i++)
        if (xyz2[i]!=xyz2[i])
        {
            badgeom = 1;
            break;
        }

    if (badgeom)
    {
        printf(" ERROR: nan detected in add_align, exiting! \n");
        //print_xyz_gen(numOfAtoms1,atomicNames1,xyz1); Mina
        //print_xyz_gen(numOfAtoms2,atomicNames2,xyz2); Mina
        exit(1);
    }


    //delete [] atomicNames1a;
    delete [] atomicNames2a;
    //delete [] xyz1a;
    delete [] xyz2a;
    delete [] v2;
    delete [] averagedBondVector1;
    delete [] averagedBondVector2;
//    delete [] bonded1;
    delete [] bonded2;
//    delete [] bonded3;
    ic1.freemem();
    ic2.freemem();
    //delete [] atomicNamesCombined;
    //delete [] atomicNumbersCombined;
    icp.freemem();

    return;
}

/*
int Align::add_align_v(int nadd1, int* add1, int wtm, double* aprv)
{
    if (inited==0)
    {
        printf(" WARNING: in add_align, not init'd \n");
        return 0;
    }
    //revert to initial geometries
    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyz1s[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyz2s[i];

#if 1
    printf("   in add_align_v, adding:");
    for (int i=0;i<nadd1;i++)
        printf(" %i-%i",add1[2*i+0]+1,add1[2*i+1]+1);
    printf("\n");
    //printf(" in add_align numOfAtoms1: %i numOfAtoms2: %i \n",numOfAtoms1,numOfAtoms2);
#endif

    ///get v1's for all add atoms
    double* averagedBondVector1 = new double[3*nadd1];
    for (int i=0;i<3*nadd1;i++) averagedBondVector1[i] = 0.;
    double* averagedBondVector2 = new double[3*nadd1];
    for (int i=0;i<3*nadd1;i++) averagedBondVector2[i] = 0.;
    double* c1 = new double[6]; //center points

    int* bonded1 = new int[8];
    int* bonded2 = new int[8];
    int* bonded3 = new int[8];

    ICoord ic1,ic2;
    ic1.init(numOfAtoms1,atomicNames1,anumbers1,xyz1);
    ic1.make_frags();
    ic1.bond_frags();
    //ic1.ic_create_nobonds();
    //ic1.print_bonds();
    ic2.init(numOfAtoms2,atomicNames2,anumbers2,xyz2);
    ic2.make_frags();
    ic2.bond_frags();
    //ic2.ic_create_nobonds();
    //ic2.print_bonds();

    int npf1 = 0;
    int npf2 = 0;
    int nvf1 = 0;
    int nvf2 = 0;
    int found = 0;
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }

        int same_frag = check_frag(atom1,atom2);
        //printf(" atom1: %i atom2: %i same_frag: %i \n",atom1+1,atom2+1,same_frag);
        atom2 -= numOfAtoms1;
        if (!same_frag)
        {
            found = 1;
            int nbondsatom1, numOfBondedToAtom2;
            nbondsatom1 = get_bonds(atom1,ic1,bonded1);
            numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);

            //first atom's vector
            for (int j=0;j<nbondsatom1;j++)
            {
                averagedBondVector1[3*nvf1+0] += xyz1[3*atom1+0] - xyz1[3*bonded1[j]+0];
                averagedBondVector1[3*nvf1+1] += xyz1[3*atom1+1] - xyz1[3*bonded1[j]+1];
                averagedBondVector1[3*nvf1+2] += xyz1[3*atom1+2] - xyz1[3*bonded1[j]+2];
            }
#if 0
            if (nbondsatom1==1 && ic1.coordn[bonded1[0]]==3)
            {
                int atom1c = bonded1[0];
                int nbondsatom1c = get_bonds(atom1c,ic1,bonded3)
                    double imptorv = ic1.torsion_val(bonded3[0],atom1,bonded3[1],bonded3[2]);
                //printf("  imptorv: %4.2f \n",imptorv);
                if (fabs(imptorv)>175.)
                {
                    //printf("  using oop v1 \n");
                    planar_cross(&averagedBondVector1[3*nvf1],a3,bonded3,xyz1);
                }
            }
#endif
            if (nbondsatom1==2)
            {
                double anglev = ic1.angle_val(bonded1[0],atom1,bonded1[1]);
                if (anglev>175.)
                {
                    //printf("  using ool v1 \n");
                    linear_right(&averagedBondVector1[3*nvf1],atom1,bonded1,xyz1);
                }
            }
            if (nbondsatom1==3)
            { 
                double imptorv = ic1.torsion_val(bonded1[0],atom1,bonded1[1],bonded1[2]);
                //printf("  imptorv: %4.2f \n",imptorv);
                if (fabs(imptorv)>175.)
                {
                    //printf("  using oop v1 \n");
                    planar_cross(&averagedBondVector1[3*nvf1],atom1,bonded1,xyz1);
                    npf1++;
                }
            }
            if (npf1>1)
            {
                align_v1(nvf1+1,averagedBondVector1);
            }
            double n1 = norm(&averagedBondVector1[3*nvf1],3);
            for (int j=0;j<3;j++)
                averagedBondVector1[3*nvf1+j] = averagedBondVector1[3*nvf1+j] / n1;
#if 0
            printf(" averagedBondVector1(atom1): \n");
            for (int j=0;j<3;j++)
                printf(" %4.3f",averagedBondVector1[3*nvf1+j]);
            printf("\n");
#endif

            nvf1++;


            //second atom's vector
            for (int j=0;j<numOfBondedToAtom2;j++)
            {
                averagedBondVector2[3*nvf2+0] += xyz2[3*atom2+0] - xyz2[3*bonded2[j]+0];
                averagedBondVector2[3*nvf2+1] += xyz2[3*atom2+1] - xyz2[3*bonded2[j]+1];
                averagedBondVector2[3*nvf2+2] += xyz2[3*atom2+2] - xyz2[3*bonded2[j]+2];
            }
            if (numOfBondedToAtom2==2)
            {
                double anglev = ic2.angle_val(bonded2[0],atom2,bonded2[1]);
                if (anglev>175.)
                {
                    //printf("  using ool v1 (linear angle found) \n");
                    linear_right(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz2);
                }
            }
            if (numOfBondedToAtom2==3)
            {
                double imptorv = ic2.torsion_val(bonded2[0],atom2,bonded2[1],bonded2[2]);
                //printf("  imptorv: %4.2f \n",imptorv);
                if (fabs(imptorv)>175.)
                {
                    //printf("  using oop v1 (planar atom found) \n");
                    planar_cross(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz2);
                    npf2++;
                }
            }
            if (npf2>1)
            {
                align_v1(nvf2+1,averagedBondVector2);
            }
            n1 = norm(&averagedBondVector2[3*nvf2],3);
            for (int j=0;j<3;j++)
                averagedBondVector2[3*nvf2+j] = averagedBondVector2[3*nvf2+j] / n1;
#if 0
            printf(" averagedBondVector2(atom2): \n");
            for (int j=0;j<3;j++)
                printf(" %4.3f",averagedBondVector2[3*nvf2+j]);
            printf("\n");
#endif

            nvf2++;

        } // if !same_frag
    } //loop i over nadd

    if (!found)
    {
        //printf("  couldn't find anything to align \n");
        return 0;
    }

#if 0
    //CPMZ NOT IMPLEMENTED
    if (npf1>0)
        point_out(averagedBondVector1,numOfAtoms1,xyz1);
    if (npf2>0)
        point_out(averagedBondVector2,numOfAtoms2,xyz2);
#endif


    //averaging over v1 to get v2
    double* v2 = new double[6];
    for (int i=0;i<6;i++) v2[i] = 0.;
    for (int i=0;i<nadd1;i++)
    {
        for (int j=0;j<3;j++)
            v2[j]   += averagedBondVector1[3*i+j];
        for (int j=0;j<3;j++)
            v2[3+j] += averagedBondVector2[3*i+j];
    }

    double n1 = norm(v2,3);
    double n2 = norm(&v2[3],3);
    if (n1<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        int nbondsatom1 = get_bonds(atom1,ic1,bonded1);
        linear_right(v2,atom1,bonded1,xyz1);
        n1 = norm(v2,3);
        //    printf(" new norm: %8.6f \n",n1);
    }
    if (n2<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        atom2 -= numOfAtoms1;
        int numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);
        linear_right(&v2[3],atom2,bonded2,xyz2);
        n2 = norm(&v2[3],3);
        //    printf(" new norm: %8.6f \n",n2);
    }
    for (int i=0;i<3;i++)
        v2[i] = v2[i] / n1;
    for (int i=3;i<6;i++)
        v2[i] = v2[i] / n2;

#if 0
    printf("     v2(1):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[j]);
    printf("\n");
    printf("     v2(2):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[3+j]);
    printf("\n");
#endif

    //get central location for each
    int cfound = 0;
    for (int i=0;i<6;i++) c1[i] = 0.;
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        int same_frag = check_frag(atom1,atom2);

        atom2 -= numOfAtoms1;
        if (!same_frag)
        {
            //printf(" c1 for %i-%i \n",atom1+1,atom2+1);
            cfound++;
            for (int j=0;j<3;j++)
                c1[j]   += xyz1[3*atom1+j];
            for (int j=0;j<3;j++)
                c1[3+j] += xyz2[3*atom2+j];
        }
    } //loop i over nadd
    if (cfound>1)
        for (int j=0;j<6;j++)
            c1[j] = c1[j] / cfound;

#if 0
    printf("  c1:");
    for (int j=0;j<6;j++)
        printf(" %4.3f",c1[j]);
    printf("\n");
#endif


    //find TM fragment, replace c1 and v terms with wtm center and aprv
    int f1 = 0;
    for (int i=0;i<numOfAtoms1;i++)
        if (ic1.isTM(i))
            f1 = 1;
    int f2 = 0;
    for (int i=0;i<numOfAtoms2;i++)
        if (ic2.isTM(i))
            f2 = 1;

    int TMaddo = 0; //is TM adding to other reactant?
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }

        int same_frag = check_frag(atom1,atom2);
        //printf(" atom1: %i atom2: %i same_frag: %i \n",atom1+1,atom2+1,same_frag);
        atom2 -= numOfAtoms1;
        if (!same_frag)
        {
            if (ic1.isTM(atom1))
            {
                TMaddo++;
                f1++;
            }
            if (ic2.isTM(atom2))
            {
                TMaddo++;
                f2++;
            }
        }
    }

    //printf("    fragment with TM: %i %i add across TM-other: %i \n",f1,f2,TMaddo);
    if (f1==f2 && TMaddo)
    {
        printf("   ERROR: neither or both fragments have TMs! \n");
        exit(1);
    }
    else if (f1>0 && TMaddo)
    {
        c1[0] = xyz1[3*wtm+0];
        c1[1] = xyz1[3*wtm+1];
        c1[2] = xyz1[3*wtm+2];
        v2[0] = aprv[0];
        v2[1] = aprv[1];
        v2[2] = aprv[2];
    }
    else if (f2>0 && TMaddo)
    {
        wtm  -= numOfAtoms1;
        c1[3] = xyz2[3*wtm+0];
        c1[4] = xyz2[3*wtm+1];
        c1[5] = xyz2[3*wtm+2];
        v2[3] = aprv[0];
        v2[4] = aprv[1];
        v2[5] = aprv[2];
    }
    else
        printf("   not using TM-specific vector \n");

#if 0
    printf("  TM c1:");
    for (int j=0;j<6;j++)
        printf(" %4.3f",c1[j]);
    printf("\n");
#endif
#if 0
    printf("  TM v2(1):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[j]);
    printf("\n");
    printf("  TM v2(2):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[3+j]);
    printf("\n");
#endif

    //move geometry to center
    double* xyz1a = new double[3*(numOfAtoms1+2)];
    double* xyz2a = new double[3*(numOfAtoms2+2)];
    string* atomicNames1a = new string[numOfAtoms1+2];
    string* atomicNames2a = new string[numOfAtoms2+2];
    for (int i=0;i<numOfAtoms1;i++) atomicNames1a[i] = atomicNames1[i];
    atomicNames1a[numOfAtoms1] = "X";
    atomicNames1a[numOfAtoms1+1] = "X";
    for (int i=0;i<numOfAtoms2;i++) atomicNames2a[i] = atomicNames2[i];
    atomicNames2a[numOfAtoms2] = "X";
    atomicNames2a[numOfAtoms2+1] = "X";

    for (int i=0;i<numOfAtoms1;i++)
    {
        for (int j=0;j<3;j++)
            xyz1a[3*i+j] = xyz1[3*i+j] - c1[j];
    }
    xyz1a[3*(numOfAtoms1)+0] = 0.; xyz1a[3*(numOfAtoms1)+1] = 0.;
    xyz1a[3*(numOfAtoms1)+2] = 0.; xyz1a[3*(numOfAtoms1+1)+0] = v2[0];
    xyz1a[3*(numOfAtoms1+1)+1] = v2[1]; xyz1a[3*(numOfAtoms1+1)+2] = v2[2];
    //print_xyz_gen(numOfAtoms1+2,atomicNames1a,xyz1a);

    for (int i=0;i<numOfAtoms2;i++)
    {
        for (int j=0;j<3;j++)
            xyz2a[3*i+j] = xyz2[3*i+j] - c1[3+j];
    }
    xyz2a[3*(numOfAtoms2)+0] = 0.; xyz2a[3*(numOfAtoms2)+1] = 0.;
    xyz2a[3*(numOfAtoms2)+2] = 0.; xyz2a[3*(numOfAtoms2+1)+0] = v2[3];
    xyz2a[3*(numOfAtoms2+1)+1] = v2[4]; xyz2a[3*(numOfAtoms2+1)+2] = v2[5];
    //print_xyz_gen(numOfAtoms2+2,atomicNames2a,xyz2a);


    //align v2 along x/-x direction
    int t1 = numOfAtoms1;
    int t2 = numOfAtoms1+1;
    align_to_z(numOfAtoms1+2,t1,t2,xyz1a,atomicNames1a,1,10.);
    t1 = numOfAtoms2;
    t2 = numOfAtoms2+1;
    align_to_z(numOfAtoms2+2,t1,t2,xyz2a,atomicNames2a,-1,10.);

    //print_xyz_gen(numOfAtoms1+2,atomicNames1a,xyz1a);
    //print_xyz_gen(numOfAtoms2+2,atomicNames2a,xyz2a);

    //face each other at X Angstroms
    double X = 8.;
    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyz1a[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyz2a[i];
    for (int i=0;i<numOfAtoms2;i++)
        xyz2[3*i+0] += X;


    //prepare for vdw opt
    string* atomicNamesCombined = new string[numOfAtoms1+numOfAtoms2];
    int* atomicNumbersCombined = new int[numOfAtoms1+numOfAtoms2];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNamesCombined[i] = atomicNames1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNamesCombined[numOfAtoms1+i] = atomicNames2[i];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNumbersCombined[i] = anumbers1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNumbersCombined[numOfAtoms1+i] = anumbers2[i];
    //full pair geometry, before vdw opt
    for (int i=0;i<3*numOfAtoms1;i++)
        xyzAugmented[i] = xyz1[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyzAugmented[3*numOfAtoms1+i] = xyz2[i];

    ICoord icp;
    icp.init(numOfAtoms1+numOfAtoms2,atomicNamesCombined,atomicNumbersCombined,xyzAugmented);

    //print_xyz();

    //rotate about central axis
    if (nadd1>1)
    {
        int atom1 = add1[0];
        int atom2 = add1[1];
        int tmp1 = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp1; }
        int a3 = add1[2];
        int a4 = add1[3];
        int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }

        int same_frag = check_frag(atom1,atom2) + check_frag(a3,a4);
        if (!same_frag && atom1!=a3 && atom2!=a4)
        {
            double torv = icp.torsion_val(atom1,a3,a4,atom2);
            //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,atom1+1,a3+1,a4+1,atom2+1);
            rotate_around_x(numOfAtoms1,numOfAtoms2,-torv,icp.coords);
        }
    }

    double* vx = new double[3];
    vx[1] = vx[2] = 0.;
    vx[0] = 1.0;
    vdw_vector_opt(numOfAtoms1,numOfAtoms2,vx,icp);
    delete [] vx;

    for (int i=0;i<3*(numOfAtoms1+numOfAtoms2);i++)
        xyzAugmented[i] = icp.coords[i];

    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyzAugmented[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyzAugmented[3*numOfAtoms1+i];

    //print_xyz_gen(numOfAtoms1+numOfAtoms2,icp.atomicNames,xyzAugmented);

    int badgeom = 0;
    for (int i=0;i<3*numOfAtoms1;i++)
        if (xyz1[i]!=xyz1[i])
        {
            badgeom = 1;
            break;
        }
    for (int i=0;i<3*numOfAtoms2;i++)
        if (xyz2[i]!=xyz2[i])
        {
            badgeom = 1;
            break;
        }

    if (badgeom)
    {
        printf(" ERROR: nan detected in add_align, exiting! \n");
        //print_xyz_gen(numOfAtoms1,atomicNames1,xyz1); Mina
        //print_xyz_gen(numOfAtoms2,atomicNames2,xyz2); Mina
        exit(1);
    }

    delete [] atomicNames1a;
    delete [] atomicNames2a;
    delete [] xyz1a;
    delete [] xyz2a;

    delete [] v2;
    delete [] averagedBondVector1;
    delete [] averagedBondVector2;
    delete [] bonded1;
    delete [] bonded2;
    delete [] bonded3;

    ic1.freemem();
    ic2.freemem();

    delete [] atomicNamesCombined;
    delete [] atomicNumbersCombined;
    icp.freemem();


    return TMaddo;
}

void Align::shuttle_align(int nadd1, int* add1)
{
    if (inited!=3)
    {
        printf(" in shuttle_align, add_third hasn't been called \n");
        exit(1);
    }

#if 0
    printf("\n in shuttle_align, adding:");
    for (int i=0;i<nadd1;i++)
        printf(" %i-%i",add1[2*i+0]+1,add1[2*i+1]+1);
    printf("\n");
    //printf(" in shuttle_align numOfAtoms1: %i numOfAtoms2: %i \n",numOfAtoms1,numOfAtoms2);
#endif
#if 0
    printf(" previous geom: \n");
    print_xyz();
#endif

    int natom12 = numOfAtoms1 + numOfAtoms2;

    ///get v1's for all add atoms
    double* averagedBondVector1 = new double[3*nadd1];
    for (int i=0;i<3*nadd1;i++) averagedBondVector1[i] = 0.;
    double* averagedBondVector2 = new double[3*nadd1];
    for (int i=0;i<3*nadd1;i++) averagedBondVector2[i] = 0.;
    double* c1 = new double[6]; //center points

    string* atomicNames12 = new string[natom12+1];
    int* anumbers12 = new int[natom12+1];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNames12[i] = atomicNames1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNames12[numOfAtoms1+i] = atomicNames2[i];
    for (int i=0;i<numOfAtoms1;i++)
        anumbers12[i] = anumbers1[i];
    for (int i=0;i<numOfAtoms2;i++)
        anumbers12[numOfAtoms1+i] = anumbers2[i];

    int* bonded1 = new int[8];
    int* bonded2 = new int[8];
    int* bonded3 = new int[8];

    ICoord ic1,ic2;
    ic1.init(natom12,atomicNames12,anumbers12,xyzAugmented);
    ic1.make_frags();
    ic1.bond_frags();
    //ic1.ic_create_nobonds();
    ic2.init(numOfAtoms3,atomicNames3,anumbers3,xyz3);
    ic2.make_frags();
    ic2.bond_frags();
    //ic2.ic_create_nobonds();

    int npf1 = 0;
    int npf2 = 0;
    int nvf1 = 0;
    int nvf2 = 0;
    int found = 0;
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }

        int same_frag = check_frag_3(atom1,atom2);
        //printf(" atom1: %i atom2: %i same_frag: %i \n",atom1+1,atom2+1,same_frag);
        atom2 -= natom12;
        if (!same_frag)
        {
            found = 1;
            int nbondsatom1, numOfBondedToAtom2;
            nbondsatom1 = get_bonds(atom1,ic1,bonded1);
            numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);

            //first atom's vector
            for (int j=0;j<nbondsatom1;j++)
            {
                averagedBondVector1[3*nvf1+0] += xyzAugmented[3*atom1+0] - xyzAugmented[3*bonded1[j]+0];
                averagedBondVector1[3*nvf1+1] += xyzAugmented[3*atom1+1] - xyzAugmented[3*bonded1[j]+1];
                averagedBondVector1[3*nvf1+2] += xyzAugmented[3*atom1+2] - xyzAugmented[3*bonded1[j]+2];
            }
            if (nbondsatom1==2)
            {
                double anglev = ic1.angle_val(bonded1[0],atom1,bonded1[1]);
                if (anglev>175.)
                {
                    //printf("  using ool v1 \n");
                    linear_right(&averagedBondVector1[3*nvf1],atom1,bonded1,xyzAugmented);
                }
            }
            if (nbondsatom1==3)
            {          
                double imptorv = ic1.torsion_val(bonded1[0],atom1,bonded1[1],bonded1[2]);
                //printf("  imptorv: %4.2f \n",imptorv);
                if (fabs(imptorv)>175.)
                {
                    //printf("  using oop v1 \n");
                    planar_cross(&averagedBondVector1[3*nvf1],atom1,bonded1,xyzAugmented);
                    npf1++;
                }
            }
            if (npf1>1)
            {
                align_v1(nvf1+1,averagedBondVector1);
            }
            double n1 = norm(&averagedBondVector1[3*nvf1],3);
            for (int j=0;j<3;j++)
                averagedBondVector1[3*nvf1+j] = averagedBondVector1[3*nvf1+j] / n1;
#if 0
            printf(" averagedBondVector1(atom1): \n");
            for (int j=0;j<3;j++)
                printf(" %4.3f",averagedBondVector1[3*nvf1+j]);
            printf("\n");
#endif

            nvf1++;


            //second atom's vector
            for (int j=0;j<numOfBondedToAtom2;j++)
            {
                averagedBondVector2[3*nvf2+0] += xyz3[3*atom2+0] - xyz3[3*bonded2[j]+0];
                averagedBondVector2[3*nvf2+1] += xyz3[3*atom2+1] - xyz3[3*bonded2[j]+1];
                averagedBondVector2[3*nvf2+2] += xyz3[3*atom2+2] - xyz3[3*bonded2[j]+2];
            }
            if (numOfBondedToAtom2==2)
            {
                double anglev = ic2.angle_val(bonded2[0],atom2,bonded2[1]);
                if (anglev>175.)
                {
                    //printf("  using ool v1 (linear angle found) \n");
                    linear_right(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz3);
                }
            }
            if (numOfBondedToAtom2==3)
            {
                double imptorv = ic2.torsion_val(bonded2[0],atom2,bonded2[1],bonded2[2]);
                //printf("  imptorv: %4.2f \n",imptorv);
                if (fabs(imptorv)>175.)
                {
                    //printf("  using oop v1 (planar atom found) \n");
                    planar_cross(&averagedBondVector2[3*nvf2],atom2,bonded2,xyz3);
                    npf2++;
                }
            }
            if (npf2>1)
            {
                align_v1(nvf2+1,averagedBondVector2);
            }
            n1 = norm(&averagedBondVector2[3*nvf2],3);
            for (int j=0;j<3;j++)
                averagedBondVector2[3*nvf2+j] = averagedBondVector2[3*nvf2+j] / n1;
#if 0
            printf(" averagedBondVector2(atom2): \n");
            for (int j=0;j<3;j++)
                printf(" %4.3f",averagedBondVector2[3*nvf2+j]);
            printf("\n");
#endif

            nvf2++;

        } // if !same_frag
    } //loop i over nadd

    if (!found)
    {
        //printf("  couldn't find anything to align \n");
        return;
    }

#if 0
    //CPMZ NOT IMPLEMENTED
    if (npf1>0)
        point_out(averagedBondVector1,numOfAtoms1,xyz1);
    if (npf2>0)
        point_out(averagedBondVector2,numOfAtoms2,xyz2);
#endif


    //averaging over v1 to get v2
    double* v2 = new double[6];
    for (int i=0;i<6;i++) v2[i] = 0.;
    for (int i=0;i<nadd1;i++)
    {
        for (int j=0;j<3;j++)
            v2[j]   += averagedBondVector1[3*i+j];
        for (int j=0;j<3;j++)
            v2[3+j] += averagedBondVector2[3*i+j];
    }

    double n1 = norm(v2,3);
    double n2 = norm(&v2[3],3);
    if (n1<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        int nbondsatom1 = get_bonds(atom1,ic1,bonded1);
        linear_right(v2,atom1,bonded1,xyzAugmented);
        n1 = norm(v2,3);
        //    printf(" new norm: %8.6f \n",n1);
    }
    if (n2<0.000001)
    {
        int atom1 = add1[2*(nadd1-1)+0];
        int atom2 = add1[2*(nadd1-1)+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        atom2 -= natom12;
        int numOfBondedToAtom2 = get_bonds(atom2,ic2,bonded2);
        linear_right(&v2[3],atom2,bonded2,xyz3);
        n2 = norm(&v2[3],3);
        //    printf(" new norm: %8.6f \n",n2);
    }
    for (int i=0;i<3;i++)
        v2[i] = v2[i] / n1;
    for (int i=3;i<6;i++)
        v2[i] = v2[i] / n2;

#if 0
    printf("  v2(1):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[j]);
    printf("\n");
    printf("  v2(2):");
    for (int j=0;j<3;j++)
        printf(" %4.3f",v2[3+j]);
    printf("\n");
#endif

    //get central location for each
    int cfound = 0;
    for (int i=0;i<6;i++) c1[i] = 0.;
    for (int i=0;i<nadd1;i++)
    {
        int atom1 = add1[2*i+0];
        int atom2 = add1[2*i+1];
        int tmp = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp; }
        int same_frag = check_frag_3(atom1,atom2);

        atom2 -= natom12;
        if (!same_frag)
        {
            //printf(" c1 for %i-%i \n",atom1+1,atom2+1);
            cfound++;
            for (int j=0;j<3;j++)
                c1[j]   += xyzAugmented[3*atom1+j];
            for (int j=0;j<3;j++)
                c1[3+j] += xyz3[3*atom2+j];
        }
    } //loop i over nadd
    if (cfound>1)
        for (int j=0;j<6;j++)
            c1[j] = c1[j] / cfound;

#if 0
    printf("  c1:");
    for (int j=0;j<6;j++)
        printf(" %4.3f",c1[j]);
    printf("\n");
#endif


    //move geometry to center
    double* xyz12a = new double[3*(natom12+2)];
    double* xyz3a = new double[3*(numOfAtoms3+2)];
    string* atomicNames12a = new string[natom12+2];
    string* atomicNames3a = new string[numOfAtoms3+2];
    for (int i=0;i<numOfAtoms1;i++) atomicNames12a[i] = atomicNames1[i];
    for (int i=0;i<numOfAtoms2;i++) atomicNames12a[numOfAtoms1+i] = atomicNames2[i];
    atomicNames12a[natom12] = "X";
    atomicNames12a[natom12+1] = "X";
    for (int i=0;i<numOfAtoms3;i++) atomicNames3a[i] = atomicNames3[i];
    atomicNames3a[numOfAtoms3] = "X";
    atomicNames3a[numOfAtoms3+1] = "X";

    for (int i=0;i<numOfAtoms1;i++)
    {
        for (int j=0;j<3;j++)
            xyz12a[3*i+j] = xyz1[3*i+j] - c1[j];
    }
    for (int i=0;i<numOfAtoms2;i++)
    {
        for (int j=0;j<3;j++)
            xyz12a[3*numOfAtoms1+3*i+j] = xyz2[3*i+j] - c1[j];
    }
    xyz12a[3*(natom12)+0] = 0.; xyz12a[3*(natom12)+1] = 0.;
    xyz12a[3*(natom12)+2] = 0.; xyz12a[3*(natom12+1)+0] = v2[0];
    xyz12a[3*(natom12+1)+1] = v2[1]; xyz12a[3*(natom12+1)+2] = v2[2];
    //print_xyz_gen(natom12+2,atomicNames12a,xyz12a);

    for (int i=0;i<numOfAtoms3;i++)
    {
        for (int j=0;j<3;j++)
            xyz3a[3*i+j] = xyz3[3*i+j] - c1[3+j];
    }
    xyz3a[3*(numOfAtoms3)+0] = 0.; xyz3a[3*(numOfAtoms3)+1] = 0.;
    xyz3a[3*(numOfAtoms3)+2] = 0.; xyz3a[3*(numOfAtoms3+1)+0] = v2[3];
    xyz3a[3*(numOfAtoms3+1)+1] = v2[4]; xyz3a[3*(numOfAtoms3+1)+2] = v2[5];
    //print_xyz_gen(numOfAtoms3+2,atomicNames3a,xyz3a);

    //align v2 along x/-x direction
    int t1 = natom12;
    int t2 = natom12+1;
    align_to_z(natom12+2,t1,t2,xyz12a,atomicNames12a,1,10.);
    t1 = numOfAtoms3;
    t2 = numOfAtoms3+1;
    align_to_z(numOfAtoms3+2,t1,t2,xyz3a,atomicNames3a,-1,10.);

    //print_xyz_gen(natom12+2,atomicNames12a,xyz12a);
    //print_xyz_gen(numOfAtoms3+2,atomicNames3a,xyz3a);

    //face each other at X Angstroms
    double X = 8.;
    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyz12a[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyz12a[3*numOfAtoms1+i];
    for (int i=0;i<3*numOfAtoms3;i++)
        xyz3[i] = xyz3a[i];
    for (int i=0;i<numOfAtoms3;i++)
        xyz3[3*i+0] += X;


    //prepare for vdw opt
    string* atomicNamesCombined = new string[natom12+numOfAtoms3];
    int* atomicNumbersCombined = new int[natom12+numOfAtoms3];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNamesCombined[i] = atomicNames1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNamesCombined[numOfAtoms1+i] = atomicNames2[i];
    for (int i=0;i<numOfAtoms3;i++)
        atomicNamesCombined[natom12+i] = atomicNames3[i];
    for (int i=0;i<numOfAtoms1;i++)
        atomicNumbersCombined[i] = anumbers1[i];
    for (int i=0;i<numOfAtoms2;i++)
        atomicNumbersCombined[numOfAtoms1+i] = anumbers2[i];
    for (int i=0;i<numOfAtoms3;i++)
        atomicNumbersCombined[natom12+i] = anumbers3[i];
    //full pair geometry, before vdw opt
    for (int i=0;i<3*numOfAtoms1;i++)
        xyzAugmented3[i] = xyz1[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyzAugmented3[3*numOfAtoms1+i] = xyz2[i];
    for (int i=0;i<3*numOfAtoms3;i++)
        xyzAugmented3[3*(natom12)+i] = xyz3[i];


    ICoord icp;
    icp.init(natom12+numOfAtoms3,atomicNamesCombined,atomicNumbersCombined,xyzAugmented3);

    //print_xyz();

    //rotate about central axis
    if (nadd1>1)
    {
        int atom1 = add1[0];
        int atom2 = add1[1];
        int tmp1 = atom1; if (atom1>atom2) { atom1 = atom2; atom2 = tmp1; }
        int a3 = add1[2];
        int a4 = add1[3];
        int tmp2 = a3; if (a3>a4) { a3 = a4; a4 = tmp2; }

        int same_frag = check_frag_3(atom1,atom2) + check_frag_3(a3,a4);
        if (!same_frag && atom1!=a3 && atom2!=a4)
        {
            double torv = icp.torsion_val(atom1,a3,a4,atom2);
            //printf("  central torv: %4.2f (%i %i %i %i) \n",torv,atom1+1,a3+1,a4+1,atom2+1);
            rotate_around_x(natom12,numOfAtoms3,-torv,icp.coords);
        }
    }

    double* vx = new double[3];
    vx[1] = vx[2] = 0.;
    vx[0] = 1.0;
    vdw_vector_opt(natom12,numOfAtoms3,vx,icp);
    delete [] vx;

    for (int i=0;i<3*(natom12+numOfAtoms3);i++)
        xyzAugmented3[i] = icp.coords[i];

    for (int i=0;i<3*numOfAtoms1;i++)
        xyz1[i] = xyzAugmented3[i];
    for (int i=0;i<3*numOfAtoms2;i++)
        xyz2[i] = xyzAugmented3[3*numOfAtoms1+i];
    for (int i=0;i<3*numOfAtoms3;i++)
        xyz3[i] = xyzAugmented3[3*(natom12)+i];

    //print_xyz_gen(natom12,icp.atomicNames,xyzAugmented);

    int badgeom = 0;
    for (int i=0;i<3*numOfAtoms1;i++)
        if (xyz1[i]!=xyz1[i])
        {
            badgeom = 1;
            break;
        }
    for (int i=0;i<3*numOfAtoms2;i++)
        if (xyz2[i]!=xyz2[i])
        {
            badgeom = 1;
            break;
        }
    for (int i=0;i<3*numOfAtoms3;i++)
        if (xyz3[i]!=xyz3[i])
        {
            badgeom = 1;
            break;
        }

    if (badgeom)
    {
        printf(" ERROR: nan detected in shuttle_align, exiting! \n");
        //print_xyz_gen(numOfAtoms1,atomicNames1,xyz1); Mina
        //print_xyz_gen(numOfAtoms2,atomicNames2,xyz2); Mina
        //print_xyz_gen(numOfAtoms3,atomicNames3,xyz3); Mina
        exit(1);
    }

    delete [] atomicNames12;
    delete [] anumbers12;

    delete [] atomicNames12a;
    delete [] atomicNames3a;
    delete [] xyz12a;
    delete [] xyz3a;

    delete [] v2;
    delete [] averagedBondVector1;
    delete [] averagedBondVector2;
    delete [] bonded1;
    delete [] bonded2;
    delete [] bonded3;

    ic1.freemem();
    ic2.freemem();

    delete [] atomicNamesCombined;
    delete [] atomicNumbersCombined;
    icp.freemem();


    return;
}
*/

void Align::print_xyz()
{
    printf(" %2i \n\n",numOfAtoms1+numOfAtoms2);
    for (int i=0;i<numOfAtoms1;i++)
        printf(" %s %9.6f %9.6f %9.6f \n",atomicNames1[i].c_str(),xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
    for (int i=0;i<numOfAtoms2;i++)
        printf(" %s %9.6f %9.6f %9.6f \n",atomicNames2[i].c_str(),xyz2[3*i+0],xyz2[3*i+1],xyz2[3*i+2]);
    printf("\n");

    return;
}

/*
void Align::print_xyz_3()
{
    if (inited!=3)
    {
        printf(" WARNING: inited is not 3 in print_xyz_3 \n");
        return;
    }

    printf(" %2i \n\n",numOfAtoms1+numOfAtoms2+numOfAtoms3);
    for (int i=0;i<numOfAtoms1;i++)
        printf(" %s %9.6f %9.6f %9.6f \n",atomicNames1[i].c_str(),xyz1[3*i+0],xyz1[3*i+1],xyz1[3*i+2]);
    for (int i=0;i<numOfAtoms2;i++)
        printf(" %s %9.6f %9.6f %9.6f \n",atomicNames2[i].c_str(),xyz2[3*i+0],xyz2[3*i+1],xyz2[3*i+2]);
    for (int i=0;i<numOfAtoms3;i++)
        printf(" %s %9.6f %9.6f %9.6f \n",atomicNames3[i].c_str(),xyz3[3*i+0],xyz3[3*i+1],xyz3[3*i+2]);
    printf("\n");

    return;
}

void Align::freemem()
{
    //printf(" Align free: %i %i \n",numOfAtoms1,numOfAtoms2);

    if (inited==1 || inited==3)
    {
        if (atomicNames1!=NULL) delete [] atomicNames1;
        if (atomicNames2!=NULL) delete [] atomicNames2;
        if (anumbers1!=NULL) delete [] anumbers1;
        if (anumbers2!=NULL) delete [] anumbers2;
        if (xyz1!=NULL) delete [] xyz1;
        if (xyz2!=NULL) delete [] xyz2;
        if (xyzAugmented!=NULL) delete [] xyzAugmented;
    }
    if (inited==3)
    {
        if (atomicNames3!=NULL) delete [] atomicNames3;
        if (anumbers3!=NULL) delete [] anumbers3;
        if (xyz3!=NULL) delete [] xyz3;
        if (xyzAugmented3!=NULL) delete [] xyzAugmented3;
    }

    inited = 0;

    return;
}

*/

double Align::norm(double* x, int size)
{
  double val = 0.;  
  for (int i=0;i<size;i++)
    val += x[i]*x[i];
  val = sqrt(val);

  return val; 
}

/*void Align::get_rotation_matrix(double** rotMat, double* thetas)
{
//  double x=thetas[0]; double y=thetas[1]; double z=thetas[2];
  rotMat[0][0] = cos(y)*cos(z);
  rotMat[0][1] = -cos(y)*sin(z);
  rotMat[0][2] = sin(y);
  rotMat[1][0] = sin(x)*sin(y)*cos(z)+cos(x)*sin(z);
  rotMat[1][1] = -sin(x)*sin(y)*sin(z)+cos(x)*cos(z);
  rotMat[1][2] = -sin(x)*cos(y);
  rotMat[2][0] = -cos(x)*sin(y)*cos(z)+sin(x)*sin(z);
  rotMat[2][1] = cos(x)*sin(y)*sin(z)+sin(x)*cos(z);
  rotMat[2][2] = cos(x)*cos(y);//

  double x=thetas[0]; double y=thetas[1]; double z=thetas[2];
  rotMat[0][0] = cos(y)*cos(z);
  rotMat[0][1] = cos(z)*sin(x)*sin(y)-cos(x)*sin(z);
  rotMat[0][2] = cos(x)*cos(z)*sin(y)+sin(x)*sin(z);
  rotMat[1][0] = cos(y)*sin(z);
  rotMat[1][1] = cos(x)*cos(z)+sin(x)*sin(y)*sin(z);
  rotMat[1][2] = cos(x)*sin(y)*sin(z)-cos(z)*sin(x);
  rotMat[2][0] = -sin(y);
  rotMat[2][1] = cos(y)*sin(x);
  rotMat[2][2] = cos(x)*cos(y);
}*/

void Align::print_xyz_gen(int natoms, string* anames, double* coords)
{
   
 printf(" %i \n",natoms);
 printf("\n");
 for (int i=0;i<natoms;i++)
 {
     cout << "  " << anames[i];
     printf(" %f %f %f \n",coords[3*i+0],coords[3*i+1],coords[3*i+2]);
 }
// printf("\n");
   
}

bool Align::writeToFile(std::string &outFile)
{
    bool success = false;
    std::ofstream ofs;
    ofs.open(outFile.c_str());

    ofs << std::to_string(numOfAtoms1+numOfAtoms2) << "\n";
    ofs << "\n";
    ofs << std::fixed << std::setprecision(15);
    for (unsigned int i=0; i<numOfAtoms1; i++)
    {
        ofs << atomicNames1[i] << "         " << xyz1[3*i+0] << "           " 
            << xyz1[3*i+1] << "            " << xyz1[3*i+2] << "\n";
    }
    for (unsigned int i=0; i<numOfAtoms2; i++)
    {
        ofs << atomicNames2[i] << "         " << xyz2[3*i+0] << "           " 
            << xyz2[3*i+1] << "            " << xyz2[3*i+2] << "\n";
    }

    ofs.close();
    success = true;
    return (success);
}

//void Align::moveToOrigin(int nadd1, int* add1, int atomIndex2)
void Align::moveToOrigin(int atom2, int atomIndex2, int i, int nadd1)
{
    // move atom2 to the origin of Cartesians, then move its attached atoms
//    for (int i=0; i<nadd1; i++) 
//    { 
        //int atom1 = add1[2*i+0];
//        int atom2 = add1[2*i+1];
        double adsorbateCenter[3] = {};

        // TODO: with the assumption that the second atom in add move ALWAYS comes from the second fragment
        // == order is preserved
//        atom2 -= numOfAtoms1;
        adsorbateCenter[0] = xyz2[3*atom2+0];
        adsorbateCenter[1] = xyz2[3*atom2+1];
        adsorbateCenter[2] = xyz2[3*atom2+2];

        // assuming we never will have a case with more tham two adsorbate molecules/fragments
        //if (i == 0)
        if (i == 0)
        {    
            for (int i=0;i<numOfAtoms2;i++)
            {    
                for (int j=0;j<3;j++)
                {    
                    //xyz2Displaced[3*i+j] = xyz2[3*i+j] + displacement[j];
                    xyz2[3*i+j] -= adsorbateCenter[j];
                }    
            }    
        }    
        else if (atomIndex2 > 0 && nadd1 > 0) 
        {    
            for (int i=atomIndex2;i<numOfAtoms2;i++)
            {    
                for (int j=0;j<3;j++)
                {    
                    xyz2[3*i+j] -= adsorbateCenter[j];
                }    
            }    
        }
//    }
}

void Align::moveToBindingSite(int atom1, int atom2, int atomIndex2, int i, int nadd1)
{
    double surfaceCenter[3], adsorbateCenter[3], displacement[3] = {0.0, 0.0, 0.0};
    surfaceCenter[0] = xyz1[3*atom1+0];
    surfaceCenter[1] = xyz1[3*atom1+1];
    surfaceCenter[2] = xyz1[3*atom1+2];

    // TODO: with the assumption that the second atom in add move ALWAYS comes from the second fragment
    // == order is preserved
    adsorbateCenter[0] = xyz2[3*atom2+0];
    adsorbateCenter[1] = xyz2[3*atom2+1];
    adsorbateCenter[2] = xyz2[3*atom2+2];

    const double DELTA_Z = 0.0;
    displacement[0] = surfaceCenter[0] - adsorbateCenter[0];
    displacement[1] = surfaceCenter[1] - adsorbateCenter[1];
    displacement[2] = surfaceCenter[2] - adsorbateCenter[2] + DELTA_Z;
    // move to binding site
    // assuming we never will have a case with more tham two adsorbate molecules/fragments
    if (i == 0)
    {
        for (int i=0;i<numOfAtoms2;i++)
        {
            for (int j=0;j<3;j++)
            {
                //xyz2Displaced[3*i+j] = xyz2[3*i+j] + displacement[j];
                xyz2[3*i+j] += displacement[j];
            }
        }
    }
    else if (atomIndex2 > 0 && nadd1 > 0)
    {
        int index = atomIndex2;
        for (int i=index;i<numOfAtoms2;i++)
        {
            for (int j=0;j<3;j++)
            {
                xyz2[3*i+j] += displacement[j];
            }
        }
    }
}
