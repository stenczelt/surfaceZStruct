// Internal coordinate class for surface reactions
#include "icoord.h"
#include "utils.h"
using namespace std;

#define MAX_FRAG_DIST 12.0

int ICoord::init(string xyzfile){


// printf(" xyzfile: %s \n",xyzfile);
 printf("\n");
 cout << " xyzfile: " << xyzfile << endl;
 structure_read(xyzfile);
 
 //print_xyz();

 alloc_mem();
 // printf(" done allocating memory\n");

 int done = ic_create();

 printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}



// initialize by feeding in xyz coordinates
int ICoord::init(int nat, string* anam, int* anum, double* xyz){

// printf(" initializing icoord via xyz structure \n");
 natoms = nat;
// printf(" natoms: %i \n",nat);
// for (int i=0;i<natoms;i++)
//    printf(" %1.3f %1.3f %1.3f \n",xyz[3*i+0],xyz[3*i+1],xyz[3*i+2]);

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

 alloc_mem();

 int done = ic_create();

// printf(" initializing MM parameters \n");
 mm_init();

 //printf("\n\n");

 return 1;
}

// initialize memory only
int ICoord::alloc(int size){

 natoms = size;

//otherwise allocated in structure_read
 anumbers = new int[1+natoms];
 amasses = new double[1+natoms];
 anames = new string[1+natoms];
 coords = new double[natoms*3];
 coordsts = new double[natoms*3];
 coords0 = new double[natoms*3];

 alloc_mem();

 return 1;
}


// initialize by feeding in xyz coordinates
int ICoord::reset(double* xyz){

// printf(" resetting icoord via xyz structure \n");

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

// initialize by feeding in xyz coordinates
int ICoord::reset(int nat, string* anam, int* anum, double* xyz){

// printf(" resetting icoord via xyz structure \n");
// printf(" in reset natoms: %i natoms_new: %i \n",natoms,nat);
 natoms = nat;

 for (int i=0;i<natoms;i++)
   anumbers[i] = anum[i];
 for (int i=0;i<natoms;i++)
   anames[i] = anam[i];

 for (int i=0;i<3*natoms;i++)
   coords[i]=xyz[i];
 for (int i=0;i<3*natoms;i++)
   coords0[i]=xyz[i];

// printf("\n");
// print_xyz();

// int done = ic_create();

// printf(" initializing MM parameters \n");
// mm_init();

 //printf("\n\n");

 return 1;
}

void ICoord::update_ic(){

  update_bonds();
  update_angles();
  update_torsion();
  update_imptor();
  update_nonbond();

  return;
} 
 
void ICoord::create_xyz()
{
  printf ("xyz_create not implemented\n");
  int* adone = new int[natoms];
  for (int i=0;i<natoms;i++) adone[i]=0;

  delete [] adone;

  return;
}


int ICoord::ic_create()
{
  make_bonds();

  if (isOpt)
  {
    printf(" isOpt: %i \n",isOpt);
    make_frags();
    bond_frags();
  }

  coord_num(); // counts # surrounding species
  if (use_xyz) get_xyzic();

  make_angles();
  make_torsions();
  make_imptor();

  if (isOpt)
    linear_ties();

  n_nonbond = make_nonbond(); //anything not connected by bond or angle

  update_ic();

  return 0;
}

int ICoord::ic_create_nobonds()
{
  if (isOpt)
  {
    printf(" isOpt: %i \n",isOpt);
    make_frags();
    bond_frags();
  }

  coord_num(); // counts # surrounding species
  make_angles();
  make_torsions();
  make_imptor();

  if (isOpt)
    linear_ties();

  n_nonbond = make_nonbond(); //anything not connected by bond or angle

  update_ic();

  return 0;
}

void ICoord::linear_ties() 
{
  printf("\n skipping linear_ties() \n");
  return;
 //Note: this works well for TM type multiple linear at same center

  update_ic();

  int maxsize = 0;
  for (int i=0;i<nangles;i++)
  if (anglev[i]>160.)
    maxsize++;

 //list of linear angles
  int* blist = new int[3*maxsize+1];
  for (int i=0;i<3*maxsize;i++) blist[i] = -1;

  int n = 0;
  for (int i=0;i<nangles;i++)
  if (anglev[i]>160.)
  {
    blist[3*n] = angles[i][0];
    blist[3*n+1] = angles[i][2];
    blist[3*n+2] = angles[i][1];
    printf(" linear angle (%i of %i): %i %i %i (%4.2f) \n",n+1,maxsize,blist[3*n],blist[3*n+2],blist[3*n+1],anglev[i]);
    n++;
  }

 //atoms attached to linear angle atoms
  int** clist = new int*[n];
  for (int i=0;i<n;i++)
    clist[i] = new int[24];
  for (int i=0;i<n;i++) 
  for (int j=0;j<24;j++)
    clist[i][j] = -1;

  int* m = new int[n];
  for (int i=0;i<n;i++) m[i] = 0;

  for (int i=0;i<n;i++)
  for (int k=0;k<nbonds;k++)
  {
    if (bonds[k][0]==blist[3*i] && bonds[k][1]!=blist[3*i+2])
      clist[i][m[i]++] = bonds[k][1];
    else if (bonds[k][1]==blist[3*i] && bonds[k][0]!=blist[3*i+2])
      clist[i][m[i]++] = bonds[k][0];

    if (bonds[k][0]==blist[3*i+1] && bonds[k][1]!=blist[3*i+2])
      clist[i][m[i]++] = bonds[k][1];
    else if (bonds[k][1]==blist[3*i+1] && bonds[k][0]!=blist[3*i+2])
      clist[i][m[i]++] = bonds[k][0];
  }

 //cross linking (torsion type)
  for (int i=0;i<n;i++)
  {
    printf(" linear angle: %i %i %i \n",blist[3*i+0],blist[3*i+2],blist[3*i+1]);
 
    int a1 = blist[3*i+0];
    int a2 = blist[3*i+1];
    if (!bond_exists(a1,a2))
    {
      printf(" adding bond via linear_ties(t): %i %i \n",a1,a2);
      bonds[nbonds][0] = a1;
      bonds[nbonds++][1] = a2;
    } //if !bond_exists

    int found;
    int b1,b2;
    for (int j=0;j<m[i];j++)
    for (int k=0;k<j;k++) 
    {
      int b1 = clist[i][j];
      int b2 = clist[i][k];
      found = 0;

      for (int l=0;l<nangles;l++)
      {
        if (b1==angles[l][0] && b2==angles[l][2])
          found = 1;
        else if (b2==angles[l][0] && b1==angles[l][2])
          found = 1;
      }
      if (!found)
      {
        int c1,c2;
        if (bond_exists(b1,a1))
          c1 = b1;
        if (bond_exists(b2,a1))
          c1 = b2;
        if (bond_exists(b1,a2))
          c2 = b1;
        if (bond_exists(b2,a2))
          c2 = b2;
        printf(" adding torsion via linear_ties(t): %i %i %i %i \n",c1,a1,a2,c2);
        torsions[ntor][0] = c1;
        torsions[ntor][1] = a1;
        torsions[ntor][2] = a2;
        torsions[ntor][3] = c2;
        ntor++;
      }
    }
      
  } //loop i over clist 

#if 0
 //cross linking (TM-type)
  for (int i=0;i<n;i++)
  for (int j=0;j<i;j++)
  {
    int a1,a2;
    double mclose = 1000.;
    int found = 0;
    for (int k=0;k<m[i];k++)
    for (int l=0;l<m[j];l++)
    if (!bond_exists(clist[i][k],clist[j][l]) && clist[i][k]!=clist[j][l])
    {
      double close = distance(clist[i][k],clist[j][l]);
      if (close<mclose)
      {
        a1 = clist[i][k];
        a2 = clist[j][l];
        mclose = close;
        found = 1;
      }
    }
    if (found)
    {
      found = 0;
      for (int j=0;j<nangles;j++)
      {
        if (a1==angles[j][0] && a2==angles[j][2])
          found = 1;
        else if (a2==angles[j][0] && a1==angles[j][2])
          found = 1;
      }
      if (!found)
      {
        printf(" adding bond via linear_ties(TM): %i %i \n",a1,a2);
        bonds[nbonds][0] = a1;
        bonds[nbonds++][1] = a2;
      }
    }
  }
#endif


#if 0
 //old cross linking
  for (int i=0;i<n;i++)
  for (int j=0;j<i;j++)
  for (int k=0;k<m[i];k++)
  for (int l=0;l<m[j];l++)
  if (!bond_exists(clist[i][k],clist[j][l]) && clist[i][k]!=clist[j][l])
  {
    printf(" adding bond via linear_ties: %i %i \n",clist[i][k],clist[j][l]);
    bonds[nbonds][0] = clist[i][k];
    bonds[nbonds++][1] = clist[j][l];
  }
#endif

#if 0
  for (int i=0;i<n;i++)
  for (int j=0;j<i;j++)
  if (!bond_exists(blist[j],blist[k]))
  {
    printf(" adding bond via linear_ties: %i %i \n",blist[j],blist[k]);
    bonds[nbonds][0] = blist[j];
    bonds[nbonds++][1] = blist[k];
  }
#endif
//  printf(" linear ties about to deallocate \n"); fflush(stdout);

  delete [] m; 
//problem deallocating these
  for (int i=0;i<n;i++)
    delete [] clist[i];
  delete [] clist;

  delete [] blist;

  printf(" linear ties complete \n"); fflush(stdout);

  return;
}

void ICoord::tm_ties()
{
  int maxsize = 1;
  for (int i=0;i<natoms;i++)
  if (coordn[i]>=4 && isTM(i))
    maxsize += coordn[i];

  int* blist = new int[maxsize];
  int nfound = 0;

  for (int i=0;i<natoms;i++)
  if (coordn[i]>=4 && isTM(i))
  {
    //add new bonds via connections
    for (int j=0;j<nbonds;j++)
    {
      if (bonds[j][0]==i)
        blist[nfound++] = bonds[j][1];
      else if (bonds[j][1]==i)
        blist[nfound++] = bonds[j][0];
    }
    for (int j=0;j<nfound;j++)
    for (int k=0;k<j;k++)
    if (!bond_exists(blist[j],blist[k]))
    {
     //adding linear angle connections only 
     //didn't help, neither did all connections between TM-L's
      if (angle_val(blist[j],i,blist[k])>160.)
      {
        printf(" adding bond via tm_ties: %i %i \n",blist[j],blist[k]);
        bonds[nbonds][0] = blist[j];
        bonds[nbonds++][1] = blist[k];
      }
    }
  }

  delete [] blist;

  return;
}

void ICoord::make_frags() 
{
  //printf(" in make_frags() \n");

  frags = new int[natoms];
  nfrags = 0;
  int* atominfrag = new int[natoms];
  for (int i=0;i<natoms;i++) atominfrag[i] = 0;

  for (int i=0;i<natoms;i++) frags[i] = -1;

  for (int i=0;i<nbonds;i++)
  {
    atominfrag[bonds[i][0]] = 1;
    atominfrag[bonds[i][1]] = 1;

    if ( frags[bonds[i][0]] == -1 && frags[bonds[i][1]] == -1 )
    {
      frags[bonds[i][0]] = nfrags;
      frags[bonds[i][1]] = nfrags;
      nfrags++;
    }
    else if ( frags[bonds[i][0]] == -1 && frags[bonds[i][1]] > -1 )
    {
      frags[bonds[i][0]] = frags[bonds[i][1]];
    }
    else if ( frags[bonds[i][0]] > -1 && frags[bonds[i][1]] == -1 )
    {
      frags[bonds[i][1]] = frags[bonds[i][0]];
    }
    else if ( frags[bonds[i][0]] > -1 && frags[bonds[i][1]] > -1 )
    {
      //both frags assigned already
      if (frags[bonds[i][0]] != frags[bonds[i][1]])
      {
        //printf(" WARNING, need to merge: %i %i \n",frags[bonds[i][0]],frags[bonds[i][1]]);
        //mergelist[2*nmerge+0] = frags[bonds[i][0]];
        //mergelist[2*nmerge+1] = frags[bonds[i][1]];
        int f1 = min(frags[bonds[i][0]],frags[bonds[i][1]]);
        int f2 = max(frags[bonds[i][0]],frags[bonds[i][1]]);
        printf(" merging: %i/%i \n",f1,f2);
        for (int j=0;j<natoms;j++)
        if (frags[j]==f2)
          frags[j] = f1;
        if (f2==nfrags-1)
          nfrags--;
      }
    }
  } //loop i over nbonds

  int foundxyz = 0;
  for (int i=0;i<natoms;i++)
  if (!atominfrag[i])
  {
    frags[i] = nfrags;
    foundxyz = 1;
  }
  if (foundxyz) nfrags++;

  for (int i=0;i<natoms;i++)
  if (frags[i]==-1)
    frags[i] = nfrags++;

  for (int i=0;i<natoms;i++)
    printf(" atom[%i] frag: %i \n",i+1,frags[i]+1);

  return;
}


void ICoord::bond_frags_xyz()
{
  printf(" in bond_frags() \n");
  if (nfrags<2) return;

  int found = 0;
  int found2 = 0;

  int a1,a2;
  int b1,b2;
  double mclose;
  double mclose2;
  for (int n1=0;n1<nfrags;n1++)
  for (int n2=0;n2<n1;n2++)
  {
    if (natoms<150)
      printf(" connecting frag %i to frag %i: ",n1+1,n2+1);

    found = 0;
    found2 = 0;
    double close = 0.;
    mclose = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      if (close<mclose)
      {
        mclose = close;
        a1 = i;
        a2 = j;
        found = 1;
      }
    }

   //now connect a1 to next nearest
    mclose2 = 1000.;
    for (int i=0;i<natoms;i++)
    if (frags[i]==n1 && i!=a1 && i!=a2)
    {
      close = distance(i,a1);
      if (close<mclose2)
      {
        mclose2 = close;
        b1 = a1;
        b2 = i;
        found2 = 1;
      }
    } 

    if (found && !bond_exists(a1,a2))
    {
      printf(" bond pair1 added : %i %i ",a1+1,a2+1);
      bonds[nbonds][0] = a1;
      bonds[nbonds][1] = a2;
      bondd[nbonds] = mclose;
      nbonds++;
    } // if found

    if (found2 && !bond_exists(b1,b2))
    {
      printf(" bond pair2 added : %i %i ",b1+1,b2+1);
      bonds[nbonds][0] = b1;
      bonds[nbonds][1] = b2;
      bondd[nbonds] = mclose2;
      nbonds++;
    } // if found2

  } //loop n2<n1 over nfrags

  return;
}


void ICoord::bond_frags() 
{
  printf(" in bond_frags() \n");
  if (nfrags<2) return;

  int found = 0;
  int found2 = 0;
  int found3 = 0;

  int a1,a2;
  int b1,b2;
  int c1,c2;
  double mclose;
  double mclose2;
  double mclose3;
  for (int n1=0;n1<nfrags;n1++)
  for (int n2=0;n2<n1;n2++)
  {
    if (natoms<50)
      printf(" connecting frag %i to frag %i: ",n1,n2);

    found = 0;
    found2 = 0;
    found3 = 0;
    double close = 0.;
    mclose = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      if (close<mclose && close < MAX_FRAG_DIST)
      {
        mclose = close;
        a1 = i;
        a2 = j;
        found = 1;
      }
    }
   
   //connect second pair, heavies or H-bond only, away from 1st pair
    b1 = -1;
    b2 = -1;
    mclose2 = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      double dia1 = distance(i,a1);
      double dja1 = distance(j,a1);
      double dia2 = distance(i,a2);
      double dja2 = distance(j,a2);
      double dist21 = (dia1+dja1)/2.;
      double dist22 = (dia2+dja2)/2.;
      if (anumbers[i] > 1 || anumbers[j] > 1)
      if (dist21 > 4.5 && dist22 > 4.5) //standard
//      if (dia1 > 4.5 && dja1 > 4.5 && dia2 > 4.5 && dja2 > 4.5) //possible change
      if (close<mclose2 && close < MAX_FRAG_DIST)
      {
        mclose2 = close;
        b1 = i;
        b2 = j;
        found2 = 1;
      }
    }

   //connect third pair, heavies or H-bond only, away from 1st pair
    c1 = -1;
    c2 = -1;
    mclose3 = 1000.;
    for (int i=0;i<natoms;i++)
    for (int j=0;j<natoms;j++)
    if (frags[i]==n1 && frags[j]==n2)
    {
      close = distance(i,j);
      double dia1 = distance(i,a1);
      double dja1 = distance(j,a1);
      double dia2 = distance(i,a2);
      double dja2 = distance(j,a2);
      double dib1 = distance(i,b1);
      double djb1 = distance(j,b1);
      double dib2 = distance(i,b2);
      double djb2 = distance(j,b2);
      double dist31 = (dia1+dja1)/2.;
      double dist32 = (dia2+dja2)/2.;
      double dist33 = (dib1+djb1)/2.;
      double dist34 = (dib2+djb2)/2.;
      if (anumbers[i] > 1 || anumbers[j] > 1)
      if (dist31 > 4.5 && dist32 > 4.5 && dist33 > 4.5 && dist34 > 4.5) //standard
//      if (dia1 > 4.5 && dja1 > 4.5 && dia2 > 4.5 && dja2 > 4.5)
//      if (dib1 > 4.5 && djb1 > 4.5 && dib2 > 4.5 && djb2 > 4.5)
      if (close<mclose3 && close < MAX_FRAG_DIST)
      {
        mclose3 = close;
        c1 = i;
        c2 = j;
        found3 = 1;
      }
    }

    if (found && !bond_exists(a1,a2))
    {
      printf(" bond pair1 added : %i %i ",a1,a2);
      bonds[nbonds][0] = a1;
      bonds[nbonds][1] = a2;
      bondd[nbonds] = mclose;
      nbonds++;
    } // if found

    if (found2 && !bond_exists(b1,b2))
    {
      printf(" bond pair2 added : %i %i ",b1,b2);
      bonds[nbonds][0] = b1;
      bonds[nbonds][1] = b2;
      bondd[nbonds] = mclose2;
      nbonds++;
    } // if found

    if (found3 && !bond_exists(c1,c2))
    {
      printf(" bond pair3 added : %i %i ",c1,c2);
      bonds[nbonds][0] = c1;
      bonds[nbonds][1] = c2;
      bondd[nbonds] = mclose3;
      nbonds++;
    } // if found

    if (found || natoms<50)
      printf("\n");
  }//loop n over nfrags



  return;
}

void ICoord::connect_1_coord_mg()
{
  for (int i=0;i<natoms;i++)
  if (!xyzic[i] && coordn[i]<2 && frozen[i])
  {
    printf("   fixing stranded single frozen atom \n");
    xyzic[i] = 1;
    nxyzic += 3;
  }

  for (int i=0;i<natoms;i++)
  if (!isTM(i) && coordn[i]<2 && !frozen[i] && anumbers[i]!=1 && anumbers[i]!=-1) 
  {
    double d1 = 1000.;
    int wa = -1;
    for (int j=0;j<natoms;j++)
    if (i!=j)
    {
      double d0 = distance(i,j);
      if (d0<d1 && !bond_exists(i,j))
      {
        d1 = d0;
        wa = j;
      }
    }
    if (wa>-1)
    {
      printf("   adding bond to correct for low coordination: %i %i \n",i+1,wa+1);
      bonds[nbonds][0] = i;
      bonds[nbonds][1] = wa; 
      bondd[nbonds] = d1;
      nbonds++;
      coord_num();
    }
  }

  return;
}

//CPMZ tune me
void ICoord::get_xyzic()
{
  printf("  in get_xyzic() \n");

  nxyzic = 0;
  for (int i=0;i<natoms;i++)
  {
      std::cout << "i: " << i << std::endl;
      xyzic[i] = 0;
  }

  for (int i=0;i<nbonds;i++)
  {
    int a1 = bonds[i][0];
    int a2 = bonds[i][1];
    int frzpair = 0;
    if (frozen!=NULL)
    if (frozen[a1] && frozen[a2])
      frzpair = 1;

//    if (coordn[a1] + coordn[a2] > 8 || coordn[a1] > 6 || coordn[a2] > 6)
//    if (coordn[a1]>4 && coordn[a2]>4)
//    if ((coordn[a1]>4 && anumbers[a2]>20) || (coordn[a2]>4 && anumbers[a1]>20))
//    if (isTM(a1) && isTM(a2) && (coordn[a1]>3 || coordn[a2]>3))
    if ((isTM(a1) && isTM(a2)) || frzpair)
    {
      //printf(" high coordn: %i/%i for atoms %i/%i \n",coordn[a1],coordn[a2],a1+1,a2+1);
      //printf(" metal-metal bond: %i-%i \n",a1+1,a2+1);
      for (int j=i;j<nbonds-1;j++)
      {
        bonds[j][0] = bonds[j+1][0];
        bonds[j][1] = bonds[j+1][1];
        bondd[j] = bondd[j+1];
      }
      nbonds--;
      i--;
    }
  }

  coord_num();

  for (int i=0;i<natoms;i++)
  if ((coordn[i]<2 && isTM(i)) || coordn[i]<1) //was 1
//  if (coordn[i]<1)
  {
    //printf(" low coordn (%i), setting xyz: %i \n",coordn[i],i);
    xyzic[i] = 1;
    nxyzic += 3;
  }

  if (nxyzic && isOpt)
  {
    connect_1_coord_mg(); //fix for oxides
    make_frags();
    bond_frags_xyz();
  }

#if 1
  if (isOpt || 0)
  {
    printf("\n XYZ vs. IC geometry \n");
    printf(" %i \n\n",natoms);
    for (int i=0;i<natoms;i++)
    if (xyzic[i])
      printf(" X %7.5f %7.5f %7.5f \n",coords[3*i+0],coords[3*i+1],coords[3*i+2]);
    else
      printf(" %s %7.5f %7.5f %7.5f \n",anames[i].c_str(),coords[3*i+0],coords[3*i+1],coords[3*i+2]);
    printf("\n");
  }
#endif

  return;
}


#if 0
int ICoord::ic_create_nobonds()
{
  coord_num(); // counts # surrounding species
  make_angles();
  make_torsions();
  make_imptor_nobonds();
  n_nonbond = make_nonbond(); //anything not connected by bond or angle
  return 0;
}
#endif

void ICoord::update_bonds(){  
  for (int i=0;i<nbonds;i++)
    bondd[i] = distance(bonds[i][0],bonds[i][1]);
  return;
}

void ICoord::update_angles(){
  for (int i=0;i<nangles;i++)
    anglev[i] = angle_val(angles[i][0],angles[i][1],angles[i][2]);
  return;
}

void ICoord::update_torsion(){
  for (int i=0;i<ntor;i++)
    torv[i]=torsion_val(torsions[i][0],torsions[i][1],torsions[i][2],torsions[i][3]);
  return;
}

void ICoord::update_imptor(){
  for (int i=0;i<nimptor;i++)
    imptorv[i]=torsion_val(imptor[i][0],imptor[i][1],imptor[i][2],imptor[i][3]);
  return;
}

void ICoord::update_nonbond(){
  for (int i=0;i<n_nonbond;i++)
    nonbondd[i] = distance(nonbond[i][0],nonbond[i][1]);
  return;
}

void ICoord::make_bonds()
{
//  printf("\n\n DEBUG: no bonds \n");
//  return;

  //printf(" in make_bonds, natoms: %i\n",natoms);
  double MAX_BOND_DIST; 
  nbonds = 0;
  nxyzic = 0;
  for (int i=0;i<natoms;i++)
    for (int j=0;j<i;j++)
    {
       MAX_BOND_DIST = (getR(i) + getR(j))/2;
       if (farBond>1.0) MAX_BOND_DIST *= farBond;
       double d = distance(i,j);
       if (d<MAX_BOND_DIST)
       {
          //printf(" found bond: %2i %2i dist: %f \n",i+1,j+1,d);
          bonds[nbonds][0]=i;
          bonds[nbonds][1]=j;
          bondd[nbonds]=d;
          nbonds++;
       }
    }

}

void ICoord::coord_num()
{ 
  for (int i=0;i<natoms;i++)
    coordn[i] = 0;
  for (int i=0;i<nbonds;i++)
  {
    coordn[bonds[i][0]]++;
    coordn[bonds[i][1]]++;
  }
}

void ICoord::make_angles()
{
  //include all consecutive connections 
  nangles=0;
  for (int i=0;i<nbonds;i++)
  {
     for (int j=0;j<i;j++)
     {
        if (bonds[i][0]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][0]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][0];
          angles[nangles][0]=bonds[i][1];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][0])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][1];
          nangles++;
        }
        else if (bonds[i][1]==bonds[j][1])
        {
          angles[nangles][1]=bonds[i][1];
          angles[nangles][0]=bonds[i][0];
          angles[nangles][2]=bonds[j][0];
          nangles++;
        }
        if (nangles>0)
          anglev[nangles-1]=angle_val(angles[nangles-1][0],angles[nangles-1][1],angles[nangles-1][2]);
     } //loop j
  } //loop i


  return;
}

void ICoord::make_torsions()
{
  int a1,b1,c1,a2,b2,c2;
  bool found;

  ntor = 0;

//  return;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       b1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       b2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,b1,c1,a2,b2,c2); fflush(stdout);

       if (b1==c2 && b2==c1)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (b1==a2 && b2==c1)
       {
          torsions[ntor][0]=a1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       else if (b1==c2 && b2==a1)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=a2;
          ntor++; found=true;
       }
       else if (b1==a2 && b2==a1)
       {
          torsions[ntor][0]=c1;
          torsions[ntor][1]=b1;
          torsions[ntor][2]=b2;
          torsions[ntor][3]=c2;
          ntor++; found=true;
       }
       if (found && torsions[ntor-1][0] == torsions[ntor-1][2]) { found = false; ntor--; }
       //if (found && bond_exists(torsions[ntor-1][0],torsions[ntor-1][3]))
       if (found && torsions[ntor-1][0] == torsions[ntor-1][3]) { found = false; ntor--; }
       if (found)
       {
         //printf(" made tor: %i %i %i %i \n",torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
         torv[ntor-1]=torsion_val(torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
       }
       //if (found && fabs(torv[ntor-1])>180 ) ntor--;
    }
  } 

  return;
}

void ICoord::make_imptor()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  nimptor = 0;
  double imptorvt;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2)
       {
         if (a1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
         else if (a1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
       } // if m1==m2
       if (found)
       {
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
         if (fabs(imptorvt) > 12.0 && fabs(imptorvt - 180) > 12.0 ) { found = false; nimptor--; }
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

#if 0
  //attempted patch for 4 atom system
  if (natoms==4)
  {
    torsions[ntor][0] = imptor[nimptor-1][0];
    torsions[ntor][1] = imptor[nimptor-1][1];
    torsions[ntor][2] = imptor[nimptor-1][2];
    torsions[ntor][3] = imptor[nimptor-1][3];
    ntor++;
    torv[ntor-1]=torsion_val(torsions[ntor-1][0],torsions[ntor-1][1],torsions[ntor-1][2],torsions[ntor-1][3]);
  }
#endif

  return;
}

void ICoord::make_imptor_nobonds()
{
  int a1,m1,c1,a2,m2,c2;
  bool found;
  nimptor = 0;
  double imptorvt;

  for (int i=0;i<nangles;i++)
  {
    for (int j=0;j<i;j++)
    {
       found = false;
       a1=angles[i][0];
       m1=angles[i][1];
       c1=angles[i][2];
       a2=angles[j][0];
       m2=angles[j][1];
       c2=angles[j][2];

       //printf(" angle1: %i %i %i angle2: %i %i %i \n",a1,m1,c1,a2,m2,c2);

       if (m1==m2)
       {
         if (a1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
         else if (a1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==c2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=a2;
           nimptor++; found=true;
         }
         else if (c1==a2)
         {
           imptor[nimptor][0]=c1;
           imptor[nimptor][1]=a1;
           imptor[nimptor][2]=m1;
           imptor[nimptor][3]=c2;
           nimptor++; found=true;
         }
       } // if m1==m2
       if (found)
       {
//FIX ME
//the following only works when center is 3 coordinate
         for (int k=0;k<nimptor-1;k++)
           if (imptor[k][2] == m1)
             { found = false; nimptor--; }
       }
       if (found)
       {
         imptorvt = torsion_val(imptor[nimptor-1][0],imptor[nimptor-1][1],imptor[nimptor-1][2],imptor[nimptor-1][3]);
//         printf(" imptorv[%i]: %1.4f \n",nimptor,imptorvt);
//       make all 3 centered atoms planar?
//         printf(" atom: %i has coordn %i \n",imptor[nimptor-1][2],coordn[imptor[nimptor-1][2]]);
         if (coordn[imptor[nimptor-1][2]]!=3)
         {
           found = false;
           nimptor--;
         }
       }
       if (found) imptorv[nimptor-1] = imptorvt;
    }
  } 

  return;
}

double ICoord::torsion_val(int i, int j, int k, int l)
{
  double tval = -999;

  double x1 = coords[3*j+0] - coords[3*i+0];
  double y1 = coords[3*j+1] - coords[3*i+1];
  double z1 = coords[3*j+2] - coords[3*i+2];
  double x2 = coords[3*k+0] - coords[3*j+0];
  double y2 = coords[3*k+1] - coords[3*j+1];
  double z2 = coords[3*k+2] - coords[3*j+2];
  
  double ux1 = y1*z2-z1*y2;
  double uy1 = z1*x2-x1*z2;
  double uz1 = x1*y2-y1*x2;

  double x3 = coords[3*l+0] - coords[3*k+0];
  double y3 = coords[3*l+1] - coords[3*k+1];
  double z3 = coords[3*l+2] - coords[3*k+2];

  double ux2 = z3*y2 - y3*z2;
  double uy2 = x3*z2 - z3*x2;
  double uz2 = y3*x2 - x3*y2;

  double u = (ux1*ux1+uy1*uy1+uz1*uz1)*(ux2*ux2+uy2*uy2+uz2*uz2);

  if (u!=0.0)
  {
     double a = (ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u);
     if (a>1) a=1; else if (a<-1) a=-1;
     tval = acos(a);
     if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
         uz1*(ux2*y2-uy2*x2) < 0.0) tval *=-1;
  }
  else
    tval = 0.;

  if (tval>3.14159) tval-=2*3.14159;
  if (tval<-3.14159) tval+=2*3.14159;

  return tval * 180/3.14159;
}

double ICoord::angle_val(int i, int j, int k)
{
   double D1 = distance(i,j);
   double D2 = distance(j,k);
   double D3 = distance(i,k);
   
   double cos = ( D1*D1 + D2*D2 - D3*D3 ) / ( 2*D1*D2);
 
   if (cos > 1) cos = 1;
   if (cos < -1) cos = -1;

   //printf(" atoms: %i %i %i cos is: %f \n",i,j,k,cos);
 
   return acos(cos) * 180/3.14159;
}

int ICoord::make_nonbond(){

  int n = 0;
  for (int i=0;i<natoms;i++)
  {
    for (int j=0;j<i;j++)
    {
      bool found = false;
      for (int k=0;k<nbonds;k++)
      {
         if (found) break;
         if ((bonds[k][0]==i && bonds[k][1]==j) ||
             (bonds[k][0]==j && bonds[k][1]==i)) found = true;
      }
      //printf(" checking for pair: %i %i \n",i,j);
      for (int k=0;k<nangles;k++)
      {
        if (found) break;
        //printf(" angle %i bonds: %i %i %i \n",k,angles[k][0],angles[k][1],angles[k][2]);
        if (angles[k][0]==i)
        {
           if (angles[k][1]==j) found = true;
           else if (angles[k][2]==j) found = true;
        }
        else if (angles[k][1]==i)
        {
           if (angles[k][0]==j) found = true;
           else if (angles[k][2]==j) found = true;
        }
        else if (angles[k][2]==i)
        {
           if (angles[k][0]==j) found = true;
           else if (angles[k][1]==j) found = true;
        }
      } // loop k over angles
      if (!found)
      {
        //printf(" not found\n");
        nonbondd[n] = distance(i,j);
        nonbond[n][0] = i;
        nonbond[n][1] = j;
        n++;
      }
    }
  }
  //printf(" n_nonbond: %i \n",n);

  return n;
}


int ICoord::isTM(int a) {

//may later be extended to all 5+ coord types
  int anum;
  if (a>-2)
    anum = anumbers[a];
  else
    return 0;

  int TM = 0;
  if (anum > 1000)
    TM = 2;
  else if (anum > 20)
  {
    if (anum < 31 || anum==-1)
      TM = 1;
    else if (38 < anum && anum < 49)
      TM = 1;
    else if (71 < anum && anum < 81)
      TM = 1;
  }
     
  return TM;
}


double ICoord::getR(int i){

  double value;
 
  int an = anumbers[i];
  if      (an==0) value = 1.0;
  else if (an==1) value = 1.3;
  else if (an==3) value = 2.65; //PT
  else if (an==4) value = 2.0; //PT
  else if (an==5) value = 1.75;
  else if (an==6) value = 1.65;
  else if (an==7) value = 1.65;
  else if (an==8) value = 1.65;
  else if (an==9) value = 1.6;
  else if (an==11) value = 3.3; //PT
  else if (an==12) value = 3.1;
  else if (an==13) value = 2.6;
  else if (an==14) value = 2.6;
  else if (an==15) value = 2.5;
  else if (an==16) value = 2.45;
  else if (an==17) value = 2.1;
  else if (an==19) value = 4.0; //PT
  else if (an==20) value = 3.5; //PT
  else if (an==21) value = 4.5;
  else if (an==22) value = 4.2;
  else if (an==23) value = 4.0;
  else if (an==24) value = 3.5;
  else if (an==25) value = 3.4;
  else if (an==26) value = 3.3;
  else if (an==27) value = 3.0;
  else if (an==28) value = 3.0;
  else if (an==29) value = 3.0;
  else if (an==35) value = 2.7;
  else if (an==44) value = 3.2;
  else if (an==45) value = 3.15;
  else if (an==46) value = 3.15;
  else if (an==53) value = 2.8; //iodine
  else if (an==73) value = 3.4;
  else if (an==74) value = 3.4;
  else if (an==77) value = 3.35;
  else if (an==78) value = 3.35;
  else if (an==79) value = 3.35;
  else if (an==-1) value = 0.; 
  else 
  {
    printf(" Need to add atomic number %i to getR! (atom %2i) \n",an,i+1);
    //int aa = anumbers[-2];
    exit(1);
  }
//  else value = 0.;

  return value;
}

double ICoord::distance(int i, int j)
{
  //printf("in distance: %i %i\n",i+1,j+1);
  return sqrt((coords[3*i+0]-coords[3*j+0])*(coords[3*i+0]-coords[3*j+0])+
              (coords[3*i+1]-coords[3*j+1])*(coords[3*i+1]-coords[3*j+1])+
              (coords[3*i+2]-coords[3*j+2])*(coords[3*i+2]-coords[3*j+2])); 
}

int ICoord::bond_exists(int b1, int b2) {

   int found = 0;
   if (bond_num(b1,b2)>-1)
     found = 1;
   return found;
}

int ICoord::bond_num(int b1, int b2) {

   int found = -1;

   for (int k1=0;k1<nbonds;k1++)
     if ( (bonds[k1][0] == b1 && bonds[k1][1] == b2)
       || (bonds[k1][1] == b1 && bonds[k1][0] == b2))
     {
       found = k1;
       break;
     }

   return found;
}

int ICoord::angle_num(int b1, int b2, int b3)
{

   int found = -1;

   for (int k1=0;k1<nangles;k1++)
   if (angles[k1][1] == b2)
   {
     if (angles[k1][0] == b1 && angles[k1][2]==b3)
     {
       found = k1;
       break;
     }
     else if (angles[k1][2] == b1 && angles[k1][0]==b3)
     {
       found = k1;
       break;
     }
   }

   return found;
}

int ICoord::tor_num(int b1, int b2, int b3, int b4)
{
   int found = -1;

   for (int k1=0;k1<ntor;k1++)
   {
     if (torsions[k1][1] == b2 && torsions[k1][2] == b3)
     {
       if (torsions[k1][0] == b1 && torsions[k1][3] == b4)
       {
         found = k1;
         break;
       }
     }
     if (torsions[k1][2] == b2 && torsions[k1][1] == b3)
     {
       if (torsions[k1][3] == b1 && torsions[k1][0] == b4)
       {
         found = k1;
         break;
       }
     }
   }

   return found;
}

int ICoord::hpair(int a1, int a2) {
  if (anumbers[a1]==1 && anumbers[a2]==1)
    return 1;
  else
    return 0;
}

int ICoord::h2count() {

  int count = 0;
  for (int i=0;i<nbonds;i++)
  {
    if (anumbers[bonds[i][0]]==1 && anumbers[bonds[i][1]]==1)
      count++;
  }

  return count;
}

double ICoord::close_bond(int i)
{
  double A = 0.2;
  double val = 0.;
  int a1 = bonds[i][0];
  int a2 = bonds[i][1];
  double d = distance(a1,a2);
  double dr = (getR(a1) + getR(a2))/2;
//  val = exp(-(d*d-dr*dr));
  val = exp(-A*(d-dr));
  if (val > 1) val = 1;
 // printf(" bond: %i %i d: %1.1f dr: %1.1f val: %1.1f \n",a1,a2,d,dr,val);

  return val;
}

double ICoord::close_angle(int i)
{
  double A = 0.1;
  double val1 = 0.;
  double val2 = 0.;
  int a1 = angles[i][0];
  int a2 = angles[i][1];
  int a3 = angles[i][2];
  double d1 = distance(a1,a2);
  double d2 = distance(a2,a3);
  double dr1 = (getR(a1) + getR(a2))/2;
  double dr2 = (getR(a2) + getR(a3))/2;
  val1 = exp(-A*(d1-dr1));
  val2 = exp(-A*(d2-dr2));
  if (val1 > 1) val1 = 1;
  if (val2 > 1) val2 = 1;

  return val1*val2;
}

double ICoord::close_tor(int i)
{
  double A = 0.01;
  double val1 = 0.;
  double val2 = 0.;
  double val3 = 0.;
  int a1 = torsions[i][0];
  int a2 = torsions[i][1];
  int a3 = torsions[i][2];
  int a4 = torsions[i][3];
  double d1 = distance(a1,a2);
  double d2 = distance(a2,a3);
  double d3 = distance(a3,a4);
  double dr1 = (getR(a1) + getR(a2))/2;
  double dr2 = (getR(a2) + getR(a3))/2;
  double dr3 = (getR(a3) + getR(a4))/2;
  val1 = exp(-A*(d1-dr1));
  val2 = exp(-A*(d2-dr2));
  val3 = exp(-A*(d3-dr3));
  if (val1 > 1) val1 = 1;
  if (val2 > 1) val2 = 1;
  if (val3 > 1) val3 = 1;

  return val1*val2*val3;
}


void ICoord::structure_read(string xyzfile){ 
   
 // cout <<" Reading and initializing string coordinates" << endl;
 // cout <<"  -Opening structure file" << endl;
  
  ifstream infile;
  infile.open(xyzfile.c_str());
  if (!infile){
    cout << "!!!!Error opening xyz file!!!!" << endl;
    exit(-1);
  } 
  
 // cout <<"  -reading file..." << endl;
  
  string line;
  bool success=true;
  success=static_cast<bool>(getline(infile, line));
  if (success){
    int length=StringTools::cleanstring(line);
    natoms=atoi(line.c_str());
  }
  cout <<"  natoms: " << natoms << endl;
  
  success=static_cast<bool>(getline(infile, line));
//  if (success){  
//    comment=line;
//  }
  
  anumbers = new int[1+natoms];
  amasses = new double[1+natoms];
  anames = new string[1+natoms];
    
  //cout <<"  -Reading the atomic names...";
  for (int i=0;i<natoms;i++){
    success=static_cast<bool>(getline(infile, line));
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    anames[i]=tok_line[0];
    anumbers[i]=PTable::atom_number(anames[i]);
    amasses[i]=PTable::atom_mass(anumbers[i]);
  }
  
  infile.close();
  
//  V_profile = new double[1+nnmax];
//  S = new double[1+nnmax];
  
  coords = new double[natoms*3];
  coords0 = new double[natoms*3];
   
  //cout <<"  -Reading coordinates...";
 // cout << "Opening the xyz file" << endl;
  infile.open(xyzfile.c_str());
  fflush(stdout);
 // cout << "xyzfile opened" << endl;
  fflush(stdout);
  
  
//  for (int i=1;i<=2;i++){
    success=static_cast<bool>(getline(infile, line));
    success=static_cast<bool>(getline(infile, line));
    for (int j=0;j<natoms;j++){
      success=static_cast<bool>(getline(infile, line));
      int length=StringTools::cleanstring(line);
      vector<string> tok_line = StringTools::tokenize(line, " \t");
      coords[3*j+0]=atof(tok_line[1].c_str());
      coords[3*j+1]=atof(tok_line[2].c_str());
      coords[3*j+2]=atof(tok_line[3].c_str());
    
    }
//  }
  
  for (int i=0;i<3*natoms;i++)
     coords0[i] = coords[i];
   
 // cout << " done" << endl;
  infile.close();
  
 // cout << "Finished reading information from structure file" << endl;
}   

void ICoord::write_ic(string filename)
{
  ofstream icfile;
  icfile.open(filename.c_str());
  icfile.setf(ios::fixed);
  icfile.setf(ios::left);
  icfile << setprecision(6);

  icfile << " natoms: " << natoms << endl;

  icfile << " bonds: " << nbonds << endl;
  for (int i=0;i<nbonds;i++)
    icfile << " " << bonds[i][0] << " " << bonds[i][1] << endl;

  icfile << " angles: " << nangles << endl;
  for (int i=0;i<nangles;i++)
    icfile << " " << angles[i][0] << " " << angles[i][1] << " " << angles[i][2] << endl;

  icfile << " torsions: " << ntor << endl;
  for (int i=0;i<ntor;i++)
    icfile << " " << torsions[i][0] << " " << torsions[i][1] << " " << torsions[i][2] << " " << torsions[i][3] << endl;

  icfile.close();

  return;
}


int ICoord::read_ics(string filename)
{
  ifstream infile;
  infile.open(filename.c_str());
  if (!infile){
    printf(" Error: couldn't open icfile: %s \n",filename.c_str());
    exit(-1);
  } 

  string line;
  bool success=true;
  int type = 1;
  success=static_cast<bool>(getline(infile, line));
  while (!infile.eof())
  {
    success=static_cast<bool>(getline(infile, line));
    //cout << "RR0: " << line << endl;
    int length=StringTools::cleanstring(line);
    vector<string> tok_line = StringTools::tokenize(line, " \t");
    if (tok_line.size()>0)
    {
      if (type==1)
      {
        nbonds = atoi(tok_line[1].c_str());
#if 0
        bonds = new int*[nbonds];
        for (int j=0;j<nbonds;j++)
          bonds[j] = new int[2];
#endif

        for (int i=0;i<nbonds;i++)
        {
          success=static_cast<bool>(getline(infile, line));
          //cout << "RR: " << line << endl;
          length=StringTools::cleanstring(line);
          tok_line = StringTools::tokenize(line, " \t");
          bonds[i][0] = atoi(tok_line[0].c_str());
          bonds[i][1] = atoi(tok_line[1].c_str());
        }
        //printf(" found %i bonds \n",nbonds);
        type++;
      } //read bonds
      else if (type==2)
      {
        //printf(" reading angles \n"); fflush(stdout);
        nangles = atoi(tok_line[1].c_str());
#if 0
        angles = new int*[nangles];
        for (int j=0;j<nangles;j++)
          angles[j] = new int[3];
#endif

        for (int i=0;i<nangles;i++)
        {
          success=static_cast<bool>(getline(infile, line));
          //cout << "RR: " << line << endl;
          length=StringTools::cleanstring(line);
          tok_line = StringTools::tokenize(line, " \t");
          angles[i][0] = atoi(tok_line[0].c_str());
          angles[i][1] = atoi(tok_line[1].c_str());
          angles[i][2] = atoi(tok_line[2].c_str());
        }
        //printf(" found %i nangles \n",nangles);
        type++;
      } //read angles
      else if (type==3)
      {
        //printf(" reading torsions \n"); fflush(stdout);
        ntor = atoi(tok_line[1].c_str());
#if 0
        torsions = new int*[ntor];
        for (int j=0;j<ntor;j++)
          torsions[j] = new int[4];
#endif

        for (int i=0;i<ntor;i++)
        {
          success=static_cast<bool>(getline(infile, line));
          //cout << "RR: " << line << endl;
          length=StringTools::cleanstring(line);
          tok_line = StringTools::tokenize(line, " \t");
          torsions[i][0] = atoi(tok_line[0].c_str());
          torsions[i][1] = atoi(tok_line[1].c_str());
          torsions[i][2] = atoi(tok_line[2].c_str());
          torsions[i][3] = atoi(tok_line[3].c_str());
        }
        //printf(" found %i torsions \n",ntor);
        type++;
      } //read torsion
      if (type>3) break;
    }
  } //while !eof

#if 0
  printf(" bonds: \n");
  for (int i=0;i<nbonds;i++)
    printf(" %i %i \n",bonds[i][0],bonds[i][1]);
  printf(" angles: \n");
  for (int i=0;i<nangles;i++)
    printf(" %i %i %i \n",angles[i][0],angles[i][1],angles[i][2]);
  printf(" torsions: \n");
  for (int i=0;i<ntor;i++)
    printf(" %i %i %i %i \n",torsions[i][0],torsions[i][1],torsions[i][2],torsions[i][3]);
#endif

  
  infile.close();

  return nbonds + nangles + ntor;
}

//double ICoord::mm_energy()
double ICoord::mm_energy(int atom1, int atom2)
{
  double Energy = 0.;
  //Energy += vdw_energy_all();
  Energy += vdw_energy_all(atom1, atom2);
  //others not yet implemented

  return Energy;
}

//double ICoord::vdw_energy_all()
double ICoord::vdw_energy_all(int atom1, int atom2)
{
  double E = 0;
  //std::cout << "vdw_energy_all\n" << n_nonbond << "     " << std::endl;
  //for (int i=0;i<n_nonbond;i++)
  //  E += vdw_energy_1(nonbond[i][0],nonbond[i][1]);
  E = vdw_energy_1(atom1, atom2);

  return E;
}

double ICoord::vdw_energy_1(int i, int j)
{
  double R = ffR[i] + ffR[j];
  std::cout << "test R:  " << ffR[i] << "   " << ffR[j] << std::endl;
  double eps = sqrt( ffeps[i] * ffeps[j] );

  double r = distance(i,j);

//  printf(" R: %1.4f r: %1.4f \n",R,r);
  double Rr = R / r;
  double Rr6 = Rr*Rr*Rr;
  Rr6 = Rr6*Rr6;

  double E = eps * ( Rr6*Rr6 - 2*Rr6 );

  return E;
}

void ICoord::freemem(){

 delete [] grad;
    
 for (int i=0;i<max_bonds;i++)
   delete [] bonds[i];
 delete [] bonds;
 delete [] bondd;

 int max_bonds_tm = 25;
 /*
 for (int i=0;i<max_bonds_tm;i++)
   delete [] bondstm[i];
 delete [] bondstm;
 for (int i=0;i<max_bonds_tm;i++)
   delete [] bondstm2[i];
 delete [] bondstm2;
 Mina */

 for (int i=0;i<max_angles;i++)
   delete [] angles[i];
 delete [] angles; 
 delete [] anglev;

 for (int i=0;i<max_torsions;i++)
   delete [] torsions[i];
 delete [] torsions;
 delete [] torv;

 for (int i=0;i<max_imptor;i++)
   delete [] imptor[i];
 delete [] imptor;
 delete [] imptorv;

 for (int i=0;i<max_nonbond;i++)
   delete [] nonbond[i];
 delete [] nonbond;
 delete [] nonbondd;

 delete [] coordn;

 delete [] ffR;
 delete [] ffeps;

  delete [] anumbers;
  delete [] amasses;
  delete [] anames;
  delete [] coords0;
  delete [] coordsts;
  delete [] coords;


  return;

}

void ICoord::alloc_mem(){
 farBond = 1.0;
 use_xyz = 0;

 frozen = NULL;
 xyzic = new int[natoms];

// pgrad = new double[3*natoms];
 grad = new double[3*natoms];
// for (int i=0;i<3*natoms;i++) pgrad[i]=0;
 for (int i=0;i<3*natoms;i++) grad[i]=0;

 nbonds = 0;
 max_bonds=natoms*natoms+8;
 bonds = new int*[max_bonds];
 for (int i=0;i<max_bonds;i++)
   bonds[i]=new int[2];
 bondd = new double[max_bonds];

 nangles = 0;
 max_angles=natoms*24;
 angles = new int*[max_angles];
 for (int i=0;i<max_angles;i++)
   angles[i]=new int[3];
 anglev = new double[max_angles];

 ntor = 0;
 max_torsions=natoms*60;
 torsions = new int*[max_torsions];
 for (int i=0;i<max_torsions;i++)
   torsions[i] = new int[4];
 torv = new double[max_torsions];

 nimptor = 0;
 max_imptor=natoms*10;
 imptor = new int*[max_imptor];
 for (int i=0;i<max_imptor;i++)
   imptor[i] = new int[4];
 imptorv = new double[max_imptor];


 n_nonbond = 0;
 max_nonbond = natoms*natoms;
 nonbond = new int*[max_nonbond];
 for (int i=0;i<max_nonbond;i++)
   nonbond[i]=new int[2];
 nonbondd = new double[max_nonbond];

 coordn = new int[natoms];
 for (int i=0;i<natoms;i++)
    coordn[i]=0;

 ffR = new double[natoms];
 ffeps = new double[natoms];

 return;
}
