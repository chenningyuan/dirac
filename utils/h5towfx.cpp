/*******************************************
--------------------------------------------
--------------------------------------------
---      DIRAC H5 to WFX INTERFACE       ---
--------------------------------------------
--------------------------------------------
-- Developed by: Dr. M. Rodriguez-Mayorga --
--             30.10.2022                 --
--     email: marm3.14@gmail.com          --
--------------------------------------------
--------------------------------------------
*******************************************/
#include<iostream>
#include<algorithm>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include<vector>
#include<complex>
#include"hdf5.h"

#define zero 0.0e0
#define one  1.0e0
#define two  2.0e0
#define ten  1.0e1
#define tol6 1.0e-6
#define tol8 1.0e-8
#define ang2au 1.8897259886

using namespace std;

// Functions
void read_int(hid_t file_id,string path,string var, vector<int> *to_store);
void read_dou(hid_t file_id,string path,string var, vector<double> *to_store);
void asign(vector<int> *s2atom_map,vector<double> *Shell_coord_L,vector<double> *Shell_coord_S);
bool coor_atoms_shells(vector<double> *Shell_coord_L,vector<double> *Atom_pos);
void h5_read();
void sort_basis();
void clean_shell2aos();
void print_wfx();

// Global variables
complex<double>czero=(zero,zero);
int ith,Ncenters,Nprimitives,Nbasis,Nbasis_L,Nbasis_S,Nshell,Nshell_L,Nshell_S,NMOs,NMOs_LS,NMOs_occ;
struct Shell2AOs
{
 int styp,atom,nprim,naos,paired2;
 int **index_min2max; 
 double *Coef,*Expon,Coord[3];
};
struct NewCenter
{
 int index;
 double Coord[3];
};
Shell2AOs *shell2aos;
double Quaternion_coef[4];
double *OCCs,**Prim2AO_Coef;
complex<double> **AO2MO_Coef,**Prim2MO_Coef;
string h5_name;
vector<int>Prim2Center_map;
vector<int>shell_types;
vector<double>Center_x;
vector<double>Center_y;
vector<double>Center_z;
vector<double>prim_exponents;
vector<double>MOsLS_occ;
vector<double>Nu_charge;
vector<double>One_Prim2MO_Coef_RE;
vector<double>One_Prim2MO_Coef_IM;


// Main  
int main(int argc, char *argv[])
{
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"---      DIRAC H5 to WFX INTERFACE       ---"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"-- Developed by: Dr. M. Rodriguez-Mayorga --"<<endl;
 cout<<"--             25.10.2022                 --"<<endl;
 cout<<"--     email: marm3.14@gmail.com          --"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<endl;
 if(argc!=2)
 {
  cout<<endl;
  cout<<"Please, include the h5 file as an argument"<<endl;
  cout<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<"--          Normal termination            --"<<endl;
  cout<<"--------------------------------------------"<<endl;
  return -1;
 } 
 bool repeated_prims;
 int ishell,ishell1,iprim,iprim1,iaos,iaos1,imos,imos1; // imos for Scalar MOs
 int naos;
 string aux(argv[1]);
 h5_name=aux;
 h5_read();
 // Read Dirac output
 // Find repeated shells ("pairing")
 for(ishell=0;ishell<Nshell-1;ishell++)
 {
  for(ishell1=ishell+1;ishell1<Nshell;ishell1++)
  {
   repeated_prims=false;
   if(shell2aos[ishell1].styp==shell2aos[ishell].styp && shell2aos[ishell1].nprim==shell2aos[ishell].nprim)
   {
    if(shell2aos[ishell1].Coord[0]==shell2aos[ishell].Coord[0] && 
       shell2aos[ishell1].Coord[1]==shell2aos[ishell].Coord[1] && 
       shell2aos[ishell1].Coord[2]==shell2aos[ishell].Coord[2])
    {
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      if(shell2aos[ishell1].Expon[iprim]==shell2aos[ishell].Expon[iprim])
      {
       repeated_prims=true;
      }
      else
      {
       repeated_prims=false;
       iprim=shell2aos[ishell].nprim;
      }
     }
    }
   }
   // Save pairing info
   if(repeated_prims && shell2aos[ishell1].paired2==-1)
   {
    shell2aos[ishell1].paired2=ishell;
   }
  }
 }
 // Print .basis file with unique primitives information
 sort_basis();
 cout<<"Number of Primitives : "<<setw(12)<<Nprimitives<<endl;
 // Find the total number of AOs and allocate indexing arrays
 naos=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  naos=naos+shell2aos[ishell].naos;
  shell2aos[ishell].index_min2max=new int*[shell2aos[ishell].naos];
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   shell2aos[ishell].index_min2max[iaos]=new int[shell2aos[ishell].nprim];
  }
 }
 cout<<"Number of AOs (calc) : "<<setw(12)<<naos<<endl;
 // Allocate the primitive coefs. matrix (used to build AOs from primitives)
 // NOTE: IT IS BUILT TRANSPOSED TO ACCELERATE THE MATRIX MULTIPLICATION
 Prim2AO_Coef=new double*[Nprimitives];
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  Prim2AO_Coef[iprim]=new double[Nprimitives];
  for(iaos=0;iaos<naos;iaos++)
  {
   Prim2AO_Coef[iprim][iaos]=zero;
  }
 }
 // Initialize the indexing for the "unpaired" shells
 iprim1=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2==-1)
  {
   for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
   {
    for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
    {
     shell2aos[ishell].index_min2max[iaos][iprim]=iprim1;
     iprim1++;
    }
   }
  }
 }
 // Copy indexing from "unpaired" shells to the "paired" ones
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2!=-1)
  {
   for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
   {
    for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
    {
     shell2aos[ishell].index_min2max[iaos][iprim]=shell2aos[shell2aos[ishell].paired2].index_min2max[iaos][iprim];
    }
   }
  } 
 }
 // Organize Prim 2 AO Coefs matrix using index_min2max
 cout<<endl;
 cout<<"AO to Primitives map (per shell):"<<endl;
 iaos1=1;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   cout<<"AO"<<setw(4)<<iaos1<<":";
   for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
   {
    cout<<"\t"<<shell2aos[ishell].index_min2max[iaos][iprim]+1;
   }
   cout<<endl;
   iaos1++;
  }
 }
 // Fill the Prim2AO_Coef matrix
 iaos1=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  for(iaos=0;iaos<shell2aos[ishell].naos;iaos++)
  {
   for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
   {
    Prim2AO_Coef[shell2aos[ishell].index_min2max[iaos][iprim]][iaos1]=shell2aos[ishell].Coef[iprim];
   }
   iaos1++;
  }
 }
 // Deallocate arrays that are not used anymore
 clean_shell2aos();
 // Compute Primitives to MO coefficients
 cout<<endl;
 cout<<"Building the matrix Prim2MO_Coefs[MO][Primitives] ((2x4MO) x 2(L+S)) = AO2MO_Coef[MO][AO]*Prim2AO_Coef[AO][Primitives]"<<endl;
 Prim2MO_Coef=new complex<double>*[NMOs_LS]; 
 for(imos=0;imos<NMOs_LS;imos++)
 {
  Prim2MO_Coef[imos]=new complex<double>[Nprimitives];
  for(iprim=0;iprim<Nprimitives;iprim++)
  {
   Prim2MO_Coef[imos][iprim]=czero;
   for(iaos=0;iaos<Nbasis;iaos++)
   {
    // Recall that Prim2AO_Coef is built transposed
    Prim2MO_Coef[imos][iprim]=Prim2MO_Coef[imos][iprim]+AO2MO_Coef[imos][iaos]*Prim2AO_Coef[iprim][iaos];
   }
  }
 }
 cout<<"Matrix Prim2MO_Coefs[MO][Primitives] for bar and unbar orbs. built."<<endl;
 cout<<endl;
 // Deallocate arrays related to AOs
 for(imos=0;imos<NMOs_LS;imos++)
 {
  delete [] AO2MO_Coef[imos];AO2MO_Coef[imos]=NULL;
 }
 delete[] AO2MO_Coef;AO2MO_Coef=NULL;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  delete[] Prim2AO_Coef[iprim];Prim2AO_Coef[iprim]=NULL;
 }
 delete[] Prim2AO_Coef;Prim2AO_Coef=NULL;
 // Print the Prim2MO_Coef matrix (coefficients are rows). 
 // WARNING! Below we may overwrite the positronic states (initial ones) with the occ. electronic states (i.e. put them on top of the list)!
 imos1=0;
 for(imos=0;imos<NMOs_LS;imos++)
 {
  if(OCCs[imos]>tol8) // Overwrite positronic MO coefs with electronic ones and print the electronic ones
  {
   for(iprim=0;iprim<Nprimitives;iprim++)
   {
    Prim2MO_Coef[imos1][iprim]=Prim2MO_Coef[imos][iprim];
   }
   MOsLS_occ.push_back(OCCs[imos]);
   imos1++;
  }
 }
 NMOs_occ=imos1;
 cout<<"Num. of occ MO (Scalar): "<<setw(12)<<NMOs_occ<<endl;
 // Print a WFX file.
 print_wfx();
 delete[] OCCs;OCCs=NULL;
 // Deallocate array Primitive to MO coefs (rows).
 for(imos=0;imos<NMOs_LS;imos++)
 {
  delete [] Prim2MO_Coef[imos];Prim2MO_Coef[imos]=NULL;
 }
 delete[] Prim2MO_Coef;Prim2MO_Coef=NULL;
 cout<<"                                            "<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"--          Normal termination            --"<<endl;
 cout<<"--------------------------------------------"<<endl;
 return 0;
}

void h5_read()
{
 bool atom_shell_coord,electronic;
 int i,j,k,nbasis_tot,coefs_quat;
 int imos,imos1,imos2,imos3,imos4,imos5,imos6,imos7,imos8;
 double **MO2AO_coefs_quat,occ;
 hid_t file_id, dset_id, dspace_id; /* identifiers */
 herr_t status;
 vector<int> Nbasis_L,Shell_type_L,Nprim_shell_L,Nshells_L;
 vector<double> Shell_coord_L,Prim_exponents_L,Prim_coefs_L;
 vector<int> Nbasis_S,Shell_type_S,Nprim_shell_S,Nshells_S;
 vector<double> Shell_coord_S,Prim_exponents_S,Prim_coefs_S;
 vector<int> Nbasis_MO,s2atom_map;
 vector<double> MO_occs,MO_coefs,Nuc_charge,Atom_pos;
 string path;

 file_id = H5Fopen(h5_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
 // Reading charges
  path="/input/molecule/";
  // Retrieve nuclear charges
  read_dou(file_id,path,"nuc_charge",&Nuc_charge);
  for(i=0;i<Nuc_charge.size();i++)
  {
   Nu_charge.push_back(Nuc_charge[i]);
  }
  // Retrieve atomic positions
  read_dou(file_id,path,"geometry",&Atom_pos);
 // Reading Large component basis
  path="/input/aobasis/1/";
  // Retrieve Num. AO Large
  read_int(file_id,path,"n_ao",&Nbasis_L);
  // Retrieve Num. large shells
  read_int(file_id,path,"n_shells",&Nshells_L);
  Nshell_L=Nshells_L[0];
  // Retrieve Shell_type_L
  read_int(file_id,path,"orbmom",&Shell_type_L);
  // Retrieve Nprim per shell L
  read_int(file_id,path,"n_prim",&Nprim_shell_L);
  // double 1D vector
  // Retrieve Shell_coord_L
  read_dou(file_id,path,"center",&Shell_coord_L);
  // Retrieve Prim_exponents_L
  read_dou(file_id,path,"exponents",&Prim_exponents_L);
  // Retrieve Prim_coefs_L
  read_dou(file_id,path,"contractions",&Prim_coefs_L);
 // Reading Small component basis
  path="/input/aobasis/2/";
  // Retrieve Num. AO Small
  read_int(file_id,path,"n_ao",&Nbasis_S);
  Nbasis=Nbasis_L[0]+Nbasis_S[0];
  // Retrieve Num. small shells
  read_int(file_id,path,"n_shells",&Nshells_S);
  Nshell_S=Nshells_S[0];
  Nshell=Nshells_L[0]+Nshells_S[0];
  shell2aos=new Shell2AOs[Nshell];
  for(i=0;i<Nshell_L;i++)
  {
   shell2aos[i].paired2=-1;
  }
  // Retrieve Shell_type_S
  read_int(file_id,path,"orbmom",&Shell_type_S);
  for(i=0;i<Nshell_L;i++)
  {
   shell2aos[i].styp=Shell_type_L[i]-1;
   if(shell2aos[i].styp==0)
   {
    shell2aos[i].naos=1;
   }
   else if(shell2aos[i].styp==1)
   {
    shell2aos[i].naos=3;
   }
   else if(shell2aos[i].styp==2)
   {
    shell2aos[i].naos=6;
   }
   else if(shell2aos[i].styp==3)
   {
    shell2aos[i].naos=10;
   }
   else if(shell2aos[i].styp==4)
   {
    shell2aos[i].naos=15;
   }
   else
   {
    cout<<"Warning! Shell type not supported "<<setw(5)<<shell2aos[i].styp<<endl;
   }
  }
  for(i=Nshell_L;i<Nshell;i++)
  {
   shell2aos[i].styp=Shell_type_S[i-Nshell_L]-1;
   if(shell2aos[i].styp==0)
   {
    shell2aos[i].naos=1;
   }
   else if(shell2aos[i].styp==1)
   {
    shell2aos[i].naos=3;
   }
   else if(shell2aos[i].styp==2)
   {
    shell2aos[i].naos=6;
   }
   else if(shell2aos[i].styp==3)
   {
    shell2aos[i].naos=10;
   }
   else if(shell2aos[i].styp==4)
   {
    shell2aos[i].naos=15;
   }
   else
   {
    cout<<"Warning! Shell type not supported "<<setw(5)<<shell2aos[i].styp<<endl;
   }
  }
  // Retrieve Nprim per shell S
  read_dou(file_id,path,"center",&Shell_coord_S);
  asign(&s2atom_map,&Shell_coord_L,&Shell_coord_S);
  for(i=0;i<s2atom_map.size();i++)
  {
   shell2aos[i].atom=s2atom_map[i];
  }
  // Retrieve Shell_coord_S
  read_int(file_id,path,"n_prim",&Nprim_shell_S);
  for(i=0;i<Nprim_shell_L.size();i++)
  {
   shell2aos[i].nprim=Nprim_shell_L[i];
   shell2aos[i].Coef=new double[shell2aos[i].nprim];
   shell2aos[i].Expon=new double[shell2aos[i].nprim];
  }
  for(i=Nprim_shell_L.size();i<Nprim_shell_S.size()+Nprim_shell_L.size();i++)
  {
   shell2aos[i].nprim=Nprim_shell_S[i-Nprim_shell_L.size()];
   shell2aos[i].Coef=new double[shell2aos[i].nprim];
   shell2aos[i].Expon=new double[shell2aos[i].nprim];
  }
  // Retrieve Prim_coefs_S
  read_dou(file_id,path,"contractions",&Prim_coefs_S);
  k=0;
  for(i=0;i<Nshell;i++)
  {
   if(i==Nshell_L){k=0;}
   for(j=0;j<shell2aos[i].nprim;j++)
   {
    if(i<Nshell_L)
    {
     shell2aos[i].Coef[j]=Prim_coefs_L[k];
    }
    else
    {
     shell2aos[i].Coef[j]=Prim_coefs_S[k];
    }
    k++;
   }
  }
  // Retrieve Prim_exponents_S
  read_dou(file_id,path,"exponents",&Prim_exponents_S);
  k=0;
  for(i=0;i<Nshell;i++)
  {
   if(i==Nshell_L){k=0;}
   for(j=0;j<shell2aos[i].nprim;j++)
   {
    if(i<Nshell_L)
    {
     shell2aos[i].Expon[j]=Prim_exponents_L[k];
    }
    else
    {
     shell2aos[i].Expon[j]=Prim_exponents_S[k];
    }
    k++;
   }
  }
  // Print shell cartes coordinates
  k=0;
  for(i=0;i<Nshell;i++)
  {
   if(i==Nshell_L){k=0;}
   for(j=0;j<3;j++)
   {
    if(i<Nshell_L)
    {
     shell2aos[i].Coord[j]=Shell_coord_L[k];
    }
    else
    {
     shell2aos[i].Coord[j]=Shell_coord_S[k];
    }
    k++;
   }
  }
 // Read results information
  path="/result/wavefunctions/scf/mobasis/";
  // Retrieve MO_occs
  read_dou(file_id,path,"occupations",&MO_occs);
  OCCs=new double[8*MO_occs.size()];
  for(imos=0;imos<8*MO_occs.size();imos++)
  {
   OCCs[imos]=-ten;
  }
  electronic=false;
  for(imos=0;imos<MO_occs.size();imos++)
  {
   occ=MO_occs[imos];
   if(abs(occ)>tol8)
   {
    electronic=true;
    imos1=8*imos;
    imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
    OCCs[imos1]=occ;OCCs[imos2]=occ;OCCs[imos3]=occ;OCCs[imos4]=occ;OCCs[imos5]=occ;OCCs[imos6]=occ;OCCs[imos7]=occ;OCCs[imos8]=occ;
   }
   if(electronic && abs(occ)<=tol8)
   {
    occ=zero;
    imos1=8*imos;
    imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
    OCCs[imos1]=occ;OCCs[imos2]=occ;OCCs[imos3]=occ;OCCs[imos4]=occ;OCCs[imos5]=occ;OCCs[imos6]=occ;OCCs[imos7]=occ;OCCs[imos8]=occ;
   }
  }
  // Retrieve Num. MO
  read_int(file_id,path,"n_mo",&Nbasis_MO);
  NMOs=Nbasis_MO[0];
  NMOs_LS=NMOs*8;
  // Retrieve MO_coefs
  read_dou(file_id,path,"orbitals",&MO_coefs);
  coefs_quat=MO_coefs.size()/4;
  nbasis_tot=coefs_quat/Nbasis_MO[0];
  if(nbasis_tot!=Nbasis)
  {
   cout<<endl;
   cout<<"Warning! The calculated size of basis Nbasis_L+Nbasis_S does not much the stored in the file."<<endl;
   cout<<endl;
  }
  MO2AO_coefs_quat=new double*[coefs_quat];
  for(i=0;i<coefs_quat;i++)
  {
   MO2AO_coefs_quat[i]=new double[4];
  } 
  k=0;
  for(i=0;i<4;i++)
  {
   for(j=0;j<coefs_quat;j++)
   {
    MO2AO_coefs_quat[j][i]=MO_coefs[k];
    k++;
   }
  }
  AO2MO_Coef=new complex<double>*[NMOs_LS];
  for(imos=0;imos<NMOs_LS;imos++)
  {
   AO2MO_Coef[imos]=new complex<double>[nbasis_tot];
   for(j=0;j<Nbasis;j++)
   {
    AO2MO_Coef[imos][j]=czero;
   }
  }
  for(imos=0;imos<Nbasis_MO[0];imos++)
  {
   imos1=imos*8;imos2=imos1+1;imos3=imos2+1;imos4=imos3+1;imos5=imos4+1;imos6=imos5+1;imos7=imos6+1;imos8=imos7+1;
   k=0;
   for(j=0;j<nbasis_tot;j++)
   {
    Quaternion_coef[0]=MO2AO_coefs_quat[j+imos*nbasis_tot][0];
    Quaternion_coef[1]=MO2AO_coefs_quat[j+imos*nbasis_tot][1];
    Quaternion_coef[2]=MO2AO_coefs_quat[j+imos*nbasis_tot][2];
    Quaternion_coef[3]=MO2AO_coefs_quat[j+imos*nbasis_tot][3];
    complex<double>ztmp4c1( Quaternion_coef[0], Quaternion_coef[1]); 
    complex<double>ztmp4c2(-Quaternion_coef[2], Quaternion_coef[3]); 
    complex<double>ztmp4c3( Quaternion_coef[2], Quaternion_coef[3]); 
    complex<double>ztmp4c4( Quaternion_coef[0],-Quaternion_coef[1]); 
    if(k<Nbasis_L[0])
    {
     AO2MO_Coef[imos1][j]=ztmp4c1;  // unbar alpha L
     AO2MO_Coef[imos2][j]=ztmp4c2;  // unbar beta L
     AO2MO_Coef[imos5][j]=ztmp4c3;  // bar alpha L
     AO2MO_Coef[imos6][j]=ztmp4c4;  // bar beta L
    }
    else
    {
     AO2MO_Coef[imos3][j]=ztmp4c1;  // unbar alpha S
     AO2MO_Coef[imos4][j]=ztmp4c2;  // unbar beta S
     AO2MO_Coef[imos7][j]=ztmp4c3;  // bar alpha S
     AO2MO_Coef[imos8][j]=ztmp4c4;  // bar beta S
    }
    k++;
   }
  }
  // Complex conjugate the coefficients because now they are stored as rows
  for(imos=0;imos<NMOs_LS;imos++)
  {
   for(j=0;j<Nbasis;j++)
   { 
    AO2MO_Coef[imos][j]=conj(AO2MO_Coef[imos][j]);
   }
  }
  for(i=0;i<coefs_quat;i++)
  {
   delete[] MO2AO_coefs_quat[i];MO2AO_coefs_quat[i]=NULL;
  } 
  delete[] MO2AO_coefs_quat;MO2AO_coefs_quat=NULL;
 // Close the file
 status = H5Fclose(file_id);
 atom_shell_coord=coor_atoms_shells(&Shell_coord_L,&Atom_pos);
 if(!atom_shell_coord)
 {
  cout<<"Warning! Atomic coordinates and shell coordinates do not coincide"<<endl;
 }
}

void read_int(hid_t file_id, string path, string var, vector<int> *to_store)
{
 int idims,ndims,res_sz;
 hid_t dset_id, dspace_id;
 herr_t status;
 string string_dset=path+var; 
 dset_id = H5Dopen(file_id,string_dset.c_str(),H5P_DEFAULT);
 dspace_id = H5Dget_space(dset_id);
 ndims = H5Sget_simple_extent_ndims(dspace_id);
 hsize_t res_dims[ndims];
 status = H5Sget_simple_extent_dims(dspace_id, res_dims, NULL);
 res_sz = 1;
 for(idims=0;idims<ndims;idims++)
 {
  res_sz *= res_dims[idims];
 }
 to_store[0].resize(res_sz);
 status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, to_store[0].data());
 status = H5Dclose(dset_id);
}

void read_dou(hid_t file_id, string path, string var, vector<double> *to_store)
{
 int idims,ndims,res_sz;
 hid_t dset_id, dspace_id;
 herr_t status;
 string string_dset=path+var; 
 dset_id = H5Dopen(file_id,string_dset.c_str(),H5P_DEFAULT);
 dspace_id = H5Dget_space(dset_id);
 ndims = H5Sget_simple_extent_ndims(dspace_id);
 hsize_t res_dims[ndims];
 status = H5Sget_simple_extent_dims(dspace_id, res_dims, NULL);
 res_sz = 1;
 for(idims=0;idims<ndims;idims++)
 {
  res_sz *= res_dims[idims];
 }
 to_store[0].resize(res_sz);
 status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, to_store[0].data());
 status = H5Dclose(dset_id);
}

void asign(vector<int> *s2atom_map,vector<double> *Shell_coord_L,vector<double> *Shell_coord_S)
{
 int i,j,k,ind,nshells=(Shell_coord_L[0].size()+Shell_coord_S[0].size())/3;
 bool newcent;
 double diff;
 vector<NewCenter>newcenter;
 NewCenter centertemp;
 ind=1;
 centertemp.index=ind;
 centertemp.Coord[0]=Shell_coord_L[0].at(0);
 centertemp.Coord[1]=Shell_coord_L[0].at(1);
 centertemp.Coord[2]=Shell_coord_L[0].at(2);
 newcenter.push_back(centertemp);
 ind++;
 for(i=1;i<Shell_coord_L[0].size()/3;i++)
 {
  newcent=true;
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_L[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    newcent=false;j=newcenter.size();
   }
  }
  if(newcent)
  {
   centertemp.index=ind;
   for(k=0;k<3;k++)
   {
    centertemp.Coord[k]=Shell_coord_L[0].at(3*i+k);
   }
   newcenter.push_back(centertemp);
   ind++;
  } 
 }
 for(i=0;i<Shell_coord_S[0].size()/3;i++)
 {
  newcent=true;
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_S[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    newcent=false;j=newcenter.size();
   }
  }
  if(newcent)
  {
   cout<<"Warning! New center found for Small component w.r.t. Large component"<<endl;
   centertemp.index=ind;
   for(k=0;k<3;k++)
   {
    centertemp.Coord[k]=Shell_coord_S[0].at(3*i+k);
   }
   newcenter.push_back(centertemp);
   ind++;
  } 
 }
 for(i=0;i<Shell_coord_L[0].size()/3;i++)
 {
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_L[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    s2atom_map[0].push_back(newcenter[j].index);
    j=newcenter.size();
   }
  }
 }
 for(i=0;i<Shell_coord_S[0].size()/3;i++)
 {
  for(j=0;j<newcenter.size();j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(newcenter[j].Coord[k]-Shell_coord_S[0].at(3*i+k),two);
   }
   if(sqrt(diff)<tol8)
   {
    s2atom_map[0].push_back(newcenter[j].index);
    j=newcenter.size();
   }
  }
 }
}

bool coor_atoms_shells(vector<double> *Shell_coord_L,vector<double> *Atom_pos)
{
 bool atshell_coord=true,found=true;
 int i,j,k;
 double diff;
 for(i=0;i<Atom_pos[0].size()/3;i++)
 {
  found=false;
  for(j=0;j<Shell_coord_L[0].size()/3;j++)
  {
   diff=zero;
   for(k=0;k<3;k++)
   {
    diff=diff+pow(ang2au*Atom_pos[0].at(3*i+k)-Shell_coord_L[0].at(3*j+k),two);
   } 
   if(sqrt(diff)<tol6)
   {
    found=true;
    j=Shell_coord_L[0].size();
   }
  }
  if(!found) 
  {
   atshell_coord=false;
   i=Atom_pos[0].size();
  }
 }  
 return atshell_coord;
}

// Function used to print the basis file and check how many unique primitives we used in the calc.
void sort_basis()
{
 bool newcenter;
 int ishell,iprim,icenter;
 vector<double>Coord_x;
 vector<double>Coord_y;
 vector<double>Coord_z;
 Nprimitives=0;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  if(shell2aos[ishell].paired2==-1)
  {
   switch(shell2aos[ishell].styp){
    case 0:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim;
     // 0 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(1);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case 1:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*3;
     // 1 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(2);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(3);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(4);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case 2:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*6;
     // 2 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(5);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(8);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(9);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(6);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(10);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(7);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case 3:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*10;
     // 3 0 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(11);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 1 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(14);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 0 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(15);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 2 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(17);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(20);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(18);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 3 0 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(12);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 1 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(16);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 2 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(19);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 3 
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(13);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    case 4:
     Nprimitives=Nprimitives+shell2aos[ishell].nprim*15;
     // 4 0 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(21);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 3 1 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(24);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 3 0 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(25);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 2 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(30);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 1 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(33);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 2 0 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(31);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 3 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(26);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 2 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(34);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 1 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(35);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 1 0 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(28);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 4 0
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(22);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 3 1
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(27);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 2 2
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(32);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 1 3
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(29);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     // 0 0 4
     for(iprim=0;iprim<shell2aos[ishell].nprim;iprim++)
     {
      shell_types.push_back(23);
      prim_exponents.push_back(shell2aos[ishell].Expon[iprim]);
      Coord_x.push_back(shell2aos[ishell].Coord[0]);Coord_y.push_back(shell2aos[ishell].Coord[1]);Coord_z.push_back(shell2aos[ishell].Coord[2]);
     }
     break;
    default:
     cout<<"Warning! Shell-type not supported. "<<setw(5)<<shell2aos[ishell].styp<<endl;
     break;
   }
  }
 }
 // Save coords. of the unique centers (Atoms) 
 Ncenters=1;
 Center_x.push_back(Coord_x[0]);Center_y.push_back(Coord_y[0]);Center_z.push_back(Coord_z[0]);
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  newcenter=true;
  for(icenter=0;icenter<Ncenters;icenter++)
  {
   if(Coord_x[iprim]==Center_x[icenter] && Coord_y[iprim]==Center_y[icenter] && Coord_z[iprim]==Center_z[icenter])
   {
    newcenter=false;
   }
  }
  if(newcenter)
  {
   Center_x.push_back(Coord_x[iprim]);Center_y.push_back(Coord_y[iprim]);Center_z.push_back(Coord_z[iprim]);
   Ncenters++;
  }
 }
 // Save Primitive to Centers map 
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  for(icenter=0;icenter<Ncenters;icenter++)
  {
   if(Coord_x[iprim]==Center_x[icenter] && Coord_y[iprim]==Center_y[icenter] && Coord_z[iprim]==Center_z[icenter])
   {
    Prim2Center_map.push_back(icenter+1);
   }
  }
 }
}

void clean_shell2aos()
{
 int ishell,iaos;
 for(ishell=0;ishell<Nshell;ishell++)
 {
  delete[] shell2aos[ishell].Coef; shell2aos[ishell].Coef=NULL;
  delete[] shell2aos[ishell].Expon; shell2aos[ishell].Expon=NULL;
  for(iaos=0;shell2aos[ishell].naos<1;iaos++)
  {
   delete[] shell2aos[ishell].index_min2max[iaos]; shell2aos[ishell].index_min2max[iaos]=NULL;
  }
  delete[] shell2aos[ishell].index_min2max;shell2aos[ishell].index_min2max=NULL;
 }
 delete[] shell2aos;shell2aos=NULL;
}

// Currently, only available for ATOMS
void print_wfx()
{
 int iprim,imos,icenter;
 string line;
 ofstream real_wfx((h5_name.substr(0,h5_name.length()-3)+"_RE.wfx").c_str());
 ofstream imag_wfx((h5_name.substr(0,h5_name.length()-3)+"_IM.wfx").c_str());
 real_wfx<<setprecision(15)<<fixed<<scientific;
 imag_wfx<<setprecision(15)<<fixed<<scientific;
 line="<Title>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="DIRAC H5 to WFX file";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="</Title>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Keywords>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="GTO";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="</Keywords>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 real_wfx<<Ncenters<<endl; 
 imag_wfx<<Ncenters<<endl; 
 line="</Number of Nuclei>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Number of Occupied Molecular Orbitals>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 real_wfx<<NMOs_occ<<endl; 
 imag_wfx<<NMOs_occ<<endl;
 line="</Number of Occupied Molecular Orbitals>";
 real_wfx<<line<<endl; 
 imag_wfx<<line<<endl; 
 line="<Electronic Spin Multiplicity>"; // Fixed for 'faked' close shell
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<1<<endl; 
 imag_wfx<<1<<endl; 
 line="</Electronic Spin Multiplicity>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Atomic Numbers>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(icenter=0;icenter<Nu_charge.size();icenter++)
 {
  real_wfx<<(int)Nu_charge.at(icenter)<<endl; 
  imag_wfx<<(int)Nu_charge.at(icenter)<<endl; 
 }
 line="</Atomic numbers>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Charges>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(icenter=0;icenter<Nu_charge.size();icenter++)
 {
  real_wfx<<Nu_charge.at(icenter)<<endl; 
  imag_wfx<<Nu_charge.at(icenter)<<endl; 
 }
 line="</Nuclear Charges>"; // Faked
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Nuclear Cartesian Coordinates>"; // Only C1 symmetry in DIRAC 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(icenter=0;icenter<Ncenters;icenter++)
 {
  real_wfx<<setw(24)<<Center_x[icenter]<<endl; // x
  imag_wfx<<setw(24)<<Center_x[icenter]<<endl; 
  real_wfx<<setw(24)<<Center_y[icenter]<<endl; // y
  imag_wfx<<setw(24)<<Center_y[icenter]<<endl; 
  real_wfx<<setw(24)<<Center_z[icenter]<<endl; // z
  imag_wfx<<setw(24)<<Center_z[icenter]<<endl; 
 }
 line="</Nuclear Cartesian Coordinates>"; // Only C1 symmetry in DIRAC 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Number of Primitives>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 real_wfx<<Nprimitives<<endl;
 imag_wfx<<Nprimitives<<endl;
 line="</Number of Primitives>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Centers>"; // Only atomic systems 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<Prim2Center_map[iprim]<<endl;
  imag_wfx<<Prim2Center_map[iprim]<<endl;
 }
 line="</Primitive Centers>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<shell_types[iprim]<<endl;
  imag_wfx<<shell_types[iprim]<<endl;
 }
 line="</Primitive Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Primitive Exponents>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(iprim=0;iprim<Nprimitives;iprim++)
 {
  real_wfx<<setw(24)<<prim_exponents[iprim]<<endl;
  imag_wfx<<setw(24)<<prim_exponents[iprim]<<endl;
 }
 line="</Primitive Exponents>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Occupation Numbers>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(imos=0;imos<NMOs_occ;imos++)
 {
  real_wfx<<setw(24)<<MOsLS_occ[imos]<<endl;
  imag_wfx<<setw(24)<<MOsLS_occ[imos]<<endl;
 }
 line="</Molecular Orbital Occupation Numbers>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Spin Types>"; // Faked spin
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(imos=0;imos<NMOs_occ;imos++)
 {
  line=" Alpha";
  real_wfx<<line<<endl;
  imag_wfx<<line<<endl;
 }
 line="</Molecular Orbital Spin Types>";
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 line="<Molecular Orbital Primitive Coefficients>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 for(imos=0;imos<NMOs_occ;imos++)
 {
  line="</MO Number>";
  real_wfx<<line<<endl;
  imag_wfx<<line<<endl;
  for(iprim=0;iprim<Nprimitives;iprim++)
  {
   real_wfx<<setw(24)<<Prim2MO_Coef[imos][iprim].real()<<endl;
   imag_wfx<<setw(24)<<Prim2MO_Coef[imos][iprim].imag()<<endl;
  }
 }
 line="</Molecular Orbital Primitive Coefficients>"; 
 real_wfx<<line<<endl;
 imag_wfx<<line<<endl;
 imag_wfx.close(); 
 real_wfx.close(); 
}
