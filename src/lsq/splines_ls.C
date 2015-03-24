#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "functions.h"
using namespace std;

int main(int argc, char* argv[])
{
  if ( argc != 2 ) {
    cout << "Usage: splines_ls <number of frames>\n" ;
    exit(1) ;
  }
  const int nframes=atoi(argv[1]);

  ////Job parameters:

  const int nlayers=1;//supercells adjacent to central cell.
  //1 is enough because of 8-Ang. limit to spline potential; 
  //electrostatic energy doesn't depend on this number.
  //setting to 2 gives identical result.

  // If true, subtract overcoordination forces.
  const bool ifsubtract_coord = false ;  

  // If true, fit coulomb parameters.
  const bool fit_coul = true ;

  // If true, fit overcoordination parameters.
  const bool fit_pover = true ;

  // If true, calculate chebyshev polynomial fit
  const bool if_cheby = true ;

  // Order of Chebyshev polynomial if used.
  const int cheby_order = 12 ;

  ////Read input file input.xyzf:
  ////xyz format in Angs, three box dimensions in the info line.
  ifstream fileread("input.xyzf");

  std::vector<int> Nat;
  Nat.reserve(nframes);
  std::vector<double*> Latcons;
  Latcons.reserve(nframes);

  std::vector<string*> Lb;
  Lb.reserve(nframes);

  std::vector<double**> Coord;
  Coord.reserve(nframes);
  std::vector<double**> Force;
  Force.reserve(nframes);

  std::vector<double***> Als;//least squares A matrix
  Als.reserve(nframes);

  std::vector<double**> Coul_oo;
  Coul_oo.reserve(nframes);
  std::vector<double**> Coul_oh;
  Coul_oh.reserve(nframes);
  std::vector<double**> Coul_hh;
  Coul_hh.reserve(nframes);

  std::vector<double**> Pover ;
  Pover.reserve(nframes);

  //define spline parameters:
  const double smin=0.75 ;
  const double smax=6.0;//spline fitting from 0-8 Angstroms.
  // const double sdelta=0.05;//0.025; //grid spacing
  const double sdelta=0.05 ; //0.025; //grid spacing
  int snum ;

  if ( if_cheby == false ) {
      snum=(1+int((smax-smin)/sdelta))*2*3;//2 is for p0/m0/p1/m1, and 3 is for oo/oh/hh.
  } else {
    snum = cheby_order * 3 ;
  }
  
  std::vector<double*> Mass;
  Mass.reserve(nframes);

  int tempnat ;
  double *tempLatcons;
  string *tempLb;

  double **tempCoord;
  double **tempForce;
  double ***tempA;
 
  double *tempMass;

  double **tempCoul_oo;
  double **tempCoul_oh;
  double **tempCoul_hh;
  double **tempPover ;

  //Each MD snapshot is pushed back in a std::vector.
  //Coordinates, forces etc. held in pointer arrays.
  for(int N=0;N<nframes;N++)
    {
      fileread >> tempnat;
      tempLatcons=new double[3];
      fileread >> tempLatcons[0] >> tempLatcons[1] >> tempLatcons[2];

      tempLb=new string[tempnat];

      tempCoord=new double *[tempnat];
      tempForce=new double *[tempnat];

      for(int a=0;a<tempnat;a++)
	{
	  tempCoord[a]=new double[3];
	  tempForce[a]=new double[3];
	}
      
      for(int a=0;a<tempnat;a++)
	{
	  fileread>>tempLb[a];
	  for(int c=0;c<3;c++)
	    fileread>>tempCoord[a][c];
	  for(int c=0;c<3;c++)
	    fileread>>tempForce[a][c];
	}


      tempA=new double **[tempnat];
      for(int a=0;a<tempnat;a++)
	{
	  tempA[a]=new double *[snum];
	  for(int n=0;n<snum;n++)
	    {
	      tempA[a][n]=new double[3];
	      for(int c=0;c<3;c++)
		tempA[a][n][c]=0;
	    }
	}

      tempMass=new double[tempnat];
      for(int a=0;a<tempnat;a++)
	{
	  if(tempLb[a]=="O")
	    tempMass[a]=15.9994;
	  else if(tempLb[a]=="H")
	    tempMass[a]=1.00794;
	  else
	    {
	      cerr<<"Atom type "<<tempLb[a]<<" not supported."<<endl;
	      exit(1);
	    }
	}

      tempCoul_oo=new double *[tempnat];
      tempCoul_oh=new double *[tempnat];
      tempCoul_hh=new double *[tempnat];
      tempPover = new double* [tempnat] ;
      for(int a=0;a<tempnat;a++)
	{
	  tempCoul_oo[a]=new double[3];
	  tempCoul_oh[a]=new double[3];
	  tempCoul_hh[a]=new double[3];
	  tempPover[a] = new double [3] ;

	  for(int c=0;c<3;c++)
	    {
	      tempCoul_oo[a][c]=0.0;
	      tempCoul_oh[a][c]=0.0;
	      tempCoul_hh[a][c]=0.0;

	      tempPover[a][c] = 0.0 ;
	    }
	}
      
      Nat.push_back(tempnat);
      Latcons.push_back(tempLatcons);
      Lb.push_back(tempLb);
      Coord.push_back(tempCoord);
      Force.push_back(tempForce);
      Als.push_back(tempA);
      Mass.push_back(tempMass);
      Coul_oo.push_back(tempCoul_oo);
      Coul_oh.push_back(tempCoul_oh);
      Coul_hh.push_back(tempCoul_hh);
      Pover.push_back(tempPover) ;
    }


  ////Wrap coordinates:
  for(int N=0;N<nframes;N++)
    {
      for(int a1=0;a1<Nat[N];a1++)
	for(int c=0;c<3;c++)
	  Coord[N][a1][c] -= floor(Coord[N][a1][c]/Latcons[N][c])*Latcons[N][c];
    }
  ////Convert forces from Hartree/bohr to kcal/Angs (Stillinger's  units).
  for(int N=0;N<nframes;N++)
    for(int a1=0;a1<Nat[N];a1++)
      for(int c=0;c<3;c++)
	Force[N][a1][c] *= 627.50960803*1.889725989;

 
  //Generate A matrix, b vector, for least-squares:

  // cout.precision(10);
  cout.precision(10);

  ////ZCalc_Deriv is the function that calculates elements of A.
  for(int N=0;N<nframes;N++) 
    {
      ZCalc_Deriv(Coord[N],Lb[N],Latcons[N],nlayers,Nat[N],Als[N],
		  smin,smax,sdelta,snum,Coul_oo[N],Coul_oh[N],Coul_hh[N],
		  if_cheby);

      // Subtract over-coordination forces from force to be output.
      if ( ifsubtract_coord ) 
	{
	  SubtractCoordForces(Coord[N], Force[N], Lb[N], Latcons[N], 
			      nlayers, Nat[N], false, Pover[N]) ;
	}
      if ( fit_pover ) 
	// Fit the overcoordination parameter.
	{
	  SubtractCoordForces(Coord[N], Force[N], Lb[N], Latcons[N], 
			      nlayers, Nat[N], true, Pover[N]) ;
	}
    }


  ////Print out least squares matrix:
  //Ax=b.
  //A printed to "A.txt", b printed to "b.txt"
  //use these files for python SVD routine.
  ofstream fileA("A.txt");
  ofstream fileb("b.txt");

  //  Reduced precision to 6 for code testing.
  fileA.precision(6);
  fileb.precision(6);

  for(int N=0;N<nframes;N++)
    for(int a=0;a<Nat[N];a++)
      {
	//print Afile:
	for(int c=0;c<3;c++)
	  {
	    //Afile
	    for(int n=0;n<snum;n++) 
	      {
		//	    if ( n % 2 == 1 ) 
		//{
		//Als[N][a][n][c] *= 100.0 ;
		//}
		fileA<<Als[N][a][n][c]<<"   ";
	      }
	  
	    if ( fit_coul ) 
	      {
		fileA<<Coul_oo[N][a][c]<<"   "<<Coul_oh[N][a][c]<<"   "<<Coul_hh[N][a][c];
	      }
	    if ( fit_pover ) 
	      {
		fileA << " " << Pover[N][a][c]  ;
	      }
	    fileA<<endl;
	  }


	  ////print bfile:
	  for(int c=0;c<3;c++)
	    {
	      fileb<<Force[N][a][c]<<endl;
	    }
      }

  //Afile charge conservation (q00 + 4 qOH + qHH = 0 )
  if ( fit_coul ) {
    for(int n=0;n<snum;n++)
      fileA<<"0.0 ";

    fileA<<"1000.0 4000.0 4000.0"; // 0.0"<<endl;
    if ( fit_pover ) 
      {
	fileA << " 0.0 " ;
      }
    fileA << endl ;

    ////bfile charge conservation
    fileb<<"0.0"<<endl;


    //Afile charge consistency (qOO - 4 qHH = 0)
    for(int n=0;n<snum;n++)
      fileA<<"0.0 ";
    fileA<<"1000.0 0000.0 -4000.0" ;
    ////bfile charge consistency
    fileb<<"0.0"<<endl;

    if ( fit_pover ) 
      {
	fileA << " 0.0 " ;
      }
    fileA << endl ;


  } 

  //cout<<endl<<endl;
  ofstream header("params.header") ;
  if ( ! if_cheby ) {
    header << smin << "\n" << smax << "\n" << sdelta << "\n" ;
  } else {
    header << smin << "\n" << smax << "\n" << snum/3 << "\n" ;
  }
  
  return 0;    
}


