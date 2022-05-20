#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define BUFSIZE 1024  /** Size of buffer for IO **/
#define MAXATOM 256  /** Maximum number of atoms **/
#define RBOND 2.0     /** Maximum bond length in Angstroms **/
#define MAXELE 10     /** Maximum number of elements **/

int not_in_molecule(int isearch , int mol_list[MAXATOM][MAXATOM],
		    int nmolecule) ;
void find_molecules(double *x, double *y, double *z, 
  		    double a, double b, double c,
		    int mol_list[MAXATOM][MAXATOM], 
		    int natom, int *nmolecule,
		    int *type, char *element,
		    int bond_list[MAXATOM][MAXATOM]) ;
int is_bonded(int j, int k, double *x, double *y, double *z,
	      double a, double b, double c,
	      int *type, char *element) ;
void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a,  double b, double c,
		   int mol_list[MAXATOM][MAXATOM],
		   int natom,
		   int nmolecule,
		   int *type,
		   char *element,
		   int bond_list[MAXATOM][MAXATOM]) ;
void wrap_atoms(double *x, double *y, double *z, 
  		    double a, double b, double c,
		    int bond_list[MAXATOM][MAXATOM], 
		int natom)  ;
void wrap_atom(int j, int k, double *x, double *y, double *z,
	      double a, double b, double c) ;
void wrap_com(double *x, double *y, double *z, 
	      double a, double b, double c,
	      int mol_list[MAXATOM][MAXATOM], int nmolecule,
	      int natom) ;

void center_all(double *x, double *y, double *z, 
	      double a, double b, double c,
	      int mol_list[MAXATOM][MAXATOM], int nmolecule) ;
double r2bond(char ea, char eb) ;
int nint(double num) ;
void print_molecule(int mol_list[MAXATOM][MAXATOM],
		    int bond_list[MAXATOM][MAXATOM],
		    int *type, char *element, int j) ;
main() {
  char buf[BUFSIZE] ;
  double x[MAXATOM], y[MAXATOM], z[MAXATOM] ;
  double a, b, c, junk1, junk2 ;
  int nmolecule, do_wrap_atoms, do_wrap_molecules, print_copies ;
  int natom, i, j, k, type[MAXATOM], mol_list[MAXATOM][MAXATOM],
    bond_list[MAXATOM][MAXATOM] ;
  char element[MAXELE] ;
  char pbc ;
  const double bohr = 0.529177 ;


  do_wrap_atoms = 1 ;  /** DO NOT wrap atoms **/
  do_wrap_molecules = 1 ;  /** DO NOT wrap molecules 
			       (requires do_wrap_atoms = 1 **/

  print_copies = 1 ;   /** Can choose 1-5 **/


  while ( fgets(buf,BUFSIZ,stdin) ) {
    if ( buf[0] != '#' ) {
      break ;
    }
  }
  sscanf(buf,"%d", &natom) ;
  fgets(buf,BUFSIZ,stdin) ; 
  for ( j = 0 ; buf[j] != 0 ; j++ ) {
    if ( buf[j] == '\'' ) {
      buf[j] = ' ' ;
    }
  }
  /** Parse the input line, get the element types **/
  k = 0 ;
  for ( j = 0 ; buf[j] ; j++ ) {
    if ( isalpha(buf[j]) ) {
      element[k] = buf[j] ;
      k++ ;
    }
  }

  for ( j = 0 ; j < natom ; j++ ) {
    scanf("%d %d %lf %lf %lf", &i, &type[j], &x[j], &y[j], &z[j]) ;
  }
  fgets(buf,BUFSIZ,stdin) ; 
  fgets(buf,BUFSIZ,stdin) ;

  scanf("%lf %lf %lf", &a, &junk1, &junk2) ;
  scanf("%lf %lf %lf", &junk1, &b, &junk2) ;
  scanf("%lf %lf %lf", &junk1, &junk2, &c) ;

  if ( a <= 0.0 || b <= 0.0 || c <= 0.0 ) {
    printf("Bad cell size read in\n") ;
    exit(1) ;
  }
  printf("a = %f, b = %f, c = %f\n", a, b, c) ;

  for ( j = 0 ; j < MAXATOM ; j++ ) {
    for ( i = 0 ; i < MAXATOM ; i++ ) {
      mol_list[i][j] = - 1 ;
      bond_list[i][j] = 0 ;
    }
  }
  nmolecule = 0 ;
  /** Each atom is assigned to a molecule **/
  find_molecules(x,y,z,a,b,c,mol_list,natom,&nmolecule,
		 type,element,bond_list) ;

  printf("The number of molecules found = %d\n", nmolecule) ;

  /** Wrap each atom to make whole molecules **/
  if ( do_wrap_atoms ) {
    wrap_atoms(x,y,z,a,b,c,bond_list,natom); 
  }

  /** Wrap the center of mass of each molecule to be inside the box
   **/
  if ( do_wrap_molecules ) {
    wrap_com(x,y,z,a,b,c,mol_list,nmolecule, natom);  
  }
  for ( j = 0 ; j < nmolecule ; j++ ) {
    print_molecule(mol_list,bond_list,type,element,j) ;
  }
}


void find_molecules(double *x, double *y, double *z, 
  		    double a, double b, double c,
		    int mol_list[MAXATOM][MAXATOM], 
		    int natom, int *nmolecule,
		    int *type, char *element, int bond_list[MAXATOM][MAXATOM]) 
{
  int j, k ;
  for ( j = 0 ; j < natom ; j++ ) {
    if ( not_in_molecule(j,mol_list,*nmolecule) ) {
      mol_list[*nmolecule][0] = j ; 
      find_molecule(j,x,y,z,a,b,c,mol_list,natom,*nmolecule,
		    type, element,bond_list) ;
#if(0)
      for ( k = 0 ; mol_list[*nmolecule][k] >= 0 ; k++ ) {
	printf("mol:%d %d %d\n", *nmolecule, k, mol_list[*nmolecule][k]) ;
      }
#endif
      (*nmolecule)++ ;
    }
  }
}



void wrap_atoms(double *x, double *y, double *z, 
  		    double a, double b, double c,
		    int bond_list[MAXATOM][MAXATOM], 
		    int natom) 
{
  int j, k ;
  for ( j = 0 ; j < natom ; j++ ) {
    for ( k = j+1 ; k < natom ; k++ ) {
      if ( bond_list[j][k] == 1 ) {
	wrap_atom(j, k, x,y,z,a,b,c) ;
      }
    }
  }
}

int not_in_molecule(int isearch , int mol_list[MAXATOM][MAXATOM],
		    int nmolecule)
     /** Return 1 if atom isearch is not in the given molecule, 0
	 otherwise **/
{
  int j, k ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      if ( mol_list[j][k] == isearch ) {
	return ( 0 ) ;
      }
    }
  }
  return (1) ;
}

void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a,  double b, double c,
		   int mol_list[MAXATOM][MAXATOM],
		   int natom,
		   int nmolecule,
		   int *type,
		   char *element,
		   int bond_list[MAXATOM][MAXATOM]) 
/** Recursively search for all atoms bonded to j, then all atoms 
    bonded to atoms bonded to j, etc. **/
{
  int k, l ;
  for ( k = 0 ; k < natom ; k++ ) {
    if ( is_bonded(j,k,x,y,z,a,b,c,type,element) ) {
      bond_list[j][k] = 1 ;
      for ( l = 0 ; mol_list[nmolecule][l] >= 0 ; l++ ) {
	if ( mol_list[nmolecule][l] == k ) break ;
      }
      if ( mol_list[nmolecule][l] < 0 ) {
	/** We found a new entry **/
	mol_list[nmolecule][l] = k ;
	find_molecule(k,x,y,z,a,b,c,mol_list,natom,nmolecule,
		      type,element,bond_list) ;
      }
    }
  }
}


int is_bonded(int j, int k, double *x, double *y, double *z,
	      double a, double b, double c,
	      int *type, char *element)
     /** Return 1 if j and k are bonded, 0 if not **/
     /** x, y, and z contain the atomic coordinates. **/
     /** a, b, and c are the rectangular box dimensions **/
{
  double r2 ;
  double rx, ry, rz, dx, dy, dz ;
  
  if ( j == k ) return 0 ;

  rx = x[j] - x[k] ;
  dx = a * nint(rx/a) ;
  rx -= dx ;
  ry = y[j] - y[k] ;
  dy = b * nint(ry/b) ;
  ry -= dy ;
  rz = z[j] - z[k] ;
  dz = c * nint(rz/c) ;
  rz -= dz ;

  r2 = rx * rx + ry * ry + rz * rz ;

#if(0)
  printf("Testing Atoms %d (%c) and %d (%c) for bonding: r = %f\n", 
	 j, element[type[j]-1], k, element[type[k]-1], sqrt(r2) ) ;  
#endif
  if ( r2 < r2bond(element[type[j]-1], element[type[k]-1]) ) {
#if(1)
#endif
    return ( 1) ;
  }
  return(0) ;
}

  
  
void wrap_atom(int j, int k, double *x, double *y, double *z,
	      double a, double b, double c)
     /** Wrap atom k so to be as close as possible to atom j **/
{
  double r2 ;
  double rx, ry, rz, dx, dy, dz ;
  
  rx = x[j] - x[k] ;
  dx = a * nint(rx/a) ;
  x[k] += dx ;

  ry = y[j] - y[k] ;
  dy = b * nint(ry/b) ;
  y[k] += dy ;

  rz = z[j] - z[k] ;
  dz = c * nint(rz/c) ;
  z[k] += dz ;
}



void wrap_com(double *x, double *y, double *z, 
	      double a, double b, double c,
	      int mol_list[MAXATOM][MAXATOM], int nmolecule, int natom)
     /** Wrap each molecule so that it's center of mass (ignoring
	 atomic weight) is inside the box **/
{
  double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0 ;

  int j , k ;
  
  rx0 = 0.0 ;
  ry0 = 0.0 ;
  rz0 = 0.0 ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    rx = 0.0 ;
    ry = 0.0 ;
    rz = 0.0 ;
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      rx += x[mol_list[j][k]] ;
      ry += y[mol_list[j][k]] ;
      rz += z[mol_list[j][k]] ;

      rx0 += x[mol_list[j][k]] ;
      ry0 += y[mol_list[j][k]] ;
      rz0 += z[mol_list[j][k]] ;

    }
    rx /= k ;
    ry /= k ;
    rz /= k ;
    dx = a * (int) (rx / a) ;
    dy = b * (int) (ry / b) ;
    dz = c * (int) (rz / c) ;
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
  }
  dx = rx0 / natom ;
  dy = ry0 / natom ;
  dz = rz0 / natom ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
  }
}

    
double r2bond(char ea, char eb)
     /** Return the square of the bond cutoff distance for 
	 element pair ea and eb.  This saves having to take the
     square root. **/
{
  int i, j ;
  static double rbond[MAXELE][MAXELE] ;
  static int called_before = 0 ;

  if ( ! called_before ) {
    /** Indices: H = 0, C = 1, N = 2, O = 3 **/
    rbond[0][0] = 0.75  ;
    rbond[1][0] = 1.70 ;
    rbond[2][0] = 1.50 ;
    rbond[3][0] = 1.50 ;
    
    rbond[1][1] = 1.70 ;
    rbond[2][1] = 2.20 ;
    rbond[3][1] = 1.80 ;
    
    rbond[2][2] = 1.60 ;
    rbond[3][2] = 1.80 ;
    
    rbond[3][3] = 1.40 ;

    for ( i = 0 ; i < MAXELE ; i++ ) {
      for ( j = i + 1 ; j < MAXELE ; j++ ) {
	rbond[i][j] = rbond[j][i] ;
      }
    }
    called_before = 1 ;
  }

  if ( ea == 'H' ) {
    i = 0 ;
  } else if ( ea == 'C' ) {
    i = 1 ;
  } else if ( ea == 'N' ) {
    i = 2 ;
  } else if ( ea == 'O' ) {
    i = 3 ; 
  } else {
    printf("Error: element %d is unknown\n", ea ) ;
    exit(1) ;
  }

  if ( eb == 'H' ) {
    j = 0 ;
  } else if ( eb == 'C' ) {
    j = 1 ;
  } else if ( eb == 'N' ) {
    j = 2 ;
  } else if ( eb == 'O' ) {
    j = 3 ; 
  } else {
    printf("Error: element %d is unknown\n", eb ) ;
    exit(1) ;
  }

  return(rbond[i][j] * rbond[i][j]) ;
}

  
 
int nint(double num)
/* returns the nearest int to the double num */
{
 
    if (  fabs( num - (int ) num ) <  0.5)
        return( (int ) num ) ;
 
    else if ( num > 0 )
         return( (int ) num + 1 ) ;
    else
         return( (int ) num - 1 ) ;
}


void print_molecule(int mol_list[MAXATOM][MAXATOM],
		    int bond_list[MAXATOM][MAXATOM],
		    int *type, char *element, int mol) 
{
  int k, j, jj, kk  ;
  int ele_count[MAXELE], bond_count[MAXELE][MAXELE] ;

  /** Calculate the overall stoichiometry **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    ele_count[k] = 0 ;
  }
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    ele_count[type[mol_list[mol][k]]-1]++ ;
  }
  for ( k = 0 ; k < MAXELE ; k++ ) {
    if ( ele_count[k] > 0 ) {
      printf("%c%d ", element[k], ele_count[k]) ;
    }
  }

  /** Calculate the number of bonds **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      bond_count[j][k] = 0 ;
    }
  }
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    kk = mol_list[mol][k] ;
    for ( j = 0 ; mol_list[mol][j] >= 0 ; j++ ) {
      jj = mol_list[mol][j] ;
      if ( bond_list[kk][jj] == 1 && jj > kk ) {
	if ( type[kk] > type[jj] ) {
	  bond_count[type[kk]-1][type[jj]-1]++ ;
	} else {
	  bond_count[type[jj]-1][type[kk]-1]++ ;
	}
      }
    }
  }

  for ( k = 0 ; k < MAXELE ; k++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      if ( bond_count[j][k] > 0 ) {
	printf("%d(%c-%c) ", bond_count[j][k],
	       element[j], element[k]) ;
      }
    }
  }
  printf("\n") ;
}
    
