/********************************************************************

A simple molecular analyzer for DFTB simulations.

usage: molanal.new <gen file> <chr file>

The chr (charge) file is optional.
The gen (coordinate) file is mandatory.

Non-orthorhombic simulation cells are supported.

Output:  Molecules found are printed out with total charges,
         if supplied.
         An xyz file suitable for jmol is created in the
	 file molanal.xyz.

Bond distances are now read from the file bonds.dat.

Larry Fried, 5/30/2005.

A bond lifetime criterion is now implemented.  Fixed a bug in 
dipole calculations when wrap_com was not called.

Larry Fried, 11/10/06

Implemented improved algorithms for large systems and OPENMP.

Larry Fried, circa 2012

Fixed a bug in OPENMP when maxbond was exceeded.

Larry Fried, 06/12/13

Fixed a bug in history mechanism (history variables searched incorrectly),

Fixed a bug in history mechanism (history of a,b,c not saved).

Implemented a new counter-based bond lifetime algorithm from L. Koziol.

Larry Fried, 07/23/13

*********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <omp.h>
#include <dirent.h>

/***   These parameters could be changed, depending on the application. */

#define MAXWRAP 2             /* Maximum number of lattice vectors to search over.
		                 Use a larger number for a highly non-orthogonal cell, 
                                 or for long chain molecules. */
#define MAXELE 6              /* Maximum number of elements.  Only change when adding a new 
                                 element. */

#define MIN_BOND_LENGTH 0.2   /* Minimum allowed bond length */

#define MAXCOORD  8           /* Maximum usual coordination of an atom */

#define BUFSIZE 1024          /* Size of buffer for IO */

#undef SPEW                   /* Define for extra output. */

#define HAVE_NEARBYINT        /* If defined, use C99 nearbyint 
                                 builtin function */

/***  The below parameters should not be changed ***/

#define NBOND_IDX 1           /* Index to store number of bonds for an atom */

#define DIM_IDX 0             /* Index to store dimension for bond array for an atom */

#define BOND_IDX 2            /* Index to storing bonding list */


#ifdef HAVE_NEARBYINT
#define nint nearbyint
#endif

void find_molecules(double *x, double *y, double *z, 
  		    double a[3], double b[3], double c[3],
		    int **mol_list, 
		    int natom, int *nmolecule,
		    int *type, char **element, 
                    int **bond_list,
                    int **bond_duration,
                    int bond_time[MAXELE][MAXELE],
		    const int **neighbor, int orthogonal) ;
int is_bonded(int j, int k, int *type, int bond_time[MAXELE][MAXELE],
	      int **bond_list) ;
void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a[3],  double b[3], double c[3],
		   int **mol_list,
		   int natom,
		   int nmolecule,
		   int *type,
		   char **element,
		   int **bond_list,
		   int **bond_duration,
                   int bond_time[MAXELE][MAXELE],
		   const int **neighbor,
		   int orthogonal,
		   int *assigned) ;
double wrap_atom(int j, int k, double *x, double *y, double *z,
		 double a[3], double b[3], double c[3], int maxwrap,
		 int do_wrap, int orthogonal) ;
void wrap_molecules(double *x, double *y, double *z, 
		    double a[3], double b[3], double c[3],
		    int **mol_list, int nmolecule,
		    int **bond_list, int natom, int orthogonal );
void wrap_com(double *x, double *y, double *z, 
	      double a[3], double b[3], double c[3],
	      int **mol_list, int nmolecule, int natom,
	      double *rcomx, double *rcomy, 
	      double *rcomz, int *type,
	      char *element[MAXELE], int orthogonal) ;
int read_charge_file(FILE *fchr, double *q) ;
int read_gen_file(FILE *fgen, char *element[MAXELE], 
		  char tmp_ele[MAXELE][MAXELE],
		  int *maxele,
		  int *type, int *maxtype,
		  double *x, double *y, double *z,
		  double *a,  double *b,  double *c) ;

#ifndef HAVE_NEARBYINT
double nint(double num) ;
#endif

void print_molecule(int **mol_list,
		    int **bond_list,
		    int *type, char **element, int mol, 
		    int read_charge,
		    double *q,
		    double *x, double *y, double *z,
		    double *rcomx, double *rcomy,
		    double *rcomz, int *printed_atom,
		    int dump_structures)  ;
void inversebox(double a[3], double b[3], double c[3],
		double invbox[3][3])  ;
static void wrap_in_box(double a[3], double b[3], double c[3],
			double invbox[3][3],
			int natom, double *x, double *y, 
			double *z)  ;
void wrap_pairs(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int **mol_list, int nmolecule,
		int **bond_list, int natom, int j, int k,
		int *wrapped, int orthogonal) ;
int ele_index(char *ea);
double atm_mass(char *ea) ;
void read_bonds(double rbond[MAXELE][MAXELE], int *xyz_copies, double *time_step,
		int *dump_structures, int bond_time[MAXELE][MAXELE]) ;
double boxvol(double a[3], double b[3], double c[3]) ;
int get_next_neighbor(int j, int *k, int *nid, const int **neighbor) ;
void find_neighbors(double *x,double *y, double *z, int *type,
		    char **element,
		    double a[3], double b[3], double c[3],
		    int **neighbor,int natom, int orthogonal,
		    int **bond_list) ;
int check_orthogonal_box(double a[3],double b[3], double c[3]) ;
void free_neighbor(int **neighbor, int natom) ;
void find_neighbors_slow(double *x,double *y, double *z, int *type,
			 char **element,
			 double a[3], double b[3], double c[3],
			 int **neighbor,int natom, int orthogonal,
			 int **bond_duration,double minimum_distance) ;
void find_neighbors_fast(double *x,double *y, double *z, int *type,
			 double a[3], double b[3], double c[3],
			 int **neighbor,int natom, int orthogonal,
			 int **bond_duration, double minium_distance) ;
int comp_int(const void *pi, const void *pj) ;
void sort_neighbor(int **neighbor, int natom);
double **dmatrix(int m, int n) ;
int **imatrix(int m, int n) ;
double *dvector(int m) ;
int *ivector(int m) ;
void allocate_arrays(double **x, double **y, double **z, double **q, 
		     double **rcomx, double **rcomy, double **rcomz, 
		     int **type0, int **type, int **printed_atom, 
		     int ***mol_list, int ***bond_list, int ***bond_duration,
		     int ***neighbor, int natom) ;
void free_arrays(double *x, double *y, double *z, double *q, 
		 double *rcomx, double *rcomy, double *rcomz, 
		 int *type0, int *type, int *printed_atom, 
		 int **mol_list, int **bond_list, 
		 int **neighbor, int natom, int nmolecule) ;
void add_bond(int j, int k, int **bond_list);
void remove_bond(int j, int k, int **bond_list) ;
int get_bond_list(int **bond_list, int j, int k) ;
void initialize_bonds(int **mol_list, int **bond_list, int natom, int *nmolecule) ;
void free_dmatrix(double **mat) ;
void free_imatrix(int **mat) ;
char *name_molecule(int **mol_list,
		    int **bond_list,
		    int *type, char **element, int mol) ;
void print_com(FILE *fcom, double *rcomx, double *rcomy, double *rcomz, char **mol_name, 
	       int nmolecule, int frame) ;
void com_molecule(int **mol_list, int *type, const double *x, 
		  const double *y, const double *z,
		  char *element[MAXELE],
		  int j, double *rxx, double *ryy, double *rzz);
void check_frame_input(int read_charge, int natom, int *natom0, 
		       int natom1, int frame, int maxele, int *maxele0,
		       int *type, int *type0) ;
int skip_frame(int frame, int skip1, int skip2, int block, int offset,
	       int *block_start);
void read_skip(int *skip1, int *skip2, int *block, int *offset);
int get_bond_duration(int **bond_list, int j, int k) ;
void reset_bonds(int **bond_list, int natom) ;

double rbond[MAXELE][MAXELE] ;       /* Bond distances indexed by internal element index */
double rbond2[MAXELE][MAXELE] ;      /* Squared bond distances indexed by internal element index */
double rbond_type[MAXELE][MAXELE] ;  /* Bond distances indexed by input element type */
int bond_index[MAXELE] ;             /* Index used in storing bond times and distances */
int bond_time[MAXELE][MAXELE] ;      /* Lifetime criterion used for bonding */
int sparse_bonds = 1 ;               /* Whether to use sparse storage for bond array */
int max_molecule_size = 2 ;          /* Maximum molecule size.  Dynamically updated as needed. */

int main(int argc, char** argv) 
{
  char buf[BUFSIZE] ;
  double *x, *y, *z ;                    /* Atomic coordinates of the latest frame */
  double *q ;                            /* Atomic charge of the latest frame */
  double *rcomx, *rcomy, *rcomz ;        /* Molecule center of mass coordinates */
  double a[3], b[3], c[3] ;              /* Unit cell vectors */
  double invbox[3][3] ;                  /* Inverse of unit cell vectors */
  int nmolecule = 0 ;                    /* The number of molecules */
  int natom ;                            /* The number of atoms */
  int natom0 = 0 ;                       /* The initial number of atoms */
  int i, j ;
  int **mol_list ;                       /* The list of atoms in a molecule */
  int **bond_list ;                      /* The list of atoms currently considered bonded according
                                            to the bond lifetime.  */
  int **bond_duration ;                  /* The list of bond durations */
  int *type ;                            /* Element type of each atom */
  int *type0 ;                           /* Element type of each atom */
  
  char *element[MAXELE] ;
  char tmp_ele[MAXELE][MAXELE] ;
  char **mol_name ;
  
  int **neighbor ; /* A neighbor list */
  FILE *fout ;
  FILE *fgen ;            /* The gen (coordinate) file. */
  FILE *fcom ;            /* The center of mass coordinate file. */
  FILE *fkl  ;				// Output file for our analysis (RKL)
  FILE *fchr = NULL ;     /* The chr (charge) file */
  int read_charge = 0 ;   /* Whether charges should be read. */
  int natom1 ;            /* Number of atoms in the chr file */
  double dx, dy, dz ;
  int xyz_copies = 0 ;  /* Print this many copies in each direction for the xyz file. */
  int frame = 1 ;
  int ix, iy, iz ;
  int *printed_atom ;         /* Whether an atom was printed out. */
  double vbox ;               /* The volume of the box. */
  double time_step ;          /* The time step between saved frames. */
  int dump_structures = 0 ;   /* Whether to dump structures. */
  int bond_time_max ;
  int orthogonal = 0 ;        /* Whether the unit cell is orthogonal */
  int maxtype ;
  int maxele, maxele0 = 0 ;
  int skip1, skip2 ;
  int block ;
  int offset ;
  int block_start = 1 ;

  /* Defaults: do not skip frames */
  skip1 = 1 ;
  skip2 = 1 ;
  block = 1 ;
  offset = 0 ;

  printf("Molecule analyzer Revision: $Revision: 130 $\n") ;
  printf("Program last changed on $Date: 2013-07-24 15:00:00 -0700 (Wed, 24 Jul 2013) $\n") ;
  printf("Program last changed by $Author: fried1 $\n") ;


  read_bonds(rbond, &xyz_copies, &time_step, &dump_structures, bond_time) ;

  read_skip(&skip1, &skip2, &block, &offset) ;

  printf("The time between frames read = %11.7e seconds\n", time_step) ;
  printf("Every %d frames will be skipped\n", skip1) ;
  printf("Every %d blocks of length %d will be skipped\n",
	 skip2, block) ;
  printf("Block skip offset = %d\n", offset) ;

  if ( dump_structures ) {
    DIR *pdir ;

    pdir = opendir("molecules") ;
    if ( pdir == NULL ) {
      printf("Creating a new molecules directory\n") ;
      system("mkdir molecules") ;
      
    } else {
      closedir(pdir) ;
      printf("Writing into existing molecules directory\n") ;
    }
  }

  bond_time_max = -1 ;
  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      bond_time[i][j] = bond_time[i][j] / skip1 ;
      if (bond_time[i][j] > bond_time_max) {
        bond_time_max = bond_time[i][j] ;
      }
    }
  }
  if ( bond_time_max * skip1 > block && skip2 > 1 ) {
    printf("Error:  block length must be larger than the maximum bond duration\n") ;
    exit(1) ;
  }

  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      /* Squared bond length */
      rbond2[i][j] = rbond[i][j] * rbond[i][j] ;
    }
  }
  
  printf ("max. bond time cutoff (time steps) = %d\n", bond_time_max);

  if ( argc >= 2 ) {
    printf ("Trajectory file: %s\n",argv[1]);
    fgen = fopen(argv[1],"r") ;
    if ( fgen == NULL ) {
      printf("Error: could not open file %s\n", argv[1]) ;
      exit(1) ;
      
    }
    /* Peak ahead to get the number of atoms for allocations */
    while ( fgets(buf,BUFSIZ,fgen) ) {
      if ( buf[0] != '#' ) {
	break ;
      }
    }
    if ( feof(fgen) || ferror(fgen) ) {
      printf("Error: could not read number of atoms\n") ;
      exit(1) ;
    }
    sscanf(buf,"%d", &natom) ;
    fclose(fgen) ;
    fgen = fopen(argv[1],"r") ;
  } else {
    printf("Error: The input file must be specified\n") ;
    exit(1) ;
  }
  if ( argc == 3 ) {
    read_charge = 1 ;
    printf ("Charge file: %s\n",argv[2]);
    fchr = fopen(argv[2],"r") ;
    if ( fchr == NULL ) {
      printf("Error: could not open file %s\n", argv[2]) ;
      exit(1) ;
    } 
  }
  fout = fopen("molanal.xyz" , "w") ;
  fcom = fopen("molanal.com" , "w") ;
  fkl  = fopen("kl-molan.dat", "w") ;

  allocate_arrays(&x, &y, &z, &q, &rcomx, &rcomy, &rcomz, 
		  &type0, &type, &printed_atom,
		  &mol_list, &bond_list, &bond_duration, &neighbor, 
		  natom) ;
  
  while ( ! feof(fgen) || ferror(fgen) ) {
    /** Defensive programming. */

    for ( j = 0 ; j < natom ; j++ ) {
      x[j] = 0.0 ;
      y[j] = 0.0 ;
      z[j] = 0.0 ;
      q[j] = 0.0 ;
      rcomx[j] = -100.0 ;
      rcomy[0] = -100.0 ;
      rcomz[j] = -100.0 ;
      type[j] = 0 ;
    }
    for ( j = 0 ; j < MAXELE ; j++ ) {
      element[j] = NULL ;
    }

    if ( read_charge ) {
      natom1 = read_charge_file(fchr, q) ;
    } else {
      natom1 = 0 ;
    }

    natom = read_gen_file(fgen, element, tmp_ele, &maxele,
			  type, &maxtype, 
			  x, y, z,
			  a, b, c) ;
    
    if ( frame > 1 && natom <= 0 ) {
      break ;
    }

    check_frame_input(read_charge, natom, &natom0, natom1,
		      frame, maxele, &maxele0, type, type0) ;
    
    if ( skip_frame(frame, skip1, skip2, block, offset,&block_start) ) {
      frame++ ;
      continue ;
    }

    for ( j = 0 ; j < maxtype ; j++ ) {
      bond_index[j] = ele_index(element[j]) ;
    }

    /** Wrap all coordinates back into the primitive unit cell. */
    orthogonal = check_orthogonal_box(a,b,c) ;

#if(0)
    if ( orthogonal ) {
      printf("Orthogonal unit cell was detected\n") ;
    }
    else {
      printf("A non-orthogonal unit cell was detected\n") ;
    }
#endif
    
    inversebox(a,b,c,invbox) ;
    wrap_in_box(a,b,c,invbox,natom,x,y,z) ;
    vbox = boxvol(a,b,c) ;

    if ( frame == block_start ) {
      /* Starting a new block.
	 Clear position history from the past block */
      reset_bonds(bond_duration, natom) ;
    }

    find_neighbors(x,y,z,type,element,a,b,c,neighbor,natom,orthogonal,
		   bond_duration) ;
    
    initialize_bonds(mol_list, bond_list, natom, &nmolecule) ;

    /** Each atom is assigned to a molecule **/
    find_molecules(x,y,z,a,b,c,mol_list,natom,&nmolecule,
		   type,element,bond_list,bond_duration,bond_time,
		   (const int**) neighbor, 
		   orthogonal) ;

    if ( frame - block_start >= bond_time_max * skip1 ) {
      printf("Beginning frame %d\n", frame) ;
      printf("The box volume = %11.7e\n", vbox) ;
      printf("The number of molecules found = %d\n", nmolecule) ;
    }

    /** Wrap each atom to make whole molecules **/
    wrap_molecules(x, y, z, a, b, c, mol_list, nmolecule,bond_list, natom,
		   orthogonal) ;
    
    /** Wrap the center of mass of each molecule to be inside the box
     **/
    wrap_com(x,y,z,a,b,c,mol_list,nmolecule, natom,
	     rcomx, rcomy, rcomz,type,element,orthogonal);  

#ifdef SPEW
    for ( j = 0 ; j < natom ; j++ ) {
      printf("New X: %12.6e %12.6e %12.6e\n",
	     x[j], y[j], z[j]) ;
    }
#endif

    /* Find a descriptive name for each molecule */
    mol_name = malloc( nmolecule * sizeof(char**) ) ;
    for ( j = 0 ; j < nmolecule ; j++ ) {
      mol_name[j] = name_molecule(mol_list, bond_list, type, element, j) ;
    }
    
    if ( frame - block_start >= bond_time_max * skip1 ) {
      for ( j = 0 ; j < natom ; j++ ) {
	printed_atom[j] = 0 ;
      }
      for ( j = 0 ; j < nmolecule  ; j++ ) {
	printf("Beginning molecule %d\n",j+1) ;
	print_molecule(mol_list,bond_list,type,element,j,read_charge,q,
		       x,y,z,rcomx,rcomy,rcomz,printed_atom,dump_structures) ;
	printf("Ending molecule %d\n",j+1) ;
      }
      for ( j = 0 ; j < natom ; j++ ) {
	if ( printed_atom[j] == 0 ) {
	  printf("Error: did not print out atom %d\n", j) ;
	  exit(1) ;
	}
      }
      fprintf(fout,"%d\n\n", natom * (2*xyz_copies+1) * (2*xyz_copies+1) * (2*xyz_copies+1)) ;
      for ( ix = -xyz_copies; ix <= xyz_copies; ix++ ) {
	for ( iy = -xyz_copies; iy <= xyz_copies; iy++ ) {
	  for ( iz = -xyz_copies; iz <= xyz_copies; iz++ ) {
	    dx = ix * a[0] + iy * b[0] + iz * c[0] ;
	    dy = ix * a[1] + iy * b[1] + iz * c[1] ;
	    dz = ix * a[2] + iy * b[2] + iz * c[2] ;
	    for ( j = 0 ; j < natom ; j++ ) {
	      fprintf(fout,"%s %f %f %f\n", element[type[j]-1], 
		      dx + x[j], dy + y[j], dz + z[j]) ;
	    }
	  }
	}
      }
      printf("Ending frame %d\n", frame) ;

      print_com(fcom, rcomx, rcomy, rcomz, mol_name, nmolecule, frame) ;
	  
    double mass;
    int    atmidx = 1;	
  	int    atmaccess;
	
	fprintf(fkl, "%f %f %f \n", a[0], b[1], c[2]);
	
  	for (int kli=0; kli<nmolecule; kli++)					// Loop over molecules?
  	{
  		for (int klj=0; mol_list[kli][klj] >= 0; klj++ )	// Loop over atoms in molecules?
  		{
	      atmaccess = mol_list[kli][klj] ;
		
		    if (      strcmp(element[type[atmaccess]-1], "C") ==0 ) 
				mass = 12.011;
			else if ( strcmp(element[type[atmaccess]-1], "H") ==0 )
				mass = 1.0079;
			else if ( strcmp(element[type[atmaccess]-1], "O") ==0 ) 
				mass = 15.999;
			else if ( strcmp(element[type[atmaccess]-1], "N") ==0 ) 
				mass = 14.007;
				
  			fprintf(fkl,"%i %i %s %i %s %f %f %f %f \n",	
  			atmidx++,
			atmaccess,
  			element[type[atmaccess]-1],
  			kli,
  			mol_name[kli],
  		 	x[atmaccess], y[atmaccess], z[atmaccess],
  			mass);
  		}
  	}
    }

    /* Free molecule names */
    for ( j = 0 ; j < nmolecule ; j++ ) {
      free(mol_name[j]) ;
    }
    free(mol_name) ;

    if ( feof(fgen) || ferror(fgen) ) {
      break ;
    }
    frame++ ;

    free_neighbor(neighbor, natom) ;

  }

  free_arrays(x, y, z, q, rcomx, rcomy, rcomz, type0, type, printed_atom,
	      mol_list, bond_list, neighbor, natom, nmolecule) ;

  fclose(fout) ;
  return(0) ;
}


void find_molecules(double *x, double *y, double *z, 
  		    double a[3], double b[3], double c[3],
		    int **mol_list, 
		    int natom, int *nmolecule,
		    int *type, char **element, 
                    int **bond_list,
                    int **bond_duration,
                    int bond_time[MAXELE][MAXELE],
		    const int **neighbor, int orthogonal) 
{
  int j, k, l ;
  int nmol ;
  int *assigned ;  /* Which molecule is assigned to the given atom */
  
  assigned = ivector(natom) ;
  for ( j = 0 ; j < natom ; j++ ) {
    assigned[j] = -1 ;
  }

  /* This loop is hard to write as parallel OMP */
  for ( j = 0 ; j < natom ; j++ ) {
    if ( assigned[j] == -1 ) {
      mol_list[*nmolecule] = calloc(max_molecule_size, sizeof(int)) ;
      for ( l = 1 ; l < max_molecule_size ; l++ ) {
	mol_list[*nmolecule][l] = -1 ;
      }
      mol_list[*nmolecule][0] = j ; 
      assigned[j] = *nmolecule ;
      find_molecule(j,x,y,z,a,b,c,mol_list,natom,*nmolecule,
		    type, element,bond_list,bond_duration,
		    bond_time,
		    neighbor,orthogonal,assigned) ;
#ifdef SPEW
 {  
   int k ;
      for ( k = 0 ; mol_list[*nmolecule][k] >= 0 ; k++ ) {
	printf("mol:%d %d %d\n", *nmolecule, k, mol_list[*nmolecule][k]) ;
      }
 }
#endif
      (*nmolecule)++ ;
    }
  }
  /* Simple sort so that the atom list is unique. */

  nmol = *nmolecule ;
#ifdef OPENMP
#pragma omp parallel for default(none) \
  shared(mol_list, nmol) \
  private(j,k)
#endif
  for ( j = 0 ; j < nmol ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      ;
    }
    qsort(mol_list[j],k,sizeof(int),comp_int) ;
  }

  free(assigned) ;
}


void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a[3],  double b[3], double c[3],
		   int **mol_list,
		   int natom,
		   int nmolecule,
		   int *type,
		   char **element,
		   int **bond_list,
		   int **bond_duration,
                   int bond_time[MAXELE][MAXELE],
		   const int **neighbor,
		   int orthogonal,
		   int *assigned) 
/** Recursively search for all atoms bonded to j, then all atoms 
    bonded to atoms bonded to j, etc. **/
{
  int k, l, m, n ;
  int nid = 0 ;
  
  k = -1 ;
  while ( get_next_neighbor(j,&k, &nid, neighbor) ) {
    if ( is_bonded(j,k,type,bond_time,bond_duration) ) {

#ifdef SPEW
      printf("FIND_MOLECULE: ") ;
#endif
      add_bond(j,k,bond_list) ;

      for ( l = 0 ; mol_list[nmolecule][l] >= 0 ; l++ ) {
	if ( mol_list[nmolecule][l] == k ) break ;
      }
      if ( mol_list[nmolecule][l] < 0 ) {
	/** We found a new entry **/
	mol_list[nmolecule][l] = k ;
	assigned[k] = nmolecule ;
	if ( l + 1 == max_molecule_size ) {
	  for ( m = 0 ; m <= nmolecule ; m++ ) {
	    mol_list[m] = realloc(mol_list[m], 2 * max_molecule_size
				  * sizeof(int) ) ;
	    for ( n = max_molecule_size ; n < 2 * max_molecule_size ; n++ ) {
	      mol_list[m][n] = -1 ;
	    }
	  }
	  max_molecule_size *= 2 ;
	}
	find_molecule(k,x,y,z,a,b,c,mol_list,natom,nmolecule,
		      type,element,bond_list,bond_duration,
		      bond_time,neighbor,orthogonal,assigned) ;
      }
    }
  }
}


int is_bonded(int j, int k, int *type, int bond_time[MAXELE][MAXELE],
	      int **bond_list)
     /** Return 1 if j and k are bonded, 0 if not **/
     /** x, y, and z contain the atomic coordinates. **/
     /** a, b, and c are the rectangular box dimensions **/
{
  int btime , iele, jele, bdur1 ;

  if ( j == k ) return 0 ;

  iele = bond_index[type[j]-1];
  jele = bond_index[type[k]-1];
  btime = bond_time[iele][jele] ;  

  if ( bond_time[iele][jele] != bond_time[jele][iele] ) {
    printf("Error: inconsistent bond time\n") ;
    exit(1) ;
  }
  if ( iele > MAXELE || iele < 0 ) {
    printf("Error: bad bond time index\n") ;
    exit(1) ;
  }
  if ( jele > MAXELE || jele < 0 ) {
    printf("Error: bad bond time index\n") ;
    exit(1) ;
  }
  if ( btime < 0 ) {
    printf("Error: bad bond time\n") ;
    printf("Between: \n");
    printf("%d\n",iele);
    printf("%d\n",jele);
    printf("with btime: \n");
    printf("%d\n",btime);
   exit(1) ;
  }

  bdur1 = get_bond_duration(bond_list,j,k)  ;
  if ( bdur1 < 0 ) {
    printf("Bad bond duration 1\n") ;
    exit(1) ;
  }

#ifdef SPEW
  int bdur2 = get_bond_duration(bond_list,k,j)  ;
  if ( bdur2 < 0 ) {
    printf("Bad bond duration 2\n") ;
    exit(1) ;
  }
  if ( bdur1 != bdur2 ) {
    printf("Inconsistent bond duration\n") ;
    exit(1) ;
  } 
#endif

  if ( bdur1 > btime ) {
#ifdef SPEW
    printf("Atoms %d and %d are bonded\n", j, k) ;
#endif
    return(1) ;
  } else {
#ifdef SPEW
    printf("Atoms %d and %d are not bonded\n", j, k) ;
#endif
    return(0) ;
  }
}

double wrap_atom(int j, int k, double *x, double *y, double *z,
		 double a[3], double b[3], double c[3], int maxwrap,
		 int do_wrap, int orthogonal)
     /** Wrap atom k so to be as close as possible to atom j if do_wrap is 
	 non-zero.
	 Returns the closest distance**2 found.  **/
{
  int l, m, n ;
  double x1[3] ;
  double bestdist = 0.0 , dist ;
  int bestl = 0 , bestm = 0 , bestn = 0 ;

  bestdist = 1.0e20 ;

  if ( ! orthogonal ) {
    for ( l = -maxwrap ; l <= maxwrap ; l++ ) {
      for ( m = -maxwrap ; m <= maxwrap ; m++ ) {
	for ( n = -maxwrap ; n <= maxwrap ; n++ ) {

	  x1[0] = x[j] + l * a[0] + m * b[0] + n * c[0] ;
	  x1[1] = y[j] + l * a[1] + m * b[1] + n * c[1] ;
	  x1[2] = z[j] + l * a[2] + m * b[2] + n * c[2] ;
	
	  dist = (x1[0]-x[k]) * (x1[0] - x[k]) +
	    (x1[1] - y[k]) * (x1[1] - y[k]) +
	    (x1[2] - z[k]) * (x1[2] - z[k]) ;
	  if ( dist < bestdist ) {
	    bestdist = dist ;
	    bestl = l ;
	    bestm = m ;
	    bestn = n ;
	  }
	}
      }
    }
    if ( do_wrap ) {
#if(0)
      printf("bestdist = %f\n", bestdist) ;
      printf("bestl = %d, bestm = %d, bestn = %d\n", bestl, bestm, bestn) ;
      printf("dx = %f\n", bestl * a[0] + bestm * b[0] + bestn * c[0] ) ;
      printf("dy = %f\n", bestl * a[1] + bestm * b[1] + bestn * c[1] ) ;
      printf("dz = %f\n", bestl * a[2] + bestm * b[2] + bestn * c[2] ) ;
#endif
      x[k] -= bestl * a[0] + bestm * b[0] + bestn * c[0] ;
      y[k] -= bestl * a[1] + bestm * b[1] + bestn * c[1] ;
      z[k] -= bestl * a[2] + bestm * b[2] + bestn * c[2] ;
    }
  } else {
    /* Orthogonal unit cell */
    double dx1, dy1, dz1 ;
    double dxb, dyb, dzb ;
    
    /* Use the minimum image convention */
    dx1 = x[j] - x[k] ;
    dxb = a[0] * nint(dx1/a[0]) ;
    dx1 = dx1 - dxb ;
    
    dy1 = y[j] - y[k] ;
    dyb = b[1] * nint(dy1/b[1]) ;
    dy1 = dy1 - dyb ;
    
    dz1 = z[j] - z[k] ;
    dzb = c[2] * nint(dz1/c[2]) ;
    dz1 = dz1 - dzb ;
    
    bestdist = dx1 * dx1 + dy1 * dy1 + dz1 * dz1 ;

    if ( do_wrap ) {
      x[k] += dxb ;
      y[k] += dyb ;
      z[k] += dzb ;
    }
  }
  
  return bestdist ;
}



void wrap_com(double *x, double *y, double *z, 
	      double a[3], double b[3], double c[3],
	      int **mol_list, int nmolecule, int natom,
	      double *rcomx, double *rcomy, 
	      double *rcomz, int *type,
	      char *element[MAXELE], int orthogonal)
     /** Wrap each molecule so that it's center of mass (ignoring
	 atomic weight) is inside the box **/
{
  double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0 ;
  double mass ;
  double x1[2], y1[2], z1[2] ;
  double mass_cell ;
  const double eps = 1.0e-08 ;
  
  int j , k, kk ;
  
  rx0 = 0.0 ;
  ry0 = 0.0 ;
  rz0 = 0.0 ;
  mass_cell = 0.0 ;

  /** Calculate the center of the unit cell, with crystal coordinates
      (0.5, 0.5, 0.5) */
  x1[1] = 0.5 * a[0] + 0.5 * b[0] + 0.5 * c[0] ;
  y1[1] = 0.5 * a[1] + 0.5 * b[1] + 0.5 * c[1] ;
  z1[1] = 0.5 * a[2] + 0.5 * b[2] + 0.5 * c[2] ;
    
#if(0)
#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(rx,ry,rz,mtot,k,kk,mass,	\
          dx,dy,dz,rcomx,rcomy,rcomz) \
  shared(mol_list,element,type,nmolecule,x,y,z,a,b,c, \
         mass_cell,rx0,ry0,rz0,orthogonal,x1,y1,z1)				
#endif
#endif
  
  for ( j = 0 ; j < nmolecule ; j++ ) {

    com_molecule(mol_list, type, x, y, z,element,
		 j, &rx, &ry, &rz) ;
    
    rcomx[j] = rx ;
    rcomy[j] = ry ;
    rcomz[j] = rz ;

    x1[0] = rx ;
    y1[0] = ry ;
    z1[0] = rz ;

    /* Find the wrapped coordinates closest to the center of the unit cell. */
    wrap_atom(1,0,x1,y1,z1,a,b,c,MAXWRAP,1, orthogonal) ;

    dx = rx - x1[0] ;
    dy = ry - y1[0] ;
    dz = rz - z1[0] ;
#if(1)
    rcomx[j] -= dx ;
    rcomy[j] -= dy ;
    rcomz[j] -= dz ;

    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      mass = atm_mass(element[type[mol_list[j][k]]-1]) ;
      /*
#ifdef OPENMP
#pragma omp critical
#endif
      */
      {
	rx0 += x[mol_list[j][k]] * mass ;
	ry0 += y[mol_list[j][k]] * mass ;
	rz0 += z[mol_list[j][k]] * mass ;
	mass_cell += mass ;
      }
    }
#endif
  }
  dx = rx0 / mass_cell ;
  dy = ry0 / mass_cell ;
  dz = rz0 / mass_cell ;
#if(0)
  /* Set cell center of mass to zero */
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
    rcomx[j] -= dx ;
    rcomy[j] -= dy ;
    rcomz[j] -= dz ;
  }
#else
  /* Do not change cell center of mass */
  for ( j = 0 ; j < nmolecule ; j++ ) {
    if ( orthogonal ) {
      if ( rcomx[j] > a[0] || rcomx[j] < -eps ) {
	printf("Error: x com out of box: %f %f\n", rcomx[j], a[0]) ;
	exit(1) ;
      }
      else if ( rcomy[j] > b[1] || rcomy[j] < -eps ) {
	printf("Error: y com out of box: %f %f\n", rcomy[j], b[1]) ;
	exit(1) ;
      }  
      else if ( rcomz[j] > c[2] || rcomz[j] < -eps ) {
	printf("Error: z com out of box: %f %f\n", rcomz[j], c[2]) ;
	exit(1) ;
      }  
    }
  }
#endif

  kk= 0 ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      kk++ ;
    }
  }
  if ( kk != natom ) {
    printf("Error: wrong number of atoms in the molecule list\n") ;
    exit(1) ;
  }

  /* Sanity check: re-compute center of mass and compare to stored value */
  for ( j = 0 ; j < nmolecule ; j++ ) {
    com_molecule(mol_list, type, x, y, z,element,
		 j, &rx, &ry, &rz) ;
    if ( fabs(rx-rcomx[j]) > eps ) {
      printf("Error: inconsistent X com for atom %d\n", j) ;
      printf("RX = %21.13e RCOMX = %12.13e\n",
	     rx, rcomx[j]) ;
      exit(1) ;
    }
    else if ( fabs(ry - rcomy[j]) > eps ) {
      printf("Error: inconsistent Y com for atom %d\n", j) ;
      printf("RY = %21.13e RCOMY = %12.13e\n",
	     ry, rcomy[j]) ;
      exit(1) ;
    } 
    else if ( fabs(rz - rcomz[j]) > eps ) {
      printf("Error: inconsistent Z com for atom %d\n:", j) ;
      printf("RZ = %21.13e RCOMZ = %12.13e\n",
	     rz, rcomz[j]) ;
      exit(1) ;
    } 
  }
}

void wrap_molecules(double *x, double *y, double *z, 
		    double a[3], double b[3], double c[3],
		    int **mol_list, int nmolecule,
		    int **bond_list, int natom, int orthogonal )
     /** Wrap the atoms of each molecule to be as close as possible
	 to bonded partners. **/
{
  /* double rx, ry, rz ;  */
  /* double dx, dy, dz ; */
  /* double rx0, ry0, rz0, mtot ; */
  /* double mass ; */
  /* double x1[2], y1[2], z1[2] ;*/

  int *wrapped  ;
  int j , k ;
  
  wrapped = ivector(natom) ;
  
  for ( j = 0 ; j < natom ; j++ ) {
    wrapped[j] = 0 ;
  }
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      if ( wrapped[mol_list[j][k]] == 0 ) {
	wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,
		   bond_list,natom,j,k, wrapped,orthogonal) ;
      }
    }
  }

  free(wrapped) ;
  
}

void wrap_pairs(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int **mol_list, int nmolecule,
		int **bond_list, int natom, int j, int k,
		int *wrapped, int orthogonal)
     /** Wrap atoms bonded to the kth atom of the jth molecule. */
{  

  /* double rx, ry, rz, dx, dy, dz ; */
  /* double rx0, ry0, rz0, mtot ; */
  /* double mass ; */
  /* double x1[2], y1[2], z1[2] ; */
  int kk, l, ll ;

  kk = mol_list[j][k] ;
  for ( l = 0 ; mol_list[j][l] >= 0 ; l++ ) {
    ll = mol_list[j][l] ;
    if ( kk >= natom || ll >= natom ) {
      printf("Error in wrap_molecules: bad index\n") ;
      exit(1) ;
    }
    if ( kk != ll && get_bond_list(bond_list,kk,ll) && wrapped[ll] == 0 ) {
      wrap_atom(kk,ll,x,y,z,a,b,c,MAXWRAP,1, orthogonal) ;
      wrapped[ll] = 1 ;
      wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,bond_list,natom,j,l,wrapped,
		 orthogonal) ;
    }
  }
}

#ifndef HAVE_NEARBYINT
double nint(double num)
/* returns the nearest int to the double num */
{
  int result ;
  
  if (  fabs( num - (int ) num ) <  0.5) 
    {
    result =  (int) num  ;
    }
  else if ( num > 0 )
    {
    result = (int) num + 1 ;
    }
  else
    result = (int) num - 1  ;

  return((double) result) ;
}
#endif


void print_molecule(int **mol_list,
		    int **bond_list,
		    int *type, char **element, int mol, 
		    int read_charge,
		    double *q,
		    double *x, double *y, double *z,
		    double *rcomx, double *rcomy,
		    double *rcomz, int *printed_atom,
		    int dump_structures) 
{
  int k, j, jj, kk  ;
  int ele_count[MAXELE], bond_count[MAXELE][MAXELE] ;
  double charge = 0.0 ;
  double dipx, dipy, dipz ;
  double dipole, q0;
  char *ele ;
  char molecule_name[BUFSIZE] ;
  char buf[BUFSIZE] ;
  int natm_mol ;
  double radius ;
  FILE *fmol ;

  memset(molecule_name, 0, BUFSIZE) ;
  memset(buf, 0, BUFSIZE) ;

  if ( dump_structures ) {
    strcat(molecule_name, "molecules/") ;
  }

  /** Calculate the overall stoichiometry **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    ele_count[k] = 0 ; 
  }
  dipx = 0.0 ;
  dipy = 0.0 ;
  dipz = 0.0 ;
  /**
  printf("COM: %f %f %f\n", rcomx[mol], rcomy[mol], rcomz[mol]) ;
  **/
  natm_mol = 0 ;
  radius = 0.0 ;
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    
    kk = mol_list[mol][k] ;

    radius += (rcomx[mol] - x[kk]) * (rcomx[mol] - x[kk])
      + (rcomy[mol] - y[kk]) * (rcomy[mol] - y[kk])
      + (rcomz[mol] - z[kk]) * (rcomz[mol] - z[kk]) ;
    /**
    printf("POS: %d %f %f %f\n", kk, x[kk], y[kk], z[kk]) ;
    **/
    q0 = - q[kk] ;
    ele = element[type[mol_list[mol][k]]-1] ;
    if ( strcmp(ele, "C") ==0 ) {
      q0 += 4 ;
    } else if ( strcmp(ele, "N") ==0 ) {
      q0 += 5 ;
    } else if ( strcmp(ele, "O") ==0) {
      q0 += 6 ;
    } else if ( strcmp(ele, "F") ==0) {
      q0 += 7 ;
    } else if ( strcmp(ele, "H") ==0) {
      q0 += 1 ;
    } else if ( strcmp(ele, "Cl") ==0) {
      q0 += 17 ;
    } else {
      printf("Error: charge for element %s is not known\n", ele) ;
      exit(1) ;
    }
    charge += q0 ;
    {
      double ddx, ddy, ddz ;
      ddx = q0 * (x[kk]-rcomx[mol]);
      ddy = q0 * (y[kk]-rcomy[mol]) ;
      ddz = q0 * (z[kk]-rcomz[mol]) ;
      /* printf("z = %e, rcomz = %e q = %e\n", z[kk], rcomz[mol], q0) ; */
      dipx += ddx ;
      dipy += ddy ;
      dipz += ddz ;
      /*      printf("The dipole contribution for atom %d = %e %e %e\n", k, ddx, ddy, ddz) ; */
    }
    ele_count[type[mol_list[mol][k]]-1]++ ;
    natm_mol ++ ;
  }
  if ( radius > 1.0e-14 ) {
    radius = sqrt(radius/natm_mol) ;
  } 
  /* printf("The number of atoms in the molecule = %d\n", natm_mol) ;  */

  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    kk = mol_list[mol][k] ;
    ele = element[type[mol_list[mol][k]]-1] ;
    if ( strcmp(ele, "C") ==0 ) {
      q0 += 4 ;
    } else if ( strcmp(ele, "N") ==0 ) {
      q0 += 5 ;
    } else if ( strcmp(ele, "O") ==0 ) {
      q0 += 6 ;
    } else if ( strcmp(ele, "F") ==0 ) {
      q0 += 7 ;
    } else if ( strcmp(ele, "H") ==0 ) {
      q0 += 1 ;
    } else if ( strcmp(ele, "Cl") ==0 ) {
      q0 += 17 ;
    } else {
      printf("Error: charge for element %s is not known\n", ele) ;
      exit(1) ;
    }
  }
  printf("Name: ") ;
  for ( k = 0 ; k < MAXELE ; k++ ) {
    if ( ele_count[k] > 0 ) {
      printf("%s%d ", element[k], ele_count[k]) ;
      sprintf(buf, "%s%d-", element[k], ele_count[k]) ;
      strncat(molecule_name, buf, BUFSIZE) ;
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
      if ( get_bond_list(bond_list,kk,jj) && jj > kk ) {
#if(0)	
	printf("Bonded: %d-%s %d-%s \n", kk, element[type[kk]-1],
	       jj, element[type[jj]-1]) ;
#endif
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
	printf("%d(%s-%s) ", bond_count[j][k],
	       element[j], element[k]) ;
	sprintf(buf,"-%d(%s-%s)", bond_count[j][k],
	       element[j], element[k]) ;
	strncat(molecule_name, buf, BUFSIZE) ;
      }
    }
  }
  printf("\n") ;


  
  for ( k = 0 ; mol_list[mol][k] >= 0 ; ) {
    printf("   Atom list: ") ;
    for ( j = 0 ; j < 15 ; j++ ) {
      kk = mol_list[mol][k] ;
      printf("%3d ", kk ) ;
      if ( printed_atom[kk] == 1 ) {
	printf("Error: attempt to print out an atom twice\n") ;
	exit(1) ;
      }
      else {
	printed_atom[kk] = 1 ;
      }
      k++ ;
      if ( mol_list[mol][k] < 0 ) {
	break ;
      }
    }
    printf("\n") ;
  }

  if ( dump_structures ) {
    strncat(molecule_name, ".xyz", BUFSIZE) ;
    fmol = fopen(molecule_name, "a") ;
    if ( fmol == NULL ) {
      printf("Error: could not open file %s\n", molecule_name) ;
      exit(1) ;
    }
    fprintf(fmol,"%d\n\n", k) ;
    for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
      j = mol_list[mol][k] ;
      ele = element[type[mol_list[mol][k]]-1] ;
      fprintf(fmol,"%s %11.4f %11.4f %11.4f\n", ele, x[j], y[j], z[j]) ;
    }
    fclose(fmol) ;
  }

  /** 4.8 is a unit conversion factor. */
  dipole = 4.8 * sqrt(dipx * dipx + dipy * dipy + dipz * dipz) ;
  if ( read_charge ) {
    printf("   Charge = %f e\n", charge) ;
    printf("   Dipole moment = %f D\n", dipole) ;
  }
  printf("RMS radius = %12.5e\n", radius ) ;

}


char *name_molecule(int **mol_list,
		    int **bond_list,
		    int *type, char **element, int mol)
/** Creates a new name for the molecule given by the index mol. */
{
  int k, j, jj, kk  ;
  int ele_count[MAXELE], bond_count[MAXELE][MAXELE] ;
  char *ele ;
  char molecule_name[BUFSIZE] ;
  char buf[BUFSIZE] ;
  int natm_mol ;
  int started ;
  
  memset(molecule_name, 0, BUFSIZE) ;
  memset(buf, 0, BUFSIZE) ;

  /** Calculate the overall stoichiometry **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    ele_count[k] = 0 ; 
  }
  natm_mol = 0 ;
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    
    kk = mol_list[mol][k] ;
    ele = element[type[mol_list[mol][k]]-1] ;
    ele_count[type[mol_list[mol][k]]-1]++ ;
    natm_mol ++ ;
  }

  started = 0 ;
  for ( k = 0 ; k < MAXELE ; k++ ) {
    if ( ele_count[k] > 0 ) {
      if ( started ) {
	sprintf(buf, "-%s%d", element[k], ele_count[k]) ;
      }
      else {
	sprintf(buf, "%s%d", element[k], ele_count[k]) ;
      }
      started = 1 ;
      strncat(molecule_name, buf, BUFSIZE) ;
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
      if ( get_bond_list(bond_list,kk,jj) && jj > kk ) {
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
	sprintf(buf,"-%d(%s-%s)", bond_count[j][k],
	       element[j], element[k]) ;
	strncat(molecule_name, buf, BUFSIZE) ;
      }
    }
  }
  
  return( strdup(molecule_name) ) ;
}
    

static void wrap_in_box(double a[3], double b[3], double c[3],
			double invbox[3][3],
			int natom, double *x, double *y, 
			double *z) 
     /** This function wraps all the atoms into the primitive simulation
	 box. **/
{
  int j ;
  double ca, cb, cc ;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(j,ca,cb,cc) \
  shared(invbox,x,y,z,a,b,c,natom)
#endif
  for ( j = 0 ; j < natom ; j++ ) {
    ca = invbox[0][0] * x[j] + invbox[1][0] * y[j] + invbox[2][0] * z[j] ;
    cb = invbox[0][1] * x[j] + invbox[1][1] * y[j] + invbox[2][1] * z[j] ;
    cc = invbox[0][2] * x[j] + invbox[1][2] * y[j] + invbox[2][2] * z[j] ;

    ca -= nint(ca) ;
    cb -= nint(cb) ;
    cc -= nint(cc) ;

    x[j] = ca * a[0] + cb * b[0] + cc * c[0] ;
    y[j] = ca * a[1] + cb * b[1] + cc * c[1] ;
    z[j] = ca * a[2] + cb * b[2] + cc * c[2] ;

  }
}

void inversebox(double a[3], double b[3], double c[3],
  double invbox[3][3]) 
     /** Calculate the inverse box maxtrix, used for finding
	 crystal coordinates. */
{
  double xhlp ;
  xhlp=-a[2]*b[1]*c[0]+ a[1]*b[2]*c[0] ;
  xhlp=xhlp+a[2]*b[0]*c[1] ;
  xhlp=xhlp-a[0]*b[2]*c[1] ;
  xhlp=xhlp-a[1]*b[0]*c[2] ;
  xhlp=xhlp+a[0]*b[1]*c[2] ;

  invbox[0][0]=(-b[2]*c[1]+b[1]*c[2])/xhlp ;
  invbox[1][0]=(b[2]*c[0]-b[0]*c[2])/xhlp ;
  invbox[2][0]=(-b[1]*c[0]+b[0]*c[1])/xhlp ;

  invbox[0][1]=(a[2]*c[1]-a[1]*c[2])/xhlp ;
  invbox[1][1]=(-a[2]*c[0]+a[0]*c[2])/xhlp ;
  invbox[2][1]=(a[1]*c[0]-a[0]*c[1])/xhlp ;

  invbox[0][2]=(-a[2]*b[1]+a[1]*b[2])/xhlp ;
  invbox[1][2]=(a[2]*b[0]-a[0]*b[2])/xhlp ;
  invbox[2][2]=(-a[1]*b[0]+a[0]*b[1])/xhlp ;
}


double atm_mass(char *ea)
     /* Return the atomic mass of the given element. */
{
  if ( strcmp(ea, "H") ==0 ) {
    return(1.0) ;
  } else if ( strcmp(ea, "C") ==0 ) {
    return(12.0) ;
  } else if ( strcmp(ea, "N") ==0 ) {
    return(14.0) ;
  } else if ( strcmp(ea, "O") ==0 ) {
    return(16.0) ;
  } else if ( strcmp(ea, "F") ==0 ) {
    return(19.0) ;
  } else if ( strcmp(ea, "Cl") ==0 ) {
    return(35.5) ;
  } else {
    printf("Error: element %s is unknown\n", ea ) ;
    exit(1) ;
  }
  return(1.0) ;
}

void read_bonds(double rbond[MAXELE][MAXELE], int *xyz_copies, double *time_step,
		int *dump_structures, int bond_time[MAXELE][MAXELE])
     /** Read in bond distances and other parameters.
	 rbond - The maximum bond distance.
	 xyz_copies - The number of copies in each direction for the xyz output file.
	 time_step - The time between each frame stored.
	 dump_structures - Whether to write out individual molecule structure files.
     **/
{
  FILE *fbond, *ftime = NULL ;
  char *ele1, *ele2 ;
  char buf1[BUFSIZ], buf2[BUFSIZ] ;
  double r ;
  int i, j, k, btime ;
  int duration ;
  
  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      /** Put in a negative value so that unspecified bonds can be caught. */
      rbond[i][j] = -1.0 ;
      bond_time[i][j] = -1 ;
    }
  }

  fbond = fopen("bonds.dat","r") ;
  if ( fbond == NULL ) {
    printf("Error: could not open bonds.dat\n") ;
    exit(1) ;
  }
  fscanf(fbond, "%d", &duration) ;
  if ( duration >= 0) {
    printf("ALL bonds require a lifetime of %d steps\n", duration)  ;
  } else {
    ftime = fopen("bond_times.dat","r") ;
    if (ftime == NULL) {
      printf ("Error: could not open bond_times.dat\n");
      exit(1);
    }
  }
  fscanf(fbond, "%d", xyz_copies) ;
  printf("%d replicas of the simulation cell will be propagated in each direction\n",
	 *xyz_copies) ;

  fscanf(fbond, "%lf\n", time_step) ;
  if ( *time_step < 0.0 || *time_step > 1.0e-10 ) {
    printf("The time step was given as %21.14e seconds.  This seems strange\n",
	   *time_step) ;
    exit(1) ;
  }
  fscanf(fbond, "%d", dump_structures) ;
  if ( *dump_structures != 0 ) {
    printf("Individual molecule xyz files will be created in the molecules directory.\n") ;
  } else {
    printf("Individual molecule xyz file will not be created.\n") ;
  }
  printf("Bond distances:\n") ;

  do {
    if ( fscanf(fbond, "%s %s %lf", buf1, buf2, &r) != 3 ) {
      break ;
    }
    ele1 = buf1 ;
    ele2 = buf2 ;
    i = ele_index(ele1) ;
    j = ele_index(ele2) ;
    if ( j < i ) {
      k = i ;
      i = j ;
      j = k ;
    }
    printf("%s %s %f\n", ele1, ele2, r) ;
    if ( rbond[i][j] > 0.0 || rbond[j][i] > 0.0 ) {
      printf("Error: Bond distance for %s %s has already been specified\n",
	     ele1, ele2) ;
      exit(1) ;
    }
    rbond[i][j] = r ;
    if ( duration >=0 ) {
      bond_time[i][j] = duration ;
    }
    
  } while ( ! feof(fbond) && ! ferror(fbond)) ;
  fclose(fbond) ;


  printf("Bond time cutoffs:\n") ;
  if ( duration >= 0) {
    for (i = 0; i < MAXELE; i++) {
      for (j = 0; j < MAXELE; j++) {
        if (i == 0) {
          ele1 = "H";
        } else if ( i == 1) {
          ele1 = "C";
        } else if ( i == 2) {
          ele1 = "N";
        } else if ( i == 3) {
          ele1 = "O";
        } else if ( i == 4) {
          ele1 = "F";
        } else if ( i == 5) {
          ele1 = "Cl";
        }
        if (j == 0) {
          ele2 = "H";
        } else if ( j == 1) {
          ele2 = "C";
        } else if ( j == 2) {
          ele2 = "N";
        } else if ( j == 3) {
          ele2 = "O";
        } else if ( j == 4) {
          ele2 = "F";
        } else if ( j == 5) {
          ele2 = "Cl";
        }
        if (bond_time[i][j] != -1) {
          printf("%s %s %d\n", ele1, ele2, bond_time[i][j]) ;
        }
      }
    }
  } else {
    do {
      if ( fscanf(ftime, "%s %s %d", buf1, buf2, &btime) != 3 ) {
        break ;
      }
      ele1 = buf1 ;
      ele2 = buf2 ;
      i = ele_index(ele1) ;
      j = ele_index(ele2) ;
      if ( j < i ) {
        k = i ;
        i = j ;
        j = k ;
      }
      printf("%s %s %d\n", ele1, ele2, btime) ;
      if ( bond_time[i][j] > 0.0 || bond_time[j][i] > 0.0 ) {
        printf("Error: Bond time cutoff for %s %s has already been specified\n",
	       ele1, ele2) ;
        exit(1) ;
      }
      bond_time[i][j] = btime ;
    
    } while ( ! feof(ftime) && ! ferror(ftime)) ;
    fclose(ftime) ;
  }

  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = i + 1 ; j < MAXELE ; j++ ) {
      rbond[j][i] = rbond[i][j] ;
      bond_time[j][i] = bond_time[i][j] ;
    }
  }
}

int ele_index(char *ea)
     /** Returns the index of the given element */
{
  int i ;

  if ( ea == NULL ) {
    printf("Error: null element name\n") ;
    exit(1) ;
  }

  if ( strcmp(ea, "H") == 0 ) {
    i = 0 ;
  } else if ( strcmp(ea, "C") == 0 ) {
    i = 1 ;
  } else if (strcmp(ea, "N") == 0 ) {
    i = 2 ;
  } else if (strcmp(ea, "O") == 0 ) {
    i = 3 ; 
  } else if (strcmp(ea, "F") == 0 ) {
    i = 4 ; 
  } else if (strcmp(ea, "Cl") == 0 ) {
    i = 5 ; 
  } else {
    printf("Error: element %s is unknown\n", ea ) ;
    exit(1) ;
  }
  return i ;
}

double boxvol(double a[3], double b[3], double c[3])
     /* Calculates the volume of the simulation box, given the
	box cell vectors a, b, and c. */
{

  double bv ;
  double xhlp ;
  double ax, ay, az, bx, by, bz, cx, cy, cz ;

  ax = a[0] ;
  ay = a[1] ;
  az = a[2] ;
  bx = b[0] ;
  by = b[1] ;
  bz = b[2] ;
  cx = c[0] ;
  cy = c[1] ;
  cz = c[2] ;

  xhlp=-az*by*cx+ ay*bz*cx ;
  xhlp=xhlp+az*bx*cy ;
  xhlp=xhlp-ax*bz*cy ;
  xhlp=xhlp-ay*bx*cz ;
  xhlp=xhlp+ax*by*cz ;

  bv = fabs(xhlp) ;
  return bv;
}

int get_next_neighbor(int j, int *k, int *nid, const int **neighbor)
/** Determine the next neighbor of atom j and store in *k. 
    
    Returns 0 when all neighbors have been found.
**/
{
  
  if ( *k < 0 ) {
    *nid = 0 ;
  }
  else {
    (*nid)++ ;
  }
  if ( neighbor[j][*nid] >= 0 ) {
    *k = neighbor[j][*nid] ;
    return(1) ;
  }
  else {
    return(0) ;
  }
}

  

void free_neighbor(int **neighbor, int natom) 
/* Frees the neighbor list */
{
  int j ;
  
  for ( j = 0 ; j < natom ; j++ ) {
    if ( neighbor[j] != NULL ) {
      free(neighbor[j]) ;
    }
  }
}


void find_neighbors(double *x,double *y, double *z, int *type,
		    char **element,
		    double a[3], double b[3], double c[3],
		    int **neighbor,int natom, int orthogonal,
		    int **bond_duration) 
/** Find the neighbor list for a given set of atoms */
{
  double minimum_distance = MIN_BOND_LENGTH ;

  if ( orthogonal && natom > 100 ) {
    /* Use the order N algorithm */
    find_neighbors_fast(x,y,z,type,a,b,c,neighbor,natom,orthogonal,bond_duration,
			minimum_distance) ;
  }
  else {
    /* Use the order N^2 algorithm */
    find_neighbors_slow(x,y,z,type,element, a,b,c,neighbor,natom,orthogonal,bond_duration,
			minimum_distance) ;
  }

  /* Sort neighbor list so results don't depend on algorithm */
  sort_neighbor(neighbor, natom) ;

}

     
void find_neighbors_slow(double *x,double *y, double *z, int *type,
			 char **element,
			 double a[3], double b[3], double c[3],
			 int **neighbor,int natom, int orthogonal,
			 int **bond_duration,double minimum_distance) 
/** Find the neighbor list for a given set of atoms */
/** Order N^2 algorithm */
{
  int init_maxbond = 10 ;
  int *maxbond ;
  double r2 = 0.0 ;
  int j, k, l ;
  double r2ele[MAXELE][MAXELE] ;

  maxbond = ivector(natom) ;

  for ( j = 0 ; j < MAXELE ; j++ ) {
    int ij ; 
    if ( element[j] == NULL ) {
      continue ;
    }
    ij = ele_index(element[j]) ;
    for ( k = 0 ; k < MAXELE ; k++ ) {
      int ik ;
      if ( element[k] == NULL ) {
	continue ;
      }
      ik = ele_index(element[k]) ;
      r2ele[j][k] = rbond[ij][ik] * rbond[ij][ik] ;
    }
  }

#ifdef OPENMP
#pragma omp parallel default(none)                       \
  private(j, k, l, r2)					 \
  shared(neighbor, maxbond, a, b, c, x, y, z, orthogonal,\
         r2ele,natom,type,init_maxbond,bond_duration,    \
         minimum_distance)
#endif
  {
    
#ifdef OPENMP
#pragma omp for
#endif
    for ( j = 0 ; j < natom ; j++ ) {
      maxbond[j] = init_maxbond ;
      neighbor[j] = malloc(maxbond[j] * sizeof(int)) ;
    
      l = 0 ;
      for ( k = 0 ; k < natom ; k++ ) {
	if ( k == j ) continue ;
      
	r2 = wrap_atom(j,k,x,y,z,a,b,c,MAXWRAP,0, orthogonal) ;
#ifdef SPEW
        printf("Testing Atoms %d and %d for bonding: r = %f\n", 
               j, k, sqrt(r2) ) ;  
#endif

	if ( r2 < minimum_distance ) {
	  printf("Error: atoms %d and %d were closer than %f\n",
		 j, k, minimum_distance) ;
	} else if ( r2 < r2ele[type[j]-1][type[k]-1] ) {
	  neighbor[j][l] = k ;
	  add_bond(j, k, bond_duration) ;
#ifdef SPEW
	  printf("Atoms %d and %d are instantaneously bonded\n",j,k) ;
#endif

	  l++ ;
	  if ( l == maxbond[j] ) {
	    maxbond[j] *= 2 ;
	    printf("Warning:  maximum number of bonded atoms increased to %d\n",
		   maxbond[j]) ;
	    neighbor[j] = realloc(neighbor[j],maxbond[j] * sizeof(int*)) ;
	  }
	} else {
	  remove_bond(j, k, bond_duration) ;
	}
      }
      neighbor[j][l] = -1 ;
    }
  }
  free(maxbond) ;
}

  
void find_neighbors_fast(double *x,double *y, double *z, int *type,
			 double a[3], double b[3], double c[3],
			 int **neighbor,int natom, int orthogonal,
			 int **bond_duration, double minimum_distance) 
/** Order N algorithm to find the neighbor list for a given set of atoms */
{
  int init_maxbond = 10 ;  /* Maximum number of bonds on an atom */
  int *maxbond ;
  double r2 = 0.0 ;
  int j, k, l ;

  int nx, ny, nz ;
  double bond_max = -1.0 ; /* Maximum bond distance */

  int **box_list ;
  int nbox ;
  int *box_end ;
  int *box_idx ;   /* The index of the box corresponding to a particular atom */
  
  int box_alloc ;  /* Allocated number of atoms in a box */
  int ix, iy, iz ;
  int idx ;
  int ixf, iyf, izf ;
  int ixx, iyy, izz ;
  int idx2 ;
  int m ;
  double r2val[MAXELE][MAXELE] ;
  
  if ( ! orthogonal ) {
    printf("Error: fast neighbor algorithm only works for an orthorhombic box\n") ;
    exit(1) ;
  }

  maxbond = ivector(natom) ;

  for ( k = 0 ; k < MAXELE ; k++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      if ( rbond[j][k] > bond_max ) {
	bond_max = rbond[j][k] ;
      }
      r2val[j][k] = rbond[j][k] * rbond[j][k] ;
    }
  }
  
  nx = floor(a[0] / bond_max) ;
  ny = floor(b[1] / bond_max) ;
  nz = floor(c[2] / bond_max) ;
  
  if ( nx <= 0 || ny <= 0 || nz <= 0 ) {
    printf("Error: box size smaller than bond size\n") ;
    exit(1) ;
  }
  nbox = nx * ny * nz ;
  
  box_list = malloc(nbox * sizeof(int*) ) ;
  box_end = malloc(nbox * sizeof(int) ) ;
  box_idx = malloc(natom * sizeof(int)) ;
  
  box_alloc = 3 * natom / nbox + 1 ;
  
  for ( j = 0 ; j < nbox ; j++ ) {
    box_list[j] = malloc( box_alloc * sizeof(int) ) ;
    box_end[j] = 0 ;
    for ( k = 0 ; k < box_alloc ; k++ ) {
      box_list[j][k] = -1 ;
    }
  }
  
  for ( j = 0 ; j < natom ; j++ ) {
    ix = x[j] / a[0] ;
    iy = y[j] / b[1] ;
    iz = z[j] / c[2] ;
    
    if ( ix < 0 || iy < 0 || iz < 0 ) {
      printf("Error: negative box index\n") ;
      exit(1) ;
    }
    if ( ix >= nx || iy >= ny || iz >= nz ) {
      printf("Error: box index was too large\n") ;
      exit(1) ;
    }

    idx = ix + iy * nx + iz * ny * nx ;
    if ( idx >= nbox ) {
      printf("Error: overran box index\n") ;
      exit(1) ;
    }    
    
    /* Place the atom index in the appropriate box. */
    box_list[idx][box_end[idx]] = j ;
    box_end[idx]++ ;
    if ( box_end[idx] >= box_alloc ) {
      box_alloc *= 2 ;
      for ( k = 0 ; k < nbox ; k++ ) {
	box_list[k] = realloc(box_list[k], box_alloc * sizeof(int) ) ;
      }
    }
    box_idx[j] = idx ;
  }
  
#ifdef OPENMP
#pragma omp parallel default(none)                               \
  private(j, k, l, idx, iz, iy, ix, izz, izf, iyy, iyf, ixx, ixf,\
          idx2, m, r2)                                           \
  shared(neighbor, maxbond, box_idx, box_list, natom, nx, ny, nz,\
         a, b, c, x, y, z, orthogonal, type, box_end,            \
	 rbond2,bond_index,init_maxbond,bond_duration,           \
         minimum_distance)
#endif
  {

#ifdef OPENMP    
#pragma omp for
#endif
    for ( j = 0 ; j < natom ; j++ ) {
      maxbond[j] = init_maxbond ;
      neighbor[j] = malloc(maxbond[j] * sizeof(int)) ;

      /** Search for neighbors across boxes */
      idx = box_idx[j] ;
      l = 0 ;
    
      iz = idx / (ny * nx) ;
      iy = (idx - iz * ny * nx) / ny ;
      ix = idx - iz * ny * nx - iy * nx ;
    
      /* Loop over nearest neighbor boxes, applying periodic boundaries where
	 appropriate */
      for ( izz = iz - 1 ; izz <= iz + 1 ; izz++ ) {
	if ( izz < 0 ) {
	  izf = izz + nz ;
	}
	else {
	  izf = izz ;
	}
	for ( iyy = iy - 1 ; iyy <= iy + 1 ; iyy++ ) {
	  if ( iyy < 0 ) {
	    iyf = iyy + ny ;
	  }
	  else {
	    iyf = iyy ;
	  }
	  for ( ixx = ix - 1 ; ixx <= ix + 1 ; ixx++ ) {
	    if ( ixx < 0 ) {
	      ixf = ixx + nx ;
	    }
	    else {
	      ixf = ixx ;
	    }

	    /* ixf, iyf, izf are indices "fixed" for periodic boundaries */

	    idx2 = ixf + iyf * nx + izf * ny * nx ;

	    for ( m = 0 ; m < box_end[idx2] ; m++ ) {
	      /* Search over all atoms in the given box. */
	      k = box_list[idx2][m] ;
	      if ( k < 0 ) {
		printf("Error: negative atom index found\n") ;
		exit(1) ;
	      } else if ( k >= natom ) {
		printf("Error: atom index was bigger than natom\n") ;
		exit(1) ;
	      }
	      if ( k == j ) continue ;
	      r2 = wrap_atom(j,k,x,y,z,a,b,c,MAXWRAP,0, orthogonal) ;
	      if ( r2 < minimum_distance ) {
		printf("Error: atoms %d and %d were closer than %f\n",
		       j, k, minimum_distance) ;
	      } else if ( r2 < rbond2[bond_index[type[j]-1]][bond_index[type[k]-1]] ) {
		/* if ( r2 < r2val[type[j]-1][type[k]-1] ) { */
		neighbor[j][l] = k ;
		add_bond(j, k, bond_duration) ;
		l++ ;
		if ( l == maxbond[j] ) {
		  maxbond[j] *= 2 ;
		  printf("Warning:  maximum number of bonded atoms increased to %d\n",
			 maxbond[j]) ;
		  neighbor[j] = realloc(neighbor[j],maxbond[j] * sizeof(int) ) ;
		}
	      } else {
		remove_bond(j, k, bond_duration) ;
	      }
	    }
	  }
	}
      }
      neighbor[j][l] = -1 ;
    }
  } /* End of parallel region */
  
  free(maxbond) ;
  free(box_end) ;
  free(box_idx) ; 
  for ( j = 0 ; j < nbox ; j++ ) {
    free(box_list[j]) ;
  }
  free(box_list) ;
}


void print_com(FILE *fcom, double *rcomx, double *rcomy, double *rcomz, char **mol_name, 
	       int nmolecule, int frame) 
/* Print out the center of mass coordinates of each molecule */
{
  int j ;
  
  fprintf(fcom,"COM Coordinates for frame %d:\n", frame) ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    fprintf(fcom, "%8.4e %8.4e %8.4e %s\n",
	    rcomx[j], rcomy[j], rcomz[j], mol_name[j]) ;
  }
}

int check_orthogonal_box(double a[3],double b[3], double c[3]) 
/* Check for box orthogonality.  */
{
  double eps = 1.0e-12 ;
  /* double ax[3], by[3], cz[3] ; */
  
  if ( fabs(a[1]) > eps || fabs(a[2]) > eps ) {
    return(0) ;
  }
  if ( fabs(b[0]) > eps || fabs(b[2]) > eps ) {
    return(0) ;
  }
  if ( fabs(c[0]) > eps || fabs(c[1]) > eps ) {
    return(0) ;
  }
  return(1) ;
}
  
  
  
void sort_neighbor(int **neighbor, int natom)
/* Sort the neighbor list in ascending order */
{
  int j, l ;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(j,l) \
  shared(natom,neighbor)
#endif
  for ( j = 0 ; j < natom ; j++ ) {
    /* Count up number of neighbors for qsort */
    for ( l = 0 ; neighbor[j][l] >= 0 ; l++ ) {
      ;
    }
    /* Sort neighbors */
    qsort(neighbor[j], l, sizeof(int), comp_int) ;
  }
}


int comp_int(const void *pi, const void *pj)
/* Compare two integers.  Used in the qsort routine */
{
  int *pii = (int*) pi ;
  int *pji = (int*) pj ;
  
  int i = *pii ;
  int j = *pji ;
  
  if ( i > j ) {
    return(1) ;
  }
  else if ( i < j ) {
    return(-1) ;
  }
  else {
    return(0) ;
  }
}

double **dmatrix(int m, int n)
/* Allocates a double precision matrix of dimension m by n */
{
  double **result ;
  double *buf ;
  int j ;
  
  buf = calloc(m*n, sizeof(double)) ;
  result = calloc(m,sizeof(double*)) ;
  
  for ( j = 0 ; j < m ; j++ ) {
    result[j] = buf + j*n ;
  }
  return(result) ;
}


double *dvector(int m)
/** Allocates a double precision vector of dimension m */
{
  double *result ;
  
  result = calloc(m, sizeof(double)) ;
  
  return(result) ;
}

int *ivector(int m)
/** Allocates an int vector of dimension m */
{
  int *result ;
  
  result = calloc(m, sizeof(int)) ;
  
  return(result) ;
}

int **imatrix(int m, int n)
/* Allocates an int matrix of dimension m by n */
{
  int **result ;
  int *buf ;
  int j ;
  
  buf = calloc(m*n, sizeof(int)) ;
  result = calloc(m,sizeof(int*)) ;
  
  for ( j = 0 ; j < m ; j++ ) {
    result[j] = buf + j*n ;
  }
  return(result) ;
}


void allocate_arrays(double **x, double **y, double **z, double **q, 
		     double **rcomx, double **rcomy, double **rcomz, 
		     int **type0, int **type, int **printed_atom, 
		     int ***mol_list, int ***bond_list, int ***bond_duration,
		     int ***neighbor, int natom)
/* Allocate arrays dynamically.  Initialize q */
{
  int i, j ;

  *x = dvector(natom) ;
  *y = dvector(natom) ;
  *z = dvector(natom) ;
  *q = dvector(natom) ;
  *rcomx = dvector(natom) ;
  *rcomy = dvector(natom) ;
  *rcomz = dvector(natom) ;

  *type = ivector(natom) ;
  *type0 = ivector(natom) ;
  *printed_atom = ivector(natom) ;

  /* Lazy allocation of molecular list to save memory */
  *mol_list = calloc(natom, sizeof(int*)) ;
    
  if ( sparse_bonds ) {
    /* Use sparse storage of bonds to avoid order N**2 memory requirement */
    *bond_list = calloc(natom, sizeof(int*)) ;
    for ( j = 0 ; j < natom ; j++ ) {
      (*bond_list)[j] = calloc(2*MAXCOORD + BOND_IDX, sizeof(int)) ;
      (*bond_list)[j][DIM_IDX] = MAXCOORD ;
      (*bond_list)[j][BOND_IDX] = 0 ;
    }

    *bond_duration = calloc(natom, sizeof(int*)) ;
    for ( j = 0 ; j < natom ; j++ ) {
      (*bond_duration)[j] = calloc(2*MAXCOORD + BOND_IDX, sizeof(int)) ;
      (*bond_duration)[j][DIM_IDX] = MAXCOORD ;
      (*bond_duration)[j][BOND_IDX] = 0 ;
    }

  }
  else {
    *bond_list = imatrix(natom,natom) ;
    *bond_duration = imatrix(natom,natom) ;
  }

  *neighbor = calloc(natom, sizeof(int*)) ;

  for ( i = 0 ; i < natom ; i++ ) {
    (*q)[i] = 0.0 ;
  }


}

void initialize_bonds(int **mol_list, int **bond_list, int natom, int *nmolecule) 
/** Initialize the bond and molecule arrays */
{
  int i, j ;
  int ndim ;

  for ( j = 0 ; j < *nmolecule ; j++ ) {
    free(mol_list[j]) ;
  }
  *nmolecule = 0 ;

  /*
  for ( j = 0 ; j < natom ; j++ ) {
    for ( i = 0 ; i < natom ; i++ ) {
      mol_list[i][j] = - 1 ;
    }
  }
  */

  if ( sparse_bonds ) {
    for ( j = 0 ; j < natom ; j++ ) {
      bond_list[j][NBOND_IDX] = 0 ;
      ndim = bond_list[j][DIM_IDX] ;
      
      for ( i = 0 ; i < ndim ; i++ ) {
	/* Defensive programming */
	bond_list[j][BOND_IDX+i] = -1 ;
      }
    }
  }
  else {
      for ( j = 0 ; j < natom ; j++ ) {
	for ( i = 0 ; i < natom ; i++ ) {
          bond_list[i][j] = 0 ;
	}
      }
  }
}

void add_bond(int j, int k, int **bond_list)
/* Add a bond between atom j and atom k */
/* Storage scheme:
   bond_list[j][NBOND_IDX] : number of bonds
   bond_list[j][DIM_IDX] :   number of allocated bonds.
   bond_list[j][BOND_IDX+2*k] : Index of the kth atom bonded to j
   bond_list[j][BOND_IDX+2*k+1] : Duration of the bond.
*/
{
  if ( sparse_bonds ) {

    int nbond ;  /* Number of bonds found so far for atom j */
    int ndim ;   /* Dimensioned maximum number of bonds */
    int l ;

    nbond = bond_list[j][NBOND_IDX] ;
    ndim  = bond_list[j][DIM_IDX] ;

    /* Check for an old bond */
    for ( l = 0 ; l < nbond ; l++ ) {
      if ( bond_list[j][2*l+BOND_IDX] == k ) {
	
	/* Increase bond duration by one */
	bond_list[j][2*l+BOND_IDX+1]++ ;
#ifdef SPEW
    printf("Add_bond: increased bond duration between %d and %d to %d\n", 
	   j, k, bond_list[j][2*l+BOND_IDX+1]) ;
#endif
	return ;
      }
    }

    /* Add a new bond */
    if ( nbond >= ndim ) {
      bond_list[j] = realloc(bond_list[j], 
			     (4 * ndim + BOND_IDX) * sizeof(int*) ) ;
      bond_list[j][DIM_IDX] = 2 * ndim ;
    }
    bond_list[j][2*nbond+BOND_IDX] = k ;
    bond_list[j][2*nbond+BOND_IDX+1] = 1 ;
    bond_list[j][NBOND_IDX]++ ;
#ifdef SPEW
    printf("Add_bond: adding a new bond between %d and %d\n", j, k) ;
#endif
  }
  else {
    bond_list[j][k]++ ;
  }
}


void reset_bonds(int **bond_list, int natom)
{
  int j, l ;
  int nbond ;
  int ndim ;

  for ( j = 0 ; j < natom ; j++ ) {
    nbond = bond_list[j][NBOND_IDX] ;
    ndim  = bond_list[j][DIM_IDX] ;

    if ( nbond > ndim ) {
      printf("Error: inconsistent nbond found\n") ;
      exit(1) ;
    }
    for ( l = 0 ; l < nbond ; l++ ) {
      bond_list[j][2*l+BOND_IDX] = -1 ;
      bond_list[j][2*l+BOND_IDX+1] = -1 ;
    }
    bond_list[j][NBOND_IDX] = 0 ;
  }
}

void remove_bond(int j, int k, int **bond_list)
/* Remove a bond between atom j and atom k */
{

  if ( sparse_bonds ) {

    int nbond ;  /* Number of bonds found so far for atom j */
    int ndim ;   /* Dimensioned maximum number of bonds */
    int l, m ;

    nbond = bond_list[j][NBOND_IDX] ;
    ndim  = bond_list[j][DIM_IDX] ;

    if ( nbond > ndim ) {
      printf("Error: inconsistent nbond found\n") ;
      exit(1) ;
    }
    /* Check for an old bond */
    for ( l = 0 ; l < nbond ; l++ ) {
      if ( bond_list[j][2*l+BOND_IDX] == k ) {
	/* The bond was found */
	break ;
      }
    }
    if ( l < nbond ) {
#ifdef SPEW
      printf("Remove_bond:  Found bond to remove between %d and %d\n",
	     j, k) ;
#endif
      for ( m = l + 1 ; m < nbond ; m++ ) {
	/* Copy index of bonded atom */
	bond_list[j][2*(m-1)+BOND_IDX] = bond_list[j][2*m+BOND_IDX] ;

	/* Copy duration of bond */
	bond_list[j][2*(m-1)+BOND_IDX+1] = bond_list[j][2*m+BOND_IDX+1] ;
      }
      bond_list[j][2*(nbond-1)+BOND_IDX] = -1 ;
      bond_list[j][2*(nbond-1)+BOND_IDX+1] = -1 ;
      /* Reduce total number of bonds to atom j */
      nbond-- ;
      bond_list[j][NBOND_IDX] = nbond ;
    } else {
#ifdef SPEW
      printf("Remove_bond:  Did not find bond to remove between %d and %d\n",
	     j, k) ;
#endif
    }
  } else {
    bond_list[j][k] = 0 ;
  }
}


int get_bond_list(int **bond_list, int j, int k)
/* Return 1 if atom j and k are bonded, 0 otherwise */
{

  if ( sparse_bonds ) {
    int l ;
    int nbond ;  /* Number of bonds found so far for atom j */

    nbond = bond_list[j][NBOND_IDX] ;
    for ( l = 0 ; l < nbond ; l++ ) {
      if ( bond_list[j][2*l+BOND_IDX] == k ) {
	return(1) ;
      } else if ( bond_list[j][2*l+BOND_IDX] < 0 ) {
	printf("Error in get_bond_list: bad array value\n") ;
	exit(1) ;
      }
    }
    return(0);
  }
  else {
    if ( bond_list[j][k] > 0 ) {
      return(1) ;
    } else {
      return(0) ;
    }
  }
}

int get_bond_duration(int **bond_list, int j, int k)
/* Return the number of frames that atom j and k have been continuously bonded, 0 otherwise */
{

  if ( sparse_bonds ) {
    int l ;
    int nbond ;  /* Number of bonds found so far for atom j */

    nbond = bond_list[j][NBOND_IDX] ;
    for ( l = 0 ; l < nbond ; l++ ) {
      if ( bond_list[j][2*l+BOND_IDX] == k ) {
	if ( bond_list[j][2*l+BOND_IDX+1] <= 0 ) {
	  printf("Error in get_bond_duration: bad array value\n") ;
	}
#ifdef SPEW
	printf("Duration of bond between %d and %d is %d\n",
	       j, k, bond_list[j][2*l+BOND_IDX+1]) ;
#endif
	return(bond_list[j][2*l+BOND_IDX+1]) ;
      } else if ( bond_list[j][2*l+BOND_IDX] < 0 ) {
	printf("Error in get_bond_duration: bad array value\n") ;
	exit(1) ;
      }
    }
    return(0);
  }
  else {
    return( bond_list[j][k] ) ;
  }
}



void free_arrays(double *x, double *y, double *z, double *q, 
		 double *rcomx, double *rcomy, double *rcomz, 
		 int *type, int *type0, int *printed_atom, 
		 int **mol_list, int **bond_list, 
		 int **neighbor, int natom, int nmolecule) 
/* Free arrays allocated at the beginning of program execution */
{
  int j ;
  
  free(x) ;
  free(y) ;
  free(z) ;
  free(q) ;
  free(rcomx) ;
  free(rcomy) ;
  free(rcomz) ;

  free(type) ;
  free(type0) ;
  free(printed_atom) ;

  /* Lazy allocation of molecular list to save memory */
  
  for ( j = 0 ; j < nmolecule ; j++ ) {
    free(mol_list[j]) ;
  }
  free(mol_list) ;
    
  if ( sparse_bonds ) {
    /* Use sparse storage of bonds to avoid order N**2 memory requirement */
    for ( j = 0 ; j < natom ; j++ ) {
      free(bond_list[j]) ;
    }
    free(bond_list) ;

  }
  else {
    free_imatrix(bond_list) ;
  }
  free(neighbor) ;
}

void free_dmatrix(double **mat)
{
  free(mat[0]) ;
  free(mat) ;
}

void free_imatrix(int **mat)
{
  free(mat[0]) ;
  free(mat) ;
}

void com_molecule(int **mol_list, int *type, const double *x, 
		  const double *y, const double *z,
		  char *element[MAXELE],
		  int j, double *rxx, double *ryy, double *rzz)
/* Calculate the center of mass of the jth molecule. */
{
  double rx, ry, rz ;
  double mtot ;
  double mass ;
  int k, kk ;
  
  
  rx = 0.0 ;
  ry = 0.0 ;
  rz = 0.0 ;
  mtot = 0.0 ;
  for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
    kk = mol_list[j][k] ;
    mass = atm_mass(element[type[kk]-1]) ;
    mtot += mass ;
    rx += x[kk] * mass ;
    ry += y[kk] * mass ;
    rz += z[kk] * mass ;
    /**
       printf("XX = %f YY = %f ZZ = %f Mass = %f type = %d\n", x[kk],
       y[kk], z[kk], mass, type[mol_list[j][k]] ) ;
       printf("rx = %f rx = %f rz = %f\n", rx, ry, rz) ;
    **/
  }
  /* rx, ry, rz is the molecular center of mass. */
  rx /= mtot ;
  ry /= mtot ;
  rz /= mtot ;

  *rxx = rx ;
  *ryy = ry ;
  *rzz = rz ;
}

int read_charge_file(FILE *fchr, double *q)
/* Read charges if provided.  
   Returns number of atoms in the charge file.
 */
{
  int i, j ;
  int natom1 ;
  char buf[BUFSIZE] ;

  for ( i = 0 ; i < 4 ; i++ ) {
    fgets(buf,BUFSIZ,fchr) ;
    if ( feof(fchr) ) {
      printf("End of charge file reached\n") ;
      exit(0) ;
    }
    if ( buf[0] != '#' ) {
      printf("Bad format in charge file.\n") ;
      exit(1) ;
    }
  }
  fscanf(fchr,"#   %d\n", &natom1) ;
  for ( i = 0 ; i < natom1 ; i++ ) {
    fgets(buf,BUFSIZ,fchr) ;
    if ( sscanf(buf,"%d  %lf", &j, &q[i]) != 2 ) {
      printf("Error: could not read in a charge\n") ;
      exit(1) ;
    } 
  }
  if ( feof(fchr) || ferror(fchr) ) {
    printf("Error: end of CHR file found\n") ;
    exit(1) ;
  }
  return(natom1) ;
}


int read_gen_file(FILE *fgen, char *element[MAXELE], 
		  char tmp_ele[MAXELE][MAXELE],
		  int *maxele,
		  int *type, int *maxtype,
		  double *x, double *y, double *z,
		  double *a,  double *b,  double *c)
/* Read the next frame in from the DFTB-style GEN file 
   Returns the number of atoms.

*/
{
  char buf[BUFSIZE] ;
  int j, k ;
  int natom = 0 ;

  while ( fgets(buf,BUFSIZ,fgen) ) {
    if ( buf[0] != '#' ) {
      break ;
    }
  }
  if ( feof(fgen) || ferror(fgen) ) {
    return (0) ;
  }
  sscanf(buf,"%d", &natom) ;

  fgets(buf,BUFSIZ,fgen) ; 
  for ( j = 0 ; buf[j] != 0 ; j++ ) {
    if ( buf[j] == '\'' ) {
      buf[j] = ' ' ;
    }
  }
  /** Parse the input line, get the element types **/
  k = 0 ;
  for ( j = 0 ; buf[j] ; j++ ) {
    if ( isalpha((int)buf[j]) ) {
      if ( isalpha((int)buf[j+1]) ) {
	tmp_ele[k][0] = buf[j];
	tmp_ele[k][1] = buf[j+1];
	tmp_ele[k][2] = '\0';
	element[k] = tmp_ele[k];
	j++;
      } else {
	tmp_ele[k][0] = buf[j];
	tmp_ele[k][1] = '\0';
	element[k] = tmp_ele[k];
      }
      k++ ;
      if ( k >= MAXELE ) {
	printf("Error: too many elements found.  Increase MAXELE and re-compile\n") ;
	exit(1) ;
      }
    }
  }
  *maxele = k ;
  *maxtype = -1 ;
  for ( j = 0 ; j < natom ; j++ ) {
    type[j] = - 1;
  }
  for ( j = 0 ; j < natom ; j++ ) {
    int ii, tt ;
    double xx, yy, zz ;
      
    fgets(buf,BUFSIZ,fgen) ;
    sscanf(buf,"%d %d %lf %lf %lf", &ii, &tt, &xx, &yy, &zz) ;
    ii-- ;

    if ( ii >= natom || ii < 0 ) {
      printf("Bad atom idex found: %d\n", ii) ;
      exit(1) ;
    }

    type[ii] = tt ;
    if ( type[ii] > *maxtype ) {
      *maxtype = type[ii] ;
    }
    x[ii] = xx ;
    y[ii] = yy ;
    z[ii] = zz ; 
  }
  for ( j = 0 ; j < natom ; j++ ) {
    if ( type[j] == - 1 ) {
      printf("Error:  did not find an entry for atom %d\n",
	     j) ;
      exit(1) ;
    }
  }

  fgets(buf,BUFSIZ,fgen) ;
    
  fgets(buf,BUFSIZ,fgen) ;
  sscanf(buf,"%lf %lf %lf", &a[0], &a[1], &a[2]) ;
  fgets(buf,BUFSIZ,fgen) ;
  sscanf(buf,"%lf %lf %lf", &b[0], &b[1], &b[2]) ;
  fgets(buf,BUFSIZ,fgen) ;
  sscanf(buf,"%lf %lf %lf", &c[0], &c[1], &c[2]) ;
  /*    fgets(buf,BUFSIZ,fgen) ; */
  if ( a[0] <= 0.0 || b[1] <= 0.0 || c[2] <= 0.0 ) {
    printf("Bad cell size read in\n") ;
    exit(1) ;
  }

  if ( feof(fgen) || ferror(fgen) ) {
    return (0) ;
  }
  
  return(natom) ;
}

void check_frame_input(int read_charge, int natom, int *natom0, 
		       int natom1, int frame, int maxele, int *maxele0,
		       int *type, int *type0)
/* Check parameters read in for the given frame against the first
   frame for consistency */
{
  int j ;

  if ( read_charge && natom != natom1 ) {
    printf("Error: number of atoms in the chr and gen files are different\n") ;
    exit(1) ;
  }

  if ( frame == 1 ) {
    *natom0 = natom ;
    if ( maxele <= 0 ) {
      printf("Error: No elements were found\n") ;
      exit(1) ;
    }
    *maxele0 = maxele ;
    for ( j = 0 ; j < natom ; j++ ) {
      type0[j] = type[j] ;
    }
  }
  else {
    if ( natom != *natom0 ) {
      printf("Error: the number of atoms changed\n") ;
      exit(1) ;
    }
    if ( maxele != *maxele0 ) {
      printf("Erorr: the number of elements changed\n") ;
    }
    for ( j = 0 ; j < natom ; j++ ) {
      if ( type0[j] != type[j] ) {
	printf("Error:  the type of atom %d changed\n", j) ;
	exit(1) ;
      }
    }
  }
}

int skip_frame(int frame, int skip1, int skip2, int block, int offset,
	       int *block_start)
/* Returns 1 if the frame should be skipped for further processing.
   Returns 0 otherwise 
   
   skip1 : Every nth frame is skipped.  For instance, only process
           every other frame when skip1 == 2.

   skip2 : Skip every nth block of frames.  For instance, if 
           block == 1000, and skip2 = 2, skip1 = 2, then frames
           1, 3, 5, ... 1000 are read.
           Frames 1001, 1003, 1005, ... 2000 are skipped.
*/
{
  frame-- ;

  if ( offset >= skip2 ) {
    offset %= skip2 ;
  }

  if (   ( (frame/block) % skip2 == offset &&
	   ((frame-1)/block) % skip2 != offset
	 )
       ||
	 ( frame == 0 && offset == 0 ) 
       ) {
    *block_start = frame + 1;
    printf("Skip2 block starts at frame %d\n", *block_start) ;
  }
  if ( frame == 0 ) {
    /* Always read the first frame to get natom0, type0, etc. */
    return 0 ;
  } else if ( skip1 != 1 && frame % skip1 != 0 ) {
    return 1 ;
  } else if ( (frame / block) % skip2 != offset ) {
    return 1 ;
  }
  return 0 ;
}

void read_skip(int *skip1, int *skip2, int *block, int *offset)
{
  FILE *fskip ;
  
  fskip = fopen("skip.dat", "r") ;

  if ( fskip == NULL ) {
    /* printf("skip.dat not found.  All frames will be read\n") ; */
    *skip1 = 1 ;
    *skip2 = 1 ;
    *block = 1 ;
    *offset = 0 ;
  } else {
    printf("skip.dat was found.  Some frames will be skipped\n") ;
    fscanf(fskip, "%d", skip1) ;
    fscanf(fskip, "%d", skip2) ;
    fscanf(fskip, "%d", block) ;
    fscanf(fskip, "%d", offset) ;
    if ( ferror(fskip) || feof(fskip) ) {
      printf("Error: did not successfully read skip.dat\n") ;
      exit(1) ;
    }
    if ( *skip1 <= 0 ) {
      printf("Error: skip1 must be > 0\n") ;
      exit(1) ;
    }
    if ( *skip2 <= 0 ) {
      printf("Error: skip2 must be > 0\n") ;
      exit(1) ;
    }
    if ( *block <= 0 ) {
      printf("Error: block must be > 0\n") ;
      exit(1) ;
    }
    if ( *offset < 0 ) {
      printf("Error: offset must be >= 0\n") ;
      exit(1) ;
    }
    if ( *offset >= *skip2 ) {
      printf("Error: offset must be < skip2\n") ;
      exit(1) ;
    }

  }
}

