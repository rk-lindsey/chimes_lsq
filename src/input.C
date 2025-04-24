#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

#ifdef USE_MPI
	#include <mpi.h>
#endif

#include "util.h"
#include "functions.h"
#include "input.h"
#include "Cheby.h"


string InputListing::operator()(int i, int j)
// Access a particular token in the InputListing.
{
	if ( i >= list.size() )
	{
		if ( RANK == 0 )
		{
			cout << "Input was not found for line " + to_string(i) << endl ;
		}
		exit_run(1) ;
	}
	if ( j >= list[i].size() )
	{
		if ( RANK == 0 )
		{
			cout << "Not enough tokens in line." << endl ;

			cout << "Working on line: " << endl ;
			for ( int j = 0 ; j < list[i].size() ; j++ ) {
				cout << list[i][j] << " " ;
			}
			cout << endl ;
		}
		exit_run(1) ;
	}
	return list[i][j] ;
}

vector<string> &InputListing::operator()(int i)
// Access a vector of tokens for a particular line.
{
	if ( i >= list.size() )
	{
		if ( RANK == 0 )
		{
			cout << "Input was not found for line " + to_string(i) << endl ;
		}
		exit_run(1) ;
	}
	return list[i] ;
}

void InputListing::push_back(vector<string> &item) 
// Add a vector of tokens to the InputListing.
{
	list.push_back(item) ;
}

int InputListing::size()
// Return the number of lines.
{
	return(list.size() ) ;
}

int InputListing::size(int i)
// Return the number of tokens for a particular line.
{
	if ( i >= list.size() )
	{
		if ( RANK == 0 )
		{
			cout << "Input was not found for line " + to_string(i) << endl ;
		}
		exit_run(1) ;
	}
	return( list[i].size() ) ;
}

// Constructor/deconstructor

INPUT::INPUT(string FNAME){ FILENAME = FNAME; }
INPUT::~INPUT(){}

// Small helpers

bool INPUT::found_input_keyword(string KEYWORD, vector<string> & LINE)
{
	// Searching for lines like "# KEYWORD #"
	
	if(LINE.size() >= 3)
		if ((LINE[0] == "#") && (LINE[2] == "#"))
			if (LINE[1] == KEYWORD)
				return true;
	return false;
}
	
	


// For storing the entire input file in memory and checkingf or unrecognized key words

void INPUT::READ_FILE()
{
	// Open file for reading
	
	INFILE.open(FILENAME.data());
	
	if (!INFILE.is_open())
		EXIT_MSG("ERROR: Cannot open input file: ", FILENAME);
	
	if ( RANK == 0 ) 	
		cout << "Input file successfully opened." << endl;
	
	// Check for unrecognized keywords
	
	CHECK_FILE();
	
	INFILE.close();
	
	INFILE.open(FILENAME.data());
	
	if ( RANK == 0 ) 
		cout << "Input file has been rewound." << endl;	
	
	// Parse the file into lines of "words" and save everything to a vec of vecs
	
	int 		KICKOUT = 0;
	string		LINE;
	vector<string>  PARSED_LINE;
	int		N_WORDS = 0;
	
	while (true)
	{
		getline(INFILE, LINE); 	
		strip_comments(LINE) ;	
		
		N_WORDS = parse_space(LINE,PARSED_LINE);
		
		if ( N_WORDS == 0 ) 
			continue;
		
		CONTENTS.push_back(PARSED_LINE);
		
		if(LINE.find("# ENDFILE #") != string::npos)
			break;
			
		KICKOUT++;
		if(KICKOUT==1000)
			EXIT_MSG("Read 1000 lines and have not found \"# ENDFILE # \" in file: ", FILENAME);
	}
	if ( RANK == 0 ) 
		cout << "Input file contents successfully stored." << endl;

	INFILE.close();
}
void INPUT::CHECK_FILE()
{
	// Do something
	return;
}


// Wrapper functions: Process entire input file

void INPUT::PARSE_INFILE_LSQ(  JOB_CONTROL	 & CONTROLS,	       
	    		       vector<PAIRS>	 & ATOM_PAIRS, 
	    		       CLUSTER_LIST		 & TRIPS, 
	    		       CLUSTER_LIST		 & QUADS, 
	    		       map<string,int>   & PAIR_MAP,
	    		       vector<int>		 & INT_PAIR_MAP,
	    		       vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
	    		       NEIGHBORS		 & NEIGHBOR_LIST,
	    		       vector<int>		 & TMP_ATOMTYPEIDX, 
	    		       vector<string>	 & TMP_ATOMTYPE,
			       A_MAT			 & A_MATRIX)
{
	// Actually read the input file 
	
	READ_FILE();
	
	// Initialize the neighbor list
	
	NEIGHBOR_LIST.USE = true;
	
	// For assigning LSQ variables: "Control Variables" 
	
	PARSE_CONTROLS_TRJFILE(CONTROLS);	
	PARSE_CONTROLS_WRAPTRJ(CONTROLS);
	PARSE_CONTROLS_SPLITFI(CONTROLS);
	PARSE_CONTROLS_HIERARC(CONTROLS);
	PARSE_CONTROLS_NFRAMES(CONTROLS);
	PARSE_CONTROLS_NLAYERS(CONTROLS);
	PARSE_CONTROLS_FITCOUL(CONTROLS);
	PARSE_CONTROLS_FITSTRS(CONTROLS);
	PARSE_CONTROLS_FITENER(CONTROLS);
	PARSE_CONTROLS_PAIRTYP(CONTROLS);
	PARSE_CONTROLS_CHBTYPE(CONTROLS);
	PARSE_CONTROLS_CHEBYFIX(CONTROLS);
	PARSE_CONTROLS_USENEIG(CONTROLS, NEIGHBOR_LIST);
	PARSE_CONTROLS_SKIP_FRAMES(CONTROLS) ;
	
	// For assigning LSQ variables: "Topology Variables" 
	
	
	PARSE_TOPOLOGY_NATMTYP(CONTROLS, ATOM_PAIRS, TRIPS, QUADS);
	PARSE_TOPOLOGY_TYPEIDX(CONTROLS, ATOM_PAIRS, PAIR_MAP, TMP_ATOMTYPEIDX, TMP_ATOMTYPE);	
//	PARSE_TOPOLOGY_EXCLUDE(TRIPS, QUADS);
	PARSE_TOPOLOGY_EXCLUDE(ATOM_PAIRS, TRIPS, QUADS, A_MATRIX, CONTROLS, TMP_ATOMTYPEIDX, TMP_ATOMTYPE);	

	PARSE_TOPOLOGY_PAIRIDX(CONTROLS, ATOM_PAIRS, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST, TMP_ATOMTYPEIDX, TMP_ATOMTYPE);
	PARSE_TOPOLOGY_CHGCONS(CHARGE_CONSTRAINTS, PAIR_MAP);
	PARSE_TOPOLOGY_SPECMIN(TRIPS, QUADS);
	PARSE_TOPOLOGY_SPECMAX(TRIPS, QUADS);
	PARSE_TOPOLOGY_FCUTTYP(CONTROLS, ATOM_PAIRS, TRIPS, QUADS);
	
	RUN_SANITY_LSQ(CONTROLS);	

}

void INPUT::PARSE_INFILE_MD (JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)
{
	READ_FILE();
	
	// MD controls initialization
	
	PARSE_CONTROLS_INITIAL(CONTROLS, NEIGHBOR_LIST);
		
	// For self-consistent fitting type runs
	
	PARSE_CONTROLS_SLFCNST();
	
	// "General control variables"
	
	PARSE_CONTROLS_RNDSEED(CONTROLS);
	PARSE_CONTROLS_TEMPERA(CONTROLS);
	PARSE_CONTROLS_PRESSUR(CONTROLS);
	PARSE_CONTROLS_CONVCUT(CONTROLS);
	PARSE_CONTROLS_CMPRFRC(CONTROLS);
	PARSE_CONTROLS_CHCKFRC(CONTROLS);
	PARSE_CONTROLS_SUBTFRC(CONTROLS);
	PARSE_CONTROLS_TIMESTP(CONTROLS);
	PARSE_CONTROLS_N_MDSTP(CONTROLS);
	PARSE_CONTROLS_PENTHRS(CONTROLS);
	PARSE_CONTROLS_KILLLEN(CONTROLS);
	PARSE_CONTROLS_NLAYERS(CONTROLS);
	PARSE_CONTROLS_USENEIG(CONTROLS, NEIGHBOR_LIST);
	PARSE_CONTROLS_PRMFILE(CONTROLS);
	PARSE_CONTROLS_SERIAL_CHIMES(CONTROLS) ;
	PARSE_CONTROLS_CRDFILE(CONTROLS);
	PARSE_CONTROLS_CHEBYFIX(CONTROLS);
	
	// " Simulation options"
	
	PARSE_CONTROLS_VELINIT(CONTROLS);
	PARSE_CONTROLS_CONSRNT(CONTROLS);
	PARSE_CONTROLS_PRSCALC(CONTROLS);
	PARSE_CONTROLS_STRSCALC(CONTROLS);	
	PARSE_CONTROLS_WRPCRDS(CONTROLS);
	
	// "Output control"
	
	PARSE_CONTROLS_ATMENER(CONTROLS);
	PARSE_CONTROLS_FRQDFTB(CONTROLS);
	PARSE_CONTROLS_TRAJEXT(CONTROLS);
	PARSE_CONTROLS_FRQENER(CONTROLS);
	PARSE_CONTROLS_PRNTFRC(CONTROLS);
	PARSE_CONTROLS_PRNTBAD(CONTROLS);
	PARSE_CONTROLS_FRQRSTR(CONTROLS);
	PARSE_CONTROLS_PRNTVEL(CONTROLS);
	PARSE_CONTROLS_GETSTRS(CONTROLS);
	PARSE_CONTROLS_GETENER(CONTROLS);
	PARSE_CONTROLS_FORDFTB(CONTROLS);
	
	// Run MD sanity checks
	
	RUN_SANITY_MD(CONTROLS, NEIGHBOR_LIST);	
	
	
}



// For assigning LSQ variables: "Control Variables" 

void INPUT::PARSE_CONTROLS_TRJFILE(JOB_CONTROL & CONTROLS)
{

	int N_CONTENTS = CONTENTS.size();

	for (int i=0; i<N_CONTENTS; i++)
	{
		if (CONTENTS(i,1) == "TRJFILE")
		{
			// Determine if we're dealing with a single or multiple trajectory files, handle appropriately
			
			if(CONTENTS.size(i+1) == 1) // We have a single file
			{
				CONTROLS.INFILE.push_back(CONTENTS(i+1,0));
				
				CONTROLS.INFILE_FORCE_FLAGS .push_back("");
				CONTROLS.INFILE_STRESS_FLAGS.push_back("");
				CONTROLS.INFILE_ENERGY_FLAGS.push_back("");
				
				if ( RANK == 0 ) 
					cout << "	# TRJFILE #: " << CONTROLS.INFILE[0] << endl;	
				break;
			}
			else if(CONTENTS.size(i+1) == 2) // We have multiple files
			{	
				if ( RANK == 0 ) 
					cout << "	# TRJFILE #: MULTI " << CONTENTS(i+1,1) << endl;
					
				// Read the files (later, will add option to skip first N frames, and to process every M frames)
				
				ifstream MULTI;
				MULTI.open(CONTENTS(i+1,1));
				
				if (!MULTI.is_open())
					EXIT_MSG("ERROR: Cannot open MULTI file: ", CONTENTS(i+1,0));	
				
				if ( RANK == 0 ) 
					cout << "		...Will read from the following files: " << endl;
				
				
				string		LINE;
				vector<string>  PARSED_LINE;
				
				getline(MULTI, LINE);
				parse_space(LINE,PARSED_LINE);
				
				if(PARSED_LINE.size() != 1)
					EXIT_MSG("ERROR: Expected to read <nfiles>, got: ", LINE);	
				
				int NFILES = convert_int(PARSED_LINE[0],i+1);
				
				for (int j=0; j<NFILES; j++)
				{
					getline(MULTI, LINE);

					parse_space(LINE,PARSED_LINE);
					
					if( (PARSED_LINE.size() != 2) && (PARSED_LINE.size() != 3) && (PARSED_LINE.size() != 6) && (PARSED_LINE.size() != 5))
						EXIT_MSG("ERROR: Expected to read <nframes> <filename> with <temperature>, <temperature> <F_flag> <S_flag>,<E_flag>, or <F_flag> <S_flag> <E_flag>,  got: ", LINE);											
					
					CONTROLS.INFILE_FRAMES.push_back(convert_int(PARSED_LINE[0],i+1));
					CONTROLS.INFILE       .push_back(     PARSED_LINE[1]);
					
					if (PARSED_LINE.size() >= 5)
					{
						string tmp_a = PARSED_LINE[ PARSED_LINE.size()-3 ];
						string tmp_b = PARSED_LINE[ PARSED_LINE.size()-2 ];
						string tmp_c = PARSED_LINE[ PARSED_LINE.size()-1 ];
						
						if (tmp_a != "None")
							CONTROLS.INFILE_FORCE_FLAGS .push_back(tmp_a);
						else
							CONTROLS.INFILE_FORCE_FLAGS .push_back("");
							
						if (tmp_b != "None")
							CONTROLS.INFILE_STRESS_FLAGS .push_back(tmp_b);
						else
							CONTROLS.INFILE_STRESS_FLAGS .push_back("");
							
						if (tmp_c != "None")
							CONTROLS.INFILE_ENERGY_FLAGS .push_back(tmp_c);
						else
							CONTROLS.INFILE_ENERGY_FLAGS .push_back("");														
					}
					else
					{
						CONTROLS.INFILE_FORCE_FLAGS .push_back("");
						CONTROLS.INFILE_STRESS_FLAGS.push_back("");
						CONTROLS.INFILE_ENERGY_FLAGS.push_back("");				     
					}
					
					
					if ( RANK == 0 ) 
						cout << "		-   " 
						     << CONTROLS.INFILE[j]        << " with " 
						     << CONTROLS.INFILE_FRAMES[j] << " frames." << endl;
				}
				
				MULTI.close();	
				
				// Clean up any files with zero frames 
				if ( RANK == 0 ) 
					cout << "Removing the following files with zero frames:" << endl;
				
				
				// Need to update 
				
				for(int i=CONTROLS.INFILE_FRAMES.size()-1; i>=0; i--)
    				{
    					if (CONTROLS.INFILE_FRAMES[i] < 1)
					{
					
						if ( RANK == 0 ) 
							cout << "		-   " 
							     << CONTROLS.INFILE[i]        << " with " 
							     << CONTROLS.INFILE_FRAMES[i] << " frames." << endl;					
						
						CONTROLS.INFILE_FORCE_FLAGS .erase( CONTROLS.INFILE_FORCE_FLAGS .begin() + i );
						CONTROLS.INFILE_STRESS_FLAGS.erase( CONTROLS.INFILE_STRESS_FLAGS.begin() + i );
						CONTROLS.INFILE_ENERGY_FLAGS.erase( CONTROLS.INFILE_ENERGY_FLAGS.begin() + i );
						CONTROLS.INFILE             .erase( CONTROLS.INFILE		.begin() + i );
						CONTROLS.INFILE_FRAMES      .erase( CONTROLS.INFILE_FRAMES	.begin() + i );
					}
				}
				break;			
			}
			else
			{
				cout << "ERROR: Unrecognized TRAJFILE option: ";
				for (int j=0; j<CONTENTS.size(i+1); j++)
					cout << CONTENTS(i+1,j) << " ";
				cout << endl;
				exit_run(0);
			}
		}
	}
}

void INPUT::PARSE_CONTROLS_WRAPTRJ(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{	
		if (found_input_keyword("WRAPTRJ", CONTENTS(i)))
		{
			if((i+1 < N_CONTENTS) && (CONTENTS.size(i+1) == 1)) // We have a single file
			{
				CONTROLS.WRAP_COORDS = convert_bool(CONTENTS(i+1,0),i+1);
				
				if ( RANK == 0 ) 
					cout << "	# WRAPTRJ #: " << bool2str(CONTROLS.WRAP_COORDS) << endl;
				break;
			}		
		}
	}
}

void INPUT::PARSE_CONTROLS_SPLITFI(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("SPLITFI", CONTENTS(i)))
		{
			CONTROLS.SPLIT_FILES = convert_bool(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# SPLITFI #: " << bool2str(CONTROLS.SPLIT_FILES) << endl;	
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_HIERARC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("HIERARC", CONTENTS(i)))
		{
			CONTROLS.HIERARCHICAL_FIT = convert_bool(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# HIERARC #: " << bool2str(CONTROLS.HIERARCHICAL_FIT) << endl;	
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_NFRAMES(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("NFRAMES", CONTENTS(i)))
		{
			CONTROLS.NFRAMES = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 	
					cout << "	# NFRAMES #: " << CONTROLS.NFRAMES << endl;			
			
			break;
		}
	}
}
/*
void INPUT::PARSE_CONTROLS_NLAYERS(JOB_CONTROL & CONTROLS)
{
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (CONTENTS(i,1) == "NLAYERS")
		{
			CONTROLS.N_LAYERS = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 		
					cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;
			
			break;
		}
	}	
}
*/ // JUST USE THE MD VERSION... IT SHOULD BE COMPATIBLE

void INPUT::PARSE_CONTROLS_FITCOUL(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FITCOUL", CONTENTS(i)))
		{
			CONTROLS.FIT_COUL = convert_bool(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# FITCOUL #: " << bool2str(CONTROLS.FIT_COUL) << endl;				
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_FITSTRS(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FITSTRS", CONTENTS(i)))
		{
			if (CONTENTS(i+1,0)=="first"  || CONTENTS(i+1,0)=="First"  || CONTENTS(i+1,0)=="FIRST")
			{
				CONTROLS.FIT_STRESS = true;
				CONTROLS.NSTRESS    = convert_int(CONTENTS(i+1,1),i+1);
			}			
			else if (CONTENTS(i+1,0)=="all"  || CONTENTS(i+1,0)=="All"  || CONTENTS(i+1,0)=="ALL"  || CONTENTS(i+1,0) == "A" || CONTENTS(i+1,0) == "a")
			{
					CONTROLS.FIT_STRESS_ALL = true;
			}
			else if(CONTENTS(i+1,0)=="firstall"  || CONTENTS(i+1,0)=="FirstAll"  || CONTENTS(i+1,0)=="FIRSTALL")
			{
					CONTROLS.FIT_STRESS_ALL = true;
					CONTROLS.NSTRESS    = convert_int(CONTENTS(i+1,1),i+1);			
			}
			else
				CONTROLS.FIT_STRESS = convert_bool(CONTENTS(i+1,0),i+1);

			if ( RANK == 0 ) 
			{
				cout << "	# FITSTRS #: ";		
							
				if (CONTROLS.FIT_STRESS_ALL)
					cout << bool2str(CONTROLS.FIT_STRESS_ALL) << " ...will fit to all tensor components" << endl;	
				else if(CONTROLS.NSTRESS>0)
					cout << bool2str(CONTROLS.FIT_STRESS) << " ...will only fit tensors for first " << CONTROLS.NSTRESS << " frames." << endl;
				else 
					cout << bool2str(CONTROLS.FIT_STRESS) << endl;
			}
			
			break;
		}
		
	}
}
void INPUT::PARSE_CONTROLS_FITENER(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FITENER", CONTENTS(i)))		
		{
			if (CONTENTS(i+1,0)=="first"  || CONTENTS(i+1,0)=="First"  || CONTENTS(i+1,0)=="FIRST")
			{
				CONTROLS.FIT_ENER = true;
				CONTROLS.NENER    = convert_int(CONTENTS(i+1,1),i+1);
			}
			else
			{
				CONTROLS.FIT_ENER = convert_bool(CONTENTS(i+1,0),i+1);
			}
			
			if ( RANK == 0 )
			{
				cout << "	# FITENER #: " << bool2str(CONTROLS.FIT_ENER) << endl;	
				if(CONTROLS.NENER>=0)
					cout << "    			 ...will fit energies for first " << CONTROLS.NENER << " frames." << endl;
			}	
		}
	}
}

void INPUT::PARSE_CONTROLS_PAIRTYP(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PAIRTYP", CONTENTS(i)))								
		{
			FF_TYPE = CONTENTS(i+1,0);
			
			if (FF_TYPE != "CHEBYSHEV") // These are not supported:  && FF_TYPE != "LJ" && FF_TYPE != "STILLIN")
			{
				cout << endl;
				cout << "ERROR: Unrecognized pair type. Acceptable options are:" << endl;
				cout << "CHEBYSHEV"    << endl;
				exit(1);
				
			}
			
			if ( RANK == 0 )
			{
				cout << "	# PAIRTYP #: " << FF_TYPE;
				cout << " ....NOTE: Forces reported in units of kcal/(mol.A), potential energy in kcal/mol." << endl;
			}
			
			if(FF_TYPE == "CHEBYSHEV")
			{
				CONTROLS.CHEBY_ORDER = convert_int(CONTENTS(i+1,1),i+1);
				
				if ( RANK == 0 ) 
					cout << "	             " << "Will use 2-body order: " << CONTROLS.CHEBY_ORDER << endl;
				
				if(FF_TYPE == "CHEBYSHEV")
				{
					if(CONTENTS.size(i+1) >= 3)
						CONTROLS.CHEBY_3B_ORDER = convert_int(CONTENTS(i+1,2),i+1);
					else
						CONTROLS.CHEBY_3B_ORDER = 0;

					if ( RANK == 0 ) 
						cout << "	             " << "Will use 3-body order: " << CONTROLS.CHEBY_3B_ORDER << endl;
						
					if(CONTROLS.CHEBY_3B_ORDER>0)
						CONTROLS.USE_3B_CHEBY = true;
						
					if(CONTENTS.size(i+1) >= 4)
						 CONTROLS.CHEBY_4B_ORDER = convert_int(CONTENTS(i+1,3),i+1);
					else
						CONTROLS.CHEBY_4B_ORDER = 0;

					if ( RANK == 0 ) 
						cout << "	             " << "Will use 4-body order: " << CONTROLS.CHEBY_4B_ORDER << endl;
						
					if(CONTROLS.CHEBY_4B_ORDER>0)
						CONTROLS.USE_4B_CHEBY = true;
					
					if(CONTENTS.size(i+1)>4)
						TMP_CHEBY_RANGE_LOW = convert_int(CONTENTS(i+1,4),i+1);
					else
						TMP_CHEBY_RANGE_LOW = -1;

					if(TMP_CHEBY_RANGE_LOW < -1.0 || TMP_CHEBY_RANGE_LOW > +1.0 )
						EXIT_MSG("ERROR: TMP_CHEBY_RANGE_LOW must be betwee -1 and 1");	

					if(CONTENTS.size(i+1)>5)
						TMP_CHEBY_RANGE_HIGH = convert_int(CONTENTS(i+1,5),i+1);
					else
						TMP_CHEBY_RANGE_HIGH = 1;

					if((TMP_CHEBY_RANGE_HIGH < -1.0) || (TMP_CHEBY_RANGE_HIGH > +1.0 ))
						EXIT_MSG("ERROR: TMP_CHEBY_RANGE_HIGH must be betwee -1 and 1");	

					if(TMP_CHEBY_RANGE_HIGH < TMP_CHEBY_RANGE_LOW)
						EXIT_MSG("ERROR: TMP_CHEBY_RANGE_HIGH must be greater than TMP_CHEBY_RANGE_LOW");	

					if ( RANK == 0 ) 
						cout << "	             " << "Will transform Chebyshev pair distances to range " << TMP_CHEBY_RANGE_LOW << " to " << TMP_CHEBY_RANGE_HIGH << endl;
				}
			}
			
			break;
		}
	}		
}
void INPUT::PARSE_CONTROLS_CHBTYPE(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CHBTYPE", CONTENTS(i)))			
		{
			if(FF_TYPE == "CHEBYSHEV" )
			{
				CONTROLS.CHEBY_TYPE = Cheby::get_trans_type(CONTENTS(i+1,0));
			
				if ( RANK == 0 ) 
					cout << "	# CHBTYPE #: " << CONTENTS(i+1,0) << endl;	
			}
			else if ( RANK == 0 ) 
			{
				cout << "Warning: CHBTYPE given for a non-Chebyshev pair type (ignored)" ;
			}
		}
	}
}

void INPUT::PARSE_CONTROLS_CHEBYFIX(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CHEBYFIX", CONTENTS(i)))			
		{
			string fix_type = CONTENTS(i+1,0) ;
			if(fix_type == "ZERO_DERIV" )
			{
				CONTROLS.cheby_fix_type = Cheby_fix::ZERO_DERIV ;
			}
			else if ( fix_type == "CONSTANT_DERIV" )
			{
				CONTROLS.cheby_fix_type = Cheby_fix::CONSTANT_DERIV ;				
			}
			else if ( fix_type == "SMOOTH" )
			{
				CONTROLS.cheby_fix_type = Cheby_fix::SMOOTH ;
				CONTROLS.cheby_smooth_distance = convert_double(CONTENTS(i+1,1),i+1) ;
			}
			else 
			{
				EXIT_MSG("ERROR: CHEBYFIX was not recognized") ;
			}
			if ( RANK == 0 )
			{
				cout << "\t# CHEBYFIX #: " << fix_type << endl ;
				if ( fix_type == "SMOOTH" )
					cout << "\t\tCheby smoothing distance = " << CONTROLS.cheby_smooth_distance << endl ;
			}
		}
	}
}



// For assigning LSQ variables: "Topology Variables" 

void INPUT::PARSE_TOPOLOGY_EXCLUDE(vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS, A_MAT & A_MATRIX, JOB_CONTROL & CONTROLS, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE)
{
	int N_CONTENTS = CONTENTS.size();
	A_MATRIX.N_EXCL_LT_3B = 0;
	
	string SEARCH;
	
	if (CONTROLS.HIERARCHICAL_FIT)
	{
		// Search for 1-body exlusion list
		
		string         TMP_TYPE;
		vector<string> TMP_TYPES;
		
		for (int i=0; i<N_CONTENTS; i++)
		{
			if(CONTENTS.size(i) >= 4)
			{
				SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
				
				if (SEARCH == "EXCLUDE 1B INTERACTION:")
				{
					A_MATRIX.DO_EXCLUDE_1B = true;
	
					// This next block should go in a A-matrix function:
					
					if (RANK == 0)
						cout << "Number of 1-body interactions to exclude: " << stoi(CONTENTS(i,3)) << endl;
											
					for(int j=1; j<=stoi(CONTENTS(i,3)); j++) // iterates over excluded atom names
					{
						for(int k=0; k<TMP_ATOMTYPEIDX.size(); k++)
						{
							if (CONTENTS(i+j,0) == TMP_ATOMTYPE[k])
							{
								A_MATRIX.EXCLUDE_1B.push_back(TMP_ATOMTYPEIDX[k]);
								A_MATRIX.N_EXCL_LT_3B += 1;
								
								if (RANK==0)
									cout << "\tEXCLUDED: " <<  A_MATRIX.EXCLUDE_1B[A_MATRIX.EXCLUDE_1B.size()-1] << ": " << TMP_ATOMTYPE[k] << endl;
							}
						}
					}
					break; 
				}
			}
		}		
		// Search for 2-body exlusion list
			
		for (int i=0; i<N_CONTENTS; i++)
		{
			if(CONTENTS.size(i) >= 4)
			{
				SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
				
				if (SEARCH == "EXCLUDE 2B INTERACTION:")
				{
					A_MATRIX.DO_EXCLUDE_2B = true;
					
					// This next block should go in a A-matrix function:
					
					if (RANK==0)
						cout << "Number of 2-body interactions to exclude: " << stoi(CONTENTS(i,3)) << endl;
					
					for(int j=1; j<=stoi(CONTENTS(i,3)); j++) // iterates over excluded atom pair names
					{
						// Get the pair type, index of that pair, and number of parameters for that pair
						for(int k=0; k<ATOM_PAIRS.size(); k++)
						{
							if (
							((ATOM_PAIRS[k].ATM1TYP == CONTENTS(i+j,0)) && (ATOM_PAIRS[k].ATM2TYP == CONTENTS(i+j,1))) ||
							((ATOM_PAIRS[k].ATM1TYP == CONTENTS(i+j,1)) && (ATOM_PAIRS[k].ATM2TYP == CONTENTS(i+j,0)))  )
							{
								A_MATRIX.  EXCLUDE_2B.push_back(ATOM_PAIRS[k].PAIRIDX);
								A_MATRIX.N_EXCLUDE_2B.push_back(CONTROLS.CHEBY_ORDER);
								A_MATRIX.N_EXCL_LT_3B += CONTROLS.CHEBY_ORDER;
								
								if (RANK==0)
									cout << "\tEXCLUDED: " << A_MATRIX.EXCLUDE_2B[A_MATRIX.EXCLUDE_2B.size()-1] << " " << ATOM_PAIRS[k].ATM1TYP << " " << ATOM_PAIRS[k].ATM2TYP << endl;
							}
						}
					}

					 break; 
				}
			}
		}	
	}
	
	// Search for 3-body exlusion list
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "EXCLUDE 3B INTERACTION:")
			{
				 TRIPS.read_exclude(CONTENTS, i); 

				 break; 
			}
		}
	}
			
	// Search for 4-body exlusion list
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "EXCLUDE 4B INTERACTION:")
			{
				 QUADS.read_exclude(CONTENTS, i); 
				 break;
			}
		}
	}
}					
void INPUT::PARSE_TOPOLOGY_NATMTYP(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("NATMTYP", CONTENTS(i)))
		{
			CONTROLS.NATMTYP = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# NATMTYP #: " << CONTROLS.NATMTYP << endl;	
			
			// Set up pairs
	
			NPAIR = CONTROLS.NATMTYP*(CONTROLS.NATMTYP+1)/2;
			ATOM_PAIRS.resize(NPAIR);
			
			// Set the default cheby range and the fcut type
	
			if (FF_TYPE == "CHEBYSHEV")
			{
				for(int j=0; j<NPAIR; j++)
				{
					ATOM_PAIRS[j].CHEBY_RANGE_LOW   = TMP_CHEBY_RANGE_LOW;
					ATOM_PAIRS[j].CHEBY_RANGE_HIGH  = TMP_CHEBY_RANGE_HIGH;
					ATOM_PAIRS[j].FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;
					ATOM_PAIRS[j].PAIRTYP           = FF_TYPE;
				}
			}	
				
			// Set up triplets
	
			NTRIP = factorial(CONTROLS.NATMTYP+3-1)/factorial(3)/factorial(CONTROLS.NATMTYP-1);
			TRIPS.allocate(NTRIP, 3, ATOM_PAIRS);
	
			// Set up quadruplets
	
			NQUAD = factorial(CONTROLS.NATMTYP+4-1)/factorial(4)/factorial(CONTROLS.NATMTYP-1);
			QUADS.allocate(NQUAD, 4, ATOM_PAIRS);
		}
	}
}
void INPUT::PARSE_TOPOLOGY_TYPEIDX(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE)
{
	int N_CONTENTS = CONTENTS.size();
	
	string TMP_STR;
	int    TEMP_INT;
	double SUM_OF_CHARGES;
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("TYPEIDX", CONTENTS(i)))		
		{
	
			if ( RANK == 0 ) 
				cout << "	# TYPEIDX #    # ATM_TYP #    # ATMCHRG #    # ATMMASS #" << endl;
	
			// Figure out the number of non-unique pairs

			TEMP_INT = CONTROLS.NATMTYP*CONTROLS.NATMTYP;
	
			SUM_OF_CHARGES = 0;
	
			TMP_ATOMTYPEIDX.resize(CONTROLS.NATMTYP);
			TMP_ATOMTYPE   .resize(CONTROLS.NATMTYP);
			
			for(int j=0; j<CONTROLS.NATMTYP; j++)
			{
		
				// Set the first atom pair types to be of type OO, HH, CC, etc...
		
				ATOM_PAIRS[j].PAIRTYP    = FF_TYPE;
				ATOM_PAIRS[j].PAIRIDX	 = j; 
				ATOM_PAIRS[j].CHEBY_TYPE = CONTROLS.CHEBY_TYPE;
				
				TMP_ATOMTYPEIDX[j]           = convert_int(CONTENTS(i+1+j,0),i+1+j) - 1;
				ATOM_PAIRS[j].ATM1TYPE_IDX   = TMP_ATOMTYPEIDX[j];

				if ( (TMP_ATOMTYPEIDX[j] < 0) || (TMP_ATOMTYPEIDX[j] >= CONTROLS.NATMTYP ))
					EXIT_MSG("Bad atom index in TYPEIDX section: " + CONTENTS(i+1+j,0) );

				ATOM_PAIRS[j].ATM1TYP = CONTENTS(i+1+j,1);
				TMP_ATOMTYPE[j]       = ATOM_PAIRS[j].ATM1TYP;
				
				if(!CONTROLS.FIT_COUL)
				{
					ATOM_PAIRS[j].ATM1CHG = convert_double(CONTENTS(i+1+j,2),i+1+j);
				}
				else
				{
					ATOM_PAIRS[j].CHRGSGN = CONTENTS(i+1+j,2);
					ATOM_PAIRS[j].ATM1CHG = 0.0;
				}

				SUM_OF_CHARGES += abs(ATOM_PAIRS[j].ATM1CHG);

				ATOM_PAIRS[j].ATM1MAS = convert_double(CONTENTS(i+1+j,3),i+1+j);
		
				ATOM_PAIRS[j].ATM2TYP      = ATOM_PAIRS[j].ATM1TYP;
				ATOM_PAIRS[j].ATM2TYPE_IDX = ATOM_PAIRS[j].ATM1TYPE_IDX;

				ATOM_PAIRS[j].ATM2CHG = ATOM_PAIRS[j].ATM1CHG;
				ATOM_PAIRS[j].ATM2MAS = ATOM_PAIRS[j].ATM1MAS;
		
				ATOM_PAIRS[j].PRPR_NM =      ATOM_PAIRS[j].ATM1TYP;
				ATOM_PAIRS[j].PRPR_NM.append(ATOM_PAIRS[j].ATM2TYP);
				
				if ( RANK == 0 ) 
				{
					cout << " 	" << setw(15) << left << j+1 << setw(15) << left << ATOM_PAIRS[j].ATM1TYP;
					if(!CONTROLS.FIT_COUL) 
						cout << setw(15) << left << ATOM_PAIRS[j].ATM1CHG;
					else
						cout << ATOM_PAIRS[j].CHRGSGN << "		";
			
					cout << setw(15) << left << ATOM_PAIRS[j].ATM1MAS << endl;
				}
			}
			
			break;
		}		
	}
			
	if(!CONTROLS.FIT_COUL)
	{
		if (SUM_OF_CHARGES>0)
		{
			CONTROLS.USE_PARTIAL_CHARGES = true;
			CONTROLS.IF_SUBTRACT_COUL    = true;
		}
		else
			CONTROLS.USE_PARTIAL_CHARGES = false;				
	}			

	// Set up all possible unique pair types
	
	TEMP_INT = CONTROLS.NATMTYP;
	
	for(int i=0; i<CONTROLS.NATMTYP; i++)
	{
		for(int j=i+1; j<CONTROLS.NATMTYP; j++)
		{	
			ATOM_PAIRS[TEMP_INT].PAIRTYP    = ATOM_PAIRS[i].PAIRTYP;
			ATOM_PAIRS[TEMP_INT].PAIRIDX    = TEMP_INT; 
			ATOM_PAIRS[TEMP_INT].CHEBY_TYPE = ATOM_PAIRS[i].CHEBY_TYPE;
												
			ATOM_PAIRS[TEMP_INT].ATM1TYP = ATOM_PAIRS[i].ATM1TYP;
			ATOM_PAIRS[TEMP_INT].ATM1TYPE_IDX = ATOM_PAIRS[i].ATM1TYPE_IDX;

			ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
			ATOM_PAIRS[TEMP_INT].ATM2TYPE_IDX = ATOM_PAIRS[j].ATM1TYPE_IDX;

			ATOM_PAIRS[TEMP_INT].ATM1CHG = ATOM_PAIRS[i].ATM1CHG;
			ATOM_PAIRS[TEMP_INT].ATM2CHG = ATOM_PAIRS[j].ATM1CHG;
			
			ATOM_PAIRS[TEMP_INT].ATM1MAS = ATOM_PAIRS[i].ATM1MAS;
			ATOM_PAIRS[TEMP_INT].ATM2MAS = ATOM_PAIRS[j].ATM1MAS;	
			
			ATOM_PAIRS[TEMP_INT].PRPR_NM =      ATOM_PAIRS[TEMP_INT].ATM1TYP;
			ATOM_PAIRS[TEMP_INT].PRPR_NM.append(ATOM_PAIRS[TEMP_INT].ATM2TYP);											
			
			TEMP_INT++;
		}	
	}

	// Set up the maps to account for non-unique pairs
	
	for(int i=0; i<CONTROLS.NATMTYP; i++)
	{
		for(int j=0; j<CONTROLS.NATMTYP; j++)
		{
			TMP_STR = ATOM_PAIRS[i].ATM1TYP;
			TMP_STR.append(ATOM_PAIRS[j].ATM1TYP);
			
			for(int k=0;k<ATOM_PAIRS.size(); k++)
			{
				if((ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM2TYP)
				 ||(ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM2TYP))
				{
					TMP_STR = ATOM_PAIRS[i].ATM1TYP;
					TMP_STR.append(ATOM_PAIRS[j].ATM2TYP);
					PAIR_MAP.insert(make_pair(TMP_STR,k));	// Maps the true pair index, k, to the string formed by joining the two atom types.
				}
			}
		}
	}
	
	//cout << "Made the following pairs: " << endl;
	//for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
	//cout << i->second << " " << i->first << endl;			
			
	if ( RANK == 0 ) 
	{
		cout << endl;
		cout << "	The following unique pair types have been identified:" << endl;
		for(int i=0;i<NPAIR; i++)
			cout << "		" << ATOM_PAIRS[i].PAIRIDX << "  " << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[i].ATM2TYP << endl;
	}	
}
		
void INPUT::PARSE_TOPOLOGY_PAIRIDX(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE)
{
	int N_CONTENTS = CONTENTS.size();
	
	string TEMP_STR;
	int    TEMP_INT; 

	for (int c=0; c< N_CONTENTS ; c++)
	{
		if (found_input_keyword("PAIRIDX", CONTENTS(c)))		
		{
			for(int i=0; i<NPAIR; i++)
			{
				TEMP_PAIR.ATM1TYP = CONTENTS(c+1+i,1);
				TEMP_PAIR.ATM2TYP = CONTENTS(c+1+i,2);

				TEMP_STR = TEMP_PAIR.ATM1TYP;
				TEMP_STR.append(TEMP_PAIR.ATM2TYP);
				TEMP_INT = PAIR_MAP[TEMP_STR];
				
				ATOM_PAIRS[TEMP_INT].PAIRIDX = TEMP_INT;

				ATOM_PAIRS[TEMP_INT].S_MINIM = convert_double(CONTENTS(c+1+i,3),c+1+i);
				ATOM_PAIRS[TEMP_INT].S_MAXIM = convert_double(CONTENTS(c+1+i,4),c+1+i);
				ATOM_PAIRS[TEMP_INT].LAMBDA  = convert_double(CONTENTS(c+1+i,6),c+1+i);
				
				ATOM_PAIRS[TEMP_INT].MIN_FOUND_DIST = 	1.0e10;	// Set an initial minimum distance	
				
				if( ATOM_PAIRS[TEMP_INT].S_MAXIM > NEIGHBOR_LIST.MAX_CUTOFF)
				{
					NEIGHBOR_LIST.MAX_CUTOFF    = ATOM_PAIRS[TEMP_INT].S_MAXIM;
					NEIGHBOR_LIST.MAX_CUTOFF_3B = ATOM_PAIRS[TEMP_INT].S_MAXIM;
					NEIGHBOR_LIST.MAX_CUTOFF_4B = ATOM_PAIRS[TEMP_INT].S_MAXIM;
				}	
			}
			break;
		}
	}
	

	if ( RANK == 0 ) 
	{
		cout << "	# PAIRIDX #     ";
		cout << "# ATM_TY1 #     ";
		cout << "# ATM_TY1 #     ";
		cout << "# S_MINIM #     ";
		cout << "# S_MAXIM #     ";
		cout << "# MORSE_LAMBDA #" << endl;	
	}
	
	for(int i=0; i<NPAIR; i++)
	{
		if ( RANK == 0 ) 
		{
			cout << "	" 
				 << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
				 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
				 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
				 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
				 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
				 << setw(16) << left << ATOM_PAIRS[i].LAMBDA << endl;
		}
	}

	build_int_pair_map(CONTROLS.NATMTYP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX, PAIR_MAP, INT_PAIR_MAP);

	for ( int i = 0; i < ATOM_PAIRS.size(); i++ ) 
		if ( ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 
			ATOM_PAIRS[i].set_cheby_vals();
}

void INPUT::PARSE_TOPOLOGY_CHGCONS(vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, map<string,int> & PAIR_MAP)
{
	int N_CONTENTS = CONTENTS.size();
	
	string SEARCH;
	
	// Search for charge constraints
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 2)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1);
			
			if (SEARCH == "CHARGE CONSTRAINTS:")
			{
				if ( RANK == 0 ) 
					cout << endl << "	Attempting to read " << NPAIR-1 << " charge constraints...:" << endl; 
				
				CHARGE_CONSTRAINTS.resize(NPAIR-1);
				
				for(int n=0; n<NPAIR-1; n++)
				{
					// Read the atom pair types
				
					for(int j=0; j<NPAIR; j++)
					{
						CHARGE_CONSTRAINTS[n].PAIRTYPE    .push_back(         CONTENTS(i+1+n,j));
						CHARGE_CONSTRAINTS[n].PAIRTYPE_IDX.push_back(PAIR_MAP[CONTENTS(i+1+n,j)]);
					}
				
					// Read the constraints 
				
					for(int j=0; j<NPAIR; j++)
						CHARGE_CONSTRAINTS[n].CONSTRAINTS.push_back(convert_double(CONTENTS(i+1+n,NPAIR+j),i+1+n)); 
				
					// Read the associated force
					CHARGE_CONSTRAINTS[n].FORCE = convert_double(CONTENTS(i+1+n,CONTENTS.size(i+1+n)-1),i+1+n);				

		
					if ( RANK == 0 ) 
					{
						cout << "		" << n+1 << "	 ";
						for(int j=0; j<NPAIR; j++)
							cout << CHARGE_CONSTRAINTS[n].PAIRTYPE[j] << " (" << CHARGE_CONSTRAINTS[n].PAIRTYPE_IDX[j] << ") ";
						for(int j=0; j<NPAIR; j++)
							cout << CHARGE_CONSTRAINTS[n].CONSTRAINTS[j] << " ";
						cout << CHARGE_CONSTRAINTS[n].FORCE << endl;	
					}	
					
				}
			
				if ( RANK == 0 ) 
					cout << endl;
				 break;
			}
		}
	}
}
void INPUT::PARSE_TOPOLOGY_SPECMIN(CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS)
{
	int N_CONTENTS = CONTENTS.size();
	
	string SEARCH;
	
	// Search for 3-body special minimum values
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "SPECIAL 3B S_MINIM:")
			{
				 TRIPS.read_cutoff_params(CONTENTS, i, "S_MINIM");
				 break;
			}
		}
	}
			
	// Search for 4-body special minimum values
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "SPECIAL 4B S_MINIM:")
			{
				 QUADS.read_cutoff_params(CONTENTS, i, "S_MINIM");
				 break;
			}
		}
	}
}	
void INPUT::PARSE_TOPOLOGY_SPECMAX(CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS)
{
	int N_CONTENTS = CONTENTS.size();
	
	string SEARCH;
	
	// Search for 3-body special maximum values
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "SPECIAL 3B S_MAXIM:")
			{
				 TRIPS.read_cutoff_params(CONTENTS, i, "S_MAXIM");
				 break;
			}
		}
	}
			
	// Search for 4-body exlusion list
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if(CONTENTS.size(i) >= 4)
		{
			SEARCH = CONTENTS(i,0) + ' ' + CONTENTS(i,1) + ' ' + CONTENTS(i,2);
			
			if (SEARCH == "SPECIAL 4B S_MAXIM:")
			{
				 QUADS.read_cutoff_params(CONTENTS, i, "S_MAXIM");
				 break;
			}
		}
	}
}	
void INPUT::PARSE_TOPOLOGY_FCUTTYP(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FCUTTYP", CONTENTS(i)))				
		{
			CONTROLS.FCUT_LINE = join_string_vec(CONTENTS(i+1),' ');
			
			parse_fcut_input(CONTROLS.FCUT_LINE, ATOM_PAIRS, TRIPS, QUADS) ;
			
			if ( RANK == 0 ) 
				cout << "# FCUTTYP #: " << CONTROLS.FCUT_LINE << "      ... for all Chebyshev interactions" << endl ;
			
			break;
		}
	}	
}		
	

// Run LSQ sanity checks

void INPUT::RUN_SANITY_LSQ(JOB_CONTROL & CONTROLS)
{

	// Run a few checks to make sure logic is correct

	if(CONTROLS.IF_SUBTRACT_COUL && CONTROLS.FIT_COUL)
	{
		cout << "LOGIC ERROR: Problem with code logic. Both fit_coul and ifsubtract_coul cannot be true." << endl;
		cout << "             ifsubtract_coul should only be true if non-zero charges have been specified " << endl;
		cout << "             and FITCOUL set false." << endl;
		exit_run(0);
	}				

	if(CONTROLS.USE_PARTIAL_CHARGES && RANK == 0)
	{
		cout << "Special feature: " << endl;
		cout << " Will subtract contributions stemming from user-specified " << endl;
		cout << " charges before generating A-matrix" << endl << endl;	
	}	

	if(CONTROLS.WRAP_COORDS && CONTROLS.N_LAYERS > 0 )
	{
		if ( RANK == 0 ) 
			cout << "WARNING: Coordinate wrapping not supported for ghost atom use. Turning option off" << endl;
		CONTROLS.WRAP_COORDS = false;
	}

	if (CONTROLS.USE_PARTIAL_CHARGES || CONTROLS.FIT_COUL)
		CONTROLS.CALL_EWALD = true;	
		
	if((CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL) && CONTROLS.CALL_EWALD)
		EXIT_MSG("ERROR: Inclusion of stress tensors currently not compatible with use of ZCalc_Ewald_Deriv.") ;
}
		
		
// MD controls initialization

void INPUT::PARSE_CONTROLS_INITIAL(JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST) 		
{
	// Set some defaults
	
	CONTROLS.IS_LSQ           = false;
	CONTROLS.SELF_CONSIST     = false;
	CONTROLS.SUBTRACT_FORCE   = false;
	CONTROLS.WRAP_COORDS      = true;
	CONTROLS.PRINT_VELOC      = false;
	CONTROLS.NVT_CONV_CUT     = 0.10;
	CONTROLS.FREEZE_IDX_START = -1; 
	CONTROLS.FREEZE_IDX_STOP  = -1; 
	CONTROLS.SCALE_SYSTEM_BY  = 1.0;
	CONTROLS.BUILD            = false;
	CONTROLS.FIT_STRESS       = false;
	CONTROLS.FIT_ENER         = false;
	CONTROLS.CHECK_FORCE      = false;
	CONTROLS.FORDFTB          = false;
	
	CONTROLS.RESTART          = false;
	CONTROLS.INIT_VEL         = false;
	
	CONTROLS.FREQ_BACKUP     	= 100;
	CONTROLS.FREQ_UPDATE_THERMOSTAT = -1.0;
	CONTROLS.USE_HOOVER_THRMOSTAT	= false;
	CONTROLS.USE_NUMERICAL_PRESS = false ;
	CONTROLS.USE_NUMERICAL_STRESS = false ;	
	CONTROLS.REAL_REPLICATES        = 0;
	NEIGHBOR_LIST.USE               = true;
	
	CONTROLS.PRINT_BAD_CFGS         = false;
	CONTROLS.TRAJ_FORMAT			= "GEN";
	
	CONTROLS.PENALTY_THRESH = -1.0;	// Default is to not enforce a max-allowed penalty
	CONTROLS.IO_ECONS_VAL   =  0.0;
}


// For self-consistent fitting type runs

void INPUT::PARSE_CONTROLS_SLFCNST()
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("SLFCNST", CONTENTS(i)))						
		{
			if ( RANK == 0 ) 
				cout << "NOTE: Feature \"SLFCNST\" is no longer supported. Ignoring." << endl;
			
			break;
		}
	}
}

// "General control variables"

void INPUT::PARSE_CONTROLS_RNDSEED(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("RNDSEED", CONTENTS(i)))								
		{
			CONTROLS.SEED = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# RNDSEED #: " << CONTROLS.SEED << endl;
			
			break;
		}
	}	
}
void INPUT::PARSE_CONTROLS_TEMPERA(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("TEMPERA", CONTENTS(i)))										
		{
			CONTROLS.TEMPERATURE = convert_double(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# TEMPERA #: " << CONTROLS.TEMPERATURE << " K" << endl;	
			
			break;
		}
	}	
}
void INPUT::PARSE_CONTROLS_PRESSUR(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRESSUR", CONTENTS(i)))										
		{
			CONTROLS.PRESSURE = convert_double(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# PRESSUR #: " << CONTROLS.PRESSURE << " GPa" << endl;	
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_CONVCUT(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CONVCUT", CONTENTS(i)))										
		{
			CONTROLS.NVT_CONV_CUT = convert_double(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# CONVCUT #: " << CONTROLS.NVT_CONV_CUT*100 << " % of set T" << endl;	
			
			break;
		}
	}
}

void INPUT::PARSE_CONTROLS_CMPRFRC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CMPRFRC", CONTENTS(i)))										
		{
			CONTROLS.COMPARE_FORCE = convert_bool(CONTENTS(i+1,0),i+1);
			
			// // Read in the name of the file
			
			if((CONTROLS.COMPARE_FORCE) && (CONTENTS.size(i+1) == 2))
				CONTROLS.COMPARE_FILE = CONTENTS(i+1,1);
			else if(CONTROLS.COMPARE_FORCE)
				EXIT_MSG("ERROR: If # CMPRFRC # is true, a force file must be specified on the same line.");

			if ( (RANK == 0)  && (CONTROLS.COMPARE_FORCE)) 
			{
				cout << "	# CMPRFRC #: true... will only do 1 md step." << endl;
				cout << "                comparing to file: " << CONTROLS.COMPARE_FILE << endl;	

				CONTROLS.N_MD_STEPS = 1;
			}
			else if( (RANK == 0) && (!CONTROLS.COMPARE_FORCE))
			{
				cout << "	# CMPRFRC #: false" << endl;	
			}
			
			break;
		}
	}
}

void INPUT::PARSE_CONTROLS_SKIP_FRAMES(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("SKPFRMS", CONTENTS(i)))										
		{
			CONTROLS.SKIP_FRAMES = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( (RANK == 0)  && (CONTROLS.SKIP_FRAMES >= 1 )) 
			{
				cout << "	# SKPFRMS #: " << CONTROLS.SKIP_FRAMES <<
					" .. will skip through frames in parallel execution (round-robin)." << endl;
			}
			else if( (RANK == 0) && CONTROLS.SKIP_FRAMES <= 0 )
			{
				cout << "	# SKPFRMS #: "<< CONTROLS.SKIP_FRAMES <<
					" ... will process all frames in contiguous order" << endl;	
			}
			break;
		}
	}
}

void INPUT::PARSE_CONTROLS_CHCKFRC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CHECKFRC", CONTENTS(i)))
		{
			CONTROLS.CHECK_FORCE =  convert_bool(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				 cout << "\t# CHECKFRC #: " << bool2str(CONTROLS.CHECK_FORCE) << endl;
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_SUBTFRC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("SUBTFRC", CONTENTS(i)))		
		{
			CONTROLS.SUBTRACT_FORCE = convert_bool(CONTENTS(i+1,0),i+1);
			
			// Read in the name of the file
			
			if(CONTENTS.size(i+1) == 2)
				CONTROLS.COMPARE_FILE = CONTENTS(i+1,1);
			else
				EXIT_MSG("ERROR: If # SUBTFRC # is true, a subtranction file must be specified on the same line.");

			// Determine if we are subtracting off forces and/or energies
			
			for(int j=2; j<CONTENTS.size(i+1); j++)
			{
				if (CONTENTS(i+1,j) == "SUBTR_ENERGY")
					CONTROLS.FIT_ENER = true;
				else if(CONTENTS(i+1,j) == "SUBTR_STRESS")
					CONTROLS.FIT_STRESS = true;
				else
					EXIT_MSG("ERROR: Unrecognized SUBTFRC option... Allowed values are SUBTR_ENERGY and SUBTR_STRESS: ", CONTENTS(i+1,2));
			}
			
			
			if ( RANK == 0  && CONTROLS.SUBTRACT_FORCE) 
			{
				cout << "	# CMPRFRC #: true... will only do 1 md step." << endl;
				cout << "		...comparing to file: " << CONTROLS.COMPARE_FILE << endl;	
				
				if(CONTROLS.FIT_ENER)
					cout << "		...Will subtract energies." << endl;	
				if(CONTROLS.FIT_STRESS)
					cout << "		...Will subtract stress tensors." << endl;					

				CONTROLS.N_MD_STEPS = 1;
			}
			else if( (RANK == 0) && (!CONTROLS.SUBTRACT_FORCE))
			{
				cout << "	# CMPRFRC #: false" << endl;	
			}
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_TIMESTP(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("TIMESTP", CONTENTS(i)))		
		{
			CONTROLS.DELTA_T_FS = convert_double(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# TIMESTP #: " << CONTROLS.DELTA_T_FS << " fs" << endl;	
			
			CONTROLS.DELTA_T = CONTROLS.DELTA_T_FS/Tfs;	// Now convert timestep to simulation units
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_N_MDSTP(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("N_MDSTP", CONTENTS(i)))		
		{
			CONTROLS.N_MD_STEPS = convert_int(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# N_MDSTP #: " << CONTROLS.N_MD_STEPS << endl;	
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_PENTHRS(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PENTHRS", CONTENTS(i)))		
		{
			CONTROLS.PENALTY_THRESH = convert_double(CONTENTS(i+1,0),i+1);
			
			if ( RANK == 0 ) 
				cout << "	# PENTHRS #: " << CONTROLS.PENALTY_THRESH << endl;	
			
			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_KILLLEN(JOB_CONTROL & CONTROLS)
{
        int N_CONTENTS = CONTENTS.size();

        for (int i=0; i<N_CONTENTS; i++)
        {
                if (found_input_keyword("USEKILL", CONTENTS(i)))
                {
                        CONTROLS.USE_KILL_LEN = convert_bool(CONTENTS(i+1,0),i+1);

                        if ( RANK == 0 )
                                cout << "        # USEKILL #: " << CONTROLS.USE_KILL_LEN << endl;

                        break;
                }
        }
}
void INPUT::PARSE_CONTROLS_NLAYERS(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("NLAYERS", CONTENTS(i)))		
		{
			CONTROLS.N_LAYERS = convert_int(CONTENTS(i+1,0),i+1);
			
			if (RANK==0)
				cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;	
			
			if (CONTENTS.size(i+1) == 3)
			{
				if(CONTENTS(i+1,1) == "REPLICATE")
					CONTROLS.REAL_REPLICATES = convert_int(CONTENTS(i+1,2),i+1);
				else
					EXIT_MSG("ERROR: Unrecognized NLAYERS command: ",CONTENTS(i+1,1));
				
				if (RANK==0)
					cout << "	             ... Creating " << CONTROLS.REAL_REPLICATES << " real replicates before ghost atom layering." << endl;	
			}

			break;
		}
	}
}
void INPUT::PARSE_CONTROLS_USENEIG(JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("USENEIG", CONTENTS(i))	 )
		{		
			NEIGHBOR_LIST.USE = convert_bool(CONTENTS(i+1,0),i+1);
		
			if (RANK==0)
				cout << "	# USENEIG #: " << bool2str(NEIGHBOR_LIST.USE) << endl;
				
			if (CONTENTS.size(i+1) == 2 && NEIGHBOR_LIST.USE )
			{
				if (CONTENTS(i+1,1) == "SMALL")
				{
					NEIGHBOR_LIST.UPDATE_WITH_BIG = false;
					if ( RANK == 0 )
						 cout << "		Will update the neighbor list through the \"small\" method " << endl;
				}
				else
				{
					EXIT_MSG("ERROR: Unrecognized # USENEIG # option: ", CONTENTS(i+1,1));
				}
			}
		
			break;
		}
	}	
}

void INPUT::PARSE_CONTROLS_PRMFILE(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRMFILE", CONTENTS(i)))		
		{
			CONTROLS.PARAM_FILE = CONTENTS(i+1,0);
		
			if (RANK==0)
				cout << "	# PRMFILE #: " << CONTROLS.PARAM_FILE << endl;	
		
			break;
		}
	}	
}


void INPUT::PARSE_CONTROLS_SERIAL_CHIMES(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();

	int i = 0 ;
	for ( ; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("SERIAL_CHIMES", CONTENTS(i)))		
		{
			 CONTROLS.SERIAL_CHIMES = convert_bool(CONTENTS(i+1,0),i+1);

			 if ( RANK == 0 )
			 {
				 cout << "	# SERIAL_CHIMES #: " ;
				 if ( CONTROLS.SERIAL_CHIMES ) {
					 cout << "TRUE" << endl ;
				 } else {
					 cout << "FALSE" << endl ;
				 }
			 }
			 break ;					 ;
		}
	}

}

void INPUT::PARSE_CONTROLS_CRDFILE(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	bool found = false;
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CRDFILE", CONTENTS(i)))		
		{
			found = true;
			
			if(CONTENTS(i+1,0) == "CAT") // Then we'll be catenating several files together
			{
				cout << "	# CRDFILE #: Creating a cell from multiple input files: ";
				int NO_FILES = convert_int(CONTENTS(i+1,2),i+1);
				
				cout << NO_FILES  << endl;
				
				for(int j=0; j<NO_FILES; j++)
				{
					CONTROLS.COORD_FILE.push_back(CONTENTS(i+1,3+j));
					
					if(RANK == 0)
						cout << "		" << CONTROLS.COORD_FILE[3+j] << endl;
				}	
				
				if(RANK == 0)
					cout << "	Note: Assumes files have the same x and y box dimensions!" << endl;			
			}
			if(CONTENTS(i+1,0) == "SCALE") // Then coordinates/boxlengths will be scaled
			{
				CONTROLS.SCALE_SYSTEM_BY = convert_double(CONTENTS(i+1,1),i+1);

				CONTROLS.COORD_FILE.push_back(CONTENTS(i+1,2));
				
				if (RANK==0)
					cout << "	# CRDFILE #: " << CONTROLS.COORD_FILE[0] << endl;
			}
			if(CONTENTS(i+1,0) == "INITIALIZE")
			{
				CONTROLS.BUILD = true;
				
				CONTROLS.BUILD_TYPE = CONTENTS(i+1,1);
				
				if (CONTROLS.BUILD_TYPE == "MOLECULAR")
					CONTROLS.BUILD_FILE = CONTENTS(i+1,2);
					
				else if (CONTROLS.BUILD_TYPE == "ATOMIC")
					CONTROLS.BUILD_ATOM = CONTENTS(i+1,2);

				if(CONTENTS(i+1,3) != "BOXL")
					EXIT_MSG("ERROR: Expected \'BOXL <value> NMOLEC <value>.");
					
				CONTROLS.BUILD_BOXL = convert_double(CONTENTS(i+1,4),i+1);

				if(CONTENTS(i+1,5) != "NMOLEC")
					EXIT_MSG("ERROR: Expected \'BOXL <value> NMOLEC <value>.");

				CONTROLS.BUILD_NMOLEC = convert_int(CONTENTS(i+1,6),i+1);	
				
			}
			else
			{
				CONTROLS.COORD_FILE.push_back(CONTENTS(i+1,0));
				
				if (RANK==0)
					cout << "	# CRDFILE #: " << CONTROLS.COORD_FILE[0] << endl;	
			}

			break;
		}
	}
	
	if(!found)
		EXIT_MSG("ERROR: # CRDFILE # must be specified!");
}

// " Simulation options"

void INPUT::PARSE_CONTROLS_VELINIT(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	bool found = false;
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("VELINIT", CONTENTS(i)))		
		{
			found = true;
			
			if (CONTENTS(i+1,0)=="READ")
				CONTROLS.INIT_VEL = false;
			
			else if (CONTENTS(i+1,0)=="GEN")
				CONTROLS.INIT_VEL = true;
			
			else if ( CONTENTS(i+1,0) == "RESTART" ) 
				CONTROLS.RESTART = true ;
			
			else
			{
				cout << "ERROR: # VELINIT # must be specified as READ or GEN." << endl;
				exit_run(1);	
			}	
			
			if (RANK==0)
			{
				if(!CONTROLS.INIT_VEL && CONTROLS.BUILD)
					EXIT_MSG("ERROR: # VELINIT # must be \'GEN\' when # CRDFILE # \'INITIALIZE\' option is used.");
	
				if(CONTROLS.INIT_VEL)
					cout << "	# VELINIT #: GEN ... generating velocites via box Muller" << endl;	
				else if(CONTROLS.RESTART)
					cout << "	# VELINIT #: RESTART ... reading velocities from coordinate (restart) file" << endl;	
				else if(!CONTROLS.INIT_VEL && !CONTROLS.RESTART)	
					cout << "	# VELINIT #: READ ... reading velocities from coordinate file" << endl;	
			}
		
			break;
		}
	}
	
	if(!found)
		EXIT_MSG("ERROR: # VELINIT # must be specified!");	
}
void INPUT::PARSE_CONTROLS_CONSRNT(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	bool found = false;
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("CONSRNT", CONTENTS(i)))		
		{
			found = true;
			
			CONTROLS.ENSEMBLE = CONTENTS(i+1,0);
			
			// Note: don't need to do ensemble checking here - CONSTRAINT initialization takes care of that. 

			/*

			RETHINK HOW THIS SECTION IS HANDLED! ... HAVE CONTROLS JUST READ IN THE STRING
			AND HAVE CONSTRAINT DO THE PARSING/STORAGE
			
			*/
			

			if(CONTROLS.ENSEMBLE != "NVE")
			{
				if (CONTROLS.ENSEMBLE == "NVT-MTK" || CONTROLS.ENSEMBLE == "NPT-MTK")
				{
					CONTROLS.USE_HOOVER_THRMOSTAT   = true;
					CONTROLS.FREQ_UPDATE_THERMOSTAT = convert_double(CONTENTS(i+1,2),i+1);
					
					if (RANK==0)
						cout << "	# CONSRNT #: HOOVER... a Hoover time of " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " will be used." << endl;

					if(CONTROLS.ENSEMBLE == "NPT-MTK")
					{
						if(CONTENTS.size(i+1) == 4)
							CONTROLS.FREQ_UPDATE_BAROSTAT = convert_double(CONTENTS(i+1,3),i+1);
						else
							CONTROLS.FREQ_UPDATE_BAROSTAT = 1000;	
							
						if (RANK==0)
							cout << "	                   ... and barostat will use  " << CONTROLS.FREQ_UPDATE_BAROSTAT << "." << endl;	
					}
				}
				else if (   CONTROLS.ENSEMBLE == "NVT-SCALE"
								 || CONTROLS.ENSEMBLE == "NPT-BEREND"
								 || CONTROLS.ENSEMBLE == "NVT-BEREND"
								 || CONTROLS.ENSEMBLE == "NPT-BEREND-ANISO")
				{
					CONTROLS.USE_HOOVER_THRMOSTAT   = false;
					CONTROLS.FREQ_UPDATE_THERMOSTAT = convert_double(CONTENTS(i+1,1),i+1);

					if (CONTROLS.ENSEMBLE == "NPT-BEREND" || CONTROLS.ENSEMBLE == "NPT-BEREND-ANISO")
					{
						if(CONTENTS.size(i+1) == 3)
							CONTROLS.FREQ_UPDATE_BAROSTAT = convert_double(CONTENTS(i+1,2),i+1);
						else
							CONTROLS.FREQ_UPDATE_BAROSTAT = 1000;		
					}
		
					if (RANK==0 && CONTROLS.ENSEMBLE == "NVT-SCALE" ) 
						 cout << "	# CONSRNT #: " << CONTROLS.ENSEMBLE << "... Velocities will be scaled every " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " MD steps." ;
					else if ( RANK == 0 ) 
						 cout << "	# CONSRNT #: " << CONTROLS.ENSEMBLE << " Velocity scaling time = " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " fs " ;

					if ( RANK == 0 && (CONTROLS.ENSEMBLE == "NPT-BEREND" || CONTROLS.ENSEMBLE == "NPT-BEREND-ANISO") ) 
						 cout << " Volume scaling time = " << CONTROLS.FREQ_UPDATE_BAROSTAT << " fs " ;

					if ( RANK == 0 ) 
						 cout << endl ;
								
				}
				else
				{
					cout << "ERROR: Unrecognized # CONSRNT # command: " << join_string_vec(CONTENTS(i+1),' ') << endl;
					exit_run(1);	
				}	
			}
			
			for (int j=0; j<CONTENTS.size(i+1); j++)
			{
				if (CONTENTS(i+1,j) == "FREEZE")
				{
					CONTROLS.FREEZE_IDX_START = convert_int(CONTENTS(i+1,j+1),i+1);
					CONTROLS.FREEZE_IDX_STOP  = convert_int(CONTENTS(i+1,j+2),i+1);
					
					if(RANK == 0)
					{
						cout << "		...Atoms over the following range will be frozen will be frozen (indexed from 0): "; 
						cout << CONTROLS.FREEZE_IDX_START << " - "  << CONTROLS.FREEZE_IDX_STOP << endl;
					}
					
					break;
					
				}
			}
			

			break;
		}
	}
	
	if(!found)
		EXIT_MSG("ERROR: # CONSRNT # must be specified!");
}

void INPUT::PARSE_CONTROLS_PRSCALC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRSCALC", CONTENTS(i)))		
		{
			if(CONTENTS(i+1,0) == "ANALYTICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = false;
			}
			else if(CONTENTS(i+1,0) == "NUMERICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = true;
			}
			else
			{
				EXIT_MSG("ERROR: # PRSCALC # must be specified as ANALYTICAL or NUMERICAL.");
			}
			
			break;
		}
	}
	
	if(CONTROLS.USE_NUMERICAL_PRESS)
	{
		if (RANK==0)
			cout << "	# PRSCALC #: NUMERICAL" << endl;
	}
	else
	{
		if (RANK==0)
			cout << "	# PRSCALC #: ANALYTICAL" << endl;
	}
		
}

void INPUT::PARSE_CONTROLS_STRSCALC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("STRSCALC", CONTENTS(i)))		
		{
			if(CONTENTS(i+1,0) == "ANALYTICAL")
			{
				CONTROLS.USE_NUMERICAL_STRESS = false;
				if (RANK==0)
					cout << "	# STRSCALC #: ANALYTICAL" << endl;
			}
			else if(CONTENTS(i+1,0) == "NUMERICAL")
			{
				CONTROLS.USE_NUMERICAL_STRESS = true;
				if (RANK==0)
					cout << "	# STRSCALC #: NUMERICAL" << endl;
				
			}
			else
			{
				EXIT_MSG("ERROR: # STRSCALC # must be specified as ANALYTICAL or NUMERICAL.");
			}
			
			break;
		}
	}
}

void INPUT::PARSE_CONTROLS_WRPCRDS(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("WRPCRDS", CONTENTS(i)))		
		{
			CONTROLS.WRAP_COORDS = convert_bool(CONTENTS(i+1,0),i+1);		
			break;
		}
	}	
	
	if (RANK==0)
		cout << "	# WRPCRDS #: " << bool2str(CONTROLS.WRAP_COORDS) << endl;	
}

// "Output control"

void INPUT::PARSE_CONTROLS_ATMENER(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	CONTROLS.INCLUDE_ATOM_OFFSETS = true ;
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("ATMENER", CONTENTS(i)))
		{
			CONTROLS.INCLUDE_ATOM_OFFSETS = convert_bool(CONTENTS(i+1,0),i+1);
			break;
		}
	}
	
	if (RANK==0)
		cout << "	# ATMENER #: " << bool2str(CONTROLS.INCLUDE_ATOM_OFFSETS) << endl;
	
	if ( CONTROLS.INCLUDE_ATOM_OFFSETS && RANK==0)
		cout << "		... All reported energies include single atom contributions if available" << endl;
	if (!CONTROLS.INCLUDE_ATOM_OFFSETS && RANK==0)
		cout << "		... All reported energies *DO NOT* include single atom contributions" << endl;
	
}

void INPUT::PARSE_CONTROLS_FRQDFTB(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FRQDFTB", CONTENTS(i)) || found_input_keyword("FRQTRAJ", CONTENTS(i)))
		{
			CONTROLS.FREQ_DFTB_GEN = convert_int(CONTENTS(i+1,0),i+1);
			break;
		}
	}	
	
	if (RANK==0)
		cout << "	# FRQTRAJ #: " << CONTROLS.FREQ_DFTB_GEN << endl;
	
	if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0 && RANK==0)
	{
		cout << "		... printing every " << CONTROLS.FREQ_DFTB_GEN*CONTROLS.DELTA_T_FS << " fs, " << endl;
		cout << "		... printing       " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_DFTB_GEN << " frames. " << endl; 
	}
			
}
void INPUT::PARSE_CONTROLS_TRAJEXT(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("TRAJEXT", CONTENTS(i)))		
		{
			CONTROLS.TRAJ_FORMAT = CONTENTS(i+1,0);		
			break;
		}
	}	
	if (RANK==0)
		cout << "	# TRAJEXT #: " << CONTROLS.TRAJ_FORMAT << endl;		
}
void INPUT::PARSE_CONTROLS_FRQENER(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FRQENER", CONTENTS(i)))		
		{
			CONTROLS.FREQ_ENER = convert_int(CONTENTS(i+1,0),i+1);		
			break;
		}
	}	
	if (RANK==0)
		cout << "	# FRQENER #: " << CONTROLS.FREQ_ENER << endl;	
}
void INPUT::PARSE_CONTROLS_PRNTFRC(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRNTFRC", CONTENTS(i)))		
		{
			CONTROLS.PRINT_FORCE = convert_bool(CONTENTS(i+1,0),i+1);
			
			if(CONTENTS.size(i+1) >= 2 )
			{
				if(CONTENTS(i+1,1) == "FRQDFTB")
				{
					CONTROLS.FREQ_FORCE = CONTROLS.FREQ_DFTB_GEN;
				}
				else if ( CONTENTS(i+1,1) == "ENERGY_STRESS" ) {
					CONTROLS.PRINT_ENERGY_STRESS = true ;
					if ( CONTENTS.size(i+1) == 3 ) 
						CONTROLS.FREQ_FORCE = convert_int(CONTENTS(i+1,2),i+1);
				} 
				else
				{
					CONTROLS.FREQ_FORCE = convert_int(CONTENTS(i+1,1),i+1);
				}
			} else {
				 // Default: print every step if PRNTFRC requested.
				 CONTROLS.FREQ_FORCE = 1 ;
			} 
		
			break;
		}
	}	
	if (RANK==0)
	{
		cout << "	# PRNTFRC #: " << bool2str(CONTROLS.PRINT_FORCE) << endl;
		
		if(CONTROLS.PRINT_FORCE)
		{
			cout << "		... and will be printed every " << CONTROLS.FREQ_FORCE << " frames."<< endl;	
			if ( CONTROLS.PRINT_ENERGY_STRESS )
				cout << "       ... with energy/stress header before forces" << endl ;
			
			if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0)
			{
				cout << "		... printing every " << CONTROLS.FREQ_FORCE*CONTROLS.DELTA_T_FS << " fs, " << endl;
				cout << "		... printing " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_FORCE << " frames. " << endl; 
			}
		}
	}	
}
void INPUT::PARSE_CONTROLS_PRNTBAD(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRNTBAD", CONTENTS(i)))		
		{
			CONTROLS.PRINT_BAD_CFGS = convert_bool(CONTENTS(i+1,0),i+1);
			break;
		}
	}	
	if (RANK==0)
		cout << "	# PRNTBAD #: " << bool2str(CONTROLS.PRINT_BAD_CFGS) << endl;	
}
void INPUT::PARSE_CONTROLS_FRQRSTR(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FRQRSTR", CONTENTS(i)))		
		{
			CONTROLS.FREQ_BACKUP = convert_int(CONTENTS(i+1,0),i+1);
			break;
		}
	}	
	if (RANK==0)
		cout << "	# FRQRSTR #: " << CONTROLS.FREQ_BACKUP << endl;	
}
void INPUT::PARSE_CONTROLS_PRNTVEL(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("PRNTVEL", CONTENTS(i)))		
		{
			CONTROLS.PRINT_VELOC = convert_bool(CONTENTS(i+1,0),i+1);
			
			if(CONTENTS.size(i+1) == 2 )
			{
				if(CONTENTS(i+1,1) == "FRQDFTB")
					CONTROLS.FREQ_VELOC = CONTROLS.FREQ_DFTB_GEN;
				else
					CONTROLS.FREQ_VELOC = convert_int(CONTENTS(i+1,1),i+1);
			}
			
			break;
		}
	}	
	if (RANK==0)
	{
		cout << "	# PRNTFRC #: " << bool2str(CONTROLS.PRINT_FORCE) << endl;
		
		if(CONTROLS.PRINT_FORCE)
		{
			cout << "		... and will be printed every " << CONTROLS.FREQ_FORCE << " frames."<< endl;	
		
			if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0)
			{
				cout << "		... printing every " << CONTROLS.FREQ_FORCE*CONTROLS.DELTA_T_FS << " fs, " << endl;
				cout << "		... printing " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_FORCE << " frames. " << endl; 
			}
		}
	}			
}
void INPUT::PARSE_CONTROLS_GETSTRS(JOB_CONTROL & CONTROLS)
{
	// Shared with LSQ ... Used for comparing fit to force field
	PARSE_CONTROLS_FITSTRS(CONTROLS);
}

void INPUT::PARSE_CONTROLS_GETENER(JOB_CONTROL & CONTROLS)
{
	// Shared with LSQ ... Used for comparing fit to force field
	PARSE_CONTROLS_FITENER(CONTROLS);
}

void INPUT::PARSE_CONTROLS_FORDFTB(JOB_CONTROL & CONTROLS)
{
	int N_CONTENTS = CONTENTS.size();
	
	for (int i=0; i<N_CONTENTS; i++)
	{
		if (found_input_keyword("FORDFTB", CONTENTS(i)))		
		{
			CONTROLS.FORDFTB = convert_bool(CONTENTS(i+1,0),i+1);
			break;
		}
	}	
	if (RANK==0)
		cout << "	# FORDFTB #: " << bool2str(CONTROLS.FORDFTB) << endl;
}

// Run MD sanity checks

void INPUT::RUN_SANITY_MD(JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)
{
	if( (CONTROLS.N_LAYERS>0) && (!NEIGHBOR_LIST.USE) )
		if (RANK == 0)
			EXIT_MSG("Neighbor lists must be used if NLAYERS > 0");
			
	if( (CONTROLS.N_LAYERS==0) && (NEIGHBOR_LIST.USE) )
			
			cout << "WARNING: Use of neighbor lists HIGHLY reccommended when NLAYERS > 0!" << endl;
	return;
}

double INPUT::convert_double(const string &str, int idx)
// Convert a string to a double with user-friendly error checking.
// Use as a replacement for stod.
// idx specifies the first index of the CONTENTS array
{
	int pos = str.find_first_not_of(" \t\n") ;
	int posdigit = str.find_first_of("0123456789") ;
		
	if ( isalpha(str[pos]) || posdigit == string::npos )
	{
		if ( RANK == 0 )
			cout << "String found where floating point expected: " + str << endl ;
		convert_error(idx) ;
	} else if ( str[pos] == '+' || str[pos] == '-' || str[pos] == '.' || isdigit(str[pos]) )
	{
		try
		{
			return stod(str) ;
		}
		catch (...)
		{
			;
		}
	}
	if ( RANK == 0 ) {
		cout << "String found where floating point number expected: " + str << endl ;
		convert_error(idx) ;
	}
		
	// Not reached.
	return 0.0 ;
}


int INPUT::convert_int(const string &str, int idx)
// Convert a string to an integer with user-friendly error checking.
// Use as a replacement for stoi.
// idx specifies the first index of the CONTENTS array
{
	int pos = str.find_first_not_of(" \t\n") ;
	int posdigit = str.find_first_of("0123456789") ;
	if ( isalpha(str[pos]) || posdigit == string::npos ) {
		if ( RANK == 0 )
			cout << "String found where integer expected: " + str << endl ;
		convert_error(idx) ;
	} else if ( str[pos] == '+' || str[pos] == '-' || isdigit(str[pos]) ) {
		try {
			return stoi(str) ;
		}
		catch (...)
		{
			;
		}
	}
	if ( RANK == 0 ) 
		cout << "String found where integer expected: " + str << endl ;
	convert_error(idx) ;

	// Not reached.
	return 0 ;
}


bool INPUT::convert_bool(const string &str, int idx)
// Convert a string to a boolean with user-friendly error checking.
// Use as a replacement for str2bool
// idx specifies the first index of the CONTENTS array
{
	int pos = str.find_first_not_of(" \t\n") ;
	if ( str[pos] == 't' || str[pos] == 'T' || str[pos] == 'F' || str[pos] == 'f' ) {
		try {
			return str2bool(str) ;
		}
		catch ( ... ) {
			;
		}
	}
	if ( RANK == 0 ) 
		cout << "String found where true/false expected: " + str << endl ;
	convert_error(idx) ;

	// Not reached.
	return false ;
}

void INPUT::convert_error(int idx)
// Call when an error occurs in converting input to a data type.
{
		if ( RANK == 0 )
		{
			cout << "Working on input: " ;
			for ( int j = 0 ; j < CONTENTS.size(idx) ; j++ ) {
				cout << CONTENTS(idx,j) << " " ;
			}
			cout << endl ;
		}
		exit_run(1) ;
}
		
		
		
		
		
