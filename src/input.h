#ifndef _INPUT_H
#define _INPUT_H

#include <vector>
#include <string>

using namespace std;

#include "util.h"
#include "A_Matrix.h"


class InputListing 
{
	// Provides safe access to input lines and tokens.

	vector< vector< string> > list;	// Entire contents of the infile - [lines]["words" ("tokens")]

	public:
	
	string operator()(int i, int j);	// Access a particular token in the InputListing.	
	vector<string> &operator()(int i);	// Access a vector of tokens for a particular line.
	void push_back(vector<string> &item);	// Add a vector of tokens to the InputListing.
  
	int size();				// Return the number of lines.
	int size(int i);			// Return the number of tokens for a particular line.	
};

class INPUT
{
	// Reads and stores the ENTIRE input file (fm_setup.in or run_md.in)
	// Then sets input options in a specified order
	// Removes the need for ordered input
	
public: // public functions
		
		// Constructor/Deconstructor
	
	INPUT(string FNAME);
	~INPUT();
		
	// Wrapper functions: Process entire input
	
	void PARSE_INFILE_LSQ(	JOB_CONTROL               & CONTROLS, 			
							vector<PAIRS>             & ATOM_PAIRS, 
							CLUSTER_LIST              & TRIPS, 
							CLUSTER_LIST              & QUADS, 
							map<string,int>           & PAIR_MAP,
							vector<int>               & INT_PAIR_MAP,
							vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
							NEIGHBORS                 & NEIGHBOR_LIST,
							vector<int>               & TMP_ATOMTYPEIDX, 
							vector<string>            & TMP_ATOMTYPE,
							A_MAT			  & A_MATRIX);
							
	void PARSE_INFILE_MD (JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);		
		
private:	// private funcitons
	
	// General parsing helpers
	
	string		MODE;			// Are we parsing fm_setup.in or run_md.in (MODE == "LSQ" or MODE == "MD")
	string 		FILENAME;		// fm_setup.in or run_md.in
	ifstream 	INFILE;

	InputListing CONTENTS ;
	
	// LSQ parsing helpers
	
	int		NPAIR;
	int		NTRIP;
	int		NQUAD;
	PAIRS 		TEMP_PAIR;
	string 		FF_TYPE;				// Chebyshev, inverse, etc... Replaces "TEMP_TYPE" in old read_lsq_input function
	int 		TMP_CHEBY_RANGE_LOW;
	int 		TMP_CHEBY_RANGE_HIGH;
		
	bool found_input_keyword(string KEYWORD, vector<string> & LINE);
		
	// For storing the entire input file in memory and checkingf or unrecognized key words
	
	void READ_FILE();
	void CHECK_FILE();					// Check for unrecognized keywords


	// For assigning LSQ variables: "Control Variables" 
	
	void PARSE_CONTROLS_WRAPTRJ(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_TRJFILE(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_SPLITFI(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_HIERARC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_NFRAMES(JOB_CONTROL & CONTROLS);
	//void PARSE_CONTROLS_NLAYERS(JOB_CONTROL & CONTROLS); // JUST USE THE MD VERSION... IT SHOULD BE COMPATIBLE
	void PARSE_CONTROLS_FITCOUL(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FITSTRS(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FITENER(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PAIRTYP(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CHBTYPE(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CHEBYFIX(JOB_CONTROL & CONTROLS) ;
	void PARSE_CONTROLS_SERIAL_CHIMES(JOB_CONTROL & CONTROLS)	;
	void PARSE_CONTROLS_SKIP_FRAMES(JOB_CONTROL &CONTROLS) ;
	
	// For assigning LSQ variables: "Topology Variables" 
	
	void PARSE_TOPOLOGY_EXCLUDE(vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS, A_MAT & A_MATRIX, JOB_CONTROL & CONTROLS, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE);
	void PARSE_TOPOLOGY_NATMTYP(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS);
	void PARSE_TOPOLOGY_TYPEIDX(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE);
	void PARSE_TOPOLOGY_PAIRIDX(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST, vector<int> & TMP_ATOMTYPEIDX, vector<string> & TMP_ATOMTYPE);
	void PARSE_TOPOLOGY_CHGCONS(vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, map<string,int> & PAIR_MAP);
	void PARSE_TOPOLOGY_SPECMIN(CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS);
	void PARSE_TOPOLOGY_SPECMAX(CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS);
	void PARSE_TOPOLOGY_FCUTTYP(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST & TRIPS, CLUSTER_LIST & QUADS);
	
	// Run LSQ sanity checks
	
	void RUN_SANITY_LSQ(JOB_CONTROL & CONTROLS);
		
	
	// MD controls initialization
	
	void PARSE_CONTROLS_INITIAL(JOB_CONTROL & CONTROLS,NEIGHBORS & NEIGHBOR_LIST);
		
	// For self-consistent fitting type runs
	
	void PARSE_CONTROLS_SLFCNST();
	
	// "General control variables"
	
	void PARSE_CONTROLS_RNDSEED(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_TEMPERA(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PRESSUR(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CONVCUT(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CMPRFRC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CHCKFRC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_SUBTFRC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_TIMESTP(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_N_MDSTP(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PENTHRS(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_KILLLEN(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_NLAYERS(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_USENEIG(JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);
	void PARSE_CONTROLS_PRMFILE(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CRDFILE(JOB_CONTROL & CONTROLS);
	
	// " Simulation options"
	
	void PARSE_CONTROLS_VELINIT(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_CONSRNT(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PRSCALC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_STRSCALC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_WRPCRDS(JOB_CONTROL & CONTROLS);
	
	// "Output control"
	
	void PARSE_CONTROLS_ATMENER(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FRQDFTB(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_TRAJEXT(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FRQENER(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PRNTFRC(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FRQRSTR(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PRNTBAD(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_PRNTVEL(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_GETSTRS(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_GETENER(JOB_CONTROL & CONTROLS);
	void PARSE_CONTROLS_FORDFTB(JOB_CONTROL & CONTROLS);
	
	// Run MD sanity checks
	
	void RUN_SANITY_MD(JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);

	// Convert a string to a double with error checking.
	double convert_double(const string &str, int idx) ;

	// Convert a string to an int with error checking.
	int convert_int(const string &str, int idx) ;

	// Convert a string to a boolean with error checking.
	bool convert_bool(const string &str, int idx) ;
			
	// Call when an error occurs in converting input to a data type.
	void convert_error(int idx) ;

};

#endif
