#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cstdlib>   // std::atof, std::atoi
#include<cmath>
#include<algorithm> // std::sort
#include<utility>   // std::pair

using namespace std;

////////////////////////////////////////////////////////////
// Notes
////////////////////////////////////////////////////////////

/*

1. Compile with -std=c++11 - Necessary!

*/


////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////

void STR_TO_VEC_OF_STRS(string STR, vector<string> & STR_VEC) // Requires sstream library
{
	// Make sure the target vector is empty
	
	STR_VEC.empty();
	STR_VEC.clear();
	
	// Setup the "tokenizer" (parser) vars

	stringstream PARSER;
	PARSER.str("");
	PARSER.clear();	
	PARSER.str(STR);
	
	// Save "tokens" ("words") to target vector
	
	while(PARSER >> STR)
		STR_VEC.push_back(STR);	
}

int STR_TO_INT(string STR) // Requires cstdlib library
{
	return atoi(STR.c_str());
}

double STR_TO_DOUB(string STR) // Requires cstdlib library
{
	return double(atof(STR.c_str()));
}

bool IS_PAIR_LESS( const pair<int,int> & PAIR1, const pair<int,int> & PAIR2)
{
	return (PAIR1.second < PAIR2.second);
}

bool IS_PAIR_EQUAL( const pair<int,int> & PAIR1, const pair<int,int> & PAIR2)
{
	return (PAIR1.second == PAIR2.second);
}

bool IS_PAIR_VEC_EQUAL( vector<pair<int,int> > & PAIRVEC1, vector<pair<int,int> > & PAIRVEC2)
{

	if(PAIRVEC1.size() != PAIRVEC2.size())
		return false;
	
	for(int i=0;i<PAIRVEC1.size(); i++)
		if(PAIRVEC1[i].second != PAIRVEC2[i].second)
			return false;
			
	return true;
}



////////////////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////////////////

struct XYZ_LINE
{
	string TYPE;
	double X;
	double Y;
	double Z;
};

struct BONDED_PAIR
{
	int IDX1;
	int IDX2;
	int PAR1;
	int PAR2;
};

////////////////////////////////////////////////////////////
// Classes
////////////////////////////////////////////////////////////

class BOND_CRIT
{
	// Read a bond criteria file:
	// Below is a sample format for a C O system
	// Note: atom/pair order doesn't matter
	// C O			// Atom types in file
	// 3			// Number of atom pairs
	// C C 1.9		// Bond criteria (one line for each pair)
	// C O 1.8
	// O O 1.7
	
	
	public:
		
		int NATOM_TYPES;
		int NTYPE_PAIRS;
		
		vector<string>          TYPES; // Possible atom types
		vector<vector<string> > PAIR;  // [npairs][2]... [i][0] = atom type 1, [i][1] = atom type 2
		vector<double>          CRIT;  // [npairs] The corresponding cutoff distance
		
		BOND_CRIT(string CRIT_FILE);
		~BOND_CRIT();
		
		void READ_CRITERIA(string CRIT_FILE); // Called by constructor
		int GET_PAIRTYPE(string a1, string a2);
};

BOND_CRIT::BOND_CRIT(string CRIT_FILE) // Constructor
{
	READ_CRITERIA(CRIT_FILE);
}

BOND_CRIT::~BOND_CRIT(){} // Deconstructor

void BOND_CRIT::READ_CRITERIA(string CRIT_FILE)
{
	cout << endl;
	cout << "		Reading bond criteria from file: " << CRIT_FILE << endl;
	
	// Setup up the fstream
	
	ifstream INFILE;
	INFILE.open(CRIT_FILE.data());
	if (!INFILE.is_open())
	{
		cout << "ERROR: Cannot open bond criteria file: " << CRIT_FILE << endl;
		cout << "Exiting." << endl;
		exit(0);
	}
	
	// Read
	
	string         LINE_STR;
	vector<string> LINE_VEC;
	int            TEMP_INT;
	double         TEMP_DOUB; 
	
	// Figure out how many atom types there are 
	
	getline(INFILE,LINE_STR);
	STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);

	if(LINE_VEC.size() == 0)
	{
		cout << "ERROR: At least one atom type is needed." << endl;
		cout << "Exiting. " << endl;
		exit(0);
	}

	NATOM_TYPES = LINE_VEC.size();
	cout << "			Expecting atom types: " << LINE_STR << "( "  <<  NATOM_TYPES << ")" << endl;
	
	TYPES.resize(NATOM_TYPES);
	
	for (int i=0; i<NATOM_TYPES; i++)
		TYPES[i] = LINE_VEC[i];
	
	// Figure out how many pairs there are
	
	getline(INFILE,LINE_STR);
	STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);

	if(LINE_VEC.size() == 0)
	{
		cout << "ERROR: Number of pairs missing in bond criteria file: " << CRIT_FILE << endl;
		cout << "Exiting. " << endl;
		exit(0);
	}	
	
	
	TEMP_INT = STR_TO_INT(LINE_VEC[0]);
	
	NTYPE_PAIRS = NATOM_TYPES*(NATOM_TYPES+1)/2.0;
	
	if (TEMP_INT != NTYPE_PAIRS)
	{
		cout << "Calculated number of pairs doesn't match read number: " << endl;
		cout << "Calculated: " << NTYPE_PAIRS << endl;
		cout << "Read:       " << TEMP_INT << endl;
		cout << "Exiting.    " << endl;
		exit(0);
	}
	
	// Read the atom pairs and thier bond criteria

	PAIR.resize(NTYPE_PAIRS);
	CRIT.resize(NTYPE_PAIRS);
	
	for (int i=0; i<NTYPE_PAIRS; i++)
	{
		getline(INFILE,LINE_STR);
		STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);
	
		if(LINE_VEC.size() != 3)
		{
			cout << "ERROR: Expect exactly 3 items on bond criteria lines in : " << CRIT_FILE << endl;
			cout << "<atom chem 1> <atom chem 2> <bond criteria> " << endl;
			cout << "Read:    " << LINE_STR << endl;
			cout << "Exiting. " << endl;
			exit(0);
		}
		
		PAIR[i].resize(2);
		
		PAIR[i][0] =             LINE_VEC[0];
		PAIR[i][1] =             LINE_VEC[1];
		CRIT[i]    = STR_TO_DOUB(LINE_VEC[2]);
	}
	
	// Output results of read
	
	cout << "			Read the following critera: " << endl;
	
	for (int i=0; i<NTYPE_PAIRS; i++)
		cout << "				" << PAIR[i][0] << " " << PAIR[i][1] << " " << CRIT[i] << endl;
	cout << endl;
}

int BOND_CRIT::GET_PAIRTYPE(string a1, string a2)
{
	for (int i=0; i<NTYPE_PAIRS; i++)
	{
		if (  ((a1 == PAIR[i][0]) && (a2 == PAIR[i][1]))  
		   || ((a1 == PAIR[i][1]) && (a2 == PAIR[i][0]))  )
		
			   return CRIT[i];
	}
	
	// If we make it to here (without "returning"), it means we haven't found
	// a match for the two atom types. Throw an error.
	
	cout << "ERROR: Found an atom pair of unknown type: " << a1 << ", " << a2 << endl;
	cout << "Update bond criteria file accordingly." << endl;
	cout << "Exiting." << endl;
	exit(0);
}


class FRAME
{
	public:
	
	int              NATOMS;
	int              NATOMS_AND_GHOSTS;
	XYZ_LINE         BOXDIM;
	vector<XYZ_LINE> COORDS;
	vector<int>      PARENT;
	
	
	FRAME();
	~FRAME();
	
	void   READ(ifstream & TRAJ_FILE);
	void   REPLICATE();
	double GET_DIST(int ai, int aj, BOND_CRIT & CRITERIA);
};

FRAME::FRAME(){}  // Constructor
FRAME::~FRAME(){} // Deconstructor

void FRAME::READ(ifstream & TRAJ_FILE)
{
	string         LINE_STR;
	vector<string> LINE_VEC;
	int            TEMP_INT;
	double         TEMP_DOUB;
	
	// Determine number of atoms
	
	getline(TRAJ_FILE,LINE_STR);
	STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);

	NATOMS = STR_TO_INT(LINE_VEC[0]);
	COORDS.resize(NATOMS);
	PARENT.resize(NATOMS);
	
	// Determine the box lengths
	
	getline(TRAJ_FILE,LINE_STR);
	STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);

	if(LINE_VEC.size() < 3)
	{
		cout << "ERROR: at least 3 items (boxlengths)." << endl;
		cout << "Got:     " << LINE_STR << endl;
		cout << "Exiting. " << endl;
		exit(0);
	}
	
	BOXDIM.X = STR_TO_DOUB(LINE_VEC[0]);
	BOXDIM.Y = STR_TO_DOUB(LINE_VEC[1]);
	BOXDIM.Z = STR_TO_DOUB(LINE_VEC[2]);
	
	// Determine the coordinates 
	
	for (int i=0; i<NATOMS; i++)
	{
		getline(TRAJ_FILE,LINE_STR);
		STR_TO_VEC_OF_STRS(LINE_STR, LINE_VEC);
		
		if(LINE_VEC.size() < 4)
		{
			cout << "ERROR: at least 4 items (atom type and xyz coords)." << endl;
			cout << "Got:     " << LINE_STR << endl;
			cout << "Exiting. " << endl;
			exit(0);
		}
		
		COORDS[i].TYPE = LINE_VEC[0];
		COORDS[i].X    = STR_TO_DOUB(LINE_VEC[1]);
		COORDS[i].Y    = STR_TO_DOUB(LINE_VEC[2]);
		COORDS[i].Z    = STR_TO_DOUB(LINE_VEC[3]);
		PARENT[i]      = i;	
	}
	
}

void FRAME::REPLICATE()
{
	// nx = 0 image is the first NATOMS atoms
	
	XYZ_LINE TEMP;
	
	for(int i=0; i<NATOMS; i++)
	{
		for (int n1=-1; n1<=1; n1++)
		{
			for (int n2=-1; n2<=1; n2++)
			{
				for (int n3=-1; n3<=1; n3++)
				{
					if((n1 == 0)&&(n2==0)&&(n3==0))
						continue;
				
					TEMP.TYPE = COORDS[i].TYPE;
					TEMP.X    = COORDS[i].X + n1 * BOXDIM.X;
					TEMP.Y    = COORDS[i].Y + n2 * BOXDIM.Y;
					TEMP.Z    = COORDS[i].Z + n3 * BOXDIM.Z;
						
					COORDS.push_back(TEMP);
					PARENT.push_back(i);
				}
			}
		}
	}
	
	NATOMS_AND_GHOSTS = COORDS.size();
}

double FRAME::GET_DIST(int ai, int aj, BOND_CRIT & CRITERIA)
{
	// Determine the 
	
	// Get distance components 
	
	double DX = COORDS[ai].X - COORDS[aj].X;
	double DY = COORDS[ai].Y - COORDS[aj].Y;
	double DZ = COORDS[ai].Z - COORDS[aj].Z;
	
	// If the distance is within the bonding criteria, return it
	// otherwise, if MIC was used, return -1, i no MIC was needed, return 
	
	double DIST = sqrt(DX*DX + DY*DY + DZ*DZ);
	
	if (  DIST < CRITERIA.CRIT[CRITERIA.GET_PAIRTYPE(COORDS[ai].TYPE,COORDS[aj].TYPE)]  )
		return DIST;
	else
		return -1;	
}


class CLUSTER
{
	// Clusters have several attributes:
	// 0. Size (natoms)
	// 1. Desity (natoms/V(Rg), where Rg is the radius of gyration for the cluster)
	// 2. Composition (w/r/t available atom types; requires information from BOND_CRIT to compute)
	// 3. A list of atom indices
	// 4. A corresponding list of atom types
	// 5. Topology (a list of all directly bonded atoms)
	
	public:

	int              NATOMS;
	double           DENSITY;
	vector<string>   COMPOSITION_TYPES; // [possible atom types]
	vector<double>   COMPOSIION;        // Corresponding compsitions
	vector<int>      ATOM_INDICES;      // [members]
	vector<int>      PARENT_INDICES;
	vector<string>   ATOM_TYPES;        // [members]
	vector<int>      ATOM_TYPES_IDX;    // [members]
	vector<pair<int,int> > IDX_AND_PAR; // 1st element is overall atom index, second is the corresponding parent atom index
	
	CLUSTER();
	~CLUSTER();
	
	void SET_NATOMS     ();
	void SET_DENSITY    (FRAME & TRAJ_FRAME);
	void SET_COMPOSITION(BOND_CRIT & CRIT);	
	void SET_ATOM_TYPES (FRAME & TRAJ_FRAME);
};

CLUSTER::CLUSTER(){}
CLUSTER::~CLUSTER(){}

void CLUSTER::SET_NATOMS()
{
	NATOMS = IDX_AND_PAR.size();
}

void CLUSTER::SET_DENSITY(FRAME & TRAJ_FRAME)
{
	// Note: All operations performed on unwrapped coordinates
	// i.e. using indices stored in ATOM_INDICES
		
	// Calculate the cluster geometric center
	
	double GC_X = 0.0; 
	double GC_Y = 0.0;
	double GC_Z = 0.0;
	
	for (int i=0; i<NATOMS; i++)
	{
		int UNWRAPPED_IDX = IDX_AND_PAR[i].first;
		
		GC_X += TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X;
		GC_Y += TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y;
		GC_Z += TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z;
	}
	
	GC_X /= NATOMS;
	GC_Y /= NATOMS;
	GC_Z /= NATOMS;
	
	// Calculate the cluster radius of gyration
	
	double ROG = 0.0;
	
	for (int i=0; i<NATOMS; i++)
	{
		int UNWRAPPED_IDX = IDX_AND_PAR[i].first;
		
		double DX = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X - GC_X;
		double DY = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y - GC_Y;
		double DZ = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z - GC_Z;
		
		ROG += (DX*DX + DY*DY + DZ*DZ);
	}
	
	ROG = pow(ROG/NATOMS,0.5);
	
	double PI = 3.14159265359;
	
	//cout << "Calculated ROG: " << ROG << endl;
	
	if(ROG == 0.0)
		DENSITY = 0.0;
	else
		DENSITY = NATOMS / (4/3*PI*pow(ROG,3.0));
	
}

void CLUSTER::SET_COMPOSITION(BOND_CRIT & CRIT)
{
	COMPOSIION       .resize(CRIT.NATOM_TYPES);
	COMPOSITION_TYPES.resize(CRIT.NATOM_TYPES);
		
	vector<double> COUNTS(CRIT.NATOM_TYPES);
	
	
	for (int i=0; i<CRIT.NATOM_TYPES; i++)
	{
		COUNTS[i] = 0.0;
		COMPOSITION_TYPES[i] = CRIT.TYPES[i];
	}
	
	for (int i=0; i<NATOMS; i++)
	{
		for (int j=0; j<CRIT.NATOM_TYPES; j++)
		{
			if (ATOM_TYPES[i] == CRIT.TYPES[j])
			{
				ATOM_TYPES_IDX[i] = j;
				COUNTS[j] += 1.0;
			}
		}
	}
	
	for (int i=0; i<CRIT.NATOM_TYPES; i++)
		COMPOSIION[i] = COUNTS[i]/NATOMS;	
}

void CLUSTER::SET_ATOM_TYPES(FRAME & TRAJ_FRAME)
{
	ATOM_TYPES    .resize(NATOMS);
	ATOM_TYPES_IDX.resize(NATOMS); // Values get assigned in SET_COMPOSITION function
	
	for (int i=0; i<NATOMS; i++)
		ATOM_TYPES[i] = TRAJ_FRAME.COORDS[ATOM_INDICES[i]].TYPE;
		
}


////////////////////////////////////////////////////////////
// General cluster processing tools
////////////////////////////////////////////////////////////

void DO_CLUSTER(BOND_CRIT CRITERIA, FRAME & TRAJ_FRAME, vector<CLUSTER> & REMAINING_CLUSTERS)
{
	
	// Setup the distance list. Only distances within pair cutoffs are stored
	
	vector<BONDED_PAIR> BONDED_PAIRS;
	int                 N_BONDED_PAIRS;
	
	// Compute pair distances (MIC) and apply cutoffs 
	
	double      TMP_DIST;
	BONDED_PAIR TMP_BONDED_PAIR;
	
	for (int i=0; i<TRAJ_FRAME.NATOMS_AND_GHOSTS; i++)
	{
		for (int j=i+1; j<TRAJ_FRAME.NATOMS_AND_GHOSTS; j++)
		{
			if(TRAJ_FRAME.PARENT[i] == TRAJ_FRAME.PARENT[j])
				continue;
		
			TMP_DIST = TRAJ_FRAME.GET_DIST(i,j, CRITERIA);
			
			if(TMP_DIST >= 0)
				BONDED_PAIRS.push_back({i,j,TRAJ_FRAME.PARENT[i],TRAJ_FRAME.PARENT[j]});
		}
	}
	
	cout << "			...Bonded pairs established" << endl;
	
	int NBONDED_PAIRS = BONDED_PAIRS.size();
	
	cout << "				...Counted bonded pairs: " << NBONDED_PAIRS << endl;
	
	// Set up the cluster objects

	int CLU_MOM_1;			// Index of first atom in given atom's parent cluster
	int CLU_MOM_2;			// Index of first atom in given atom's parent cluster
		
	vector<CLUSTER> CLUSTERS(TRAJ_FRAME.NATOMS_AND_GHOSTS);	// Init list of clu's as set of natoms clu's w/ a single atom each	

	for (int i=0; i<TRAJ_FRAME.NATOMS_AND_GHOSTS; i++)
		CLUSTERS[i].ATOM_INDICES.push_back(i);
		
	// Determine clusters based on computed distances

	for (int i=0; i<NBONDED_PAIRS; i++)
	{
		// Overview:
		//
		// CLU_MOM can have 2 possible values: 
		//
		// (1)  CLU_MOM == BONDED_PAIRS[i].IDX1
		//         ... This means that frame  BONDED_PAIRS[i].IDX1 has not been assigned to another 
		//             cluster yet. Thus, we can add frame  BONDED_PAIRS[i].IDX2 to it.
		//
		// (2) CLU_MOM != DISTS[i].IDX1
		//         ... This means that frame  BONDED_PAIRS[i].IDX1 has been assigned to a different 
		//                cluster, so we need to add frame  BONDED_PAIRS[i].IDX2 to CLUSTERS[CLU_MOM] 
		//                instead of  BONDED_PAIRS[i].IDX1

		int A, B;

		A = BONDED_PAIRS[i].IDX1;        // Figure out which atom index the first atom is 
		B = CLUSTERS[A].ATOM_INDICES[0]; // First guess at atom 1's parent cluster

		while (A != B) // If A != B, then it has been assigned to a different cluster. Figure out which one that is
		{
			A = B;
			B = CLUSTERS[A].ATOM_INDICES[0];
		}

		CLU_MOM_1 = B; // We've now determined which cluster atom 1 belongs to

		// Repeat the process to figure out which clustrer atom 2 belongs to 

		A = BONDED_PAIRS[i].IDX2;
		B = CLUSTERS[A].ATOM_INDICES[0];// First guess at atom 2's parent cluster

		while (A != B)
		{
			A = B;
			B = CLUSTERS[A].ATOM_INDICES[0];
		}

		CLU_MOM_2 = B;					// We've now determined which cluster atom 2 belongs to

		// If atom 1 and 2 are not already in the same cluster, combine thier clusters

		if(CLU_MOM_1 != CLU_MOM_2)
		{
			CLUSTERS[CLU_MOM_1].ATOM_INDICES.insert(CLUSTERS[CLU_MOM_1].ATOM_INDICES.end(), CLUSTERS[CLU_MOM_2].ATOM_INDICES.begin(), CLUSTERS[CLU_MOM_2].ATOM_INDICES.end());
			CLUSTERS[CLU_MOM_2].ATOM_INDICES.resize(1);
			CLUSTERS[CLU_MOM_2].ATOM_INDICES[0] = CLU_MOM_1;
		}
	}

	int sanity = 0;
	int san    = 0;
	
	for(int i=0; i<CLUSTERS.size(); i++)
	{
		if (CLUSTERS[i].ATOM_INDICES[0] == i)
			sanity += CLUSTERS[i].ATOM_INDICES.size();
		san += CLUSTERS[i].ATOM_INDICES.size();
	}
	
	cout << "			...Redundant cluster list built" << endl;

	
	// Now we have a list of all clusters, even clusters comprised entirely of ghosts
	// Lets clean this up so only unique clusters with at least one real atom are left
	
	vector<bool> MARK_FOR_REMOVAL (CLUSTERS.size());
	vector<bool> DOES_PARENT_EXIST(TRAJ_FRAME.NATOMS);
	
	for(int i=0; i<TRAJ_FRAME.NATOMS; i++)
		DOES_PARENT_EXIST[i] = false;
	
	for(int i=0; i<TRAJ_FRAME.NATOMS_AND_GHOSTS; i++)
	{		
	
		//////////////////////////////
		// Initialize removal marker
		//////////////////////////////
		
		MARK_FOR_REMOVAL[i] = false;

		//////////////////////////////
		// Deletion case 0: The "Cluster" is just a pointer to a different cluster
		//////////////////////////////
		
		if(CLUSTERS[i].ATOM_INDICES[0] != i)
		{
			MARK_FOR_REMOVAL[i] = true;				
			continue;
		}
		
		//////////////////////////////
		// Deletion case 1: No real atoms exist in the cluster
		//////////////////////////////

		int MIN = CLUSTERS[i].ATOM_INDICES[0];
		
		for(int j=0; j<CLUSTERS[i].ATOM_INDICES.size(); j++)
		{
			if(CLUSTERS[i].ATOM_INDICES[j] < MIN)
				MIN = CLUSTERS[i].ATOM_INDICES[j];
		}

		if(MIN >= TRAJ_FRAME.NATOMS)
		{
			MARK_FOR_REMOVAL[i] = true;				
			continue;
		}
		
		//////////////////////////////
		// Deletion case 2: cluster is a replicate (parent atoms identical to existing cluster)
		//////////////////////////////
		
		// make a vector of pairs that map an atom index to its parent
		
		CLUSTERS[i].IDX_AND_PAR.resize(CLUSTERS[i].ATOM_INDICES.size());
		
		for(int j=0; j<CLUSTERS[i].ATOM_INDICES.size(); j++)
			CLUSTERS[i].IDX_AND_PAR[j] = make_pair(CLUSTERS[i].ATOM_INDICES[j],TRAJ_FRAME.PARENT[ CLUSTERS[i].ATOM_INDICES[j] ]); 
										
		// Sort that vector based on the parents, remove duplicates
		
		sort  (CLUSTERS[i].IDX_AND_PAR.begin(), CLUSTERS[i].IDX_AND_PAR.end(), IS_PAIR_LESS);
		
		// Remove duplicates based on parents (requires resizing)
						
		vector<pair<int,int> >::iterator it;
		it = unique(CLUSTERS[i].IDX_AND_PAR.begin(), CLUSTERS[i].IDX_AND_PAR.end(), IS_PAIR_EQUAL);
		CLUSTERS[i].IDX_AND_PAR.resize(distance(CLUSTERS[i].IDX_AND_PAR.begin(),it));
		
		for(int j=0; j<i; j++)
		{
			bool EQUIVALENT = IS_PAIR_VEC_EQUAL(CLUSTERS[i].IDX_AND_PAR, CLUSTERS[j].IDX_AND_PAR);
		
			if(EQUIVALENT && (!MARK_FOR_REMOVAL[j]))
			{
				MARK_FOR_REMOVAL[i] = true;
				break;
			}
		}					
	}
	
	sanity = 0;
	san    = 0;		
	
	//vector<CLUSTER> REMAINING_CLUSTERS;

	for(int i=0; i<TRAJ_FRAME.NATOMS_AND_GHOSTS; i++)	
	{
		if (!MARK_FOR_REMOVAL[i])
		{
			REMAINING_CLUSTERS.push_back(CLUSTERS[i]);
			sanity += CLUSTERS[i].IDX_AND_PAR.size();
			san ++;
		}
	}		
	
	cout << "			...Unique cluster list built." << endl;
	cout << "				...Counted unique clusters: " << san << endl;
	cout << "				...Total atoms in unique clusters: " << sanity << endl;
	
	// ... At this point, the only clusters left should be "real" and unique
	// the actual ATOM_INDICES will give the unwrapped coordinates, and 
	// PARENT_INDICES will give the wrapped (original) coordinates
	
	
	
}

void CALC_AND_PRINT_CLU_INFO(vector<CLUSTER> & CLUSTERS, FRAME & TRAJ_FRAME, BOND_CRIT & CRITERIA, int f)
{
	// Determine cluster characteristics:
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
		CLUSTERS[i].SET_NATOMS     ();				  
		CLUSTERS[i].SET_DENSITY    (TRAJ_FRAME);	  
		CLUSTERS[i].SET_ATOM_TYPES (TRAJ_FRAME);	  
		CLUSTERS[i].SET_COMPOSITION(CRITERIA);
	}
	
	// Print cluster characteristics
	
	int clu_count = 0;
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
		ofstream OUTFILE;
		string   OUTNAME = "frame." + to_string(f) + ".cluster." + to_string(clu_count) + ".stats";
		OUTFILE.open(OUTNAME.data());
		
		OUTFILE << "# NATOMS DENSITY COMPOSITION (";
		for (int j=0; j<CLUSTERS[i].COMPOSITION_TYPES.size(); j++)
			OUTFILE << " " << CLUSTERS[i].COMPOSITION_TYPES[j];
		OUTFILE << ")" << endl;
		
		OUTFILE << CLUSTERS[i].NATOMS  << " ";
		OUTFILE << CLUSTERS[i].DENSITY << " ";
		
		for (int j=0; j<CLUSTERS[i].COMPOSIION.size(); j++)
			OUTFILE << CLUSTERS[i].COMPOSIION[j] << " ";
		
		OUTFILE << endl;
		
		OUTFILE.close();
		clu_count++;
	}
}

void PRINT_WRAPPED_AND_UNWRAPPED_CLUSTERS(string KIND,vector<CLUSTER> & CLUSTERS, FRAME & TRAJ_FRAME, int f)
{
	// Print out cluster structures (unwrapped)
	
	int clu_count = 0;
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
			ofstream OUTFILE;
			string   OUTNAME = KIND + ".unwrapped.frame." + to_string(f) + ".cluster." + to_string(clu_count) + ".xyz";
			OUTFILE.open(OUTNAME.data());

			OUTFILE << CLUSTERS[i].NATOMS << endl;
			OUTFILE << TRAJ_FRAME.BOXDIM.X << " " <<  TRAJ_FRAME.BOXDIM.Y << " " << TRAJ_FRAME.BOXDIM.Z << endl;

			for(int j=0; j<CLUSTERS[i].NATOMS; j++)
			{
				int UNWRAPPED_IDX = CLUSTERS[i].ATOM_INDICES[j];
				
				
				OUTFILE << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE;
				OUTFILE << " " << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X;
				OUTFILE << " " << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y;
				OUTFILE << " " << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z << endl; " ";
			}
			
			OUTFILE.close();
			clu_count++;
	}	
	
	// Print out cluster structures (wrapped)
	
	clu_count = 0;
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
		ofstream OUTFILE;
		string   OUTNAME =  KIND + ".wrapped.frame." + to_string(f) + ".cluster." + to_string(clu_count) + ".xyz";
		OUTFILE.open(OUTNAME.data());

		OUTFILE << CLUSTERS[i].NATOMS << endl;
		OUTFILE << TRAJ_FRAME.BOXDIM.X << " " <<  TRAJ_FRAME.BOXDIM.Y << " " << TRAJ_FRAME.BOXDIM.Z << endl;

		for(int j=0; j<CLUSTERS[i].NATOMS; j++)
		{
			int UNWRAPPED_IDX = CLUSTERS[i].IDX_AND_PAR[j].first;
			
			XYZ_LINE TEMP_CRDS;
			TEMP_CRDS.TYPE = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE;
			TEMP_CRDS.X    = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X;
			TEMP_CRDS.Y    = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y;
			TEMP_CRDS.Z    = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z;

			TEMP_CRDS.X -= floor(TEMP_CRDS.X/TRAJ_FRAME.BOXDIM.X)*TRAJ_FRAME.BOXDIM.X;
			TEMP_CRDS.Y -= floor(TEMP_CRDS.Y/TRAJ_FRAME.BOXDIM.Y)*TRAJ_FRAME.BOXDIM.Y;
			TEMP_CRDS.Z -= floor(TEMP_CRDS.Z/TRAJ_FRAME.BOXDIM.Z)*TRAJ_FRAME.BOXDIM.Z;				    

			OUTFILE <<        TEMP_CRDS.TYPE;
			OUTFILE << " " << TEMP_CRDS.X;
			OUTFILE << " " << TEMP_CRDS.Y;
			OUTFILE << " " << TEMP_CRDS.Z << endl; " ";
		}
		
		OUTFILE.close();
		clu_count++;
	}
}

void PRINT_WRAPPED_AND_UNWRAPPED_FRAME(string KIND,vector<CLUSTER> & CLUSTERS, FRAME & TRAJ_FRAME, BOND_CRIT & CRITERIA, int f)
{
	// Print out frame with wrapped molecules
	
	string OUTNAME = KIND + ".wrapped.frame." + to_string(f) + ".lammpstrj";
	ofstream OUTFILE;
	OUTFILE.open(OUTNAME.data());
	
	OUTFILE << "ITEM: TIMESTEP" << endl;						     			 
	OUTFILE << f << endl;
	OUTFILE << "ITEM: NUMBER OF ATOMS" << endl;
	OUTFILE << TRAJ_FRAME.NATOMS << endl;
	OUTFILE << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.X << " 0.0" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.Y << " 0.0" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.Z << " 0.0" << endl;
	
	OUTFILE << "ITEM: ATOMS id mol type element xu yu zu" << endl;
	
	int ATM_COUNT = 0;
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
	
		for (int j=0; j<CLUSTERS[i].NATOMS; j++) 
		{
			int UNWRAPPED_IDX = CLUSTERS[i].IDX_AND_PAR[j].first;
			int ATOM_TYPE_IDX_IDX = 0;
			
			for(int k=0; k<CRITERIA.TYPES.size(); k++)
				if (CRITERIA.TYPES[k] == TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE)
					ATOM_TYPE_IDX_IDX = k;
					
			 XYZ_LINE TEMP_CRDS;
			 TEMP_CRDS.TYPE = TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE;
			 TEMP_CRDS.X	= TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X;
			 TEMP_CRDS.Y	= TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y;
			 TEMP_CRDS.Z	= TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z;

			 TEMP_CRDS.X -= floor(TEMP_CRDS.X/TRAJ_FRAME.BOXDIM.X)*TRAJ_FRAME.BOXDIM.X;
			 TEMP_CRDS.Y -= floor(TEMP_CRDS.Y/TRAJ_FRAME.BOXDIM.Y)*TRAJ_FRAME.BOXDIM.Y;
			 TEMP_CRDS.Z -= floor(TEMP_CRDS.Z/TRAJ_FRAME.BOXDIM.Z)*TRAJ_FRAME.BOXDIM.Z;

			OUTFILE 
			<< right << setw(20) << ATM_COUNT+1
			<< right << setw(20) << i+1 
			<< right << setw(4 ) << ATOM_TYPE_IDX_IDX+1 
			<< right << setw(4 ) << TEMP_CRDS.TYPE
			<< fixed << setprecision(5) << setw(15) << TEMP_CRDS.X
			<< fixed << setprecision(5) << setw(15) << TEMP_CRDS.Y
			<< fixed << setprecision(5) << setw(15) << TEMP_CRDS.Z
			<< endl;
			
			ATM_COUNT++;
		}
	}
	OUTFILE.close();						
	
	// Print out frame with unwrapped molecules
	
	OUTNAME = KIND + ".unwrapped.frame." + to_string(f) + ".lammpstrj";
	OUTFILE.open(OUTNAME.data());
	
	OUTFILE << "ITEM: TIMESTEP" << endl;						     			 
	OUTFILE << f << endl;
	OUTFILE << "ITEM: NUMBER OF ATOMS" << endl;
	OUTFILE << TRAJ_FRAME.NATOMS << endl;
	OUTFILE << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.X << " 0.0" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.Y << " 0.0" << endl;
	OUTFILE << "0.0 " << TRAJ_FRAME.BOXDIM.Z << " 0.0" << endl;
	
	OUTFILE << "ITEM: ATOMS id mol type element xu yu zu" << endl;
	
	ATM_COUNT = 0;
	
	for (int i=0; i<CLUSTERS.size(); i++)
	{
	
		for (int j=0; j<CLUSTERS[i].NATOMS; j++) 
		{
			int UNWRAPPED_IDX = CLUSTERS[i].IDX_AND_PAR[j].first;
			int ATOM_TYPE_IDX_IDX = 0;
			
			for(int k=0; k<CRITERIA.TYPES.size(); k++)
				if (CRITERIA.TYPES[k] == TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE)
					ATOM_TYPE_IDX_IDX = k;

			OUTFILE 
			<< right << setw(20) << ATM_COUNT+1
			<< right << setw(20) << i+1 
			<< right << setw(4 ) << ATOM_TYPE_IDX_IDX+1 
			<< right << setw(4 ) << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].TYPE
			<< fixed << setprecision(5) << setw(15) << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].X
			<< fixed << setprecision(5) << setw(15) << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Y
			<< fixed << setprecision(5) << setw(15) << TRAJ_FRAME.COORDS[UNWRAPPED_IDX].Z
			<< endl;
			
			ATM_COUNT++;
		}
	}
	OUTFILE.close();
}


////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////



int main(int argc, char** argv)	// The main idea: We're going to process the trajectory one frame at a time
{
	////////////////////////////////////////////////////////////
	// Read user inputs (traj file and nframes)
	////////////////////////////////////////////////////////////
	
	string XYZ_FILE  =            argv[1];
	int    NFRAMES   = STR_TO_INT(argv[2]);	
	
	////////////////////////////////////////////////////////////
	// Read criteria file names
	////////////////////////////////////////////////////////////	
	
	vector<string> CRIT_FILE;
	
	for (int i=3; i<argc; i++)
		CRIT_FILE.push_back(argv[i]); 
	
	if (CRIT_FILE.size() == 1)
		cout << endl << "Will only produce \"real\" clusters. " << endl;
	else
		cout << endl << "Will produce \"real\" and \"transition state\" clusters." << endl;
	
	////////////////////////////////////////////////////////////
	// Begin processing the file by iterating over frames
	////////////////////////////////////////////////////////////

	ifstream INTRAJ;
	INTRAJ.open(XYZ_FILE.data());
	if (!INTRAJ.is_open())
	{
		cout << "ERROR: Cannot open trajectory file: " << XYZ_FILE << endl;
		cout << "Exiting." << endl;
		exit(0);
	}
	
	for (int f=0; f<NFRAMES; f++)
	{
		cout << endl;
		cout << "Processing frame: " << f << endl;
		
		////////////////////////////////////////////////////////////
		// Setup/read the frame coordinates, replicate the sytsem to make ghost images
		////////////////////////////////////////////////////////////
		
		FRAME TRAJ_FRAME;
		
		TRAJ_FRAME.READ(INTRAJ);
		cout << "	Original frame has natoms: " << TRAJ_FRAME.COORDS.size() << endl;
		TRAJ_FRAME.REPLICATE();
		cout << "	Replicated frame has natoms: " << TRAJ_FRAME.COORDS.size() << endl;
		
		cout << "	...Frame coordinates read and replicated" << endl;
				
		// Do the clustering for the first set of bond constraints
		// We only ever expect two, so we'll do this explicitly.
		
		////////////////////////////////////////////////////////////
		// Process clusters based on "tight" constraints
		////////////////////////////////////////////////////////////
		
		cout << endl << "	...Building clusters based on TIGHT constraints..." << endl;
		
		vector<CLUSTER> TIGHT_CLUSTERS;
	
		BOND_CRIT TIGHT_CRITERIA(CRIT_FILE[0]);
		
		DO_CLUSTER(TIGHT_CRITERIA, TRAJ_FRAME, TIGHT_CLUSTERS);

		CALC_AND_PRINT_CLU_INFO(TIGHT_CLUSTERS, TRAJ_FRAME, TIGHT_CRITERIA, f);
		
		PRINT_WRAPPED_AND_UNWRAPPED_CLUSTERS("tight", TIGHT_CLUSTERS, TRAJ_FRAME, f);
		
		PRINT_WRAPPED_AND_UNWRAPPED_FRAME("tight", TIGHT_CLUSTERS, TRAJ_FRAME, TIGHT_CRITERIA, f);
		
		////////////////////////////////////////////////////////////
		// Process clusters based on "loose" constraints
		////////////////////////////////////////////////////////////
		
		if(CRIT_FILE.size() == 1) // Don't process this if only one set of constraints have been provided
			continue;
		
		cout << "	...Building clusters based on LOOSE constraints..." << endl;
		
		vector<CLUSTER> LOOSE_CLUSTERS;
		
		BOND_CRIT LOOSE_CRITERIA(CRIT_FILE[1]);
		
		DO_CLUSTER(LOOSE_CRITERIA, TRAJ_FRAME, LOOSE_CLUSTERS);
		
		////////////////////////////////////////////////////////////
		// Determine the "transition state" clusters as LOOSE - intersection(LOOSE,TIGHT)
		// And print out the corresponding results
		////////////////////////////////////////////////////////////
		
		bool EQUIVALENT; 
		
		vector<CLUSTER> TS_CLUSTERS;
		
		for (int i=LOOSE_CLUSTERS.size()-1; i>=0; i--)
		{
			EQUIVALENT = false; 
			
			for (int j=0; j<TIGHT_CLUSTERS.size(); j++)
			{	
				if (IS_PAIR_VEC_EQUAL(LOOSE_CLUSTERS[i].IDX_AND_PAR, TIGHT_CLUSTERS[j].IDX_AND_PAR))
				{
					EQUIVALENT = true;
					break;
				}
			}
			
			if (!EQUIVALENT)
				TS_CLUSTERS.push_back(LOOSE_CLUSTERS[i]);
		}

		CALC_AND_PRINT_CLU_INFO(TS_CLUSTERS, TRAJ_FRAME, LOOSE_CRITERIA, f);
		
		PRINT_WRAPPED_AND_UNWRAPPED_CLUSTERS("ts", TS_CLUSTERS, TRAJ_FRAME, f);

	}
}
