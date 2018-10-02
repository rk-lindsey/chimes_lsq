
// Setup with typical headers

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<cstring>
#include<string>
#include<sstream>
#include<map>
#include<algorithm> 


#ifdef USE_MPI
#include <mpi.h>
#endif 

// Utility functions.
#include "functions.h"
#include "util.h"

int parse_space(string line, vector<string>& tokens)
// Break a line up into tokens based on space separators.
// Returns the number of tokens parsed.
{
  string buf ;
  
  stringstream stream_parser ;

  // Strip comments beginining with ! or ##
	strip_comments(line) ;

  // Strip off terminal new line.
  int pos = line.find('\n') ;
  if ( pos != string::npos ) 
	 line.erase(pos, 1) ;

  stream_parser.str(line);
	 
  tokens.clear() ;

  while ( stream_parser >> buf ) 
	 tokens.push_back(buf) ;

  return(tokens.size() ) ;
}



void strip_comments(string &line)
// Remove all comments from the line.  Does NOT remove terminal new line.
{
  // Strip comments beginining with !
  int pos = line.find('!') ;
  if ( pos != string::npos ) 
	 line.erase(pos, line.length() - pos) ;

	// Strip comments beginning with ## 
	pos = line.find("##") ;
	if ( pos != string::npos ) 
		line.erase(pos, line.length()-pos) ;

}


void validate_num_args(int nargs, int required, string line)
// Validate that there are the required number of arguments in the line.  nargs is the current number of arguments (tokens).
{
  if ( nargs < required )
	 EXIT_MSG("Wrong number of input parameters: " + line) ;
}

string get_next_line(istream& str)
// Read a line and return it, with error checking.
{
  string line ;

  getline(str, line) ;
  if ( ! str.good() )
	 EXIT_MSG("Error reading line") ;
  
  return line ;
}

	 
bool is_true(string s) 
// Return true if the string has a recognized spelling of the word "TRUE".
{
  return (s=="true"  || s=="True"  || s=="TRUE"  || s == "T" || s == "t") ;
}


int factorial(int input)
{
	int result = 1;
	for(int i=input; i>0; i--)
		result *= i;
	
	return result;
}


void SORT_THREE_DESCEND(int & a, int & b, int & c)
{
	static int tmp;
	
	tmp = a;
	
	if(b>a)
	{
		tmp = a;
		a = b;
		b = tmp;
	}
	
	tmp = a;
	
	if(c>a)
	{
		tmp = a;
		a = c;
		c = tmp;
	}
	
	tmp = b;
	
	if(c>b)
	{
		tmp = b;
		b = c;
		c = tmp;
	}
} 


void exit_run(int value)
// Call this instead of exit(1) to properly terminate all MPI processes.
{
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD,value);
	#else
		exit(value);
	#endif

}

  
bool operator==(const vector<int>& lhs, const vector<int>& rhs) 
{
  if ( lhs.size() != rhs.size() ) return false ;

  for ( int i = 0 ; i < lhs.size() ; i++ ) 
  {
	 if ( lhs[i] != rhs[i] ) return false ;
  }
  return true ;
}

bool operator==(const vector<string>& lhs, const vector<string>& rhs) 
{
  if ( lhs.size() != rhs.size() ) return false ;

  for ( int i = 0 ; i < lhs.size() ; i++ ) 
  {
	 if ( lhs[i] != rhs[i] ) return false ;
  }
  return true ;
}


//////////////////////////////////////////
// Overloaded error message exit functions
//////////////////////////////////////////
 
void EXIT_MSG(string EXIT_STRING) // Single error message
{
	cout << EXIT_STRING << endl;
	exit_run(0);
}

void EXIT_MSG(string EXIT_STRING, string EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit_run(0);
}

void EXIT_MSG(string EXIT_STRING, double EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit_run(0);
}
