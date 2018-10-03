#ifndef _util_h
#define _util_h

// The factorial function
int factorial(int input) ;

// Return true if the string has a recognized spelling of the word "TRUE".
bool is_true(string s)  ; 

// Break a line up into tokens based on space separators.
// Returns the number of tokens parsed.
int parse_space(string line, vector<string>& tokens) ;

// Remove all comments from the line.  Does NOT remove terminal new line.
void strip_comments(string &line) ;

// Validate that there are the required number of arguments in the line. 
void validate_num_args(int nargs, int required, string line) ;

struct ANSI_COLORS
{
	string MAGENTA;  
	string BLUE;	
	string GREEN;	
	string PINK;	
	string RED;	
	string BOLD;	
	string UNDERLINE;
	string ENDSTYLE; 
};

static const ANSI_COLORS COUT_STYLE  = 
{
	"\033[95m",    // MAGENTA  
	"\033[94m",    // BLUE     
	"\033[92m",    // GREEN    
	"\033[93m",    // PINK     
	"\033[91m",    // RED	   
	"\033[1m ",    // BOLD     
	"\033[4m ",    // UNDERLINE
	"\033[0m ",    // ENDSTYLE 

};

void SORT_THREE_DESCEND(int & a, int & b, int & c);

void exit_run(int val);
void normal_exit()  ;

bool operator==(const vector<int>& lhs, const vector<int>& rhs)  ;
bool operator==(const vector<string>& lhs, const vector<string>& rhs)  ;

void enable_fp_exceptions();

// Single error message
void EXIT_MSG(string EXIT_STRING) ; 

// error message: var
void EXIT_MSG(string EXIT_STRING, string EXIT_VAR) ;

 // error message: var
void EXIT_MSG(string EXIT_STRING, double EXIT_VAR) ;

// Read the next input line, with error checking.
string get_next_line(istream& str) ;

#endif // defined(_util_h)
