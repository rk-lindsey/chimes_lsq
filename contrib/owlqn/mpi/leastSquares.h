#pragma once

#include <fstream>
#include <string>

#include "OWLQN.h"
#include "mpiexterns.h"

struct LeastSquaresObjective;

class LeastSquaresProblem {
	std::vector<double> Amat;
	std::vector<double> b;
	size_t mstore, m, n, mstart, mend ;
	
	friend struct LeastSquaresObjective;
	void read_split_files(const char* matFilename, const char* bFilename)	;

public:
LeastSquaresProblem(size_t m, size_t mstore, size_t n) : Amat(m * n), b(m), mstore(mstore), m(m), n(n) {
		calc_storage() ;
	}

	LeastSquaresProblem(const char* matfile, const char* bFile, bool split_files);

	void calc_storage() ;
	double A(size_t i, size_t j) const {
#ifdef DEBUG
		if ( i >= mstart && i <= mend )
			return Amat[i - mstart + mstore * j];
		else
			ErrorMsg("A matrix access violation") ;
		return 0.0 ;
#else
		return Amat[i - mstart + mstore * j];
#endif		
	}

	double& A(size_t i, size_t j) {
#ifdef DEBUG		
		if ( i >= mstart && i <= mend )
			return Amat[i - mstart + mstore * j];
		else
			ErrorMsg("A matrix access violation") ;
		return Amat[0] ;
#else
	return Amat[i - mstart + mstore * j];
#endif	
	}
	size_t NumFeats() const { return n; }
	size_t NumInstances() const { return m; }
	size_t NumInstancesStored() const { return mstore; }
};

struct LeastSquaresObjective : public DifferentiableFunction {
	const LeastSquaresProblem& problem;
	const double l2weight;

	LeastSquaresObjective(const LeastSquaresProblem& p, double l2weight = 0) : problem(p), l2weight(l2weight) { }
	void Eval_Ax(const DblVec& input, DblVec& Ax) ;
	double Eval(const DblVec& input, DblVec& gradient);
};

void skipEmptyAndComment(std::ifstream& file, std::string& s) ;

