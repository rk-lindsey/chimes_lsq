#include <vector>
#include <deque>
#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

#include "mpiexterns.h"
#include "OWLQN.h"
#include "TerminationCriterion.h"

double OptimizerState::dotProduct(const DblVec& a, const DblVec& b) {
	double result = 0;
	for (size_t i=0; i<a.size(); i++) {
		result += a[i] * b[i];
	}
	return result;
}

void OptimizerState::addMult(DblVec& a, const DblVec& b, double c) {
	for (size_t i=0; i<a.size(); i++) {
		a[i] += b[i] * c;
	}
}

void OptimizerState::add(DblVec& a, const DblVec& b) {
	for (size_t i=0; i<a.size(); i++) {
		a[i] += b[i];
	}
}

void OptimizerState::addMultInto(DblVec& a, const DblVec& b, const DblVec& c, double d) {
	for (size_t i=0; i<a.size(); i++) {
		a[i] = b[i] + c[i] * d;
	}
}

void OptimizerState::scale(DblVec& a, double b) {
	for (size_t i=0; i<a.size(); i++) {
		a[i] *= b;
	}
}

void OptimizerState::scaleInto(DblVec& a, const DblVec& b, double c) {
	for (size_t i=0; i<a.size(); i++) {
		a[i] = b[i] * c;
	}
}

void OptimizerState::MapDirByInverseHessian() {
	int count = (int)sList.size();

	if (count != 0) {
		for (int i = count - 1; i >= 0; i--) {
			alphas[i] = -dotProduct(*sList[i], dir) / roList[i];
			addMult(dir, *yList[i], alphas[i]);
		}

		const DblVec& lastY = *yList[count - 1];
		double yDotY = dotProduct(lastY, lastY);
		double scalar = roList[count - 1] / yDotY;
		scale(dir, scalar);

		for (int i = 0; i < count; i++) {
			double beta = dotProduct(*yList[i], dir) / roList[i];
			addMult(dir, *sList[i], -alphas[i] - beta);
		}
	}
}

void OptimizerState::MakeSteepestDescDir() {
	if (l1weight == 0) {
		scaleInto(dir, grad, -1);
	} else {

		for (size_t i=0; i<dim; i++) {
			if (x[i] < 0) {
				dir[i] = -grad[i] + l1weight;
			} else if (x[i] > 0) {
				dir[i] = -grad[i] - l1weight;
			} else {
				if (grad[i] < -l1weight) {
					dir[i] = -grad[i] - l1weight;
				} else if (grad[i] > l1weight) {
					dir[i] = -grad[i] + l1weight;
				} else {
					dir[i] = 0;
				}
			}
		}
	}

	steepestDescDir = dir;
}

void OptimizerState::FixDirSigns() {
	if (l1weight > 0) {
		for (size_t i = 0; i<dim; i++) {
			if (dir[i] * steepestDescDir[i] <= 0) {
				dir[i] = 0;
			}
		}
	}
}

void OptimizerState::UpdateDir() {
	MakeSteepestDescDir();
	MapDirByInverseHessian();
	FixDirSigns();

#ifdef _DEBUG
	TestDirDeriv();
#endif
}

void OptimizerState::TestDirDeriv()
{
	double dirNorm = sqrt(dotProduct(dir, dir));
	double eps = 1.05e-8 / dirNorm;
	GetNextPoint(eps);
	double val2 = EvalL1();
	double numDeriv = (val2 - value) / eps;
	double deriv = DirDeriv();
	if (!quiet && RANK == 0) {
		cout.precision(10) ;
		cout << scientific << "  Grad check: " << numDeriv << " vs. " << deriv << "  ";
	}
}

double OptimizerState::DirDeriv() const {
	if (l1weight == 0) {
		return dotProduct(dir, grad);
	} else {
		double val = 0.0;
		for (size_t i = 0; i < dim; i++) {
			if (dir[i] != 0) { 
				if (x[i] < 0) {
					val += dir[i] * (grad[i] - l1weight);
				} else if (x[i] > 0) {
					val += dir[i] * (grad[i] + l1weight);
				} else if (dir[i] < 0) {
					val += dir[i] * (grad[i] - l1weight);
				} else if (dir[i] > 0) {
					val += dir[i] * (grad[i] + l1weight);
				}
			}
		}

		return val;
	}
}

void OptimizerState::GetNextPoint(double alpha) {
	addMultInto(newX, x, dir, alpha);
	if (l1weight > 0) {
		for (size_t i=0; i<dim; i++) {
			if (x[i] * newX[i] < 0.0) {
				newX[i] = 0.0;
			}
		}
	}
#ifdef USE_MPI
	// Sync the new X so that every process has same state.
	// Compensates for possible floating point error drift.
	MPI_Bcast(newX.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
#endif
	
}

double OptimizerState::EvalL1() {
	double val = func.Eval(newX, newGrad);
	if (l1weight > 0) {
		for (size_t i=0; i<dim; i++) {
			val += fabs(newX[i]) * l1weight;
		}
	}

	return val;
}

void OptimizerState::Print(const char *filename)
{
	printVector(x, filename) ;
}

void OptimizerState::BackTrackingLineSearch() {
	double origDirDeriv = DirDeriv();
	// if a non-descent direction is chosen, the line search will break anyway, so throw here
	// The most likely reason for this is a bug in your function's gradient computation
	if (origDirDeriv >= 0) {
		cerr << "L-BFGS chose a non-descent direction: check your gradient!" << endl;
		exit(1);
	}

	double alpha = 1.0;
	double backoff = 0.5;
	if (iter == 1) {
		//alpha = 0.1;
		//backoff = 0.5;
		double normDir = sqrt(dotProduct(dir, dir));
		alpha = (1 / normDir);
		backoff = 0.1;
	}

	const double c1 = 1e-4;
	double oldValue = value;

	while (true) {
		GetNextPoint(alpha);
		value = EvalL1();

		if (value <= oldValue + c1 * origDirDeriv * alpha) break;

		if (!quiet && RANK == 0) cout << "." << flush;

		alpha *= backoff;
	}

	if (!quiet && RANK == 0) cout << endl;
}

void OptimizerState::Shift() {
	DblVec *nextS = NULL, *nextY = NULL;

	int listSize = (int)sList.size();

	if (listSize < m) {
		try {
			nextS = new vector<double>(dim);
			nextY = new vector<double>(dim);
		} catch (bad_alloc) {
			m = listSize;
			if (nextS != NULL) {
				delete nextS;
				nextS = NULL;
			}
		}
	}

	if (nextS == NULL) {
		nextS = sList.front();
		sList.pop_front();
		nextY = yList.front();
		yList.pop_front();
		roList.pop_front();
	}

	addMultInto(*nextS, newX, x, -1);
	addMultInto(*nextY, newGrad, grad, -1);
	double ro = dotProduct(*nextS, *nextY);

	sList.push_back(nextS);
	yList.push_back(nextY);
	roList.push_back(ro);

	x.swap(newX);
	grad.swap(newGrad);

	iter++;
}

void OWLQN::Minimize(DifferentiableFunction& function, const DblVec& initial, DblVec& minimum, double l1weight, double tol, int m) const {
	OptimizerState state(function, initial, m, l1weight, quiet);

	if (!quiet && RANK == 0) {
		cout << setprecision(10) << scientific << right;
		cout << endl << "Optimizing function of " << state.dim << " variables with OWL-QN parameters:" << endl;
		cout << "   l1 regularization weight: " << l1weight << "." << endl;
		cout << "   L-BFGS memory parameter (m): " << m << endl;
		cout << "   Convergence tolerance: " << tol << endl;
		cout << endl;
		cout << "Iter    n:  new_value    (conv_crit)   line_search" << endl << flush;
		cout << "Iter    0:  " << setw(16) << state.value << "  (***********) " << flush;
	}

	ostringstream str;
	termCrit->GetValue(state, str);

	while (true) {
		state.UpdateDir();
		state.BackTrackingLineSearch();

		ostringstream str;
		double termCritVal = termCrit->GetValue(state, str);
		if (!quiet && RANK == 0) {
			cout.precision(10) ;
			cout << "Iter " << setw(4) << state.iter << ":  " << setw(16) << state.value;
			cout << str.str() << flush;
		}

		if (termCritVal < tol) break;

		state.Shift();
		state.Print("restart.txt") ;
	}

	if (!quiet && RANK == 0) cout << endl;

	minimum = state.newX;
}
