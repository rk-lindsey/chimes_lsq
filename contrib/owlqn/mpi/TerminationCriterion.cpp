#include "TerminationCriterion.h"

#include "OWLQN.h"

#include <limits>
#include <iomanip>
#include <cmath>

using namespace std;

double RelativeMeanImprovementCriterion::GetValue(const OptimizerState& state, std::ostream& message) {
	double retVal = numeric_limits<double>::infinity();
	int memory = 5 ;
	if (prevVals.size() > memory) {
		double prevVal = prevVals.front();
		if (prevVals.size() == 2 * memory) prevVals.pop_front();
		double averageImprovement = (prevVal - state.GetValue()) / prevVals.size();
		double relAvgImpr = averageImprovement / fabs(state.GetValue());
		message << setprecision(4) << scientific << right;
		message << "  (" << setw(10) << relAvgImpr << ") " << flush;
		retVal = relAvgImpr;
	} else {
		message << "  (wait for " << memory << " iters) " << flush;
	}

	prevVals.push_back(state.GetValue());
	return retVal;
}
