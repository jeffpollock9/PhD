#ifndef LIKELIHOOD_HPP_
#define LIKELIHOOD_HPP_

#include "Data.hpp"
#include "PiecewiseLinearModel.hpp"
#include "../Math/Stats.hpp"

inline double normalLogLik(const PiecewiseLinearModel& model, const double sigma,
		const PiecewiseLinearData& data) {

	int n = data.size();
	double logLik = 0.0;

	for (int i = 0; i < n; ++i) {
		//TODO: temp while X is a column
		double tmp = data.atX(i)[0];
		logLik += dnorm(data.atY(i), model.at(tmp), sigma, true);
	}
	return logLik;
}

#endif /* LIKELIHOOD_HPP_ */
