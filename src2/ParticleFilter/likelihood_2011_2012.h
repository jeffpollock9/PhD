#ifndef LIKELIHOOD_2011_2012_H_
#define LIKELIHOOD_2011_2012_H_

#include "../Data/Data2011_2012.h"
#include "../Models/Model9.h"
#include <armadillo>

inline double logLikelihood9_2011_2012(const arma::rowvec &theta, const int week) {

	double c[3] = { theta[2], theta[3], theta[4] };
	double d[3] = { theta[5], theta[6], theta[7] };
	double rho[2] = { theta[8], theta[9] };

	if (c[0] > 1.0 || c[1] > 1.0 || c[2] > 1.0) {
		return -HUGE_VAL;
	}

	if (c[0] < 0.0 || c[1] < 0.0 || c[2] < 0.0) {
		return -HUGE_VAL;
	}

	if (d[0] > 1.0 || d[1] > 1.0 || d[2] > 1.0) {
		return -HUGE_VAL;
	}

	if (d[0] < 0.0 || d[1] < 0.0 || d[2] < 0.0) {
		return -HUGE_VAL;
	}

	Model9 model(theta[0], theta[1], c, d, rho);

	double logLikelihood = 0.0;

	for (int i = 0; i < DATA_LENGTH_2011_2012; ++i) {

		if (WEEK_2011_2012[i] > week) {
			break;
		}

		if (WEEK_2011_2012[i] == week) {
			//NOTE: added exp()
			model.setHomeR(exp(theta[10 + HOME_ID_2011_2012[i]]));
			model.setAwayR(exp(theta[10 + AWAY_ID_2011_2012[i]]));

			logLikelihood += model.logLikelihood(START_2011_2012[i], END_2011_2012[i],
					IS_HOME_GOAL_2011_2012[i], IS_AWAY_GOAL_2011_2012[i], HOME_SCORE_2011_2012[i],
					AWAY_SCORE_2011_2012[i]);
		}
	}

	return logLikelihood;
}

#endif /* LIKELIHOOD_2011_2012_H_ */
