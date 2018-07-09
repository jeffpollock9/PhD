#ifndef LIKELIHOOD9_SIM_H_
#define LIKELIHOOD9_SIM_H_

#include "../Data/SimulatedData.h"
#include "../Models/Model9.h"
#include <armadillo>

inline double logLikelihood9_sim(const arma::rowvec &theta, const int week) {

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

	for (int i = 0; i < DATA_LENGTH_SIM; ++i) {

		if (WEEK_SIM[i] > week) {
			break;
		}

		if (WEEK_SIM[i] == week) {
			//NOTE: added exp()
			model.setHomeR(exp(theta[10 + HOME_ID_SIM[i]]));
			model.setAwayR(exp(theta[10 + AWAY_ID_SIM[i]]));

			logLikelihood += model.logLikelihood(START_SIM[i], END_SIM[i], IS_HOME_GOAL_SIM[i],
					IS_AWAY_GOAL_SIM[i], HOME_SCORE_SIM[i], AWAY_SCORE_SIM[i]);
		}
	}

	return logLikelihood;
}

#endif /* LIKELIHOOD9_SIM_H_ */
