#ifndef LIKELIHOODBT_H_
#define LIKELIHOODBT_H_

#include "optimisation.h"
#include "../Models/ModelBT.h"

#ifdef bet
#include "../Data/BetData.h"
#else
#ifdef getGamma
#include "../Data/Data2010_2011.h"
#else
#include "../Data/Data2011_2012.h"
#endif
#endif

static const int nParametersBT = 21;
static const double parametersBT[nParametersBT] = { 0.366287, 0.624293,
		0.678609, -0.409694, -0.975585, -0.816468, 0.537722, 0.147756,
		0.0181379, -0.047004, 1.59125, 1.61518, 0.470065, -0.216199, -0.757202,
		-0.207892, -0.262299, -0.260536, 0.671319, -0.260545, -0.431824 };

inline double logLikelihoodBT(const DoubleVector &theta, const bool constrained,
		const int team, const int upToWeek, const bool printTrace) {

	double h = theta[0];
	double delta = theta[1];

	if (constrained) {
		if (delta < 0.0 || h < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelBT model = ModelBT(delta, h);

	double logLikelihood = 0.0, rSum = 0.0;

	DoubleVector R(20);
	for (int i = 0; i < 19; ++i) {
		R[i] = theta[2 + i];
		rSum += R[i];
	}
	R[19] = -rSum;

	for (int g = 0; g < 380; ++g) {

		if (MATCH_WEEK[g] > upToWeek) {
			break;
		}

		if (team < 0 || MATCH_HOME_ID[g] == team || MATCH_AWAY_ID[g] == team) {

			model.setHomeR(R[MATCH_HOME_ID[g]]);
			model.setAwayR(R[MATCH_AWAY_ID[g]]);

			logLikelihood += model.logLikelihood(MATCH_HOME_SCORE[g],
					MATCH_AWAY_SCORE[g]);

		}
	}

	if (printTrace) {
		trace(iteration, logLikelihood, theta);
	}

	return logLikelihood;
}

inline double logLikelihoodBT(const DoubleVector &theta, DoubleVector &grad,
		void *data) {

	Data *d = reinterpret_cast<Data*>(data);

	return logLikelihoodBT(theta, d->constrained, d->team, d->upToWeek,
			d->printTrace);
}

#endif /* LIKELIHOODBT_H_ */
