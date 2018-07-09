#ifndef LIKELIHOODDR_H_
#define LIKELIHOODDR_H_

#include "optimisation.h"
#include "../Models/ModelDR.h"
#ifdef bet
#include "../Data/BetData.h"
#else
#include "../Data/Data2011_2012.h"
#endif

static const int nParametersDR = 52;
static const double parametersDR[nParametersDR] = { 1.222736, 4.943017,
		8.133927, 1.138752, 1.294172, 0.909510, 0.784157, 1.354709, 1.052938,
		0.931802, 0.890484, 0.959813, 0.598448, 1.582998, 0.602961, 0.904078,
		0.882268, 1.280912, 0.820037, 0.834084, 0.838870, 1.927805, 1.855339,
		1.101487, 0.960295, 0.686063, 0.543381, 0.759551, 0.817576, 1.295968,
		0.866458, 0.759087, 0.784868, 0.831330, 1.333189, 1.285346, 0.621969,
		0.465720, 0.788005, 0.553661, 0.177578, 0.409243, 0.718917, 1.086896,
		0.999517, 0.783985, 0.615708, 0.724675, 0.563460, 0.780059, 0.975879,
		1.420491 };

inline double logLikelihoodDR(const DoubleVector &theta, const bool constrained,
		const int team, const int upToWeek, const bool printTrace) {

	double h = theta[0];
	double rho[2] = { theta[1], theta[2] };
	double lam[4] = { theta[3], theta[4], theta[5], theta[6] };
	double mu[4] = { theta[7], theta[8], theta[9], theta[10] };
	double eps[2] = { theta[11], theta[12] };

	if (constrained) {
		for (int i = 0; i < 11; ++i) {
			if (theta[i] < 0.0) {
				return -HUGE_VAL;
			}
		}

		if (h < 1.0) {
			return -HUGE_VAL;
		}
	}

	ModelDR model(h, rho, lam, mu, eps, false);

	double logLikelihood = 0.0, alphaSum = 0.0;

	DoubleVector alpha(20), beta(20);
	for (int i = 0; i < 19; ++i) {
		alpha[i] = theta[13 + i];
		alphaSum += alpha[i];

		beta[i] = theta[32 + i];
	}
	alpha[19] = 20.0 - alphaSum;
	beta[19] = theta[51];

	for (int i = 0; i < 20; ++i) {
		if (beta[i] < 0.0 || alpha[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	for (int i = 0; i < DATA_LENGTH_2011_2012; ++i) {

		if (WEEK_2011_2012[i] > upToWeek) {
			break;
		}

		if (team < 0 || HOME_ID_2011_2012[i] == team || AWAY_ID_2011_2012[i] == team) {

			model.setAlphaHome(alpha[HOME_ID_2011_2012[i]]);
			model.setAlphaAway(alpha[AWAY_ID_2011_2012[i]]);

			model.setBetaHome(beta[HOME_ID_2011_2012[i]]);
			model.setBetaAway(beta[AWAY_ID_2011_2012[i]]);

			logLikelihood += model.logLikelihood(START_2011_2012[i], END_2011_2012[i],
					IS_HOME_GOAL_2011_2012[i], IS_AWAY_GOAL_2011_2012[i], HOME_SCORE_2011_2012[i],
					AWAY_SCORE_2011_2012[i]);

		}
	}

	if (printTrace) {
		trace(iteration, logLikelihood, theta);
	}

	return logLikelihood;
}

inline double logLikelihoodDR(const DoubleVector &theta, DoubleVector &grad,
		void *data) {

	Data *d = reinterpret_cast<Data*>(data);

	return logLikelihoodDR(theta, d->constrained, d->team, d->upToWeek,
			d->printTrace);
}

inline ScoresMatrix simulateScoresMatrixFromTheta(const int &nIterations,
		const DoubleVector &theta, const int &homeTeam, const int &awayTeam) {

	double rho[2] = { theta[1], theta[2] };
	double lam[4] = { theta[3], theta[4], theta[5], theta[6] };
	double mu[4] = { theta[7], theta[8], theta[9], theta[10] };
	double eps[2] = { theta[11], theta[12] };

	ModelDR model(theta[0], rho, lam, mu, eps, false);

	double alphaSum = 0.0;

	DoubleVector alpha(20), beta(20);
	for (int i = 0; i < 19; ++i) {
		alpha[i] = theta[13 + i];
		alphaSum += alpha[i];

		beta[i] = theta[32 + i];
	}
	alpha[19] = 20.0 - alphaSum;
	beta[19] = theta[51];

	model.setAlphaHome(alpha[homeTeam]);
	model.setAlphaAway(alpha[awayTeam]);

	model.setBetaHome(beta[homeTeam]);
	model.setBetaAway(beta[awayTeam]);

	ScoreLine score;
	int dim = 30;
	ScoresMatrix scoresMatrix(dim);

	for (int i = 0; i < nIterations; ++i) {
		score = model.simulateScoreLine();
		if (score.getHomeScore() < dim && score.getAwayScore() < dim) {
			scoresMatrix.at(score.getHomeScore(), score.getAwayScore())++;
		}
	}

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			scoresMatrix.at(i, j) /= nIterations;
		}
	}

	return scoresMatrix;
}

#endif /* LIKELIHOODDR_H_ */
