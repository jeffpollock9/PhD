#ifndef LIKELIHOOD9_H_
#define LIKELIHOOD9_H_

#include "optimisation.h"
#include "../Models/Model9.h"

#ifdef bet
#include "../Data/BetData.h"
#else
#ifdef getGamma
#include "../Data/Data2010_2011.h"
#else
#include "../Data/Data2011_2012.h"
#endif
#endif

static const int nParameters9 = 30;
static const double parameters9[nParameters9] = { 0.376913, 0.109862, 0.465363,
		0.226152, 0.219615, 0.682706, 0.620086, 0.627205, 1.359511, 1.519903,
		1.133597, 0.390395, 0.159546, 0.205165, 1.042218, 0.929504, 0.686801,
		0.808856, 1.934131, 1.720620, 0.798399, 0.421770, 0.271415, 0.406772,
		0.704348, 0.554350, 1.131852, 0.530889, 0.382694, 0.000000 };

inline double logLikelihood9(const DoubleVector &theta, const bool constrained,
		const int team, const int upToWeek, const bool printTrace) {

	double c[3] = { theta[2], theta[3], theta[4] };
	double d[3] = { theta[5], theta[6], theta[7] };
	double rho[2] = { theta[8], theta[9] };

	if (constrained) {
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

		for (int i = 10; i < 30; ++i) {
			if (theta[i] < 0.0) {
				return -HUGE_VAL;
			}
		}

		if (theta[1] > theta[0]) {
			return -HUGE_VAL;
		}
	}

	Model9 model(theta[0], theta[1], c, d, rho);

	double logLikelihood = 0.0;

	for (int i = 0; i < DATA_LENGTH_2011_2012; ++i) {

		if (WEEK_2011_2012[i] > upToWeek) {
			break;
		}

		if (team < 0 || HOME_ID_2011_2012[i] == team
				|| AWAY_ID_2011_2012[i] == team) {

			model.setHomeR(theta[10 + HOME_ID_2011_2012[i]]);
			model.setAwayR(theta[10 + AWAY_ID_2011_2012[i]]);

			logLikelihood += model.logLikelihood(START_2011_2012[i],
					END_2011_2012[i], IS_HOME_GOAL_2011_2012[i],
					IS_AWAY_GOAL_2011_2012[i], HOME_SCORE_2011_2012[i],
					AWAY_SCORE_2011_2012[i]);

		}
	}

	if (printTrace) {
		trace(iteration, logLikelihood, theta);
	}

	return logLikelihood;
}

inline double logLikelihood9(const DoubleVector &theta, DoubleVector &grad,
		void *data) {

	Data *d = reinterpret_cast<Data*>(data);

	return logLikelihood9(theta, d->constrained, d->team, d->upToWeek,
			d->printTrace);
}

inline ScoresMatrix simulateScoresMatrixFromTheta(const int &nIterations,
		const DoubleVector &theta, const int &homeTeam, const int &awayTeam) {

	double c[3] = { theta[2], theta[3], theta[4] };
	double d[3] = { theta[5], theta[6], theta[7] };
	double rho[2] = { theta[8], theta[9] };

	Model9 model(theta[0], theta[1], c, d, rho);

	model.setHomeR(theta[10 + homeTeam]);
	model.setAwayR(theta[10 + awayTeam]);

	ScoreLine score;
	int dim = 30;
	ScoresMatrix scoresMatrix(dim);

	for (int i = 0; i < nIterations; ++i) {
		score = model.simulateScoreLine();

		if (score.getHomeScore() >= dim || score.getAwayScore() >= dim) {
//			std::cout << score.getHomeScore() << " - " << score.getAwayScore()
//					<< std::endl;
		} else {
			scoresMatrix.at(score.getHomeScore(), score.getAwayScore())++;}
		}

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			scoresMatrix.at(i, j) /= nIterations;
		}
	}

	return scoresMatrix;
}

#endif /* LIKELIHOOD9_H_ */
