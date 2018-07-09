#ifndef MODEL10_H_
#define MODEL10_H_

#include "optimisation.h"
#include "../Models/Model10.h"

static const int nParameters10 = 31;
static const double parameters10[nParameters10] = { 0.5, 0.240626, 0.368427,
		0.536078, 4.94576, 0.9438, 4.50368, 3.55438, 1.38198, 1.36424, 1.58777,
		1.10215, 0.466303, 0.177636, 0.100824, 0.965037, 0.96654, 0.631891,
		0.914437, 1.86194, 1.72711, 0.77002, 0.435941, 0.347615, 0.403674,
		0.776573, 0.618617, 1.18875, 0.519704, 0.365687, 0.0};

inline double logLikelihood10(const DoubleVector &theta, const bool constrained,
		const int team, const int upToWeek, const bool printTrace) {

	double c = theta[2];
	double e = theta[3];
	double u[5];
	for (int i = 0; i < 5; ++i) {
		u[i] = theta[i + 4];
	}
	double rho[2] = { theta[9], theta[10] };

	if (constrained) {
		if (c > 1.0 || c < 0.0) {
			return -HUGE_VAL;
		}

		if (e > 1.0 || e < 0.0) {
			return -HUGE_VAL;
		}

		for (int i = 11; i < 31; ++i) {
			if (theta[i] < 0.0) {
				return -HUGE_VAL;
			}
		}

	}

	Model10 model(theta[0], theta[1], c, e, u, rho);

	double logLikelihood = 0.0;

	for (int i = 0; i < DATA_LENGTH; ++i) {

		if (WEEK[i] > upToWeek) {
			break;
		}

		if (team < 0 || HOME_ID[i] == team || AWAY_ID[i] == team) {

			if (i == 36) {
				std::cout << "broke" << std::endl;
			}
			model.setHomeR(theta[11 + HOME_ID[i]]);
			model.setAwayR(theta[11 + AWAY_ID[i]]);

			logLikelihood += model.logLikelihood(START[i], END[i],
					IS_HOME_GOAL[i], IS_AWAY_GOAL[i], HOME_SCORE[i],
					AWAY_SCORE[i]);

		}
	}

	if (printTrace) {
		trace(iteration, logLikelihood, theta);
	}

	return logLikelihood;
}

inline double logLikelihood10(const DoubleVector &theta, DoubleVector &grad,
		void *data) {

	Data *d = reinterpret_cast<Data*>(data);

	return logLikelihood10(theta, d->constrained, d->team, d->upToWeek,
			d->printTrace);
}
/*

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
 std::cout << score.getHomeScore() << " - " << score.getAwayScore()
 << std::endl;
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
 */

#endif /* MODEL10_H_ */
