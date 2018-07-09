#ifndef SIMULATEBT_H_
#define SIMULATEBT_H_

#include "inferenceModelBT.h"

inline void updateHistory(DoubleMatrix &history, const int &nIterationsMCMC,
		DoubleVector &theta, const int &upToWeek, const double rGamma1,
		const double rGamma2, const bool &trace = false) {

	double logLikelihood = 0.0;

	for (int iteration = 0; iteration < nIterationsMCMC; ++iteration) {
		for (int parameter = 0; parameter < nParametersBT; ++parameter) {
			updateParameter(parameter, theta, history, iteration, logLikelihood,
					upToWeek, rGamma1, rGamma2);
		}
		if (trace && iteration % 500 == 0) {
			std::cout << iteration << " / " << nIterationsMCMC << std::endl;
		}
	}

}

inline DoubleVector getProbsFromHistory(const DoubleMatrix &history,
		const int homeTeam, const int awayTeam, const int burn,
		const int thin) {

	int nIterations = history.size();
	int divisor = (nIterations - burn) / thin;

	ModelBT model(0.0, 0.0);

	DoubleVector out(3), tmp(3);

	for (int i = burn; i < nIterations; i += thin) {
		model.setH(history[i][0]);
		model.setDelta(history[i][1]);

		model.setHomeR(history[i][homeTeam + 2]);
		model.setAwayR(history[i][awayTeam + 2]);

		tmp = model.getHomeDrawAway();

		out[0] += tmp[0];
		out[1] += tmp[1];
		out[2] += tmp[2];
	}
	out[0] /= divisor;
	out[1] /= divisor;
	out[2] /= divisor;

	return out;
}

inline void updateParameters(DoubleVector &theta, const int &upToWeek,
		const double rGamma1, const double rGamma2) {

	double logLikelihood = 0.0;

	for (int p = 0; p < nParametersBT; ++p) {
		updateParameter(p, theta, logLikelihood, upToWeek, rGamma1, rGamma2);
	}
}

inline DoubleMatrix getProbsMatrix(DoubleVector &theta, const int upToWeek,
		const int nSamplesMCMC, const int burn, const int thin,
		const double rGamma1, const double rGamma2) {

	int gameIdx;
	DoubleMatrix mat = createDoubleMatrix(10, 3);
	DoubleVector homeDrawAwayProbs(3);
	ModelBT model(0.0, 0.0);

	for (int b = 0; b < burn; ++b) {
		updateParameters(theta, upToWeek, rGamma1, rGamma2);
	}

	for (int i = 0; i < nSamplesMCMC; ++i) {
		for (int j = 0; j < thin; ++j) {
			updateParameters(theta, upToWeek, rGamma1, rGamma2);
		}

		for (int g = 0; g < 10; ++g) {
			gameIdx = upToWeek * 10 + g;

			model.setH(theta[0]);
			model.setDelta(theta[1]);

			model.setHomeR(theta[MATCH_HOME_ID[gameIdx] + 2]);
			model.setAwayR(theta[MATCH_AWAY_ID[gameIdx] + 2]);

			homeDrawAwayProbs = model.getHomeDrawAway();

			mat[g][0] += homeDrawAwayProbs[0];
			mat[g][1] += homeDrawAwayProbs[1];
			mat[g][2] += homeDrawAwayProbs[2];
		}
	}

	for (int g = 0; g < 10; ++g) {
		mat[g][0] /= nSamplesMCMC;
		mat[g][1] /= nSamplesMCMC;
		mat[g][2] /= nSamplesMCMC;
	}

	return mat;
}

#endif /* SIMULATEBT_H_ */
