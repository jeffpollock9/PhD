#ifndef SIMULATE9_H_
#define SIMULATE9_H_

#include "../Utils/ScoresMatrix.h"
#include "inferenceModel9.h"

inline void updateHistory(DoubleMatrix &history, const int &nIterationsMCMC,
		DoubleVector &theta, const int &upToWeek, const double &rGamma10,
		const double &rGamma20, const bool &trace = false) {

	double logLikelihood = 0.0;

	for (int iteration = 0; iteration < nIterationsMCMC; ++iteration) {
		for (int parameter = 0; parameter < 30; ++parameter) {
			updateParameter(parameter, theta, history, iteration, logLikelihood,
					upToWeek, rGamma10, rGamma20);
		}
		if (trace && iteration % 500 == 0) {
			std::cout << iteration << " / " << nIterationsMCMC << std::endl;
		}
	}

}
#ifdef deprecated

inline ScoresMatrix simulateScoresMatrixFromHistory(const DoubleMatrix &history,
		const int homeTeam, const int awayTeam, const int burn,
		const int thin) {

	double c[3] = { 0.0, 0.0, 0.0 };
	double d[3] = { 0.0, 0.0, 0.0 };
	double rho[2] = { 0.0, 0.0 };

	Model9 model = Model9(0.0, 0.0, c, d, rho);
	ScoreLine score;
	int dim = 50; //huge
	ScoresMatrix scoresMatrix(dim);
	int nIterations = history.size();
	int divisor = (nIterations - burn) / thin;

	for (int i = burn; i < nIterations; i += thin) {
		model.setH(history[i][0]);
		model.setA(history[i][1]);

		model.setC(history[i][2], history[i][3], history[i][4]);
		model.setD(history[i][5], history[i][6], history[i][7]);
		model.setRho(history[i][8], history[i][9]);

		model.setHomeR(history[i][homeTeam + 10]);
		model.setAwayR(history[i][awayTeam + 10]);

		score = model.simulateScoreLine();

		if (score.getHomeScore() >= 20 || score.getAwayScore() >= 20) {
			std::cout << "i: " << i << " " << score.getHomeScore() << " - "
					<< score.getAwayScore() << std::endl;
		}

		if (score.getHomeScore() >= 50 || score.getAwayScore() >= 50) {
			std::cout << "i: " << i << " " << score.getHomeScore() << " - "
					<< score.getAwayScore() << std::endl;
			--divisor;
		} else {
			scoresMatrix.at(score.getHomeScore(), score.getAwayScore())++;}

	}

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			scoresMatrix.at(i, j) /= divisor;
		}
	}

	return scoresMatrix;
}

DoubleMatrix simulateExpectedLeagueScores(const DoubleMatrix &history,
		const int burn, const int thin) {

	double c[3] = { 0.0, 0.0, 0.0 };
	double d[3] = { 0.0, 0.0, 0.0 };
	double rho[2] = { 0.0, 0.0 };

	Model9 model(0.0, 0.0, c, d, rho);
	ScoreLine score;
	int dim = 50; //huge
	ScoresMatrix scoresMatrix(dim);
	int nIterations = history.size();
	int divisor = (nIterations - burn) / thin;

	std::vector<ScoreLine> leagueScores;

	int homeScore, awayScore;

	double R[20];

	for (int i = burn; i < nIterations; i += thin) {

		model.setH(history[i][0]);
		model.setA(history[i][1]);

		model.setC(history[i][2], history[i][3], history[i][4]);
		model.setD(history[i][5], history[i][6], history[i][7]);
		model.setRho(history[i][8], history[i][9]);

		for (int j = 0; j < 20; ++j) {
			R[j] = history[i][j + 10];
		}

		leagueScores = model.simulateLeague(R);

		for (int g = 0; g < NUMBER_GAMES_LEAGUE; ++g) {

			homeScore = leagueScores[g].getHomeScore();
			awayScore = leagueScores[g].getAwayScore();

			if (homeScore < dim && awayScore < dim) {
				scoresMatrix.at(homeScore, awayScore)++;
			} else {
				divisor--;
			}
		}
	}

	DoubleMatrix out = createDoubleMatrix(dim, dim);

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			out[i][j] = scoresMatrix.at(i, j) / divisor;
		}
	}

	return out;
}

ScoresMatrix simulateLeagueScores(const DoubleVector &historyRow) {

	double c[3] = { 0.0, 0.0, 0.0 };
	double d[3] = { 0.0, 0.0, 0.0 };
	double rho[2] = { 0.0, 0.0 };

	Model9 model(0.0, 0.0, c, d, rho);
	ScoreLine score;
	int dim = 50; //huge
	ScoresMatrix scoresMatrix(dim);

	std::vector<ScoreLine> leagueScores;

	int homeScore, awayScore;

	double R[20];

	model.setH(historyRow[0]);
	model.setA(historyRow[1]);

	model.setC(historyRow[2], historyRow[3], historyRow[4]);
	model.setD(historyRow[5], historyRow[6], historyRow[7]);
	model.setRho(historyRow[8], historyRow[9]);

	for (int j = 0; j < 20; ++j) {
		R[j] = historyRow[j + 10];
	}

	leagueScores = model.simulateLeague(R);

	for (int g = 0; g < NUMBER_GAMES_LEAGUE; ++g) {

		homeScore = leagueScores[g].getHomeScore();
		awayScore = leagueScores[g].getAwayScore();

		if (homeScore < dim && awayScore < dim) {
			scoresMatrix.at(homeScore, awayScore)++;
		}
	}

	return scoresMatrix;
}
#endif

#endif /* SIMULATE9_H_ */
