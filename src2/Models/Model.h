#ifndef MODEL_H_
#define MODEL_H_

#include "../Utils/ScoreLine.h"
#include "../Utils/ScoresMatrix.h"
#include "../Utils/Stats.h"

#include <armadillo>

#define sim
#ifdef sim
	#include "../Data/SimulatedData.h"
#endif

static const int SUBDIVISIONS = 100000;
static const double EPS = 1e-8;
static const int NUMBER_GAMES_LEAGUE = 380;
static const int NUMBER_GAMES_TEAM = 38;
static const int NUMBER_OF_TEAMS = 20;

class Model {
public:
	virtual ~Model() {
	};
	virtual double getRate(const double &t, const int &homeScore,
			const int &awayScore, const bool &isHome,
			const bool &isLog = false) = 0;
	virtual double getRateIntegral(const double &start, const double &end,
			const int &homeScore, const int &awayScore, const bool &isHome);
	virtual double simulateGoalTime(const double &start, const int &homeScore,
			const int &awayScore);
	virtual double logLikelihood(const double &start, const double &end,
			const bool &isHomeGoal, const bool &isAwayGoal,
			const int &homeScore, const int &awayScore);

	void setHomeR(const double &R) {
		this->homeR = R;
	};
	void setAwayR(const double &R) {
		this->awayR = R;
	};

	void setH(const double h) {
		this->h = h;
	};
	void setA(const double a) {
		this->a = a;
	};

	void setRho(const double rho[2]) {
		this->rho[0] = 0.0;
		this->rho[1] = rho[0];
		this->rho[2] = rho[1];
	};

	void setRho(const double rho1, const double rho2) {
		this->rho[0] = 0.0;
		this->rho[1] = rho1;
		this->rho[2] = rho2;
	};


	ScoreLine simulateScoreLine();
	ScoresMatrix simulateScoresMatrix(int nGames);

	void writeScores(int homeID,int awayID, int week, std::ofstream &file);
	void simulateLeagueData();
protected:
	double h, a, homeR, awayR, rho[3];
	int homeTeam, awayTeam;
	int getRhoIdx(const double t) {
		if (t <= 44.0 / 90.0) {
			return 0;
		} else if (t <= 45.0 / 90) {
			return 1;
		} else if (t <= 89.0 / 90.0) {
			return 0;
		} else {
			return 2;
		}
	}
};

inline double Model::getRateIntegral(const double &start, const double &end,
		const int &homeScore, const int &awayScore, const bool &isHome) {

	double l = (end - start) / SUBDIVISIONS;

	double sum = (getRate(start, homeScore, awayScore, isHome, false)
			+ getRate(end, homeScore, awayScore, isHome, false)) / 2.0;

	for (int i = 1; i < SUBDIVISIONS; ++i) {
		sum += getRate(start + i * l, homeScore, awayScore, isHome, false);
	}

	return l * sum;
}

inline ScoreLine Model::simulateScoreLine() {

	int homeScore = 0, awayScore = 0;
	double homeRateAtT, awayRateAtT, t = 0.0;

	while (true) {
		t = simulateGoalTime(t, homeScore, awayScore);

		if (t > 1.0) {
			break;
		}

		homeRateAtT = getRate(t, homeScore, awayScore, true);
		awayRateAtT = getRate(t, homeScore, awayScore, false);

		if (arma::randu() < homeRateAtT / (homeRateAtT + awayRateAtT)) {
			homeScore++;
		} else {
			awayScore++;
		}
	}

	return ScoreLine(homeScore, awayScore);
}

inline double Model::simulateGoalTime(const double &start, const int &homeScore,
		const int &awayScore) {

	double T = -std::log(arma::randu());
	double integral = 0.0;
	double homeIntegral, awayIntegral;
	double tau = start;

	while (integral < T) {
		homeIntegral = (getRate(tau, homeScore, awayScore, true)
				+ getRate(tau + EPS, homeScore, awayScore, true)) / 2.0;
		awayIntegral = (getRate(tau, homeScore, awayScore, false)
				+ getRate(tau + EPS, homeScore, awayScore, false)) / 2.0;
		integral += EPS * (homeIntegral + awayIntegral);
		tau += EPS;
	}

	return tau;
}

inline ScoresMatrix Model::simulateScoresMatrix(int nGames) {

	ScoresMatrix scores;
	ScoreLine scoreLine;

	for (int game = 0; game < nGames; ++game) {
		scoreLine = simulateScoreLine();
		scores.at(scoreLine.getHomeScore(), scoreLine.getAwayScore())++;}

	int dim = scores.getDim();

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			scores.at(i, j) /= nGames;
		}
	}

	return scores;
}

inline double Model::logLikelihood(const double &start, const double &end,
		const bool &isHomeGoal, const bool &isAwayGoal, const int &homeScore,
		const int &awayScore) {

	double homeIntegral = getRateIntegral(start, end, homeScore, awayScore, true);
	double awayIntegral = getRateIntegral(start, end, homeScore, awayScore, false);

	if (isHomeGoal || isAwayGoal) {
		return getRate(end, homeScore, awayScore, isHomeGoal, true)
				- homeIntegral - awayIntegral;
	} else {
		return -homeIntegral - awayIntegral;
	}
}

inline void Model::simulateLeagueData() {

	static const double pi = 3.14159265359;

	arma::mat R(38, 20);
	arma::rowvec a = arma::randu(1, 20) * 2.0 * pi;

	for (int t = 0; t < 38; ++t) {
		arma::rowvec b = a + t / 15.0;
		R.row(t) = arma::cos(b) + 0.05 * arma::randn(1, 20);
	}
	writeCSV(R, "/home/jeff/workspace/Latex-Particle-Filter/data/R-simulated.csv");

	int homeID, awayID, week;
	std::ofstream file;
	file.open("/home/jeff/workspace/Latex-Particle-Filter/data/results-simulated.csv");
	file << "matchID,time,event" << std::endl;
	for (int g = 0; g < NUMBER_GAMES_LEAGUE; ++g) {

		homeID = MATCH_HOME_ID_SIM[g];
		awayID = MATCH_AWAY_ID_SIM[g];
		week = MATCH_WEEK_SIM[g];

		setHomeR(exp(R.at(week - 1, homeID)));
		setAwayR(exp(R.at(week - 1, awayID)));
		writeScores(homeID, awayID, week, file);
	}
	file.close();
}

inline void Model::writeScores(int homeID,int awayID, int week, std::ofstream &file) {
	int homeScore = 0, awayScore = 0;
	double homeRateAtT, awayRateAtT, t = 0.0;

	while (true) {
		t = simulateGoalTime(t, homeScore, awayScore);

		if (t > 1.0) {
			break;
		}

		homeRateAtT = getRate(t, homeScore, awayScore, true);
		awayRateAtT = getRate(t, homeScore, awayScore, false);

		if (arma::randu() < homeRateAtT / (homeRateAtT + awayRateAtT)) {
			file << homeID << "-" << awayID << "," << t << ",home" << std::endl;
			homeScore++;
		} else {
			file << homeID << "-" << awayID << "," << t << ",away" << std::endl;
			awayScore++;
		}
	}

}

#endif /* MODEL_H_ */
