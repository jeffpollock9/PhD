#ifndef LIKELIHOOD10_H_
#define LIKELIHOOD10_H_

#include "Model.h"

static const int KNOTS[8] = { 1, 2, 4, 5, 6, 17, 18, 20 };

class Model10: public Model {
public:
	Model10(double h, double a, double c[3], double d[3], double u[8], double A, double B, double rho[3]);
	double getRate(const double &t, const int &homeScore, const int &awayScore, const bool &isHome,
			const bool &isLog = false);
	double getRate(const double &t, const int &homeScore, const int &awayScore, const int week,
			const bool &isHome, const bool &isLog = false);
//	double getRateIntegral(const double &start, const double &end,
//			const int &homeScore, const int &awayScore, const bool &isHome);
//	double simulateGoalTime(const double &start, const int &homeScore,
//			const int &awayScore);

	void setC(const double c[3]);
	void setC(const double c0, const double c1, const double c2);

	void setD(const double d[3]);
	void setD(const double d0, const double d1, const double d2);

	void setAB(const double A, const double B);
	void setU(const double u[8]);

	double utility(const int position);
	//TODO: write this
	// calculate importance at the same time as utility contributions
	void calculateAlphaLeague(double &alpha, double &utilityImportance, const int homeScore,
			const int awayScore, const int team, const bool isHome, const int week);
private:
	double A, B;
	double c[3], d[3];
	double U[8];
};

inline Model10::Model10(double h, double a, double c[3], double d[3], double U[8], double A, double B,
		double rho[3]) {

	setHomeR(0.0);
	setAwayR(0.0);

	setH(h);
	setA(a);
	setRho(rho);

	setC(c);
	setD(d);
	setAB(A, B);
	setU(U);

}

inline void Model10::setC(const double c[3]) {
	this->c[0] = c[0];
	this->c[1] = c[1];
	this->c[2] = c[2];
}

inline void Model10::setC(const double c0, const double c1, const double c2) {
	this->c[0] = c0;
	this->c[1] = c1;
	this->c[2] = c2;
}

inline void Model10::setD(const double d[3]) {
	this->d[0] = d[0];
	this->d[1] = d[1];
	this->d[2] = d[2];
}

inline void Model10::setD(const double d0, const double d1, const double d2) {
	this->d[0] = d0;
	this->d[1] = d1;
	this->d[2] = d2;
}

inline void Model10::setU(const double U[8]) {
	for (int i = 0; i < 8; ++i) {
		this->U[i] = U[i];
	}
}

inline void Model10::setAB(const double A, const double B) {
	this->A = A;
	this->B = B;
}

inline double Model10::getRate(const double& t, const int& homeScore, const int& awayScore,
		const bool& isHome, const bool& isLog) {
	throw std::runtime_error("normal rate function not implemented for Model10");
}

inline double Model10::utility(const int position) {

	if (position == 1) {
		return U[0];
	}

	double mix;

	for (int i = 1; i < 8; ++i) {
		if (position <= KNOTS[i]) {
			mix = (double) (position - KNOTS[i - 1]) / (KNOTS[i] - KNOTS[i - 1]);
			return mix * U[i] + (1.0 - mix) * U[i - 1];
		}
	}

	throw std::runtime_error("utility function doesn't work");
}

inline double Model10::getRate(const double &t, const int &homeScore, const int &awayScore, const int week,
		const bool &isHome, const bool &isLog) {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = 2 - homeIdx;

	double alphaMatchHome = (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t);
	double alphaMatchAway = (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t);

	double alphaLeagueHome, alphaLeagueAway;
	double utilityImportanceHome, utilityImportanceAway;

	calculateAlphaLeague(alphaLeagueHome, utilityImportanceHome, homeScore, awayScore, homeTeam, true, week);
	calculateAlphaLeague(alphaLeagueAway, utilityImportanceAway, homeScore, awayScore, awayTeam, false, week);

	double betaHome = A * utilityImportanceHome * exp((week - 38) / B);
	double betaAway = A * utilityImportanceAway * exp((week - 38) / B);

	double alphaHome = betaHome * alphaLeagueHome + (1.0 - betaHome) * alphaMatchHome;
	double alphaAway = betaAway * alphaLeagueAway + (1.0 - betaAway) * alphaMatchAway;

	double logRate;

	if (isHome) {
		logRate = h + rho[getRhoIdx(t)] + alphaHome * homeR - (1.0 - alphaAway) * awayR;
	} else {
		logRate = a + rho[getRhoIdx(t)] + alphaAway * awayR - (1.0 - alphaHome) * homeR;
	}

	if (isLog) {
		return logRate;
	} else {
		return exp(logRate);
	}
}

#endif /* LIKELIHOOD10_H_ */
