#ifndef MODELRJ_H_
#define MODELRJ_H_

#include "Model.hpp"

class ModelRJ: public Model {
public:
	ModelRJ() {
	}

	ModelRJ(const bool allSeasons) {
		this->teamMap = TEAMS_2007_2012;
	}

	DoublePair rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
			const Team& awayTeam, const League& league, const bool isLog) const;

	DoublePair rateIntegral(const double start, const double end, const int homeScore, const int awayScore,
			const Team& homeTeam, const Team& awayTeam, const League& league) const;

	void setParametersAllSeasons(const DoubleVector& theta);

	void setC(const DoubleVector& c) {
		this->c = c;
	}

	void setD(const DoubleVector& d) {
		this->d = d;
	}

	void setAB(const double A, const double B) {
		this->A = A;
		this->B = B;
	}

	const DoubleVector getKnotX() const {
		return knotX;
	}

	void addRandomKnot(Knot &addedKnot, const double proposalKnotSigma);
	void removeRandomKnot(Knot& removedKnot);
	void moveRandomKnot(Knot& oldKnot, Knot& newKnot, const double sigma);

private:
	double A, B;
	DoubleVector c, d;
	DoubleVector knotX = { 1.0, 2.0, 4.0, 5.0, 6.0, 17.0, 18.0, 20.0 };
};

//ofstream fileAlpha;

inline DoublePair ModelRJ::rate(const double t, const int homeScore, const int awayScore,
		const Team& homeTeam, const Team& awayTeam, const League& league, const bool isLog) const {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = 2 - homeIdx;

	double alphaMatchHome = c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t;
	double alphaMatchAway = c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t;

	DoublePair alphaLeagueHomeAway;

	calculateAlphaLeague(alphaLeagueHomeAway, homeScore, awayScore, homeTeam, awayTeam, league);

	int homeMatchesPlayed = league.matchesPlayedForTeam(homeTeam);
	int awayMatchesPlayed = league.matchesPlayedForTeam(awayTeam);

	double betaHome = logistic(A + B * homeMatchesPlayed);
	double betaAway = logistic(A + B * awayMatchesPlayed);

	double alphaHome = betaHome * alphaLeagueHomeAway.first + (1.0 - betaHome) * alphaMatchHome;
	double alphaAway = betaAway * alphaLeagueHomeAway.second + (1.0 - betaAway) * alphaMatchAway;

//	if (homeTeam == Team("Manchester City") && awayTeam == Team("Queens Park Rangers")) {
//		fileAlpha << t << "," << alphaHome << "," << alphaAway << std::endl;
//	}

	double homeR = resourceOfTeam(homeTeam);
	double awayR = resourceOfTeam(awayTeam);

	double logRateHome = h + rho[getRhoIdx(t)] + alphaHome * homeR - (1.0 - alphaAway) * awayR;
	double logRateAway = a + rho[getRhoIdx(t)] + alphaAway * awayR - (1.0 - alphaHome) * homeR;

	if (isLog) {
		return DoublePair(logRateHome, logRateAway);
	} else {
		return DoublePair(exp(logRateHome), exp(logRateAway));
	}
}

inline DoublePair ModelRJ::rateIntegral(const double start, const double end, const int homeScore,
		const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const {

	if (start == end) {
		return DoublePair(0.0, 0.0);
	}

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = 2 - homeIdx;
	int rhoIdx = getRhoIdx((end + start) / 2.0);

	double alphaMatchHomeStart = c[homeIdx] + (d[homeIdx] - c[homeIdx]) * start;
	double alphaMatchHomeEnd = c[homeIdx] + (d[homeIdx] - c[homeIdx]) * end;

	double alphaMatchAwayStart = c[awayIdx] + (d[awayIdx] - c[awayIdx]) * start;
	double alphaMatchAwayEnd = c[awayIdx] + (d[awayIdx] - c[awayIdx]) * end;

	DoublePair alphaLeagueHomeAway;

	calculateAlphaLeague(alphaLeagueHomeAway, homeScore, awayScore, homeTeam, awayTeam, league);

	int homeMatchesPlayed = league.matchesPlayedForTeam(homeTeam);
	int awayMatchesPlayed = league.matchesPlayedForTeam(awayTeam);

	double betaHome = logistic(A + B * homeMatchesPlayed);
	double betaAway = logistic(A + B * awayMatchesPlayed);

	double alphaHomeStart = betaHome * alphaLeagueHomeAway.first + (1.0 - betaHome) * alphaMatchHomeStart;
	double alphaHomeEnd = betaHome * alphaLeagueHomeAway.first + (1.0 - betaHome) * alphaMatchHomeEnd;

	double alphaAwayStart = betaAway * alphaLeagueHomeAway.second + (1.0 - betaAway) * alphaMatchAwayStart;
	double alphaAwayEnd = betaAway * alphaLeagueHomeAway.second + (1.0 - betaAway) * alphaMatchAwayEnd;

	double homeR = resourceOfTeam(homeTeam);
	double awayR = resourceOfTeam(awayTeam);

	double deriv = (1.0 - betaHome) * (d[homeIdx] - c[homeIdx]) * homeR
			+ (1.0 - betaAway) * (d[awayIdx] - c[awayIdx]) * awayR;

	double homeStartRate = exp(h + rho[rhoIdx] + alphaHomeStart * homeR - (1.0 - alphaAwayStart) * awayR);
	double homeEndRate = exp(h + rho[rhoIdx] + alphaHomeEnd * homeR - (1.0 - alphaAwayEnd) * awayR);

	double awayStartRate = exp(a + rho[rhoIdx] + alphaAwayStart * awayR - (1.0 - alphaHomeStart) * homeR);
	double awayEndRate = exp(a + rho[rhoIdx] + alphaAwayEnd * awayR - (1.0 - alphaHomeEnd) * homeR);

	return DoublePair((homeEndRate - homeStartRate) / deriv, (awayEndRate - awayStartRate) / deriv);
}

inline void ModelRJ::setParametersAllSeasons(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setAB(theta[39], theta[40]);
}

inline void ModelRJ::addRandomKnot(Knot &addedKnot, const double proposalKnotSigma) {

	DoubleVector xAvailable(20);
	for (int i = 0; i < 20; ++i) {
		xAvailable[i] = i + 1;
	}

	vector<Knot> knots = utilityModel.getKnots();

	for (const Knot knot : knots) {
		xAvailable.erase(remove(xAvailable.begin(), xAvailable.end(), round(knot.x)), xAvailable.end());
	}

	int nAvailable = xAvailable.size();
	double u = arma::randu();
	double sum = 0.0;
	for (const double x : xAvailable) {
		sum += 1.0 / nAvailable;
		if (sum > u) {
			addedKnot.x = x;
			addedKnot.y = rnorm(utilityModel.at(x), proposalKnotSigma);
			utilityModel.addKnot(addedKnot);
			return;
		}
	}
}

inline void ModelRJ::removeRandomKnot(Knot& removedKnot) {
	//can't choose first or last knot for death (always first or last knot)
	int n = utilityModel.getK() - 2;
	double sum = 0.0;
	double u = arma::randu();
	for (int i = 1; i <= n; ++i) {
		sum += 1.0 / n;
		if (u < sum) {
			utilityModel.removeKnot(i, removedKnot);
			return;
		}
	}
}

inline void ModelRJ::moveRandomKnot(Knot& oldKnot, Knot& newKnot, const double sigma) {
	double u = arma::randu();
	double sum = 0.0;
	int n = utilityModel.getK() - 1; //can't remove knot at (20, 0) (always last knot)
	for (int i = 0; i < n; ++i) {
		sum += 1.0 / n;
		if (u < sum) {
			utilityModel.moveKnot(i, oldKnot, newKnot, sigma);
			return;
		}
	}
}

#endif /* MODELRJ_H_ */
