/*
#ifndef MODEL10_H_
#define MODEL10_H_

#include "Model.hpp"
#include "../Data/League.hpp"

class Model10: public Model {
public:
	Model10() {
	}
	Model10(const bool allSeasons) {
		this->teamMap = TEAMS_2007_2012;
	}
	Model10(const DoubleVector& logR, const double h, const double a, const DoubleVector& c,
			const DoubleVector& d, const DoubleVector& y, const double A, const double B, const double C,
			const DoubleVector& rho, const Season& season);

	DoublePair rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
			const Team& awayTeam, const League& league, const bool isLog) const;

	DoublePair rateIntegral(const double start, const double end, const int homeScore, const int awayScore,
			const Team& homeTeam, const Team& awayTeam, const League& league) const;

	void calculateAlphaLeague(DoublePair& alphaHomeAway, DoublePair& importanceHomeAway, const int homeScore,
			const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const;

#ifdef HAVE_NLOPT
	void findMLE(DoubleVector& theta, const Data& data, const Season& season);
	void findMLE(DoubleVector& theta, const vector<Data>& allData);
#endif

	void setParameters(const DoubleVector& theta);
	void setParametersAllSeasons(const DoubleVector& theta, const bool withKnots);

	void setC(const DoubleVector& c) {
		this->c = c;
	}

	void setD(const DoubleVector& d) {
		this->d = d;
	}

	void setABC(const double A, const double B, const double C) {
		this->A = A;
		this->B = B;
		this->C = C;
	}

private:
	double A, B, C;
	DoubleVector c, d;
};

inline Model10::Model10(const DoubleVector& logR, const double h, const double a, const DoubleVector& c,
		const DoubleVector& d, const DoubleVector& y, double A, double B, double C, const DoubleVector& rho,
		const Season& season) {

	setLogR(logR);
	setH(h);
	setA(a);
	setRho(rho);

	setC(c);
	setD(d);
	setABC(A, B, C);
	setKnots(y);
}

inline DoublePair Model10::rate(const double t, const int homeScore, const int awayScore,
		const Team& homeTeam, const Team& awayTeam, const League& league, const bool isLog) const {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = 2 - homeIdx;

	double alphaMatchHome = c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t;
	double alphaMatchAway = c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t;

	DoublePair alphaLeagueHomeAway, importanceHomeAway;

	calculateAlphaLeague(alphaLeagueHomeAway, importanceHomeAway, homeScore, awayScore, homeTeam, awayTeam,
			league);

	int homeMatchesPlayed = league.matchesPlayedForTeam(homeTeam);
	int awayMatchesPlayed = league.matchesPlayedForTeam(awayTeam);

	double betaHome = logistic(A + B * homeMatchesPlayed + C * importanceHomeAway.first);
	double betaAway = logistic(A + B * awayMatchesPlayed + C * importanceHomeAway.second);

	double alphaHome = betaHome * alphaLeagueHomeAway.first + (1.0 - betaHome) * alphaMatchHome;
	double alphaAway = betaAway * alphaLeagueHomeAway.second + (1.0 - betaAway) * alphaMatchAway;

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

inline DoublePair Model10::rateIntegral(const double start, const double end, const int homeScore,
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

	DoublePair alphaLeagueHomeAway, importanceHomeAway;

	calculateAlphaLeague(alphaLeagueHomeAway, importanceHomeAway, homeScore, awayScore, homeTeam, awayTeam,
			league);

	int homeMatchesPlayed = league.matchesPlayedForTeam(homeTeam);
	int awayMatchesPlayed = league.matchesPlayedForTeam(awayTeam);

	double betaHome = logistic(A + B * homeMatchesPlayed + C * importanceHomeAway.first);
	double betaAway = logistic(A + B * awayMatchesPlayed + C * importanceHomeAway.second);

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

inline void Model10::calculateAlphaLeague(DoublePair& alphaHomeAway, DoublePair& importanceHomeAway,
		const int homeScore, const int awayScore, const Team& homeTeam, const Team& awayTeam,
		const League& league) const {

	static const IntegerVector HOME_POS = { 1, 2, 0, 0 };
	static const IntegerVector AWAY_POS = { 0, 0, 1, 2 };

	double homeAlphaSum = 0.0, awayAlphaSum = 0.0;
	double homeImpSum = 0.0, awayImpSum = 0.0;

	double currentHomeUtility = utilityModel.at(league.positionOfTeam(homeTeam));
	double currentAwayUtility = utilityModel.at(league.positionOfTeam(awayTeam));

	for (int i = 0; i < 4; ++i) {
		League leagueCopy = league;
		leagueCopy.currentMatch(homeTeam, awayTeam, homeScore + HOME_POS[i], awayScore + AWAY_POS[i]);

		double homeDiff = utilityModel.at(leagueCopy.positionOfTeam(homeTeam)) - currentHomeUtility;
		double awayDiff = utilityModel.at(leagueCopy.positionOfTeam(awayTeam)) - currentAwayUtility;

		homeAlphaSum += homeDiff;
		awayAlphaSum += awayDiff;

		homeImpSum += abs(homeDiff);
		awayImpSum += abs(awayDiff);
	}

	alphaHomeAway.first = logistic(homeAlphaSum);
	alphaHomeAway.second = logistic(awayAlphaSum);

	importanceHomeAway.first = homeImpSum;
	importanceHomeAway.second = awayImpSum;
}

inline void Model10::setParameters(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 30);
	setLogR(logR);

	setABC(theta[30], theta[31], theta[32]);

	setKnots( { theta[33], theta[34], theta[35], theta[36], theta[37], theta[38], theta[39], 0.0 });
}

inline double wrapperModel10(const DoubleVector& theta, DoubleVector& grad, void* modelDataAndSeason) {

	//constraints here
	DoubleVector c = { theta[2], theta[3], theta[4] };
	DoubleVector d = { theta[5], theta[6], theta[7] };

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

	//B and C should be >= 0
	if (theta[31] < 0.0 || theta[32] < 0.0) {
		return -HUGE_VAL;
	}
	for (int i = 33; i <= 39; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelDataAndSeason<Model10>* m = static_cast<ModelDataAndSeason<Model10>* >(modelDataAndSeason);

	Model10 model = m->model;
	model.setParameters(theta);

	double logLik = model.logLikelihoodSingleSeason(m->data, m->season);

	trace(iteration, logLik, theta);

	return logLik;
}

#ifdef HAVE_NLOPT
inline void Model10::findMLE(DoubleVector& theta, const Data& data, const Season& season) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, 40);

	ModelDataAndSeason<Model10> modelDataAndSeason;
	modelDataAndSeason.model = *this;
	modelDataAndSeason.data = data;
	modelDataAndSeason.season = season;

	opt.set_max_objective(wrapperModel10, &modelDataAndSeason);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}
#endif

inline void Model10::setParametersAllSeasons(const DoubleVector& theta, const bool withKnots) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setABC(theta[39], theta[40], theta[41]);

	setKnots( { theta[42], theta[43], theta[44], theta[45], theta[46], theta[47], theta[48], 0.0 });
}

inline double wrapperModel10AllSeasons(const DoubleVector& theta, DoubleVector& grad,
		void* modelAndAllData) {

	//constraints here
	DoubleVector c = { theta[2], theta[3], theta[4] };
	DoubleVector d = { theta[5], theta[6], theta[7] };

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

	//B and C >= 0
	if (theta[40] < 0.0 || theta[41] < 0.0) {
		return -HUGE_VAL;
	}

	for (int i = 42; i <= 48; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelAndAllData<Model10>* m = static_cast<ModelAndAllData<Model10>* >(modelAndAllData);

	Model10 model = m->model;
	model.setParametersAllSeasons(theta, true);

	double logLik = model.logLikelihoodAllSeasons(m->allData);

	trace(iteration, logLik, theta);

	return logLik;
}

#ifdef HAVE_NLOPT
inline void Model10::findMLE(DoubleVector& theta, const vector<Data>& allData) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, 49);

	ModelAndAllData<Model10> modelAndAllData;
	modelAndAllData.model = *this;
	modelAndAllData.allData = allData;

	opt.set_max_objective(wrapperModel10AllSeasons, &modelAndAllData);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}
#endif

#endif  MODEL10_H_
*/
