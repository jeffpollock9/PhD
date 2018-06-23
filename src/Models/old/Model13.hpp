/*
#ifndef MODEL13_H_
#define MODEL13_H_

#include "Model.hpp"
#include "../Data/League.hpp"

class Model13: public Model {
public:
	Model13() {
	}
	Model13(const bool allSeasons) {
		this->teamMap = TEAMS_2007_2012;
	}
	Model13(const DoubleVector& logR, const double h, const double a, const DoubleVector& c,
			const DoubleVector& d, const DoubleVector& U, const double A, const double B,
			const DoubleVector& rho, const Season& season);

	DoublePair rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
			const Team& awayTeam, const League& league, const bool isLog) const;

	DoublePair rateIntegral(const double start, const double end, const int homeScore, const int awayScore,
			const Team& homeTeam, const Team& awayTeam, const League& league) const;

	void calculateAlphaLeague(DoublePair& alphaHomeAway, DoublePair& importanceHomeAway,
			const int homeScore, const int awayScore, const Team& homeTeam, const Team& awayTeam,
			const League& league) const;

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

	void setAB(const double A, const double B) {
		this->A = A;
		this->B = B;
	}

private:
	double A, B;
	DoubleVector c, d;
};

inline Model13::Model13(const DoubleVector& logR, const double h, const double a, const DoubleVector& c,
		const DoubleVector& d, const DoubleVector& y, double A, double B, const DoubleVector& rho,
		const Season& season) {

	setLogR(logR);
	setH(h);
	setA(a);
	setRho(rho);

	setC(c);
	setD(d);
	setAB(A, B);
	setKnots(y);
}

inline DoublePair Model13::rate(const double t, const int homeScore, const int awayScore,
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

	double betaHome = importanceHomeAway.first / (2.0 * utilityModel.at(1.0)) * logistic(A + B * homeMatchesPlayed);
	double betaAway = importanceHomeAway.second / (2.0 * utilityModel.at(1.0)) * logistic(A + B * awayMatchesPlayed);

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

inline DoublePair Model13::rateIntegral(const double start, const double end, const int homeScore,
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

	double betaHome = importanceHomeAway.first / (2.0 * utilityModel.at(1.0)) * logistic(A + B * homeMatchesPlayed);
	double betaAway = importanceHomeAway.second / (2.0 * utilityModel.at(1.0)) * logistic(A + B * awayMatchesPlayed);

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

inline void Model13::calculateAlphaLeague(DoublePair& alphaHomeAway, DoublePair& importanceHomeAway,
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

inline void Model13::setParameters(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 30);
	setLogR(logR);

	setAB(theta[30], theta[31]);

	setKnots( { theta[32], theta[33], theta[34], theta[35], theta[36], theta[37], theta[38], 0.0 });
}

inline double wrapperModel13(const DoubleVector& theta, DoubleVector& grad, void* modelDataAndSeason) {

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

	//B should be +ve
	if (theta[31] < 0.0) {
		return -HUGE_VAL;
	}

	for (int i = 32; i <= 38; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelDataAndSeason<Model13>* m = static_cast<ModelDataAndSeason<Model13>* >(modelDataAndSeason);

	Model13 model = m->model;
	model.setParameters(theta);

	double logLik = model.logLikelihoodSingleSeason(m->data, m->season);

	trace(iteration, logLik, theta);

	return logLik;
}

#ifdef HAVE_NLOPT
inline void Model13::findMLE(DoubleVector& theta, const Data& data, const Season& season) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, 39);

	ModelDataAndSeason<Model13> modelDataAndSeason;
	modelDataAndSeason.model = *this;
	modelDataAndSeason.data = data;
	modelDataAndSeason.season = season;

	opt.set_max_objective(wrapperModel13, &modelDataAndSeason);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}
#endif

inline void Model13::setParametersAllSeasons(const DoubleVector& theta, const bool withKnots) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setAB(theta[39], theta[40]);

	if (withKnots) {
		setKnots( { theta[41], theta[42], theta[43], theta[44], theta[45], theta[46], theta[47], 0.0 });
	}
}

inline double wrapperModel13AllSeasons(const DoubleVector& theta, DoubleVector& grad,
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

	//B should be +ve
	if (theta[40] < 0.0) {
		return -HUGE_VAL;
	}

	for (int i = 41; i <= 47; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelAndAllData<Model13>* m = static_cast<ModelAndAllData<Model13>* >(modelAndAllData);

	Model13 model = m->model;
	vector<Data> allData = m->allData;

	model.setParametersAllSeasons(theta, true);

	int nSeasons = ALL_SEASONS.size();
	double logLik = 0.0;

	for (int i = 0; i < nSeasons; ++i) {
		logLik += model.logLikelihoodSingleSeason(allData[i], ALL_SEASONS[i]);
	}

	trace(iteration, logLik, theta);

	return logLik;
}

#ifdef HAVE_NLOPT
inline void Model13::findMLE(DoubleVector& theta, const vector<Data>& allData) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, 48);

	ModelAndAllData<Model13> modelAndAllData;
	modelAndAllData.model = *this;
	modelAndAllData.allData = allData;

	opt.set_max_objective(wrapperModel13AllSeasons, &modelAndAllData);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}
#endif

#endif  MODEL13_H_
*/
