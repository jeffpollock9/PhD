#ifndef MODEL1_HPP_
#define MODEL1_HPP_

#include "Model.hpp"

class Model1: public Model {
public:
	Model1() {
	}
	Model1(const bool allSeasons) {
		this->teamMap = TEAMS_2007_2012;
	}

	DoublePair rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
			const Team& awayTeam, const League& league, const bool isLog) const;

	DoublePair rateIntegral(const double start, const double end, const int homeScore, const int awayScore,
			const Team& homeTeam, const Team& awayTeam, const League& league) const;

	void setParametersAllSeasons(const DoubleVector& theta);

	void findMLE(DoubleVector& theta, const vector<Data>& allData);

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

	static const int getNParameters() {
		return nParameters;
	}

	const DoubleVector getStartingValues() const {
		return startingValues;
	}

	const DoubleVector getKnotX() const {
		return knotX;
	}

private:
	static const int nParameters = 45;
	static const int nKnots = 5;
	DoubleVector knotX = { 1.0, 5.0, 17.0, 18.0, 20.0 };
	double A, B;
	DoubleVector c, d;
	DoubleVector startingValues = { 0.294087, -0.0353254, 0.355654, 0.415833, 0.492032, 0.718001, 0.556365,
			0.491126, 0.708273, 1.70319, 0.872544, 0.5799, 0.369559, 0.392573, 0.290132, 0.332518, 0.00305986,
			0.919296, -0.800548, 0.695654, 0.505194, 0.123791, 0.813829, 0.782743, 0.991527, 0.296443,
			0.462828, 0.43379, 0.363338, 0.289296, 0.239037, 0.407146, 0.390451, 0.490379, 0.699346, 0.342557,
			0.365358, 0.294507, 0.207645, -5.07017, 0.117218, 10.0, 8.0, 6.0, 1.0 };
};

inline DoublePair Model1::rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
		const Team& awayTeam, const League& league, const bool isLog) const {

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

inline DoublePair Model1::rateIntegral(const double start, const double end, const int homeScore,
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

inline void Model1::setParametersAllSeasons(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setAB(theta[39], theta[40]);

	DoubleVector knotY = { theta[41], theta[42], theta[43], theta[44], 0.0 };

	vector<Knot> knots(nKnots);

	for (int i = 0; i < nKnots; ++i) {
		knots[i] = Knot(knotX[i], knotY[i]);
	}

	setKnots(knots);
}

inline double wrapperModel1(const DoubleVector& theta, DoubleVector& grad, void* modelAndAllData) {

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

	for (int i = 41; i < 45; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelAndAllData<Model1>* m = static_cast<ModelAndAllData<Model1>*>(modelAndAllData);

	Model1 model = m->model;
	model.setParametersAllSeasons(theta);

	double logLik = model.logLikelihoodAllSeasons(m->allData);

	trace(iteration, logLik, theta);

	return logLik;
}

inline void Model1::findMLE(DoubleVector& theta, const vector<Data>& allData) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, nParameters);

	ModelAndAllData<Model1> modelAndAllData;
	modelAndAllData.model = *this;
	modelAndAllData.allData = allData;

	opt.set_max_objective(wrapperModel1, &modelAndAllData);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}

#endif /* MODEL1_HPP_ */
