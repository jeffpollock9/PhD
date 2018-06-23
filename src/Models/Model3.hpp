#ifndef MODEL3_HPP_
#define MODEL3_HPP_

#include "Model.hpp"

class Model3: public Model {
public:
	Model3() {
	}
	Model3(const bool allSeasons) {
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
	static const int nParameters = 48;
	static const int nKnots = 8;
	DoubleVector knotX =  {1.0, 2.0, 4.0, 5.0, 6.0, 17.0, 18.0, 20.0};
	double A, B;
	DoubleVector c, d;
	DoubleVector startingValues = { 0.293642, -0.0357756, 0.347895, 0.4171, 0.499926, 0.711756, 0.556388,
			0.494807, 0.70792, 1.70526, 0.875828, 0.581649, 0.375842, 0.398453, 0.296476, 0.341083,
			0.00631813, 0.924121, -0.773783, 0.701635, 0.509805, 0.131521, 0.81902, 0.786637, 0.993145,
			0.302762, 0.467028, 0.44408, 0.369224, 0.299156, 0.246674, 0.412463, 0.395222, 0.495794, 0.702215,
			0.35111, 0.373068, 0.298448, 0.213093, -5.10682, 0.11722, 11.0, 9.0, 8.0, 7.0, 7.5,
			6.0, 1.0 };
};

inline DoublePair Model3::rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
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

inline DoublePair Model3::rateIntegral(const double start, const double end, const int homeScore,
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

inline void Model3::setParametersAllSeasons(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setAB(theta[39], theta[40]);

	DoubleVector knotY = { theta[41], theta[42], theta[43], theta[44], theta[45], theta[46], theta[47], 0.0 };

	vector<Knot> knots(nKnots);

	for (int i = 0; i < nKnots; ++i) {
		knots[i] = Knot(knotX[i], knotY[i]);
	}

	setKnots(knots);
}

inline double wrapperModel3(const DoubleVector& theta, DoubleVector& grad, void* modelAndAllData) {

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

	for (int i = 41; i < 48; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelAndAllData<Model3>* m = static_cast<ModelAndAllData<Model3>*>(modelAndAllData);

	Model3 model = m->model;
	model.setParametersAllSeasons(theta);

	double logLik = model.logLikelihoodAllSeasons(m->allData);

	trace(iteration, logLik, theta);

	return logLik;
}

inline void Model3::findMLE(DoubleVector& theta, const vector<Data>& allData) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, nParameters);

	ModelAndAllData<Model3> modelAndAllData;
	modelAndAllData.model = *this;
	modelAndAllData.allData = allData;

	opt.set_max_objective(wrapperModel3, &modelAndAllData);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}

#endif /* MODEL3_HPP_ */
