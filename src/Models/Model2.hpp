#ifndef MODEL2_HPP_
#define MODEL2_HPP_

#include "Model.hpp"

class Model2: public Model {
public:
	Model2() {
	}
	Model2(const bool allSeasons) {
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
	static const int nParameters = 46;
	static const int nKnots = 6;
	DoubleVector knotX = { 1.0, 2.0, 5.0, 17.0, 18.0, 20.0 };
	double A, B;
	DoubleVector c, d;
	DoubleVector startingValues = { 0.294744, -0.0347762, 0.348593, 0.416168, 0.500336, 0.711778, 0.55549,
			0.494364, 0.707684, 1.70528, 0.874879, 0.583193, 0.373204, 0.39631, 0.295084, 0.337824,
			0.00315937, 0.921391, -0.776456, 0.698781, 0.508501, 0.129328, 0.816169, 0.784872, 0.992863,
			0.301103, 0.466745, 0.438052, 0.367768, 0.295105, 0.244711, 0.411691, 0.394704, 0.49456, 0.701899,
			0.347178, 0.370088, 0.29855, 0.210107, -5.10682, 0.114754, 11.0, 9.0, 8.0, 6.0,
			1.0 };
};

inline DoublePair Model2::rate(const double t, const int homeScore, const int awayScore, const Team& homeTeam,
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

inline DoublePair Model2::rateIntegral(const double start, const double end, const int homeScore,
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

inline void Model2::setParametersAllSeasons(const DoubleVector& theta) {

	setH(theta[0]);
	setA(theta[1]);

	setC( { theta[2], theta[3], theta[4] });
	setD( { theta[5], theta[6], theta[7] });
	setRho( { 0.0, theta[8], theta[9] });

	DoubleVector logR(theta.begin() + 10, theta.begin() + 39);
	setLogR(logR);

	setAB(theta[39], theta[40]);

	DoubleVector knotY = { theta[41], theta[42], theta[43], theta[44], theta[45], 0.0 };

	vector<Knot> knots(nKnots);

	for (int i = 0; i < nKnots; ++i) {
		knots[i] = Knot(knotX[i], knotY[i]);
	}

	setKnots(knots);
}

inline double wrapperModel2(const DoubleVector& theta, DoubleVector& grad, void* modelAndAllData) {

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

	for (int i = 41; i < 46; ++i) {
		if (theta[i] < 0.0) {
			return -HUGE_VAL;
		}
	}

	ModelAndAllData<Model2>* m = static_cast<ModelAndAllData<Model2>*>(modelAndAllData);

	Model2 model = m->model;
	model.setParametersAllSeasons(theta);

	double logLik = model.logLikelihoodAllSeasons(m->allData);

	trace(iteration, logLik, theta);

	return logLik;
}

inline void Model2::findMLE(DoubleVector& theta, const vector<Data>& allData) {

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, nParameters);

	ModelAndAllData<Model2> modelAndAllData;
	modelAndAllData.model = *this;
	modelAndAllData.allData = allData;

	opt.set_max_objective(wrapperModel2, &modelAndAllData);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	opt.optimize(theta, f);

	std::cout << "logLik: " << f << std::endl;
}

#endif /* MODEL2_HPP_ */
