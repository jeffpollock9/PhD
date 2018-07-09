#ifndef MODEL9_H_
#define MODEL9_H_

#include <ctime>
#include "Model.h"

class ModelDR: public Model {
public:
	ModelDR(const double h, const double rho[3],
			const double lam[4], const double mu[4], const double eps[2],
			const bool randomise);
	double getRate(const double &t, const int &homeScore, const int &awayScore,
			const bool &isHome, const bool &isLog = false);
	double getRateIntegral(const double &start, const double &end,
			const int &homeScore, const int &awayScore, const bool &isHome);
	double simulateGoalTime(const double &start, const int &homeScore, const int &awayScore);
	void setAlphaHome(const double alpha) {
		this->alphaHome = alpha;
	}
	void setBetaHome(const double beta) {
		this->betaHome = beta;
	}
	void setAlphaAway(const double alpha) {
		this->alphaAway = alpha;
	}
	void setBetaAway(const double beta) {
		this->betaAway = beta;
	}
	void setLam(const double lam[4]) {
		for (int i = 0; i < 4; ++i) {
			this->lam[i + 1] = lam[i];
		}
	}
	void setMu(const double mu[4]) {
		for (int i = 0; i < 4; ++i) {
			this->mu[i + 1] = mu[i];
		}
	}
	void setEps(const double eps[2]) {
		this->epsHome = eps[0];
		this->epsAway = eps[1];
	}

private:
	double lam[5], mu[5], epsHome, alphaHome, betaHome, epsAway, alphaAway,
			betaAway;
};

inline ModelDR::ModelDR(const double h, const double rho[3],
		const double lam[4], const double mu[4], const double eps[2],
		const bool randomise) {

	setH(h);
	setRho(rho);
	this->rho[0] = 1.0;

	this->lam[0] = 1.0;
	this->mu[0] = 1.0;
	setLam(lam);
	setMu(mu);

	setAlphaHome(0.0);
	setBetaHome(0.0);
	setAlphaAway(0.0);
	setBetaAway(0.0);

	setEps(eps);
}

inline double ModelDR::getRate(const double &t, const int &homeScore,
		const int &awayScore, const bool &isHome, const bool &isLog) {

	int idx;
	if (homeScore == awayScore) {
		idx = 0;
	} else if (homeScore == 1 && awayScore == 0) {
		idx = 1;
	} else if (homeScore == 0 && awayScore == 1) {
		idx = 2;
	} else if (homeScore > awayScore) {
		idx = 3;
	} else {
		idx = 4;
	}

	if (isHome) {
		if (isLog) {
			return log(h * alphaHome * betaAway * rho[getRhoIdx(t)] * lam[idx] + epsHome * t);
		} else {
			return h * alphaHome * betaAway * rho[getRhoIdx(t)] * lam[idx] + epsHome * t;
		}
	} else {
		if (isLog) {
			return log(alphaAway * betaHome * rho[getRhoIdx(t)] * mu[idx] + epsAway * t);
		} else {
			return alphaAway * betaHome * rho[getRhoIdx(t)] * mu[idx] + epsAway * t;
		}
	}
}

inline double ModelDR::getRateIntegral(const double &start, const double &end,
		const int &homeScore, const int &awayScore, const bool &isHome) {

	int idx;
	if (homeScore == awayScore) {
		idx = 0;
	} else if (homeScore == 1 && awayScore == 0) {
		idx = 1;
	} else if (homeScore == 0 && awayScore == 1) {
		idx = 2;
	} else if (homeScore > awayScore) {
		idx = 3;
	} else {
		idx = 4;
	}

	int rhoIdx = getRhoIdx((end + start) / 2.0);

	if (isHome) {
		double x = h * alphaHome * betaAway * lam[idx] * rho[rhoIdx];
		return x * (end - start) + 0.5 * epsHome * (end * end - start * start);
	} else {
		double x = alphaAway * betaHome * mu[idx] * rho[rhoIdx];
		return x * (end - start) + 0.5 * epsAway * (end * end - start * start);
	}
}

inline double ModelDR::simulateGoalTime(const double &start, const int &homeScore, const int &awayScore) {

	int idx;
	if (homeScore == awayScore) {
		idx = 0;
	} else if (homeScore == 1 && awayScore == 0) {
		idx = 1;
	} else if (homeScore == 0 && awayScore == 1) {
		idx = 2;
	} else if (homeScore > awayScore) {
		idx = 3;
	} else {
		idx = 4;
	}

	double T = -log(arma::randu());

	double k0 = h * alphaHome * betaAway * lam[idx] + alphaAway * betaHome * mu[idx];
	double k1 = rho[1] * k0;
	double k2 = rho[2] * k0;
	double e, s;
	double eps = epsHome + epsAway;

	e = 44.0 / 90.0;
	if (start < e) {
		double integralTo44 = k0 * (e - start) + 0.5 * eps * (e * e - start * start);
		if (integralTo44 > T) {
			return (-k0 + sqrt(k0 * k0 + 2.0 * eps * (k0 * start + 0.5 * eps * start * start + T))) / eps;
		}
		e = 45.0 / 90.0;
		s = 44.0 / 90.0;
		double integralTo45 = integralTo44 + k1 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo45 > T) {
			return (-k1 + sqrt(k1 * k1 + 2.0 * eps * (k1 * s + 0.5 * eps * s * s + (T - integralTo44)))) / eps;
		}
		e = 89.0 / 90.0;
		s = 45.0 / 90.0;
		double integralTo89 = integralTo45 + k0 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo89 > T) {
			return (-k0 + sqrt(k0 * k0 + 2.0 * eps * (k0 * s + 0.5 * eps * s * s + (T - integralTo45)))) / eps;
		}
		e = 90.0 / 90.0;
		s = 89.0 / 90.0;
		double integralTo90 = integralTo89 + k2 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo90 > T) {
			return (-k2 + sqrt(k2 * k2 + 2.0 * eps * (k2 * s + 0.5 * eps * s * s + (T - integralTo89)))) / eps;
		} else {
			return 100.0;
		}
	}

	e = 45.0 / 90.0;
	if (start < e) {
		double integralTo45 = k1 * (e - start) + 0.5 * eps * (e * e - start * start);
		if (integralTo45 > T) {
			return (-k1 + sqrt(k1 * k1 + 2.0 * eps * (k1 * start + 0.5 * eps * start * start + T))) / eps;
		}
		e = 89.0 / 90.0;
		s = 45.0 / 90.0;
		double integralTo89 = integralTo45 + k0 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo89 > T) {
			return (-k0 + sqrt(k0 * k0 + 2.0 * eps * (k0 * s + 0.5 * eps * s * s + (T - integralTo45)))) / eps;
		}
		e = 90.0 / 90.0;
		s = 89.0 / 90.0;
		double integralTo90 = integralTo89 + k2 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo90 > T) {
			return (-k2 + sqrt(k2 * k2 + 2.0 * eps * (k2 * s + 0.5 * eps * s * s + (T - integralTo89)))) / eps;
		} else {
			return 100.0;
		}
	}

	e = 89.0 / 90.0;
	if (start < e) {
		double integralTo89 = k0 * (e - start) + 0.5 * eps * (e * e - start * start);
		if (integralTo89 > T) {
			return (-k0 + sqrt(k0 * k0 + 2.0 * eps * (k0 * start + 0.5 * eps * start * start + T))) / eps;
		}
		e = 90.0 / 90.0;
		s = 89.0 / 90.0;
		double integralTo90 = integralTo89 + k2 * (e - s) + 0.5 * eps * (e * e - s * s);
		if (integralTo90 > T) {
			return (-k2 + sqrt(k2 * k2 + 2.0 * eps * (k2 * s + 0.5 * eps * s * s + (T - integralTo89)))) / eps;
		} else {
			return 100.0;
		}
	}

	e = 90.0 / 90.0;
	if (start < e) {
		double integralTo90 = k2 * (e - start) + 0.5 * eps * (e * e - start * start);
		if (integralTo90 > T) {
			return (-k2 + sqrt(k2 * k2 + 2.0 * eps * (k2 * start + 0.5 * eps * start * start + T))) / eps;
		} else {
			return 100.0;
		}
	}
	return 100.0;

}

#endif /* MODEL9_H_ */

