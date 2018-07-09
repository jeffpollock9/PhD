#ifndef MODEL9_H_
#define MODEL9_H_

#include "Model.h"

class Model9: public Model {
public:
	Model9(const double h, const double a, const double c[3], const double d[3], const double rho[3]);
	double getRate(const double &t, const int &homeScore, const int &awayScore,
			const bool &isHome, const bool &isLog = false);
	double getRateIntegral(const double &start, const double &end, const int &homeScore,
			const int &awayScore, const bool &isHome);
	double simulateGoalTime(const double &start, const int &homeScore, const int &awayScore);

	void setC(const double c[3]);
	void setC(const double &c0, const double &c1, const double &c2);

	void setD(const double d[3]);
	void setD(const double &d0, const double &d1, const double &d2);

private:
	double c[3], d[3];
};

inline Model9::Model9(const double h, const double a, const double c[3],
		const double d[3], const double rho[3]) {

	setHomeR(0.0);
	setAwayR(0.0);

	setH(h);
	setA(a);
	setC(c);
	setD(d);

	setRho(rho);
}

inline double Model9::getRate(const double &t, const int &homeScore, const int &awayScore,
		const bool &isHome, const bool &isLog) {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = 2 - homeIdx;

	if (isHome) {
		if (isLog) {
			return h + rho[getRhoIdx(t)] + (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t) * homeR
					- (1.0 - (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t)) * awayR;
		} else {
			return exp(h + rho[getRhoIdx(t)] + (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t) * homeR
					- (1.0 - (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t)) * awayR);
		}
	} else {
		if (isLog) {
			return a + rho[getRhoIdx(t)] + (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t) * awayR
					- (1.0 - (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t)) * homeR;
		} else {
			return exp(a + rho[getRhoIdx(t)] + (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t) * awayR
					- (1.0 - (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t)) * homeR);
		}
	}
}

inline double Model9::getRateIntegral(const double &start, const double &end, const int &homeScore,
		const int &awayScore, const bool &isHome) {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = -(homeIdx - 1) + 1;
	int rhoIdx = getRhoIdx((end + start) / 2.0);
	double z = (d[homeIdx] - c[homeIdx]) * homeR + (d[awayIdx] - c[awayIdx]) * awayR;

	if (z == 0.0) {
		return getRate((end + start) / 2.0, homeScore, awayScore, isHome) * (end - start);
	} else {
		if (isHome) {
			double x = h + rho[rhoIdx] + c[homeIdx] * homeR + (c[awayIdx] - 1.0) * awayR;
			return (exp(x + end * z) - exp(x + start * z)) / z;
		} else {
			double x = a + rho[rhoIdx] + c[awayIdx] * awayR + (c[homeIdx] - 1.0) * homeR;
			return (exp(x + end * z) - exp(x + start * z)) / z;
		}
	}

}

inline double Model9::simulateGoalTime(const double &start, const int &homeScore, const int &awayScore) {

	int homeIdx = homeScore > awayScore ? 2 : homeScore == awayScore ? 1 : 0;
	int awayIdx = -(homeIdx - 1) + 1;

	double T = -log(arma::randu());

	double z = (d[homeIdx] - c[homeIdx]) * homeR + (d[awayIdx] - c[awayIdx]) * awayR;

	double y0 = exp(h + c[homeIdx] * homeR + (c[awayIdx] - 1.0) * awayR)
			+ exp(a + c[awayIdx] * awayR + (c[homeIdx] - 1.0) * homeR);
	double y1 = y0 * exp(rho[1]);
	double y2 = y0 * exp(rho[2]);

	if (start < 44.0 / 90.0) {
		double integralTo44 = 1.0 / z * (exp(44.0 / 90.0 * z) - exp(start * z)) * y0;
		if (integralTo44 > T) {
			return 1.0 / z * log(T * z / y0 + exp(start * z));
		} else {
			double integralTo45 = integralTo44 + 1.0 / z * (exp(45.0 / 90.0 * z) - exp(44.0 /90.0 * z)) * y1;
			if (integralTo45 > T) {
				return 1.0 / z * log((T - integralTo44) * z / y1 + exp(44.0 / 90.0 * z));
			} else {
				double integralTo89 = integralTo45 + 1.0 / z * (exp(89.0 / 90.0 * z) - exp(45.0 /90.0 * z)) * y0;
				if (integralTo89 > T) {
					return 1.0 / z * log((T - integralTo45) * z / y0 + exp(45.0 / 90.0 * z));
				} else {
					double integralTo90 = integralTo89 + 1.0 / z * (exp(z) - exp(89.0 /90.0 * z)) * y2;
					if (integralTo90 > T) {
						return 1.0 / z * log((T - integralTo89) * z / y2 + exp(89.0 / 90.0 * z));
					} else {
						return 100.0;
					}
				}
			}
		}
	} else if (start < 45.0 / 90.0) {
		double integralTo45 = 1.0 / z * (exp(45.0 / 90.0 * z) - exp(start * z)) * y1;
		if (integralTo45 > T) {
			return 1.0 / z * log(T * z / y1 + exp(start * z));
		} else {
			double integralTo89 = integralTo45 + 1.0 / z * (exp(89.0 / 90.0 * z) - exp(45.0 / 90.0 * z)) * y0;
			if (integralTo89 > T) {
				return 1.0 / z * log((T - integralTo45) * z / y0 + exp(45.0 / 90.0 * z));
			} else {
				double integralTo90 = integralTo89 + 1.0 / z * (exp(z) - exp(89.0 / 90.0 * z)) * y2;
				if (integralTo90 > T) {
					return 1.0 / z * log((T - integralTo89) * z / y2 + exp(89.0 / 90.0 * z));
				} else {
					return 100.0;
				}
			}
		}
	} else if (start < 89.0 / 90.0) {
		double integralTo89 = 1.0 / z * (exp(89.0 / 90.0 * z) - exp(start * z)) * y0;
		if (integralTo89 > T) {
			return 1.0 / z * log(T * z / y0 + exp(start * z));
		} else {
			double integralTo90 = integralTo89 + 1.0 / z * (exp(z) - exp(89.0 / 90.0 * z)) * y2;
			if (integralTo90 > T) {
				return 1.0 / z * log((T - integralTo89) * z / y2 + exp(89.0 / 90.0 * z));
			} else {
				return 100.0;
			}
		}
	} else if (start < 1.0) {
		double integralTo90 =  1.0 / z * (exp(z) - exp(89.0 / 90.0 * z)) * y2;
		if (integralTo90 > T) {
			return 1.0 / z * log(T * z / y2 + exp(start * z));
		} else {
			return 100.0;
		}
	}

	return 100.0;
}

inline void Model9::setC(const double c[3]) {
	this->c[0] = c[0];
	this->c[1] = c[1];
	this->c[2] = c[2];
}

inline void Model9::setC(const double &c0, const double &c1, const double &c2) {
	this->c[0] = c0;
	this->c[1] = c1;
	this->c[2] = c2;
}

inline void Model9::setD(const double d[3]) {
	this->d[0] = d[0];
	this->d[1] = d[1];
	this->d[2] = d[2];
}

inline void Model9::setD(const double &d0, const double &d1, const double &d2) {
	this->d[0] = d0;
	this->d[1] = d1;
	this->d[2] = d2;
}

#endif /* MODEL9_H_ */
