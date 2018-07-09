#ifndef INFERENCEMODEL9_H_
#define INFERENCEMODEL9_H_

#include "../Optimisation/likelihood9.h"
#include "../Utils/Stats.h"
#include "../Utils/distance.h"

//prior parameters
static const double h0 = 0.4048495;
static const double a0 = 0.080374;
static const double cAlpha0 = 1.5;
static const double cBeta0 = 1.5;
static const double dAlpha0 = 3.0;
static const double dBeta0 = 1.0;
static const double rho10 = 1.098;
static const double rho20 = 1.504;
static const double rAlpha0 = 2.030572;
static const double rBeta0 = 3.934892;
//static const double rGamma10 = 0.10;
//static const double rGamma20 = 0.05;

static const double start9[nParameters9] = { h0, a0, 0.5, 0.5, 0.5, 0.75, 0.75,
		0.75, rho10, rho20, 0.4, 0.41053, 0.42105, 0.43158, 0.44211, 0.45263,
		0.46316, 0.47368, 0.48421, 0.49474, 0.50526, 0.51579, 0.52632, 0.53684,
		0.54737, 0.55789, 0.56842, 0.57895, 0.58947, 0.6 };

static const double sigma0 = 0.5;

//proposal sd
static const double proposalSigma_h = 0.2;
static const double proposalSigma_a = 0.2;
static const double proposalSigma_c_losing = 0.5;
static const double proposalSigma_c_drawing = 0.2;
static const double proposalSigma_c_winning = 0.4;
static const double proposalSigma_d_losing = 0.3;
static const double proposalSigma_d_drawing = 0.25;
static const double proposalSigma_d_winning = 0.3;
static const double proposalSigma_rho = 0.7;
static const double proposalSigma_r = 0.4;

IntegerVector accept(nParameters9);

//proposal distribution
inline double eps(const int &parameter) {
	if (parameter < 1) {
		return arma::randn() * proposalSigma_h;
	} else if (parameter < 2) {
		return arma::randn() * proposalSigma_a;
	} else if (parameter < 3) {
		return arma::randn() * proposalSigma_c_losing;
	} else if (parameter < 4) {
		return arma::randn() * proposalSigma_c_drawing;
	} else if (parameter < 5) {
		return arma::randn() * proposalSigma_c_winning;
	} else if (parameter < 6) {
		return arma::randn() * proposalSigma_d_losing;
	} else if (parameter < 7) {
		return arma::randn() * proposalSigma_d_drawing;
	} else if (parameter < 8) {
		return arma::randn() * proposalSigma_d_winning;
	} else if (parameter < 10) {
		return arma::randn() * proposalSigma_rho;
	} else {
		return arma::randn() * proposalSigma_r;
	}
}

inline double logPrior(const double &x, const int &parameter,
		const DoubleVector &theta, const double &rGamma10,
		const double &rGamma20) {
	if (parameter < 1) {
		return dnorm(x, h0, sigma0, true);
	} else if (parameter < 2) {
		return dnorm(x, a0, sigma0, true);
	} else if (parameter < 5) {
		return dbeta(x, cAlpha0, cBeta0, true);
	} else if (parameter < 8) {
		return dbeta(x, dAlpha0, dBeta0, true);
	} else if (parameter < 9) {
		return dnorm(x, rho10, sigma0, true);
	} else if (parameter < 10) {
		return dnorm(x, rho20, sigma0, true);
	} else {
		DoubleVector R(20);
		for (int i = 0; i < 20; ++i) {
			R[i] = theta[i + 10];
		}
		int D1, D2;
		bubbleSwaps(D1, D2, R);

		return dgamma(x, rAlpha0, rBeta0, true) - rGamma10 * D1 - rGamma20 * D2;
	}
}

inline void updateParameter(const int &parameter, DoubleVector &theta,
		DoubleMatrix &history, const int &iteration,
		double &fullLogLikelihoodAtTheta, const int &upToWeek,
		const double &rGamma10, const double &rGamma20) {

	double original = theta[parameter];
	double proposal = original + eps(parameter);

	if (parameter == 0) {
		fullLogLikelihoodAtTheta = logLikelihood9(theta, false, -1, upToWeek,
				false);
	}

	double logDenomenator;
	if (parameter < 10) {
		logDenomenator = fullLogLikelihoodAtTheta
				+ logPrior(original, parameter, theta, rGamma10, rGamma20);
	} else {
		logDenomenator = logLikelihood9(theta, false, parameter - 10, upToWeek,
				false)
				+ logPrior(original, parameter, theta, rGamma10, rGamma20);
	}

	theta[parameter] = proposal;
	double proposalLogPrior = logPrior(proposal, parameter, theta, rGamma10,
			rGamma20);
	if (proposalLogPrior > -HUGE_VAL) {

		double proposalLogLikelihood = logLikelihood9(theta, false,
				parameter - 10, upToWeek, false);
		double logNumerator = proposalLogLikelihood + proposalLogPrior;

		if (logNumerator >= logDenomenator
				|| arma::randu() < std::exp(logNumerator - logDenomenator)) {
			accept[parameter]++;
			history[iteration][parameter] = proposal;
			if (parameter < 10) {
				fullLogLikelihoodAtTheta = proposalLogLikelihood;
			}
		} else {
			history[iteration][parameter] = original;
			theta[parameter] = original;
		}
	} else {
		history[iteration][parameter] = original;
		theta[parameter] = original;
	}

}

#endif /* INFERENCEMODEL9_H_ */
