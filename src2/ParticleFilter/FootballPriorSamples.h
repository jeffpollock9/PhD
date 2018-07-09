#ifndef FOOTBALLPRIORSAMPLES_H_
#define FOOTBALLPRIORSAMPLES_H_

#include <armadillo>
#include "likelihood_2010_2011.h"
#include "../Utils/Stats.h"

//prior parameters
static const double h0 = 0.4;
static const double a0 = 0.08;
static const double cAlpha0 = 1.5;
static const double cBeta0 = 1.5;
static const double dAlpha0 = 3.0;
static const double dBeta0 = 1.0;
static const double rho10 = 1.098;
static const double rho20 = 1.504;
static const double logRMu = -0.7;

static const int nParameters9 = 30;
static const double start9[nParameters9] = { h0, a0, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, rho10, rho20, -0.7,
		-0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7, -0.7,
		-0.7, -0.7 };

static const double sigma0 = 0.5;
static const double sigmaLogR0 = 1.0;

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
static const double proposalSigma_r = 1.0;

static arma::vec accept(nParameters9);

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

inline double logPrior(const double &x, const int &parameter) {
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
		return dnorm(x, logRMu, sigmaLogR0, true);
	}
}

inline void updateParameter(const int &parameter, arma::rowvec &theta, double &logLikelihoodAtTheta,
		const int &upToWeek) {

	double original = theta[parameter];
	double proposal = original + eps(parameter);

	double logDenomenator = logLikelihoodAtTheta + logPrior(original, parameter);

	theta[parameter] = proposal;
	double proposalLogPrior = logPrior(proposal, parameter);
	if (proposalLogPrior > -HUGE_VAL) {

		double proposalLogLikelihood = logLikelihood9_2010_2011(theta, upToWeek);
		double logNumerator = proposalLogLikelihood + proposalLogPrior;

		if (logNumerator >= logDenomenator || arma::randu() < std::exp(logNumerator - logDenomenator)) {
			++accept[parameter];
			logLikelihoodAtTheta = proposalLogLikelihood;
		} else {
			theta[parameter] = original;
		}
	} else {
		theta[parameter] = original;
	}
}

inline void updateParameters(arma::rowvec &theta, double &logLikelihoodAtTheta, const int &upToWeek) {
	for (int p = 0; p < nParameters9; ++p) {
		updateParameter(p, theta, logLikelihoodAtTheta, upToWeek);
	}
}

inline arma::mat getPriorSamples(const int n, const int burn, const int thin) {

	arma::mat samples(n, nParameters9);
	arma::rowvec theta(nParameters9);

	accept.fill(0.0);

	for (int i = 0; i < nParameters9; ++i) {
		theta[i] = start9[i];
	}

	double logLikelihoodAtTheta = logLikelihood9_2010_2011(theta, 38);

	for (int i = 0; i < burn; ++i) {
		updateParameters(theta, logLikelihoodAtTheta, 38);
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < thin; ++j) {
			updateParameters(theta, logLikelihoodAtTheta, 38);
		}
		if (i % 500 == 0) {
			std::cout << "prior sample: " << i << " of " << n << std::endl;

		}
		samples.row(i) = theta;
	}
	std::cout << "|====FINISHED PRIOR SAMPLES====|" << std::endl;

	//assign 3 newly promoted teams posterior samples from 3 relegated teams
	int match[30] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 15, 13, 24,
			25, 28, 26, 27, 29, 30 };

	arma::mat y(samples.n_rows, samples.n_cols);

	for (int i = 0; i < 30; ++i) {
		y.col(i) = samples.col(match[i] - 1);
	}

	return y;
}

#endif /* FOOTBALLPRIORSAMPLES_H_ */
