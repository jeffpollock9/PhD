#ifndef INFERENCEMODELBT_H_
#define INFERENCEMODELBT_H_

#include "../Optimisation/likelihoodBT.h"
#include "../Utils/distance.h"

//prior parameters
static const double h0 = 0.35;
static const double deltaAlpha0 = 3.0;
static const double deltaBeta0 = 5.0;
static const double rMu0 = 0.0;
static const double sigma0 = 1.0;

static const double startBT[nParametersBT] = { h0, 0.6, -0.1, -0.0889, -0.0778,
		-0.0667, -0.0556, -0.0444, -0.0333, -0.0222, -0.0111, 0.134, 0.0111, 0.0222,
		0.0333, 0.0444, 0.0556, 0.0667, 0.0778, 0.0889, 0.1 };

//proposal sd
static const double proposalSigma_h = 0.5;
static const double proposalSigma_delta = 0.5;
static const double proposalSigma_r = 0.5;

//proposal distribution
static NormalDistribution proposalDist_h(0.0, proposalSigma_h);
static NormalGenerator eps_h(seed, proposalDist_h);

static NormalDistribution proposalDist_delta(0.0, proposalSigma_delta);
static NormalGenerator eps_delta(seed, proposalDist_delta);

static NormalDistribution proposalDist_r(0.0, proposalSigma_r);
static NormalGenerator eps_r(seed, proposalDist_r);

IntegerVector accept(nParametersBT, 0);

inline double eps(const int parameter) {
	if (parameter < 1) {
		return eps_h();
	} else if (parameter < 2) {
		return eps_delta();
	} else {
		return eps_r();
	}
}

inline double logPrior(const double &x, const int &parameter,
		const DoubleVector &theta, const double rGamma1, const double rGamma2) {
	if (parameter < 1) {
		return dnorm(x, h0, sigma0, true);
	} else if (parameter < 2) {
		return dgamma(x, deltaAlpha0, deltaBeta0, true);
	} else {
		DoubleVector R(20);
		double rSum = 0.0;
		for (int i = 0; i < 19; ++i) {
			R[i] = theta[2 + i];
			rSum += R[i];
		}
		R[19] = -rSum;
		int D1, D2;
		bubbleSwaps(D1, D2, R);

		return dnorm(x, rMu0, sigma0, true) - rGamma1 * D1 - rGamma2 * D2;
	}
}

inline void updateParameter(const int &parameter, DoubleVector &theta,
		DoubleMatrix &history, const int &iteration,
		double &fullLogLikelihoodAtTheta, const int &upToWeek,
		const double rGamma1, const double rGamma2) {

	double original = theta[parameter];
	double proposal = original + eps(parameter);

	if (parameter == 0) {
		fullLogLikelihoodAtTheta = logLikelihoodBT(theta, false, -1, upToWeek,
				false);
	}

	double logDenomenator;
	if (parameter < 2) {
		logDenomenator = fullLogLikelihoodAtTheta
				+ logPrior(original, parameter, theta, rGamma1, rGamma2);
	} else {
		logDenomenator = logLikelihoodBT(theta, false, parameter - 2, upToWeek,
				false) + logPrior(original, parameter, theta, rGamma1, rGamma2);
	}

	theta[parameter] = proposal;
	double proposalLogPrior = logPrior(proposal, parameter, theta, rGamma1,
			rGamma2);
	if (proposalLogPrior > -HUGE_VAL) {

		double proposalLogLikelihood = logLikelihoodBT(theta, false,
				parameter - 2, upToWeek, false);
		double logNumerator = proposalLogLikelihood + proposalLogPrior;

		if (logNumerator >= logDenomenator
				|| uniform() < std::exp(logNumerator - logDenomenator)) {
			accept[parameter]++;
			history[iteration][parameter] = proposal;
			if (parameter < 2) {
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

inline void updateParameter(const int &parameter, DoubleVector &theta,
		double &fullLogLikelihoodAtTheta, const int &upToWeek,
		const double rGamma1, const double rGamma2) {

	double original = theta[parameter];
	double proposal = original + eps(parameter);

	if (parameter == 0) {
		fullLogLikelihoodAtTheta = logLikelihoodBT(theta, false, -1, upToWeek,
				false);
	}

	double logDenomenator;
	if (parameter < 2) {
		logDenomenator = fullLogLikelihoodAtTheta
				+ logPrior(original, parameter, theta, rGamma1, rGamma2);
	} else {
		logDenomenator = logLikelihoodBT(theta, false, parameter - 2, upToWeek,
				false) + logPrior(original, parameter, theta, rGamma1, rGamma2);
	}

	theta[parameter] = proposal;
	double proposalLogPrior = logPrior(proposal, parameter, theta, rGamma1,
			rGamma2);
	if (proposalLogPrior > -HUGE_VAL) {

		double proposalLogLikelihood = logLikelihoodBT(theta, false,
				parameter - 2, upToWeek, false);
		double logNumerator = proposalLogLikelihood + proposalLogPrior;

		if (logNumerator >= logDenomenator
				|| uniform() < std::exp(logNumerator - logDenomenator)) {
			accept[parameter]++;
			if (parameter < 2) {
				fullLogLikelihoodAtTheta = proposalLogLikelihood;
			}
		} else {
			theta[parameter] = original;
		}
	} else {
		theta[parameter] = original;
	}

}

#endif /* INFERENCEMODELBT_H_ */
