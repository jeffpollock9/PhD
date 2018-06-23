#ifndef RJMCMC11_HPP_
#define RJMCMC11_HPP_

#include "../Models/ModelRJ.hpp"

static const string dataPath("/home/jp148/Data/");
//static const string dataPath("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

inline double logPriorForParameter(const int i, const double x) {

	static const double sigma0 = 0.5;
	static const double sigmaA0 = 2.0;
	static const double sigmaLogR0 = 1.0;

	static const double h0 = 0.4;
	static const double a0 = 0.08;
	static const double cAlpha0 = 1.5;
	static const double cBeta0 = 1.5;
	static const double dAlpha0 = 3.0;
	static const double dBeta0 = 1.0;
	static const double rho10 = 1.098;
	static const double rho20 = 1.504;
	static const double logR0 = -0.7;
	static const double A0 = 0.0;
	static const double BLambda = 0.5;

	if (i < 1) {
		return dnorm(x, h0, sigma0, true);
	} else if (i < 2) {
		return dnorm(x, a0, sigma0, true);
	} else if (i < 5) {
		return dbeta(x, cAlpha0, cBeta0, true);
	} else if (i < 8) {
		return dbeta(x, dAlpha0, dBeta0, true);
	} else if (i < 9) {
		return dnorm(x, rho10, sigma0, true);
	} else if (i < 10) {
		return dnorm(x, rho20, sigma0, true);
	} else if (i < 39) {
		return dnorm(x, logR0, sigmaLogR0, true);
	} else if (i < 40) {
		return dnorm(x, A0, sigmaA0, true);
	} else {
		return dexp(x, BLambda, true);
	}
}

inline double logPriorForKnot(const Knot knot) {

	if (knot.x == 20.0) {
		throw runtime_error("trying to get prior for knot 20");
	}

	static const double alphaKnot = 2.5;
	static const DoubleVector means = { 10.0, 8.0, 7.75, 7.5, 7.0, 6.5, 6.27273, 6.04545, 5.81818, 5.59091,
			5.36364, 5.13636, 4.90909, 4.68182, 4.45455, 4.22727, 4.0, 1.0, 0.5, 0.0 };

	double beta = alphaKnot / means[round(knot.x - 1)];

	return dgamma(knot.y, alphaKnot, beta, true);
}

inline void runRJMCMC() {

	//number of iterations
	static const int nIterations = 100000;

	//proposal values
	static const double proposalKnotSigma = 1.0;
	static const double proposalParameterSigma = 0.1;
	static const double proposalLogResourceSigma = 0.45;
	static const double proposalASigma = 0.7;
	static const double proposalBSigma = 0.035;

	//binomial prior on knots
	static const double priorKnotsProb = 0.05;
	static const int priorMaxKnots = 20;

	const int nSeasons = ALL_SEASONS.size();
	vector<Data> allData(nSeasons);

	for (int i = 0; i < nSeasons; ++i) {
		allData[i].readEvents(dataPath + "events-" + ALL_SEASONS[i] + ".csv");
		allData[i].readMatches(dataPath + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	ModelRJ model(allSeasons);

	//chain starting values
	DoubleVector theta = { 0.289096, -0.0400848, 0.354038, 0.410594, 0.483351, 0.730718, 0.56244, 0.495277,
			0.707502, 1.70607, 0.844082, 0.536788, 0.323115, 0.345956, 0.242618, 0.282746, -0.0585018,
			0.894378, -0.99492, 0.662983, 0.462712, 0.0614064, 0.785704, 0.751695, 0.969919, 0.245218,
			0.418171, 0.395916, 0.31563, 0.233076, 0.181737, 0.360352, 0.342446, 0.445395, 0.661706, 0.297798,
			0.315175, 0.240284, 0.15638, -4.20749, 0.100131 };

	DoubleVector xStart = { 1.0, 2.0, 4.0, 5.0, 6.0, 17.0, 18.0, 20.0 };
	DoubleVector yStart = { 11.0, 9.0, 8.0, 7.0, 7.5, 6.0, 1.0, 0.0 };

	int k = yStart.size();
	vector<Knot> knots(k);

	for (int i = 0; i < k; ++i) {
		knots[i] = Knot(xStart[i], yStart[i]);
	}

	model.setParametersAllSeasons(theta);
	model.utilityModel.setKnots(knots);

	IntegerVector kHistory(nIterations);
	arma::mat thetaHistory(nIterations, 41);

	double birthRateToPrime = -1.0, birthRateFromPrime = -1.0;
	double deathRateToPrime = -1.0, deathRateFromPrime = -1.0;
	double logLik = model.logLikelihoodAllSeasons(allData);

	ofstream kHistoryFile, knotXHistoryFile, knotYHistoryFile, thetaHistoryFile;

	kHistoryFile.open("kHistory.txt");
	knotXHistoryFile.open("knotXHistory.txt");
	knotYHistoryFile.open("knotYHistory.txt");
	thetaHistoryFile.open("thetaHistory.txt");

	for (int i = 0; i < nIterations; ++i) {

		if (k == 2) {
			birthRateToPrime = 1.0 / 2.0;
			deathRateToPrime = 0.0;
			//moveRate = 1.0 / 2.0
		} else if (k == priorMaxKnots) {
			birthRateToPrime = 0.0;
			deathRateToPrime = 1.0 / 2.0;
			//moveRate = 1.0 / 2.0
		} else {
			birthRateToPrime = 1.0 / 3.0;
			deathRateToPrime = 1.0 / 3.0;
			//moveRate <- 1.0 / 3.0
		}

		//for a move, k is constant so probability of proposing a
		//move from and to prime will be equal
		if (k == 3) {
			birthRateFromPrime = 1.0 / 2.0;
			deathRateFromPrime = 1.0 / 3.0;
		} else if (k == 19) {
			birthRateFromPrime = 1.0 / 3.0;
			deathRateFromPrime = 1.0 / 2.0;
		} else {
			birthRateFromPrime = 1.0 / 3.0;
			deathRateFromPrime = 1.0 / 3.0;
		}

		double u = arma::randu();
		if (u <= birthRateToPrime) {
			ModelRJ modelPrime = model;
			int kPrime = k + 1;
			Knot addedKnot;
			modelPrime.addRandomKnot(addedKnot, proposalKnotSigma);

			double logPriorPrime = dbinom(kPrime, priorMaxKnots, priorKnotsProb, true)
					+ logPriorForKnot(addedKnot);
			if (logPriorPrime > -HUGE_VAL) {
				double logPrior = dbinom(k, priorMaxKnots, priorKnotsProb, true);

				double logProposalToPrime = log(1.0 / (20.0 - k)) + log(birthRateToPrime)
						+ dnorm(addedKnot.y, model.utilityModel.at(addedKnot.x), proposalKnotSigma, true);
				double logProposalFromPrime = log(1.0 / (kPrime - 2)) + log(deathRateFromPrime);

				double logLikPrime = modelPrime.logLikelihoodAllSeasons(allData);

				double logNumerator = logLikPrime + logPriorPrime + logProposalFromPrime;
				double logDenomenator = logLik + logPrior + logProposalToPrime;

				if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
					logLik = logLikPrime;
					k = kPrime;
					model = modelPrime;
				}
			}
		} else if (u <= birthRateToPrime + deathRateToPrime) {
			ModelRJ modelPrime = model;
			int kPrime = k - 1;
			Knot removedKnot;
			modelPrime.removeRandomKnot(removedKnot);

			double logPriorPrime = dbinom(kPrime, priorMaxKnots, priorKnotsProb, true);
			double logPrior = dbinom(k, priorMaxKnots, priorKnotsProb, true) + logPriorForKnot(removedKnot);

			double logProposalToPrime = log(1.0 / (k - 2)) + log(deathRateToPrime);
			double logProposalFromPrime = log(1.0 / (20.0 - kPrime)) + log(birthRateFromPrime)
					+ dnorm(removedKnot.y, modelPrime.utilityModel.at(removedKnot.x), proposalKnotSigma,
							true);

			double logLikPrime = modelPrime.logLikelihoodAllSeasons(allData);

			double logNumerator = logLikPrime + logPriorPrime + logProposalFromPrime;
			double logDenomenator = logLik + logPrior + logProposalToPrime;

			if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
				logLik = logLikPrime;
				k = kPrime;
				model = modelPrime;
			}
		} else {
			ModelRJ modelPrime = model;

			Knot oldKnot, newKnot;
			modelPrime.moveRandomKnot(oldKnot, newKnot, proposalKnotSigma);

			double logPriorPrime = logPriorForKnot(newKnot);
			if (logPriorPrime > -HUGE_VAL) {
				double logPrior = logPriorForKnot(oldKnot);

				double logLikPrime = modelPrime.logLikelihoodAllSeasons(allData);

				double logNumerator = logLikPrime + logPriorPrime;
				double logDenomenator = logLik + logPrior;

				if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
					logLik = logLikPrime;
					model = modelPrime;
				}
			}
		}

		for (int j = 0; j < 41; ++j) {
			ModelRJ modelPrime = model;
			DoubleVector thetaPrime = theta;
			if (j < 10) {
				thetaPrime[j] += rnorm(0.0, proposalParameterSigma);
			} else if (j < 39) {
				thetaPrime[j] += rnorm(0.0, proposalLogResourceSigma);
			} else if (j < 40) {
				thetaPrime[j] += rnorm(0.0, proposalASigma);
			} else {
				thetaPrime[j] += rnorm(0.0, proposalBSigma);
			}

			modelPrime.setParametersAllSeasons(thetaPrime);

			double logPriorPrime = logPriorForParameter(j, thetaPrime[j]);
			if (logPriorPrime > -HUGE_VAL) {
				double logPrior = logPriorForParameter(j, theta[j]);

				double logLikPrime = modelPrime.logLikelihoodAllSeasons(allData);

				double logNumerator = logLikPrime + logPriorPrime;
				double logDenomenator = logLik + logPrior;

				if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
					logLik = logLikPrime;
					model.setParametersAllSeasons(thetaPrime);
					theta = thetaPrime;
				}
			}
		}

		//store current chain values
		vector<Knot> knots = model.utilityModel.getKnots();
		kHistoryFile << k << endl;
		for (int j = 0; j < k; ++j) {
			knotXHistoryFile << knots[j].x << ",";
			knotYHistoryFile << knots[j].y << ",";
		}
		for (int j = 0; j < 41; ++j) {
			thetaHistoryFile << theta[j];
			if (j < 40) {
				thetaHistoryFile << ",";
			}
		}
		knotXHistoryFile << endl;
		knotYHistoryFile << endl;
		thetaHistoryFile << endl;
	}

	kHistoryFile.close();
	knotXHistoryFile.close();
	knotYHistoryFile.close();
	thetaHistoryFile.close();

	std::cout << "done" << std::endl;
}

#endif /* RJMCMC_HPP_ */
