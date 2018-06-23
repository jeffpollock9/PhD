#ifndef MCMC_HPP_
#define MCMC_HPP_

#include "../Models/Model0.hpp"
#include "../Models/Model1.hpp"
#include "../Models/Model2.hpp"
#include "../Models/Model3.hpp"
#include "../Models/Model4.hpp"

//TODO: update these lines
typedef Model0 MCMCModel;
static const string modelNumber = "0";
static const string dataPath("/home/jp148/Data/");
//static const string dataPath("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

inline DoubleVector priorMeans(const Model0& model) {
									//{ 1.0,  17., 18. };
	static const DoubleVector means = { 10.0, 4.0, 1.0 };
	return means;
}

inline DoubleVector priorMeans(const Model1& model) {
									//{ 1.0,  5.0, 17., 18.0 };
	static const DoubleVector means = { 10.0, 7.0, 4.0, 1.0 };
	return means;
}

inline DoubleVector priorMeans(const Model2& model) {
									//{ 1.0,  2.0, 5.0, 17., 18.0 };
	static const DoubleVector means = { 10.0, 8.0, 7.0, 4.0, 1.0 };
	return means;
}

inline DoubleVector priorMeans(const Model3& model) {
									//{ 1.0,  2.0, 4.0, 5.0, 6.0, 17., 18.0 };
	static const DoubleVector means = { 10.0, 8.0, 7.5, 7.0, 6.5, 4.0, 1.0 };
	return means;
}

inline DoubleVector priorMeans(const Model4& model) {
									//{ 1.0,  2.0, 4.0, 5.0, 6.0, 17., 18.0 };
	static const DoubleVector means = { 10.0, 8.0, 7.5, 7.0, 6.5, 4.0, 1.0 };
	return means;
}

inline double logPriorForParameter(const int i, const double x, const MCMCModel& model) {

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
	static const double alphaKnot = 2.5;

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
	} else if (i < 41) {
		return dexp(x, BLambda, true);
	} else {
		DoubleVector means = priorMeans(model);

		int knotIdx = i - 41;
		double beta = alphaKnot / means[knotIdx];

		return dgamma(x, alphaKnot, beta, true);
	}
}

inline void runMCMC() {

	//number of iterations
	static const int nIterations = 40000;

	//proposal values
	static const double proposalKnotSigma = 1.0;
	static const double proposalParameterSigma = 0.1;
	static const double proposalLogResourceSigma = 0.45;
	static const double proposalASigma = 0.7;
	static const double proposalBSigma = 0.035;

	const int nSeasons = ALL_SEASONS.size();
	vector<Data> allData(nSeasons);

	for (int i = 0; i < nSeasons; ++i) {
		allData[i].readEvents(dataPath + "events-" + ALL_SEASONS[i] + ".csv");
		allData[i].readMatches(dataPath + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	const bool allSeasons = true;
	MCMCModel model(allSeasons);

	const int nParameters = model.getNParameters();

	DoubleVector theta = model.getStartingValues();

	model.setParametersAllSeasons(theta);

	ofstream thetaHistoryFile, logLikHistoryFile;

	string thetaString = "thetaHistory" + modelNumber + ".txt";
	string logLikString = "logLikHistory" + modelNumber + ".txt";
	thetaHistoryFile.open(thetaString);
	logLikHistoryFile.open(logLikString);

	double logLik = model.logLikelihoodAllSeasons(allData);

	for (int i = 0; i < nIterations; ++i) {

		for (int j = 0; j < nParameters; ++j) {
			MCMCModel modelPrime = model;
			DoubleVector thetaPrime = theta;
			if (j < 10) {
				thetaPrime[j] += rnorm(0.0, proposalParameterSigma);
			} else if (j < 39) {
				thetaPrime[j] += rnorm(0.0, proposalLogResourceSigma);
			} else if (j < 40) {
				thetaPrime[j] += rnorm(0.0, proposalASigma);
			} else if (j < 41) {
				thetaPrime[j] += rnorm(0.0, proposalBSigma);
			} else {
				thetaPrime[j] += rnorm(0.0, proposalKnotSigma);
			}

			modelPrime.setParametersAllSeasons(thetaPrime);

			double logPriorPrime = logPriorForParameter(j, thetaPrime[j], modelPrime);
			if (logPriorPrime > -HUGE_VAL) {
				double logPrior = logPriorForParameter(j, theta[j], model);

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
		logLikHistoryFile << logLik << endl;
		for (int j = 0; j < nParameters; ++j) {
			thetaHistoryFile << theta[j];
			if (j < nParameters - 1) {
				thetaHistoryFile << ",";
			}
		}
		thetaHistoryFile << endl;
	}

	thetaHistoryFile.close();
	logLikHistoryFile.close();

	std::cout << "done" << std::endl;
}

#endif /* MCMC_HPP_ */
