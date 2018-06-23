#ifndef RJMCMC_HPP_
#define RJMCMC_HPP_

#include "../Defn/defn.hpp"
#include "likelihood.hpp"

void runRJMCMC1() {

	string path("/home/jeff/workspace/Latex-Utility-Model/data/test/");
	PiecewiseLinearData data(path + "X.csv", path + "Y.csv", true);

	const arma::rowvec minX = arma::min(data.getX());
	const arma::rowvec maxX = arma::max(data.getX());

	static const int nIterations = 30000;

	//proposal values
	static const double proposalKnotSigma = 0.25;
	static const double proposalErrorSigma = 0.01;

	//accept proportions
	int acceptBirth = 0;
	int acceptDeath = 0;
	int acceptMove = 0;
	int acceptSigma = 0;

	//binomial prior on knots
	static const double priorKnotsProb = 0.05;
	static const int priorMaxKnots = 20;

	//gamma prior on error
	static const double priorErrorAlpha = 0.001;
	static const double priorErrorBeta = 0.001;

	IntegerVector kHistory(nIterations);
	DoubleVector sigmaErrorHistory(nIterations);
	vector<vector<Knot> > knotsHistory(nIterations);

	//chain starting values
	int k = floor(priorKnotsProb * priorMaxKnots);
	vector<Knot> knots(k);
	for (int i = 0; i < k; ++i) {
		knots[i] = data.randomKnot();
	}

//	double sigmaError = priorErrorAlpha / priorErrorBeta;
	double sigmaError = 0.3;

	PiecewiseLinearModel model(knots);

	double birthRate = -1.0;
	double deathRate = -1.0;

	for (int i = 0; i < nIterations; ++i) {

		if (i % 2500 == 0L) {
			cout << "iteration: " << i << " of " << nIterations << endl;
		}

		//firstly move the knots
		if (k == 1) {
			birthRate = 1.0 / 2.0;
			deathRate = 0.0;
			//moveRate = 1.0 / 2.0
		} else if (k == priorMaxKnots) {
			birthRate = 0.0;
			deathRate = 1.0 / 2.0;
			//moveRate = 1.0 / 2.0
		} else {
			birthRate = 1.0 / 3.0;
			deathRate = 1.0 / 3.0;
			//moveRate <- 1.0 / 3.0
		}

		double u = arma::randu();
		if (u <= birthRate) {
			//birth of a new knot
			//TODO: temp while there is only one variable in X
			double xPrime = runif(minX[0], maxX[0]);
			double yPrime = rnorm(model.at(xPrime), proposalKnotSigma);
			int kPrime = k + 1;

			PiecewiseLinearModel modelPrime = model;

			modelPrime.addKnot(Knot(xPrime, yPrime));

			//TODO: once data is normalised, will use N(0, sigma^2) prior on heights of knots
			//in addition to prior on number of knots
			double logPriorPrime = dbinom(kPrime, priorMaxKnots, priorKnotsProb, true);
			double logPrior = dbinom(k, priorMaxKnots, priorKnotsProb, true);

			double logProposalToPrime = dunif(xPrime, minX[0], maxX[0], true) + dnorm(yPrime, model.at(xPrime), proposalKnotSigma, true);
			double logProposalFromPrime = log(1.0 / kPrime);

			double logLikPrime = normalLogLik(modelPrime, sigmaError, data);
			double logLik = normalLogLik(model, sigmaError, data);

			double logNumerator = logLikPrime + logPriorPrime + logProposalFromPrime;
			double logDenomenator = logLik + logPrior + logProposalToPrime;

			if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
				k = kPrime;
				model = modelPrime;
				++acceptBirth;
			}
		} else if (u <= birthRate + deathRate) {
			//death of randomly chosen existing knot
			int kPrime = k - 1;

			PiecewiseLinearModel modelPrime = model;

			Knot removedKnot;
			modelPrime.removeRandomKnot(removedKnot);

			double logPriorPrime = dbinom(kPrime, priorMaxKnots, priorKnotsProb, true);
			double logPrior = dbinom(k, priorMaxKnots, priorKnotsProb, true);

			double logProposalToPrime = log(1.0 / k);
			double logProposalFromPrime = dunif(removedKnot.x, minX[0], maxX[0], true) + dnorm(removedKnot.y, modelPrime.at(removedKnot.x), proposalKnotSigma, true);

			double logLikPrime = normalLogLik(modelPrime, sigmaError, data);
			double logLik = normalLogLik(model, sigmaError, data);

			double logNumerator = logLikPrime + logPriorPrime + logProposalFromPrime;
			double logDenomenator = logLik + logPrior + logProposalToPrime;

			if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
				k = kPrime;
				model = modelPrime;
				++acceptDeath;
			}
		} else {
			//move of an existing knot y coords
			PiecewiseLinearModel modelPrime = model;

			Knot oldKnot, newKnot;
			modelPrime.moveRandomKnot(oldKnot, newKnot, proposalKnotSigma);

			//TODO: once data is normalised, will use N(0, sigma^2) prior on heights of knots
			double logPriorPrime = 0.0;
			double logPrior = 0.0;

			double logLikPrime = normalLogLik(modelPrime, sigmaError, data);
			double logLik = normalLogLik(model, sigmaError, data);

			//symmetric proposal dist.
			double logNumerator = logLikPrime + logPriorPrime;
			double logDenomenator = logLik + logPrior;

			if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
				model = modelPrime;
				++acceptMove;
			}
		}

		//secondly move sigmaError
		double sigmaErrorPrime = sigmaError + rnorm(0.0, proposalErrorSigma);

		double logPriorPrime = dgamma(sigmaErrorPrime, priorErrorAlpha, priorErrorBeta, true);
		double logPrior = dgamma(sigmaError, priorErrorAlpha, priorErrorBeta, true);

		double logLikPrime = normalLogLik(model, sigmaErrorPrime, data);
		double logLik = normalLogLik(model, sigmaError, data);

		//symmetric proposal dist.
		double logNumerator = logLikPrime + logPriorPrime;
		double logDenomenator = logLik + logPrior;

		if (logNumerator >= logDenomenator || arma::randu() < exp(logNumerator - logDenomenator)) {
			sigmaError = sigmaErrorPrime;
			++acceptSigma;
		}

		//store current chain values
		kHistory[i] = k;
		knotsHistory[i] = model.getKnots();
		sigmaErrorHistory[i] = sigmaError;
	}

	std::cout << "number of iterations: " << nIterations << std::endl;
	std::cout << "Births: " << acceptBirth << std::endl;
	std::cout << "Deaths: " << acceptDeath << std::endl;
	std::cout << "Moves: " << acceptMove << std::endl;
	std::cout << "Sigmas: " << acceptSigma << std::endl;

	ofstream kHistoryFile, knotXHistoryFile, knotYHistoryFile, sigmaErrorFile;

	kHistoryFile.open("kHistory.txt");
	knotXHistoryFile.open("knotXHistory.txt");
	knotYHistoryFile.open("knotYHistory.txt");
	sigmaErrorFile.open("sigmaErrorHistory.txt");

	for (int i = 0; i < nIterations; ++i) {
		kHistoryFile << kHistory[i] << endl;
		sigmaErrorFile << sigmaErrorHistory[i] << endl;
		for (int j = 0; j < kHistory[i]; ++j) {
			knotXHistoryFile << knotsHistory[i][j].x << ",";
			knotYHistoryFile << knotsHistory[i][j].y << ",";
		}
		knotXHistoryFile << endl;
		knotYHistoryFile << endl;
	}

	kHistoryFile.close();
	knotXHistoryFile.close();
	knotYHistoryFile.close();
	sigmaErrorFile.close();
}

#endif /* RJMCMC_HPP_ */
