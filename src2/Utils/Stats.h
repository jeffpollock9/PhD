#ifndef STATS_H_
#define STATS_H_

#include "ScoresMatrix.h"
#include <armadillo>

arma::mat rMultiNorm(const int n, const arma::rowvec &mu, const arma::mat &cov) {

	arma::mat M;

	try {
		M = arma::randn(n, mu.n_elem) * arma::chol(cov);
		M.each_row() += mu;
	} catch (const std::runtime_error &e) {
		std::cout << e.what() << mu << std::endl;
		std::cout << "cov\n: " << cov << std::endl;
	}

	return M;
}

arma::vec sample(const int n, const arma::vec &x, const arma::vec &weights) {

	arma::vec sample(n);
	double u, sum;

	for (int i = 0; i < n; ++i) {
		u = arma::randu();
		sum = 0.0;
		for (int j = 0; j < n; ++j) {
			sum += weights[j];
			if (u < sum) {
				sample[i] = x[j];
				break;
			}
		}
	}

	return sample;
}

arma::colvec rMultinomial(const int nSamples, const arma::colvec &weights) {

	const int n = weights.n_rows;
	arma::colvec counts(n);
	counts.fill(0);
	double u, sum;

	for (int i = 0; i < nSamples; ++i) {
		u = arma::randu();
		sum = 0.0;
		for (int j = 0; j < n; ++j) {
			sum += weights[j];
			if (u < sum) {
				++counts[j];
				break;
			}
		}
	}
	return counts;
}

arma::mat residualResampleRows(const arma::mat &x, const arma::colvec &weights) {

	const int nParticles = x.n_rows;
	arma::mat sample(x.n_rows, x.n_cols);

	arma::colvec expectedCounts = nParticles * weights;
	arma::colvec floorExpectedCounts = arma::floor(expectedCounts);

	int nResidualSamples = nParticles - arma::sum(floorExpectedCounts);

	arma::colvec residualWeights = (expectedCounts - floorExpectedCounts) / nResidualSamples;
	arma::colvec counts = floorExpectedCounts + rMultinomial(nResidualSamples, residualWeights);

	int row = 0;
	for (int i = 0; i < nParticles; ++i) {
		for (int j = 0; j < counts[i]; ++j) {
			sample.row(row) = x.row(i);
			++row;
		}
	}

	return sample;
}

arma::uvec residualResampleRowsAux(const arma::mat &x, const arma::colvec &weights) {

	const int nParticles = x.n_rows;
	arma::uvec sample(x.n_rows);

	arma::colvec expectedCounts = nParticles * weights;
	arma::colvec floorExpectedCounts = arma::floor(expectedCounts);

	int nResidualSamples = nParticles - arma::sum(floorExpectedCounts);

	arma::colvec residualWeights = (expectedCounts - floorExpectedCounts) / nResidualSamples;
	arma::colvec counts = floorExpectedCounts + rMultinomial(nResidualSamples, residualWeights);

	int row = 0;
	for (int i = 0; i < nParticles; ++i) {
		for (int j = 0; j < counts[i]; ++j) {
			sample.at(row) = i;
			++row;
		}
	}

	return sample;
}

arma::mat sampleRows(const int n, const arma::mat &x, const arma::colvec &weights) {

	arma::mat sample(x.n_rows, x.n_cols);
	double u, sum;

	for (int i = 0; i < n; ++i) {
		u = arma::randu();
		sum = 0.0;
		for (int j = 0; j < n; ++j) {
			sum += weights[j];
			if (u < sum) {
				sample.row(i) = x.row(j);
				break;
			}
		}
	}

	return sample;
}

double dnorm(const double &x, const double &mu, const double &sigma, const bool &isLog = false) {

	static const double oneOverRootTwoPi = 0.3989422804014327;
	double z = (x - mu) / sigma;

	if (isLog) {
		return log(oneOverRootTwoPi / sigma) - 0.5 * z * z;
	} else {
		return oneOverRootTwoPi / sigma * exp(-0.5 * z * z);
	}
}

double dunif_01(const double &x, const bool &isLog = false) {

	if (isLog) {
		return x >= 0.0 && x <= 1.0 ? 0.0 : -HUGE_VAL;
	} else {
		return x >= 0.0 && x <= 1.0 ? 1.0 : 0.0;
	}
}

double dexp(const double &x, const double &lambda, const bool &isLog = false) {

	if (isLog) {
		return x > 0.0 ? log(lambda) - lambda * x : -HUGE_VAL;
	} else {
		return x > 0.0 ? lambda * exp(-lambda * x) : 0.0;
	}
}

double dgamma(const double &x, const double &alpha, const double &beta, const bool &isLog = false) {

	if (isLog) {
		return x >= 0.0 ? alpha * log(beta) - lgamma(alpha) + (alpha - 1.0) * log(x) - beta * x: -HUGE_VAL;
	} else {
		return x >= 0.0 ? pow(beta, alpha) * pow(x, alpha - 1.0) * exp(-beta * x) / tgamma(alpha) : 0.0;
	}
}

double dbeta(const double &x, const double &alpha, const double &beta, const bool &isLog = false) {

	if (isLog) {
		return x >= 0.0 && x <= 1.0 ? lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
				(alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x) : -HUGE_VAL;
	} else {
		return x >= 0.0 && x <= 1.0 ? tgamma(alpha + beta) / (tgamma(alpha) * tgamma(beta)) *
				pow(x, alpha - 1.0) * pow(1.0 - x, beta - 1.0) : 0.0;
	}
}

double dbinom(const double &x, const double &alpha, const double &beta, const bool &isLog = false) {

	if (isLog) {
		return x >= 0.0 && x <= 1.0 ? lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
				(alpha - 1.0) * log(x) + (beta - 1.0) * log(1.0 - x) : -HUGE_VAL;
	} else {
		return x >= 0.0 && x <= 1.0 ? tgamma(alpha + beta) / (tgamma(alpha) * tgamma(beta)) *
				pow(x, alpha - 1.0) * pow(1.0 - x, beta - 1.0) : 0.0;
	}
}

double plogis(const double x) {
	return 1.0 / (1.0 + exp(-x));
}

double logit(const double x) {
	return log(x / (1.0 - x));
}

static const double values[15][15] = { { 27, 20, 21, 7, 3, 3, 0, 0, 0, 0, 0, 0,
		0, 0, 0 }, { 33, 45, 34, 8, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0 }, { 30, 30,
		14, 11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 21, 23, 7, 5, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0 }, { 9, 4, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 4, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 1, 1, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0 }, { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };

ScoresMatrix getObservedScores() {
	ScoresMatrix observedValues;

	for (int i = 0; i < 15; ++i) {
		for (int j = 0; j < 15; ++j) {
			observedValues.at(i, j) = values[i][j];
		}
	}

	return observedValues;
}

class ChiSqStat {
public:
	ChiSqStat(double chi, int df) {
		setChi(chi);
		setDf(df);
	}
	ChiSqStat() {
		setChi(-1);
		setDf(-1);
	}
	double getChi() const {
		return chi;
	}

	void setChi(double chi) {
		this->chi = chi;
	}

	int getDf() const {
		return df;
	}

	void setDf(int df) {
		this->df = df;
	}

private:
	double chi;
	int df;
};

ChiSqStat chiSqTest(ScoresMatrix expected, ScoresMatrix observed) {

	double sum = 0.0, expectedOther = 0.0, observedOther = 0.0, diff;
	int df = 0, dim = observed.getDim();

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			if (expected.at(i, j) >= 5.0) {
				diff = observed.at(i, j) - expected.at(i, j);
				sum += diff * diff / expected.at(i, j);
				df++;
			} else {
				expectedOther += expected.at(i, j);
				observedOther += observed.at(i, j);
			}
		}
	}

	double otherDiff = observedOther - expectedOther;
	double chi = sum + otherDiff * otherDiff / expectedOther;

	return ChiSqStat(chi, df);
}

#endif /* STATS_H_ */
