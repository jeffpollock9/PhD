#ifndef STATS_HPP_
#define STATS_HPP_

#include <armadillo>

double runif(const double a, const double b) {
	return a + (b - a) * arma::randu();
}

double rnorm(const double mu, const double sigma) {
	return mu + sigma * arma::randn();
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

double dnorm(const double x, const double &mu, const double sigma, const bool isLog = false) {

	static const double oneOverRootTwoPi = 0.3989422804014327;
	double z = (x - mu) / sigma;

	if (isLog) {
		return log(oneOverRootTwoPi / sigma) - 0.5 * z * z;
	} else {
		return oneOverRootTwoPi / sigma * exp(-0.5 * z * z);
	}
}

double dunif(const double x, const double a, const double b, const bool isLog = false) {

	if (isLog) {
		return x >= a && x <= b ? log(1.0 / (b - a)) : -HUGE_VAL;
	} else {
		return x >= a && x <= b ? 1.0 / (b - a) : 0.0;
	}
}

double dexp(const double x, const double lambda, const bool isLog = false) {

	if (isLog) {
		return x > 0.0 ? log(lambda) - lambda * x : -HUGE_VAL;
	} else {
		return x > 0.0 ? lambda * exp(-lambda * x) : 0.0;
	}
}

double dgamma(const double x, const double alpha, const double beta, const bool isLog = false) {

	if (isLog) {
		return x >= 0.0 ? alpha * log(beta) - lgamma(alpha) + (alpha - 1.0) * log(x) - beta * x : -HUGE_VAL;
	} else {
		return x >= 0.0 ? pow(beta, alpha) * pow(x, alpha - 1.0) * exp(-beta * x) / tgamma(alpha) : 0.0;
	}
}

double dgammaMeanVar(const double x, const double mu, const double sigma2, const bool isLog = false) {

	double alpha = mu * mu / sigma2;
	double beta = alpha / mu;

	if (isLog) {
		return x >= 0.0 ? alpha * log(beta) - lgamma(alpha) + (alpha - 1.0) * log(x) - beta * x : -HUGE_VAL;
	} else {
		return x >= 0.0 ? pow(beta, alpha) * pow(x, alpha - 1.0) * exp(-beta * x) / tgamma(alpha) : 0.0;
	}
}

double dbeta(const double x, const double alpha, const double beta, const bool isLog = false) {

	if (isLog) {
		return x >= 0.0 && x <= 1.0 ?
				lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) + (alpha - 1.0) * log(x)
						+ (beta - 1.0) * log(1.0 - x) :
				-HUGE_VAL;
	} else {
		return x >= 0.0 && x <= 1.0 ?
				tgamma(alpha + beta) / (tgamma(alpha) * tgamma(beta)) * pow(x, alpha - 1.0)
						* pow(1.0 - x, beta - 1.0) :
				0.0;
	}
}

int nChooseK(const int n, int k) {
    if (k > n) {
    	return 0;
    }
    if (k * 2 > n) {
    	k = n - k;
    }
    if (k == 0) {
    	return 1;
    }
    int result = n;
    for(int i = 2; i <= k; ++i ) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

double dbinom(const int x, const int n, const double p, const bool isLog = false) {

	if (p > 1.0 || p < 0.0) {
		throw runtime_error("p not in [0, 1] in dbinom");
	}

	if (isLog) {
		return x >= 0 && x <= n ? log(nChooseK(n, x)) + x * log(p) + (n - x) * log(1.0 - p) : -HUGE_VAL;
	} else {
		return x >= 0 && x <= n ? nChooseK(n, x) * pow(p, x) * pow(1.0 - p, n - x) : 0;
	}
}

double logistic(const double x) {
	return 1.0 / (1.0 + exp(-x));
}

double logit(const double x) {
	return log(x / (1.0 - x));
}

#endif /* STATS_HPP_ */
