#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <armadillo>
#include "../Utils/Stats.h"

static const double DELTA = 0.95;
static const double ALPHA = (3.0 * DELTA - 1.0) / (2.0 * DELTA);
static const double H2 = 1 - ALPHA * ALPHA;

double normalLogLikelihood(const arma::vec &y, const double mu, const double sigma) {

	int n = y.n_elem;
	double sum = 0.0;

	for (int i = 0; i < n; ++i) {
		sum += dnorm(y[i], mu, sigma, true);
	}

	return sum;
}

// example
// prior for mu is N(5, 1^2)
// 2 data points from N(1, 2^2)
// sigma = 2^2 is known
// infer mu = 1 unknown
std::pair<arma::mat, arma::vec> exampleFilter() {

	int T = 25, n = 50000;
	double sigma = 2.0, var;

	arma::mat x(n, T + 1);
	arma::vec y(2), yT(2 * T), kernelMeans(n), particles(n),
			w1(n), w2(n), aux(n);

	// prior samples at t = 0
	particles = arma::randn(n) + 5.0;
	x.col(0) = particles;

	// go through time updating particles
	for (int t = 1; t <= T; ++t) {
		// two random data points
		y = arma::randn(2) * 2.0 + 1.0;
		yT[t * 2 - 2] = y[0];
		yT[t * 2 - 1] = y[1];

		kernelMeans = ALPHA * particles + (1.0 - ALPHA) * arma::mean(particles);

		// first stage weights
		for (int i = 0; i < n; ++i) {
			w1[i] = normalLogLikelihood(y, kernelMeans[i], sigma);
		}
		w1 = exp(w1 - w1.max());
		w1 /= arma::sum(w1);

		// artificial evolution of aux particles
		aux = sample(n, kernelMeans, w1);
		var = H2 * arma::var(x.col(t-1));
		particles = arma::randn(n) * sqrt(var) + aux;

		// second stage weights
		for (int i = 0; i < n; ++i) {
			w2(i) = normalLogLikelihood(y, particles[i], sigma) -
					normalLogLikelihood(y, aux[i], sigma);
		}
		w2 = exp(w2 - w2.max());
		w2 /= arma::sum(w2);
		particles = sample(n, particles, w2);
		x.col(t) = particles;
	}

	return std::pair<arma::mat, arma::vec>(x, yT);
}


#endif /* PARTICLE_H_ */
