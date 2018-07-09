#ifndef SIMULATEDPARTICLEFILTER_H_
#define SIMULATEDPARTICLEFILTER_H_

#include <armadillo>
#include "likelihood9_sim.h"

inline arma::rowvec matchBetting(const arma::mat &particles, const int homeTeam, const int awayTeam);
inline double metric(const arma::mat &particles, const int t);

//outputs .csv files for each time point.
//each row in each .csv is a particle of theta.
inline void filter() {

	arma::rowvec sigmasq;// = arma_seq(50, 0.0001, 0.5);
	sigmasq << 0.0613122 << arma::endr;
	int nSigmas = sigmasq.n_cols;

	arma::rowvec deltas; // = arma_seq(10, 0.8, 0.99);
	deltas << 0.99 << arma::endr;
	int nDeltas = deltas.n_cols;

	int nParticles = 100000;
	arma::mat particles(nParticles, 30);

	arma::mat metrics(nSigmas + 1, nDeltas + 1);
	metrics.fill(0.0);

	double alpha, h2;

	for (int i = 0; i < nSigmas; ++i) {
		metrics.at(i + 1, 0) = sigmasq.at(i);
	}

	for (int j = 0; j < nDeltas; ++j) {
		metrics.at(0, j + 1) = deltas.at(j);
	}

	arma::mat m(nParticles, 30);
	arma::uvec aux(nParticles);
	arma::colvec w1(nParticles), w2(nParticles);

	//prior samples
	particles.col(0) = 0.4 + 0.5 * arma::randn(nParticles, 1);
	particles.col(1) = 0.1 + 0.5 * arma::randn(nParticles, 1);

	//#double c[3] = { 0.7, 0.4, 0.6 };
	//#double d[3] = { 0.9, 0.5, 0.3 };
	particles.col(2) = 0.6 + 0.2 * arma::randu(nParticles, 1);
	particles.col(3) = 0.3 + 0.2 * arma::randu(nParticles, 1);
	particles.col(4) = 0.5 + 0.2 * arma::randu(nParticles, 1);

	particles.col(5) = 0.8 + 0.2 * arma::randu(nParticles, 1);
	particles.col(6) = 0.4 + 0.2 * arma::randu(nParticles, 1);
	particles.col(7) = 0.2 + 0.2 * arma::randu(nParticles, 1);

	//rho
	particles.col(8) = 1.0 + 0.5 * arma::randn(nParticles, 1);
	particles.col(9) = 1.5 + 0.5 * arma::randn(nParticles, 1);

	//R
	//TODO: change me when new simulation is run
	double Rprior[20] = { -0.336698, -1.01648, 1.01182, 0.87005, -0.594951, 0.664923, -0.875045, 0.0949707,
			-0.404757, -0.602735, -0.88483, -0.486396, 0.670291, 1.04919, 0.999636, 0.977015, -0.88648,
			-0.714131, 0.485152, -0.345473 };

	for (int i = 10; i < 30; ++i) {
		particles.col(i) = Rprior[i - 10] + 0.5 * arma::randn(nParticles, 1);
	}

	arma::mat prior = particles;

	const bool makeMetricCSV = nSigmas > 1 || nDeltas > 1;
	const bool makeParticleCSV = !makeMetricCSV;

	for (int s = 0; s < nSigmas; ++s) {
		std::cout << "sigma: " << s + 1 << " of " << nSigmas << std::endl;
		for (int d = 0; d < nDeltas; ++d) {
			std::cout << "\tdelta: " << d + 1 << " of " << nDeltas << std::endl;

			particles = prior;

			if (makeParticleCSV) {
				writeCSV(particles, fileNameSim(0));
			}

			alpha = (3.0 * deltas[d] - 1.0) / (2.0 * deltas[d]);
			h2 = 1 - alpha * alpha;

			for (int t = 1; t <= 38; ++t) {
				if (t % 5 == 0) {
					std::cout << "\t\tweek: " << t << std::endl;
				}

				//forecast games before updating on that data
				if (makeMetricCSV) {
					metrics.at(s + 1, d + 1) += metric(particles, t);
				}

				//kernel means of first 10 parameters (static)
				for (int i = 0; i < 10; ++i) {
					m.col(i) = alpha * particles.col(i) + (1.0 - alpha) * arma::mean(particles.col(i));
				}

				//random walk has mean 0 so particles not expected to move
				for (int i = 10; i < 30; ++i) {
					m.col(i) = particles.col(i);
				}

				//first stage weights based on aux locations
				#pragma omp parallel for
				for (int i = 0; i < nParticles; ++i) {
					w1[i] = logLikelihood9_sim(m.row(i), t);
				}
				w1 = exp(w1 - w1.max());
				w1 /= arma::sum(w1);

				//re-sample aux integers based on first stage weights
				aux = residualResampleRowsAux(m, w1);

				m = m.rows(aux);

				// artificial evolution of first 10 parameters
				arma::mat temp = particles.rows(aux);
				arma::mat cov = h2 * arma::cov(temp.cols(0, 9));

				arma::rowvec mu(10);
				mu.fill(0.0);
				particles.cols(0, 9) = m.cols(0, 9) + rMultiNorm(nParticles, mu, cov);

				// evolution of team parameters
				for (int i = 10; i < 30; ++i) {
					particles.col(i) = m.col(i) + arma::randn(nParticles) * sqrt(sigmasq[s]);
				}

				//second stage weights based on newly moved particle locations
				#pragma omp parallel for
				for (int i = 0; i < nParticles; ++i) {
					w2[i] = logLikelihood9_sim(particles.row(i), t) - logLikelihood9_sim(m.row(i), t);
				}
				w2 = exp(w2 - w2.max());
				w2 /= arma::sum(w2);

				//re-sample based on second stage weights
				particles = residualResampleRows(particles, w2);
				if (makeParticleCSV) {
					writeCSV(particles, fileNameSim(t));
				}
			}
			if (makeMetricCSV) {
				std::cout << metrics << std::endl;
			}
		}
	}
	if (makeMetricCSV) {
		writeCSV(metrics, "metrics-simulated-openmp.csv");
	}
}

inline arma::rowvec matchBetting(const arma::mat &particles, const int homeTeam, const int awayTeam) {

	double c[3] = { 0.0, 0.0, 0.0 };
	double d[3] = { 0.0, 0.0, 0.0 };
	double rho[2] = { 0.0, 0.0 };

	Model9 model(0.0, 0.0, c, d, rho);
	ScoreLine score;
	int nParticles = particles.n_rows;
	double home = 0.0, away = 0.0, draw = 0.0;

	for (int i = 0; i < nParticles; ++i) {

		model.setH(particles.at(i, 0));
		model.setA(particles.at(i, 1));

		model.setC(particles.at(i, 2), particles.at(i, 3), particles.at(i, 4));
		model.setD(particles.at(i, 5), particles.at(i, 6), particles.at(i, 7));
		model.setRho(particles.at(i, 8), particles.at(i, 9));

		model.setHomeR(exp(particles.at(i, homeTeam + 10)));
		model.setAwayR(exp(particles.at(i, awayTeam + 10)));

		score = model.simulateScoreLine();

		if (score.getHomeScore() > score.getAwayScore()) {
			++home;
		} else if (score.getHomeScore() == score.getAwayScore()) {
			++draw;
		} else {
			++away;
		}
	}

	home /= nParticles;
	draw /= nParticles;
	away /= nParticles;

	arma::rowvec out(3);
	out.at(0) = home;
	out.at(1) = draw;
	out.at(2) = away;

	return out;
}

inline double metric(const arma::mat &particles, const int t) {

	int gameIdx;
	double metric = 0.0;
	arma::rowvec probs;

	for (int game = 0; game < 10; ++game) {
		gameIdx = (t - 1) * 10 + game;

		probs = matchBetting(particles, MATCH_HOME_ID_SIM[gameIdx], MATCH_AWAY_ID_SIM[gameIdx]);

		if (MATCH_HOME_SCORE_SIM[gameIdx] > MATCH_AWAY_SCORE_SIM[gameIdx]) {
			metric += log(probs[0]);
		} else if (MATCH_HOME_SCORE_SIM[gameIdx] == MATCH_AWAY_SCORE_SIM[gameIdx]) {
			metric += log(probs[1]);
		} else {
			metric += log(probs[2]);
		}
	}

	return metric;
}

#endif /* SIMULATEDPARTICLEFILTER_H_ */
