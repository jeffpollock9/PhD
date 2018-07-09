#include "Utils/csv.h"
#include "Utils/Timer.h"

//#define check
#ifdef check
int main(int argc, char **argv) {
	double c[3] = {0.7, 0.4, 0.6};
	double d[3] = {0.9, 0.5, 0.3};
	double rho[2] = {1.0, 1.5};
	double U[8] = {5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.25, 0.01};
	double A = 0.5;
	double B = 10.0;

	Model10 model = Model10(0.4, 0.1, c, d, U, A, B, rho);

	for (int p = 1; p <= 20; ++p) {
		std::cout << "position = " << p << ", utility = " << model.utility(p) << std::endl;
	}

}
#endif

#define realfilter
#ifdef realfilter
#include "ParticleFilter/FootballParticleFilter.h"
int main(int argc, char **argv) {
//	INIT_TIMER
	filter();
//	STOP_TIMER("filter")
}
#endif

//#define simFilter
#ifdef simFilter
#include "ParticleFilter/SimulatedParticleFilter.h"
int main(int argc, char **argv) {
	filter();
}
#endif

//#define simLeague
#ifdef simLeague
#include "Models/Model9.h"
int main(int argc, char **argv) {

	arma::arma_rng::set_seed_random();

	double c[3] = {0.7, 0.4, 0.6};
	double d[3] = {0.9, 0.5, 0.3};
	double rho[2] = {1.0, 1.5};

	Model9 model(0.4, 0.1, c, d, rho);
	model.simulateLeagueData();

}
#endif

//#define getGamma
//#define bayesianBT
#ifdef bayesianBT
#include "Bayesian/simulateBT.h"
#include <ctime>
int main(int argc, char **argv) {
	clock_t startTime = clock();

	static const int nIterationsMCMC = 45000;
	static const int burn = 5000;
	static const int thin = 1;

	double rGamma1[] = {log(1), log(3.5), log(5)}; //max @ log(3.5), log(1.5), -391.849
	double rGamma2[] = {log(1), log(1.5), log(5)};
	static const int nRGamma1 = sizeof(rGamma1) / sizeof(double);
	static const int nRGamma2 = sizeof(rGamma2) / sizeof(double);

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParametersBT);
	DoubleMatrix metrics = createDoubleMatrix(nRGamma1, nRGamma2);

	int gameIdx;
	double metric;

	DoubleVector probs(3);
	for (int i = 0; i < nRGamma1; ++i) {
		for (int j = 0; j < nRGamma2; ++j) {
			std::cout << "gamma_1: " << rGamma1[i] << " gamma_2: " << rGamma2[j] << std::endl;
			metric = 0.0;
			for (int week = 0; week < NUM_WEEKS; ++week) {
				std::cout << "\t week: " << week << std::endl;
				DoubleVector theta(parametersBT, parametersBT + nParametersBT);
				updateHistory(history, nIterationsMCMC, theta, week, rGamma1[i],
						rGamma2[j]);

				for (int game = 0; game < 10; ++game) {
					gameIdx = week * 10 + game;

					probs = getProbsFromHistory(history, MATCH_HOME_ID[gameIdx], MATCH_AWAY_ID[gameIdx], burn, thin);

					if (MATCH_HOME_SCORE[gameIdx] > MATCH_AWAY_SCORE[gameIdx]) {
						metric += log(probs[2]);
					} else if (MATCH_HOME_SCORE[gameIdx] == MATCH_AWAY_SCORE[gameIdx]) {
						metric += log(probs[1]);
					} else {
						metric += log(probs[0]);
					}

				}
			}
			metrics[i][j] = metric;
		}
	}
	writeCSV(metrics, "metricsBT.csv");

	std::cout << "time: " << (double) (clock() - startTime) / CLOCKS_PER_SEC
	<< " seconds" << std::endl;
}
#endif

//#define bayesianBTprobs
#ifdef bayesianBTprobs
#include "Bayesian/simulateBT.h"
#include <ctime>
int main(int argc, char **argv) {
	clock_t startTime = clock();

	static const int nIterationsMCMC = 45000;
	static const int burn = 5000;
	static const int thin = 1;
	static const double bestGamma1 = log(3.5);
	static const double bestGamma2 = log(1.5);
	static const int startModelWeek = 0;

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParametersBT);
	DoubleMatrix modelProbs = createDoubleMatrix(10 * (NUM_WEEKS - startModelWeek), 3);

	int gameIndex, betIndex;

	DoubleVector probs(3);
	for (int week = startModelWeek; week < NUM_WEEKS; ++week) {
		std::cout << "\t week: " << week << std::endl;
		DoubleVector theta(parametersBT, parametersBT + nParametersBT);
		updateHistory(history, nIterationsMCMC, theta, week, bestGamma1,
				bestGamma2);

		for (int game = 0; game < 10; ++game) {
			gameIndex = week * 10 + game;
			betIndex = gameIndex - startModelWeek * 10;

			probs = getProbsFromHistory(history, MATCH_HOME_ID[gameIndex], MATCH_AWAY_ID[gameIndex], burn, thin);

			modelProbs[betIndex][0] = probs[2];
			modelProbs[betIndex][1] = probs[1];
			modelProbs[betIndex][2] = probs[0];

		}
	}
	writeCSV(modelProbs, "modelProbsBTBayesian.csv");

	std::cout << "time: " << (double) (clock() - startTime) / CLOCKS_PER_SEC
	<< " seconds" << std::endl;
}
#endif

//found maximum at f -358.4701804
//#define testModelBT
#ifdef testModelBT
#include "Optimisation/likelihoodBT.h"
int main(int argc, char **argv) {

	DoubleVector theta(parametersBT, parametersBT + nParametersBT);
	nlopt::opt opt(nlopt::LN_SBPLX, nParametersBT);

	Data d;
	d.constrained = true;
	d.team = -1;
	d.upToWeek = 38;
	d.printTrace = false;

	const int startModelWeek = 1;
	const int nIterations = 100000;

	double f, rSum;
	opt.set_max_objective(logLikelihoodBT, &d);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(1000000);

	DoubleMatrix modelProbs = createDoubleMatrix(10 * (NUM_WEEKS - startModelWeek), 5);
	DoubleMatrix modelParameters = createDoubleMatrix(NUM_WEEKS - startModelWeek + 1, nParametersBT);

	int gameIndex, betIndex;
	DoubleVector R(20);
	ModelBT model(0.0, 0.0);

	for (int w = startModelWeek; w <= NUM_WEEKS; ++w) {
		d.upToWeek = w;
		DoubleVector theta(parametersBT, parametersBT + nParametersBT);

		opt.optimize(theta, f);

		modelParameters[w - startModelWeek] = theta;

		std::cout << "||------------- WEEK: " << w << " -------------||" << std::endl;
		std::cout << "result: " << opt.last_optimize_result() << std::endl;
		std::cout << "log-likelihood: " << f << std::endl;
		std::cout << "theta: " << std::endl;
		for (int i = 0; i < nParametersBT; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << std::endl << std::endl;

		rSum = 0.0;
		for (int i = 0; i < 19; ++i) {
			R[i] = theta[2 + i];
			rSum += R[i];
		}
		R[19] = -rSum;

		model.setH(theta[0]);
		model.setDelta(theta[1]);

		if (w < NUM_WEEKS) {
			for (int g = 0; g < 10; ++g) {
				gameIndex = w * 10 + g;
				betIndex = gameIndex - startModelWeek * 10;

				model.setHomeR(R[MATCH_HOME_ID[gameIndex]]);
				model.setAwayR(R[MATCH_AWAY_ID[gameIndex]]);

				modelProbs[betIndex][0] = model.getHome();
				modelProbs[betIndex][1] = model.getDraw();
				modelProbs[betIndex][2] = model.getAway();

				modelProbs[betIndex][3] = MATCH_HOME_ID[gameIndex];
				modelProbs[betIndex][4] = MATCH_AWAY_ID[gameIndex];
			}
		}
	}
	writeCSV(modelProbs, "modelProbsBT.csv");
	writeCSV(modelParameters, "modelParametersBT.csv");
}
#endif

//#define testModelBT_MLE
//#define getGamma
#ifdef testModelBT_MLE
#include "Optimisation/likelihoodBT.h"
int main(int argc, char **argv) {

	DoubleVector theta(/*parametersBT, parametersBT + */nParametersBT);
	nlopt::opt opt(nlopt::LN_SBPLX, nParametersBT);

	Data d;
	d.constrained = true;
	d.team = -1;
	d.upToWeek = 38;
	d.printTrace = true;

	double f;
	opt.set_max_objective(logLikelihoodBT, &d);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(1000000);

	opt.optimize(theta, f);

	printf("result = %x\n", opt.last_optimize_result());
	printf("found maximum at f %0.10g\n\n", f);

	DoubleVector R(20);
	double rSum = 0.0;
	for (int i = 0; i < 19; ++i) {
		R[i] = theta[2 + i];
		rSum += R[i];
	}
	R[19] = -rSum;

	for (int t = 0; t < 20; ++t) {
		std::cout << TeamName[t] << ": " << R[t] << std::endl;
	}

}
#endif
//found maximum at f -518.053668
//#define DRbetThroughSeason
#ifdef DRbetThroughSeason
#include "Optimisation/likelihoodDR.h"
int main(int argc, char **argv) {

	const int startModelWeek = 5;
	const int nIterations = 100000;

	Data d;
	d.constrained = true;
	d.team = -1;
	d.printTrace = false;

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, nParametersDR);

	opt.set_max_objective(logLikelihoodDR, &d);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(1000000);

	DoubleMatrix modelProbs = createDoubleMatrix(10 * (NUM_WEEKS_2011_2012 - startModelWeek), 5);
	DoubleMatrix modelParameters = createDoubleMatrix(NUM_WEEKS_2011_2012 - startModelWeek + 1, nParametersDR);

	ScoresMatrix scores;
	int gameIndex, betIndex;

	for (int w = startModelWeek; w <= NUM_WEEKS_2011_2012; ++w) {
		d.upToWeek = w;
		DoubleVector theta(parametersDR, parametersDR + nParametersDR);

		opt.optimize(theta, f);

		modelParameters[w - startModelWeek] = theta;

		std::cout << "||------------- WEEK: " << w << " -------------||" << std::endl;
		std::cout << "result: " << opt.last_optimize_result() << std::endl;
		std::cout << "log-likelihood: " << f << std::endl;
		std::cout << "theta: " << std::endl;
		for (int i = 0; i < nParametersDR; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << std::endl << std::endl;

		if (w < NUM_WEEKS_2011_2012) {
			for (int g = 0; g < 10; ++g) {
				gameIndex = w * 10 + g;
				betIndex = gameIndex - startModelWeek * 10;

				scores = simulateScoresMatrixFromTheta(nIterations, theta,
						MATCH_HOME_ID_2011_2012[gameIndex], MATCH_AWAY_ID_2011_2012[gameIndex]);

				modelProbs[betIndex][0] = scores.getHomeWin();
				modelProbs[betIndex][1] = scores.getDraw();
				modelProbs[betIndex][2] = scores.getAwayWin();

				modelProbs[betIndex][3] = MATCH_HOME_ID_2011_2012[gameIndex];
				modelProbs[betIndex][4] = MATCH_AWAY_ID_2011_2012[gameIndex];
			}
		}
	}
	writeCSV(modelProbs, "modelProbsDR.csv");
	writeCSV(modelParameters, "modelParametersDR.csv");
}
#endif

//found maximum at f -531.2389234
//#define testModel10
#ifdef testModel10
#include "Optimisation/likelihood10.h"
int main(int argc, char **argv) {

	DoubleVector theta(parameters10, parameters10 + nParameters10);

	double l = logLikelihood10(theta, true, -1, 38, true);

	//	nlopt::opt opt(nlopt::LN_SBPLX, nParameters10);
	//
	//	Data d;
	//	d.constrained = true;
	//	d.team = -1;
	//	d.upToWeek = 38;
	//	d.printTrace = true;
	//
	//	opt.set_max_objective(logLikelihood10, &d);
	//
	//	opt.set_ftol_rel(1e-10);
	//	opt.set_maxeval(100000);
	//	double f;
	//	opt.optimize(theta, f);
	//
	//	printf("result = %x\n", opt.last_optimize_result());
	//	printf("found maximum at f %0.10g\n\n", f);
	//	printf("theta: ");
	//	for (int i = 0; i < nParameters10; ++i) {
	//		printf("%f, ", theta[i]);
	//	}

}
#endif

//#define RpriorSamples
#ifdef RpriorSamples
#include "Bayesian/inferenceModel9.h"
int main(int argc, char **argv) {

	static const double gamma_1 = 1.0;
	static const double gamma_2 = 0.5;

	int nIterations = 2000000;
	DoubleMatrix history = createDoubleMatrix(nIterations, nParameters9);
	DoubleVector theta(parameters9, parameters9 + nParameters9);
	double l = 0.0;

	for (int i = 0; i < nIterations; ++i) {
		if (i % 10000 == 0) {
			std::cout << i << std::endl;
		}
		for (int r = 10; r < nParameters9; ++r) {
			updateParameter(r, theta, history, i, l, 0, gamma_1, gamma_2);
		}
	}

	int thin = 10;
	DoubleMatrix thinnedParameters = createDoubleMatrix(nIterations / thin, 20);

	for (int i = 0; i < nIterations / thin; ++i) {
		for (int j = 0; j < 20; ++j) {
			thinnedParameters[i][j] = history[i * thin][j + 10];
		}
	}

	writeCSV(thinnedParameters, "priorSamples-1.0-0.5-thinned.csv");
}
#endif

//#define bet
#ifdef bet
#include "Bayesian/simulate9.h"
int main(int argc, char **argv) {
	const int nIterationsMCMC = 150000;
	const int burn = 20000;
	const int thin = 1;
	const double bestGamma1 = log(2.0);
	const double bestGamma2 = log(1.5);

	static const int nGames = 10;

	DoubleMatrix modelProbs = createDoubleMatrix(nGames, 3);
	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);

	ScoresMatrix scores;
	int gameIndex;
	int dataUpToWeek = 20; //TODO: UPDATE ME!!!

	DoubleVector theta(start9, start9 + nParameters9);
	updateHistory(history, nIterationsMCMC, theta, dataUpToWeek, bestGamma1,
			bestGamma2, true);

	int homeTeam[nGames] = {
		HullCity,
		CardiffCity,
		Everton,
		Fulham,
		Southampton,
		TottenhamHotspur,
		ManchesterUnited,
		NewcastleUnited,
		StokeCity,
		AstonVilla};

	int awayTeam[nGames] = {
		Chelsea,
		WestHamUnited,
		NorwichCity,
		Sunderland,
		WestBromwichAlbion,
		CrystalPalace,
		SwanseaCity,
		ManchesterCity,
		Liverpool,
		Arsenal};

	for (int g = 0; g < nGames; ++g) {

		scores = simulateScoresMatrixFromHistory(history,
				homeTeam[g], awayTeam[g], burn, thin);

		modelProbs[g][0] = scores.getHomeWin();
		modelProbs[g][1] = scores.getDraw();
		modelProbs[g][2] = scores.getAwayWin();

	}

	writeCSV(modelProbs, "csv/bet/modelProbsBayesian.csv");
}
#endif

//#define play
#ifdef play
#include "Optimisation/staticModel9.h"
#include "Data/Team.h"
int main(int argc, char **argv) {

	DoubleVector theta(parameters9, parameters9 + nParameters9);

	double c[3] = {theta[2], theta[3], theta[4]};
	double d[3] = {theta[5], theta[6], theta[7]};
	double rho[2] = {theta[8], theta[9]};

	Model9 model(theta[0], theta[1], c, d, rho, false);

	double R[20];

	for (int i = 0; i < 20; ++i) {
		R[i] = theta[i + 10];
	}

	model.setHomeR(R[ManchesterCity]);
	model.setAwayR(R[AstonVilla]);

	int n = 5000000;
	double noGoal = 0.0, homeGoal = 0.0, awayGoal = 0.0;
	double goalTime, homeRate, awayRate;

	for (int i = 0; i < n; ++i) {
		goalTime = model.simulateGoalTime(0.0, 0, 0);

		if (goalTime > 1.0) {
			++noGoal;
		} else {

			homeRate = model.getRate(goalTime, 0, 0, true, false);
			awayRate = model.getRate(goalTime, 0, 0, false, false);

			if (uniform() < homeRate / (homeRate + awayRate)) {
				++homeGoal;
			} else {
				++awayGoal;
			}
		}

	}

	noGoal /= n;
	homeGoal /= n;
	awayGoal /= n;

	std::cout << "simulation:" << std::endl;
	std::cout << "\tno goal: " << noGoal << std::endl;
	std::cout << "\thome goal: " << homeGoal << std::endl;
	std::cout << "\taway goal: " << awayGoal << std::endl;

	//analytic answer

	double homeIntegral = model.getRateIntegral(0.0 / 90.0, 44.0 / 90.0, 0, 0, true)
	+ model.getRateIntegral(44.0 / 90.0, 45.0 / 90.0, 0, 0, true)
	+ model.getRateIntegral(45.0 / 90.0, 89.0 / 90.0, 0, 0, true)
	+ model.getRateIntegral(89.0 / 90.0, 90.0 / 90.0, 0, 0, true);

	double awayIntegral = model.getRateIntegral(0.0 / 90.0, 44.0 / 90.0, 0, 0, false)
	+ model.getRateIntegral(44.0 / 90.0, 45.0 / 90.0, 0, 0, false)
	+ model.getRateIntegral(45.0 / 90.0, 89.0 / 90.0, 0, 0, false)
	+ model.getRateIntegral(89.0 / 90.0, 90.0 / 90.0, 0, 0, false);

	noGoal = exp(-(homeIntegral + awayIntegral));
	homeGoal = (1.0 - noGoal) * homeIntegral / (homeIntegral + awayIntegral);
	awayGoal = (1.0 - noGoal) * awayIntegral / (homeIntegral + awayIntegral);

	std::cout << "analytic:" << std::endl;
	std::cout << "\tno goal: " << noGoal << std::endl;
	std::cout << "\thome goal: " << homeGoal << std::endl;
	std::cout << "\taway goal: " << awayGoal << std::endl;

}
#endif

//#define bayesianBetThroughSeason
#ifdef bayesianBetThroughSeason
#include "Bayesian/simulate9.h"
#include <chrono>
int main(int argc, char **argv) {

//	INIT_TIMER

	const int nIterationsMCMC = 100000;
	const int burn = 0;
	const int thin = 1;
	const double bestGamma1 = log(2.0);
	const double bestGamma2 = log(1.5);

	DoubleMatrix modelProbs = createDoubleMatrix(10 * NUM_WEEKS_2011_2012, 5);
	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);
	DoubleVector theta(start9, start9 + nParameters9);

	ScoresMatrix scores;
	int gameIndex;

	auto start = std::chrono::high_resolution_clock::now();

	for (int w = 1; w <= NUM_WEEKS_2011_2012; ++w) {
//		std::cout << "||------------- WEEK: " << w << " -------------||" << std::endl;
		updateHistory(history, nIterationsMCMC, theta, w, bestGamma1, bestGamma2);

		std::cout << "week" << w << ","
				<< std::chrono::duration_cast<std::chrono::seconds>(
						std::chrono::high_resolution_clock::now() - start).count() << std::endl;

		/*
		 for (int g = 0; g < 10; ++g) {
		 gameIndex = w * 10 + g;

		 scores = simulateScoresMatrixFromHistory(history,
		 MATCH_HOME_ID_2011_2012[gameIndex], MATCH_AWAY_ID_2011_2012[gameIndex], burn, thin);

		 modelProbs[gameIndex][0] = scores.getHomeWin();
		 modelProbs[gameIndex][1] = scores.getDraw();
		 modelProbs[gameIndex][2] = scores.getAwayWin();

		 modelProbs[gameIndex][3] = MATCH_HOME_ID_2011_2012[gameIndex];
		 modelProbs[gameIndex][4] = MATCH_AWAY_ID_2011_2012[gameIndex];
		 }
		 */
	}
//	writeCSV(modelProbs, "modelProbsM.csv");

//	STOP_TIMER("MH updates")
}
#endif

//#define classicalBetThroughSeason
#ifdef classicalBetThroughSeason
#include "Optimisation/staticModel9.h"
int main(int argc, char **argv) {

	const int startModelWeek = 10;
	const int nIterations = 100000;

	Data d;
	d.constrained = true;
	d.team = -1;
	d.printTrace = false;

	double f;
	nlopt::opt opt(nlopt::LN_SBPLX, nParameters9);

	opt.set_max_objective(staticLogLikelihood9, &d);
	opt.set_ftol_rel(1e-12);
	opt.set_maxeval(100000);

	DoubleMatrix modelProbs = createDoubleMatrix(10 * (NUM_WEEKS - startModelWeek), 5);
	DoubleMatrix modelParameters = createDoubleMatrix(NUM_WEEKS - startModelWeek + 1, nParameters9);

	ScoresMatrix scores;
	int gameIndex, betIndex;

	for (int w = startModelWeek; w <= NUM_WEEKS; ++w) {
		d.upToWeek = w;
		DoubleVector theta(parameters9, parameters9 + nParameters9);

		opt.optimize(theta, f);

		modelParameters[w - startModelWeek] = theta;

		std::cout << "||------------- WEEK: " << w << " -------------||" << std::endl;
		std::cout << "result: " << opt.last_optimize_result() << std::endl;
		std::cout << "log-likelihood: " << f << std::endl;
		std::cout << "theta: " << std::endl;
		for (int i = 0; i < nParameters; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << std::endl << std::endl;

		if (w < NUM_WEEKS) {
			for (int g = 0; g < 10; ++g) {
				gameIndex = w * 10 + g;
				betIndex = gameIndex - startModelWeek * 10;

				scores = simulateScoresMatrixFromTheta(nIterations, theta,
						MATCH_HOME_ID[gameIndex], MATCH_AWAY_ID[gameIndex]);

				modelProbs[betIndex][0] = scores.getHomeWin();
				modelProbs[betIndex][1] = scores.getDraw();
				modelProbs[betIndex][2] = scores.getAwayWin();

				modelProbs[betIndex][3] = MATCH_HOME_ID[gameIndex];
				modelProbs[betIndex][4] = MATCH_AWAY_ID[gameIndex];
			}
		}
	}
	writeCSV(modelProbs, "modelProbs.csv");
	writeCSV(modelParameters, "modelParameters.csv");
}
#endif

//#define getGamma
//#define getGammaM
#ifdef getGammaM
#include "Bayesian/inferenceModel9.h"
#include "Bayesian/simulate.h"
int main(int argc, char **argv) {

	clock_t startTime = clock();

	static const int nIterationsMCMC = 45000;
	static const int burn = 5000;
	static const int thin = 1;

	double rGamma1[] = {log(1), log(1.5), log(2), log(2.5), log(3), log(4), log(5), log(10)};
	double rGamma2[] = {log(1), log(1.5), log(2), log(2.5), log(3), log(4), log(5), log(10)};
	//	double rGamma2[] = {log(1), log(2), log(5)}; // max @ log(3.5), log(2): -385.161 _done_

	static const int nRGamma1 = sizeof(rGamma1) / sizeof(double);
	static const int nRGamma2 = sizeof(rGamma2) / sizeof(double);

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);
	DoubleMatrix metrics = createDoubleMatrix(nRGamma1, nRGamma2);

	int gameIdx;
	double metric;
	ScoresMatrix scores;

	for (int i = 0; i < nRGamma1; ++i) {
		for (int j = 0; j < nRGamma2; ++j) {
			std::cout << "gamma_1: " << rGamma1[i] << " gamma_2: " << rGamma2[j] << std::endl;
			metric = 0.0;
			for (int week = 0; week < NUM_WEEKS; ++week) {
				std::cout << "\t week: " << week << std::endl;
				DoubleVector theta(start9, start9 + nParameters9);
				updateHistory(history, nIterationsMCMC, theta, week, rGamma1[i],
						rGamma2[j]);

				for (int game = 0; game < 10; ++game) {
					gameIdx = week * 10 + game;
					scores = simulateScoresMatrixFromHistory(history,
							MATCH_HOME_ID[gameIdx], MATCH_AWAY_ID[gameIdx],
							burn, thin);

					if (MATCH_HOME_SCORE[gameIdx] > MATCH_AWAY_SCORE[gameIdx]) {
						metric += log(scores.getHomeWin());
					} else if (MATCH_HOME_SCORE[gameIdx] == MATCH_AWAY_SCORE[gameIdx]) {
						metric += log(scores.getDraw());
					} else {
						metric += log(scores.getAwayWin());
					}

				}
			}
			metrics[i][j] = metric;
		}
	}
	writeCSV(metrics, "metricsM.csv");

	std::cout << "time: " << (double) (clock() - startTime) / CLOCKS_PER_SEC
	<< " seconds" << std::endl;

}
#endif

//#define bayesian9
#ifdef bayesian9
#include "Bayesian/simulate9.h"
#include "Bayesian/inferenceModel9.h"
int main(int argc, char **argv) {

	const int nIterationsMCMC = 100000;
	const int nWeeks = 38;
	const double bestGamma1 = log(2.0);
	const double bestGamma2 = log(1.5);

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);

	clock_t start = clock();
	DoubleVector theta(parameters9, parameters9 + nParameters9);

	updateHistory(history, nIterationsMCMC, theta, nWeeks, bestGamma1, bestGamma2, true);

	writeCSV(history, "csv/history.csv");

	std::cout << nIterationsMCMC << " iterations took "
	<< (double) (clock() - start) / CLOCKS_PER_SEC << "seconds"
	<< std::endl;
	for (int i = 0; i < nParameters9; ++i) {
		std::cout << "acceptance rate for parameter " << i << " : "
		<< (double) accept[i] / nIterationsMCMC << std::endl;
	}
}
#endif

//found maximum at f -530.1039412
//#define testModel9
#ifdef testModel9
#include "Optimisation/likelihood9.h"
int main(int argc, char **argv) {

	DoubleVector theta(nParameters9, 0.0);
	nlopt::opt opt(nlopt::LN_SBPLX, nParameters9);

	Data d;
	d.constrained = true;
	d.team = -1;
	d.upToWeek = 38;
	d.printTrace = true;

	opt.set_max_objective(logLikelihood9, &d);

	opt.set_ftol_rel(1e-10);
	opt.set_maxeval(100000);
	double f;
	opt.optimize(theta, f);

	printf("result = %x\n", opt.last_optimize_result());
	printf("found maximum at f %0.10g\n\n", f);
	printf("theta: ");
	for (int i = 0; i < nParameters9; ++i) {
		printf("%f, ", theta[i]);
	}

	std::cout << std::endl;

	for (int i = 0; i < 20; ++i) {
		std::cout << TeamName[i] << ": " << theta[i + 10] << std::endl;
	}

}
#endif

//#define makeChiSqStats
#ifdef makeChiSqStats
#include "Bayesian/simulate9.h"
#include "Utils/csv.h"
int main(int argc, char **argv) {

	const int nIterationsMCMC = 45000;
	const int burn = 5000;
	const int thin = 1;
	const int nWeeks = 38;
	const int dim = 50;
	const int nLeagues = (nIterationsMCMC - burn) / thin;
	const double bestGamma1 = log(2.0);
	const double bestGamma2 = log(1.5);

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);
	DoubleVector theta(parameters9, parameters9 + nParameters9);

	updateHistory(history, nIterationsMCMC, theta, nWeeks, bestGamma1, bestGamma2, true);

	DoubleMatrix expectedScoresMat = simulateExpectedLeagueScores(history, burn, thin);
	ScoresMatrix expectedScores(dim);

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			expectedScores.at(i, j) = expectedScoresMat[i][j];
		}
	}

	ScoresMatrix leagueScores;
	std::vector<ChiSqStat> chiSqStats(nLeagues);

	for (int l = 0; l < nLeagues; ++l) {
		leagueScores = simulateLeagueScores(history[l * thin + burn]);
		chiSqStats[l] = chiSqTest(expectedScores, leagueScores);
	}

	DoubleMatrix stats = createDoubleMatrix(nLeagues, 2);
	for (int i = 0; i < nLeagues; ++i) {
		stats[i][0] = chiSqStats[i].getDf();
		stats[i][1] = chiSqStats[i].getChi();
	}

//	writeCSV(stats, "chi-sq1.csv");

	ChiSqStat chiSqStat = chiSqTest(expectedScores, getObservedScores());
	std::cout << "chiSqStat = " << chiSqStat.getChi() << std::endl;
	std::cout << "df = " << chiSqStat.getDf() << std::endl;

}
#endif

//#define makeExpectedScores
#ifdef makeExpectedScores
#include "Bayesian/simulate9.h"
int main(int argc, char **argv) {

	const int nIterationsMCMC = 45000;
	const int burn = 5000;
	const int thin = 1;
	const int nWeeks = 38;
	const double bestGamma1 = log(2.0);
	const double bestGamma2 = log(1.5);

	DoubleMatrix history = createDoubleMatrix(nIterationsMCMC, nParameters9);
	DoubleVector theta(parameters9, parameters9 + nParameters9);

	updateHistory(history, nIterationsMCMC, theta, nWeeks, bestGamma1, bestGamma2, true);

	DoubleMatrix expectedScores = simulateExpectedLeagueScores(history, burn, thin);

	writeCSV(expectedScores, "expectedScores.csv");

}
#endif
