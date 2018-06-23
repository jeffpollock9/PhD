//#define MCMC
#ifdef MCMC
#include "Inference/MCMC.hpp"
int main(int argc, char **argv) {
	runMCMC();
}
#endif

//#define RJMCMC
#ifdef RJMCMC
#include "Inference/RJMCMC.hpp"
int main(int argc, char **argv) {
	runRJMCMC();
}
#endif

//#define RATES
//need to uncomment some hacks for this to work
#ifdef RATES
#include "Models/ModelRJ.hpp"
static const string dataPath("/home/jeff/workspace/Latex-Utility-Model/data/processed/");
int main(int argc, char **argv) {

	Data data;

	data.readEvents(dataPath + "events-" + Season("2011-2012") + ".csv");
	data.readMatches(dataPath + "matches-" + Season("2011-2012") + ".csv");

	const bool AllSeasons = true;

	ModelRJ model(AllSeasons);

	DoubleVector x(20);

	for (int i = 0; i < 20; ++i) {
		x[i] = i + 1;
	}

	DoubleVector y = { 11.4835, 8.9914, 6.6695, 4.4751, 4.1354, 5.0013, 5.7721, 6.3387, 6.716, 6.6787, 6.2198,
			5.3751, 4.3238, 3.8981, 3.4289, 2.9672, 2.429, 1.5284, 0.734, 0 };

	vector<Knot> knots(20);

	for (int i = 0; i < 20; ++i) {
		knots[i] = Knot(x[i], y[i]);
	}

	DoubleVector theta = { 0.3282, 0.0, 0.2303, 0.2829, 0.3714, 0.8104, 0.6298, 0.6193, 0.7152, 1.7222, 0.3025,
			-0.2812, -0.8699, -0.7716, -1.2268, -0.9832, -1.8452, 0.3797, -2.4408, -0.0109, -0.4727, -1.7831,
			0.2134, 0.1458, 0.5058, -1.1457, -0.5922, -0.8179, -0.8923, -1.2841, -1.2191, -0.7399, -0.7951,
			-0.6968, -0.0229, -1.0002, -0.892, -1.111, -1.5259, -3.1242, 0.0743 };

	model.setParametersAllSeasons(theta);
	model.utilityModel.setKnots(knots);

	fileAlpha.open("alpha.csv");

	model.logLikelihoodSingleSeason(data, Season("2011-2012"));

	fileAlpha.close();
}
#endif

#define DIC_CALC
#ifdef DIC_CALC
#include "Math/DIC.hpp"
#include "Models/Model0.hpp"
#include "Models/Model1.hpp"
#include "Models/Model2.hpp"
#include "Models/Model3.hpp"
int main(int argc, char **argv) {
	string logLik;
	string theta;
	DIC dic;
	bool allSeasons = true;

	logLik = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC0/logLikHistory0.txt";
	theta = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC0/thetaHistory0.txt";
	Model0 model0(allSeasons);
	dic = calculateDIC(logLik, theta, Model0::getNParameters(), model0);
	std::cout << "Model0:" << std::endl;
	std::cout << "Dbar: " << dic.value - 2.0 * dic.pD << std::endl;
	std::cout << dic << std::endl << std::endl;

	logLik = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC1/logLikHistory1.txt";
	theta = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC1/thetaHistory1.txt";
	Model1 model1(allSeasons);
	dic = calculateDIC(logLik, theta, Model1::getNParameters(), model1);
	std::cout << "Model1:" << std::endl;
	std::cout << "Dbar: " << dic.value - 2.0 * dic.pD << std::endl;
	std::cout << dic << std::endl << std::endl;

	logLik = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC2/logLikHistory2.txt";
	theta = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC2/thetaHistory2.txt";
	Model2 model2(allSeasons);
	dic = calculateDIC(logLik, theta, Model2::getNParameters(), model2);
	std::cout << "Model2:" << std::endl;
	std::cout << "Dbar: " << dic.value - 2.0 * dic.pD << std::endl;
	std::cout << dic << std::endl << std::endl;

	logLik = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC3/logLikHistory3.txt";
	theta = "/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC3/thetaHistory3.txt";
	Model3 model3(allSeasons);
	dic = calculateDIC(logLik, theta, Model3::getNParameters(), model3);
	std::cout << "Model3:" << std::endl;
	std::cout << "Dbar: " << dic.value - 2.0 * dic.pD << std::endl;
	std::cout << dic << std::endl << std::endl;

}
#endif

//#define FIND_MLE_0
#ifdef FIND_MLE_0
#include "Models/Model0.hpp"
int main(int argc, char **argv) {
	int nSeasons = ALL_SEASONS.size();
	vector<Data> allData(nSeasons);

	static const string path("/home/jeff/workspace/Latex-Utility-Model/data/processed/");
	//static const string path("/home/jp148/Data/");

	for (int i = 0; i < nSeasons; ++i) {
		allData[i].readEvents(path + "events-" + ALL_SEASONS[i] + ".csv");
		allData[i].readMatches(path + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	Model0 model(allSeasons);

	//logLik0: -2776.15 (44)
	//logLik1: -2773.54 (45)
	//logLik2: -2773.14 (46)
	//logLik3: -2770.38 (48)
	DoubleVector theta = model.getStartingValues();

	model.findMLE(theta, allData);
}
#endif

//#define FIND_MLE_14_ALL_SEASONS
#ifdef FIND_MLE_14_ALL_SEASONS
#define HAVE_NLOPT
#include "Models/Model14.hpp"
int main(int argc, char **argv) {

	int nSeasons = ALL_SEASONS.size();
	vector<Data> data(nSeasons);

	string path("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

	for (int i = 0; i < nSeasons; ++i) {
		data[i].readEvents(path + "events-" + ALL_SEASONS[i] + ".csv");
		data[i].readMatches(path + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	Model14 model(allSeasons);

	//logLik: -2770.83
	//n. par.: 60
	DoubleVector theta = {0.282976, -0.0465221, 0.356382, 0.426047, 0.510884, 0.710569, 0.552829, 0.48403,
		0.707881, 1.70352, 0.938272, 0.664458, 0.473778, 0.499912, 0.407389, 0.440319, 0.148742, 0.98324,
		-0.503488, 0.773381, 0.598093, 0.256671, 0.88053, 0.854411, 1.04912, 0.403081, 0.558772, 0.550057,
		0.467406, 0.410409, 0.356402, 0.508597, 0.493327, 0.587797, 0.776501, 0.458691, 0.47313, 0.409468,
		0.33136, -5.65107, 0.12711, 28.8914, 16.2711, 14.4736, 11.1271, 10.3455, 10.3455, 7.68901,
		7.66282, 7.36445, 7.34588, 7.34587, 7.08953, 3.89458, 3.89458, 3.56235, 3.56235, 3.00556, 1.02925,
		3.28287e-05};

	model.findMLE(theta, data);
}
#endif

//#define FIND_MLE_13_ALL_SEASONS
#include "Data/Data.hpp"
#ifdef FIND_MLE_13_ALL_SEASONS
#include "Models/Model13.hpp"
int main(int argc, char **argv) {

	int nSeasons = ALL_SEASONS.size();
	vector<Data> data(nSeasons);

	string path("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

	for (int i = 0; i < nSeasons; ++i) {
		data[i].readEvents(path + "events-" + ALL_SEASONS[i] + ".csv");
		data[i].readMatches(path + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	Model13 model13(allSeasons);

	//logLik: -2771.89
	//n. par.: 48
	DoubleVector theta = {0.33816, 0.00920651, 0.353269, 0.421273, 0.491995, 0.669566, 0.540596, 0.481584,
		0.708091, 1.70408, 0.923125, 0.644804, 0.458657, 0.473897, 0.373632, 0.421847, 0.101387, 0.976685,
		-0.501617, 0.76255, 0.583819, 0.231653, 0.872221, 0.842404, 1.04295, 0.388859, 0.539908, 0.515562,
		0.454533, 0.382279, 0.337233, 0.491711, 0.477385, 0.566313, 0.757479, 0.431869, 0.45512, 0.387122,
		0.296786, -23.7956, 0.877321, 13.2509, 13.262, 6.16238, 3.94994, 7.73081, 1.59789, 8.2504e-11};

	model13.findMLE(theta, data);
}
#endif

//#define FIND_MLE_11_ALL_SEASONS
#ifdef FIND_MLE_11_ALL_SEASONS
#define HAVE_NLOPT
#include "Models/Model11.hpp"
int main(int argc, char **argv) {

	int nSeasons = ALL_SEASONS.size();
	vector<Data> data(nSeasons);

	string path("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

	for (int i = 0; i < nSeasons; ++i) {
		data[i].readEvents(path + "events-" + ALL_SEASONS[i] + ".csv");
		data[i].readMatches(path + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	Model11 model11(allSeasons);

	//logLik: -2769.30
	//n. par.: 48
	DoubleVector theta = {0.289096, -0.0400848, 0.354038, 0.410594, 0.483351, 0.730718, 0.56244, 0.495277,
		0.707502, 1.70607, 0.844082, 0.536788, 0.323115, 0.345956, 0.242618, 0.282746, -0.0585018,
		0.894378, -0.99492, 0.662983, 0.462712, 0.0614064, 0.785704, 0.751695, 0.969919, 0.245218,
		0.418171, 0.395916, 0.31563, 0.233076, 0.181737, 0.360352, 0.342446, 0.445395, 0.661706, 0.297798,
		0.315175, 0.240284, 0.15638, -4.20749, 0.100131, 11.1791, 11.1262, 7.47513, 6.81639, 7.72237,
		6.16655, 0.990518};

	model11.findMLE(theta, data);
}
#endif

//#define FIND_MLE_10_ALL_SEASONS
#ifdef FIND_MLE_10_ALL_SEASONS
#include "Models/Model10.hpp"
int main(int argc, char **argv) {

	int nSeasons = ALL_SEASONS.size();
	vector<Data> data(nSeasons);

	string path("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

	for (int i = 0; i < nSeasons; ++i) {
		data[i].readEvents(path + "events-" + ALL_SEASONS[i] + ".csv");
		data[i].readMatches(path + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	bool allSeasons = true;
	Model10 model10(allSeasons);

	//logLik: -2769.5
	//n. par.: 49
	DoubleVector theta = {0.295659, -0.0334682, 0.341361, 0.402159, 0.478532, 0.738071, 0.565233, 0.49996,
		0.707445, 1.70679, 0.795383, 0.470199, 0.242479, 0.264987, 0.151411, 0.197619, -0.18076, 0.848129,
		-1.32747, 0.604895, 0.391207, -0.0456696, 0.734788, 0.698277, 0.927861, 0.155906, 0.342353,
		0.319455, 0.232653, 0.141447, 0.0883005, 0.281078, 0.260499, 0.371372, 0.603641, 0.211786,
		0.230629, 0.151173, 0.0569774, -4.15411, 0.100906, 1.0789e-10, 10.9501, 10.8999, 7.37369, 6.72951,
		7.65119, 6.13325, 0.9861};

	model10.findMLE(theta, data);
}
#endif

//#define FIND_MLE_12
#ifdef FIND_MLE_12
#include "Models/Model12.hpp"
int main(int argc, char **argv) {
	Data data;
	data.readEvents("/home/jeff/workspace/Latex-Utility-Model/data/processed/events-2011-2012.csv");
	data.readMatches("/home/jeff/workspace/Latex-Utility-Model/data/processed/matches-2011-2012.csv");

	Model12 model12("2011-2012");

	//logLik: -521.234
	//n. par.: 39
	DoubleVector theta = {0.225146, -0.0416059, 0.883856, 0.322791, 0.166807, 0.583482, 0.643812, 0.685659,
		1.3487, 1.55579, 0.689571, 0.0971835, -0.0913643, 0.0866448, 0.596694, 0.512515, 0.390362,
		0.444817, 0.969246, 0.910037, 0.444979, 0.205979, 0.0475617, 0.143785, 0.331957, 0.226996,
		0.639206, 0.267187, 0.0852214, -0.255749, -0.669455, 0.019745, 8.83791, 8.79551, 8.20694, 8.03372,
		8.24719, 7.17986, 4.77685};

	model12.findMLE(theta, data);
}
#endif

//#define FIND_MLE_11
#ifdef FIND_MLE_11
#include "Models/Model11.hpp"
int main(int argc, char **argv) {
	Data data;
	data.readEvents("/home/jeff/workspace/Latex-Utility-Model/data/processed/events-2011-2012.csv");
	data.readMatches("/home/jeff/workspace/Latex-Utility-Model/data/processed/matches-2011-2012.csv");

	Model11 model11("2011-2012");

	//logLik: -518.89
	//n. par.: 39
	DoubleVector theta = {0.188824, -0.0809313, 1, 0.215857, 0.00751937, 1, 0.740634, 0.66176, 1.34974,
		1.53529, 0.672781, -0.0140278, -0.223774, 0.00707319, 0.527203, 0.417801, 0.339209, 0.366976,
		0.936696, 0.88616, 0.427141, 0.156961, -0.0345712, 0.0669311, 0.213739, 0.13252, 0.598149,
		0.154338, -0.0498836, -0.312303, -0.334847, 0.0274298, 8.76051, 8.71939, 8.13026, 7.95839,
		8.08659, 6.974, 4.77627};

	model11.findMLE(theta, data);
}
#endif

//#define FIND_MLE_10
#ifdef FIND_MLE_10
#include "Models/Model10.hpp"
int main(int argc, char **argv) {
	Data data;
	data.readEvents("/home/jeff/workspace/Latex-Utility-Model/data/processed/events-2011-2012.csv");
	data.readMatches("/home/jeff/workspace/Latex-Utility-Model/data/processed/matches-2011-2012.csv");

	Model10 model10("2011-2012");

	//logLik: -518.893
	//n. par.: 40
	DoubleVector theta = {0.191267, -0.0784153, 1, 0.210589, 0.00518262, 0.999999, 0.742153, 0.66637,
		1.34958, 1.53581, 0.666079, -0.0267703, -0.239141, -0.00486614, 0.519874, 0.409885, 0.330906,
		0.358489, 0.931612, 0.880551, 0.418008, 0.145759, -0.0479787, 0.0544067, 0.203561, 0.121088,
		0.590575, 0.143919, -0.0636211, -0.330133, -0.32721, 0.0275887, 7.74555e-11, 8.7474, 8.7052,
		8.10827, 7.93604, 8.06435, 6.95077, 4.76816,};

	model10.findMLE(theta, data);
}
#endif

//#define TEST_LOGLIK
#ifdef TEST_LOGLIK
int main(int argc, char **argv) {
	Data data;
	data.readEvents("/home/jeff/workspace/Latex-Utility-Model/data/processed/events-2011-2012.csv");
	data.readMatches("/home/jeff/workspace/Latex-Utility-Model/data/processed/matches-2011-2012.csv");

	double h = 0.4, a = 0.1;
	DoubleVector logR(20, 0.0);
	DoubleVector rho = {0, 1.0, 1.5};
	Season season("2011-2012");
	double A = 0.0, B = 0.0, C = 0.0;
	DoubleVector c = {0.4, 0.2, 0.2};
	DoubleVector d = {0.8, 0.7, 0.7};
	DoubleVector U = {5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.25, 0.01};

	Model10 model(logR, h, a, c, d, U, A, B, C, rho, season);
	//RUNTIME of logLik: 0.596934 seconds average using numeric integral 10 divisions (-610.183)
	//RUNTIME of logLik: ~0.04 seconds average using analytic integral (-610.13, A = B = C = 0)

	INIT_TIMER

	double l = model.logLikelihood(data);

	std::cout << l << std::endl;

	STOP_TIMER("logLik")

	return 0;
}
#endif

