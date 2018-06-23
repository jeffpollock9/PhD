#ifndef DIC_HPP_
#define DIC_HPP_

#include "../Inference/MCMC.hpp"

struct DIC {
	double pD;
	double value;
};

ostream& operator <<(std::ostream& o, const DIC& dic) {
	o << "pD: " << dic.pD << "\nvalue: " << dic.value << endl;
	return o;
}

inline DIC calculateDIC(const string& logLikFileName, const string& thetaFileName, const int nParameters, Model& model) {

	const int nSeasons = ALL_SEASONS.size();
	vector<Data> allData(nSeasons);

	static const string dataPath("/home/jeff/workspace/Latex-Utility-Model/data/processed/");

	for (int i = 0; i < nSeasons; ++i) {
		allData[i].readEvents(dataPath + "events-" + ALL_SEASONS[i] + ".csv");
		allData[i].readMatches(dataPath + "matches-" + ALL_SEASONS[i] + ".csv");
	}

	ifstream logLikFile(logLikFileName.c_str());
	string field1;
	DoubleVector logLikVec;

	while (getline(logLikFile, field1)) {
		logLikVec.push_back(strtod(field1.c_str(), NULL));
	}

	int n = logLikVec.size();

	arma::colvec logLik(logLikVec);
	arma::mat theta(n, nParameters);

	ifstream thetaFile(thetaFileName.c_str());
	string line, field2;

	int i = 0;
	while (getline(thetaFile, line)) {
		stringstream ss(line);

		int j = 0;
		while (getline(ss, field2, ',')) {
			theta.at(i, j) = strtod(field2.c_str(), NULL);
			++j;
		}
		++i;
	}

	double Dbar = -2.0 * arma::mean(logLik);
	arma::rowvec thetaBar = arma::mean(theta);

	model.setParametersAllSeasons(DoubleVector(thetaBar.begin(), thetaBar.end()));

	double DThetaBar = -2.0 * model.logLikelihoodAllSeasons(allData);
	double pD = Dbar - DThetaBar;
	double value = Dbar + pD;

	DIC out;
	out.pD = pD;
	out.value = value;

	return out;
}



#endif /* DIC_HPP_ */
