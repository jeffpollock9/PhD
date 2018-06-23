#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <nlopt.hpp>
#include <set>

#include "../Data/Data.hpp"
#include "../Data/League.hpp"
#include "../Defn/defn.hpp"
#include "../PiecewiseLinear/PiecewiseLinearModel.hpp"

static const int SUBDIVISIONS = 1000;
static const int DIM = 50;
static const double EPS = 1e-8;

class Model {
public:
  virtual ~Model() {
  }

  virtual DoublePair rate(const double t, const int homeScore, const int awayScore,
      const Team& homeTeam, const Team& awayTeam, const League& league, const bool isLog) const = 0;

  virtual DoublePair rateIntegral(const double start, const double end, const int homeScore,
      const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const;

  virtual double simulateGoalTime(const double start, const int homeScore, const int awayScore,
      const Team& homeTeam, const Team& awayTeam, const League& league) const;

  virtual void setParametersAllSeasons(const DoubleVector& theta) = 0;

  virtual void calculateAlphaLeague(DoublePair& alphaHomeAway, const int homeScore,
      const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const;

  double logLikelihoodSingleSeason(const Data& data, const Season& season) const;
  double logLikelihoodAllSeasons(const vector<Data>& allData) const;

  ScoreLine simulateScoreLine(const Team& homeTeam, const Team& awayTeam,
      const League& league) const;

  Scores simulateScoresMatrix(const int nGames, const Team& homeTeam, const Team& awayTeam,
      const League& league) const;

  double resourceOfTeam(const Team& team) const {
    return exp(this->logR[teamMap.at(team)]);
  }

  int getRhoIdx(const double t) const {
    if (t <= 44.0 / 90.0) {
      return 0;
    } else if (t <= 45.0 / 90.0) {
      return 1;
    } else if (t <= 89.0 / 90.0) {
      return 0;
    } else {
      return 2;
    }
  }

  void setLogR(const DoubleVector& logR) {
    this->logR = logR;
  }

  void setH(const double h) {
    this->h = h;
  }

  void setA(const double a) {
    this->a = a;
  }

  void setRho(const DoubleVector& rho) {
    this->rho = rho;
  }

  void setKnots(const vector<Knot>& knots) {
    utilityModel.setKnots(knots);
  }

  PiecewiseLinearModel utilityModel;
protected:
  double h, a;
  DoubleVector logR, rho;
  map<Team, int> teamMap;
};

template<class M>
struct ModelDataAndSeason {
  M model;
  Data data;
  Season season;
};

template<class M>
struct ModelAndAllData {
  M model;
  vector<Data> allData;
};

inline DoublePair Model::rateIntegral(const double start, const double end, const int homeScore,
    const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const {

  if (start == end) {
    return DoublePair(0.0, 0.0);
  }

  double l = (end - start) / SUBDIVISIONS;

  DoublePair sum = (rate(start, homeScore, awayScore, homeTeam, awayTeam, league, false)
      + rate(end, homeScore, awayScore, homeTeam, awayTeam, league, false)) / 2.0;

  for (int i = 1; i < SUBDIVISIONS; ++i) {
    sum += rate(start + i * l, homeScore, awayScore, homeTeam, awayTeam, league, false);
  }

  return sum * l;
}

inline ScoreLine Model::simulateScoreLine(const Team& homeTeam, const Team& awayTeam,
    const League& league) const {

  int homeScore = 0, awayScore = 0;
  double t = 0.0;

  while (true) {
    t = simulateGoalTime(t, homeScore, awayScore, homeTeam, awayTeam, league);

    if (t > 1.0) {
      break;
    }

    DoublePair ratesAtT = rate(t, homeScore, awayScore, homeTeam, awayTeam, league, false);

    double homeRateAtT = ratesAtT.first;
    double awayRateAtT = ratesAtT.second;

    if (arma::randu() < homeRateAtT / (homeRateAtT + awayRateAtT)) {
      homeScore++;
    } else {
      awayScore++;
    }
  }

  return make_pair(homeScore, awayScore);
}

inline void Model::calculateAlphaLeague(DoublePair& alphaHomeAway, const int homeScore,
    const int awayScore, const Team& homeTeam, const Team& awayTeam, const League& league) const {

  static const IntegerVector HOME_POS = { 1, 2, 0, 0 };
  static const IntegerVector AWAY_POS = { 0, 0, 1, 2 };

  double homeAlphaSum = 0.0, awayAlphaSum = 0.0;

  double currentHomeUtility = utilityModel.at(league.positionOfTeam(homeTeam));
  double currentAwayUtility = utilityModel.at(league.positionOfTeam(awayTeam));

  set<int> Phome, Paway;

  for (int i = 0; i < 4; ++i) {
    League leagueCopy = league;
    leagueCopy.currentMatch(homeTeam, awayTeam, homeScore + HOME_POS[i], awayScore + AWAY_POS[i]);

    Phome.insert(leagueCopy.positionOfTeam(homeTeam));
    Paway.insert(leagueCopy.positionOfTeam(awayTeam));
  }

  for (const int p : Phome) {
    homeAlphaSum += utilityModel.at(p) - currentHomeUtility;
  }

  for (const int p : Paway) {
    awayAlphaSum += utilityModel.at(p) - currentAwayUtility;
  }

  alphaHomeAway.first = logistic(homeAlphaSum);
  alphaHomeAway.second = logistic(awayAlphaSum);
}

inline double Model::simulateGoalTime(const double start, const int homeScore, const int awayScore,
    const Team& homeTeam, const Team& awayTeam, const League& league) const {

  double tau = -log(arma::randu());
  double integral = 0.0;
  double end = start;
  DoublePair homeAwayIntegral;

  while (integral < tau) {
    homeAwayIntegral = (rate(end, homeScore, awayScore, homeTeam, awayTeam, league, false)
        + rate(end + EPS, homeScore, awayScore, homeTeam, awayTeam, league, false)) / 2.0;

    integral += EPS * (homeAwayIntegral.first + homeAwayIntegral.second);
    end += EPS;
  }

  return end;
}

inline Scores Model::simulateScoresMatrix(const int nGames, const Team& homeTeam,
    const Team& awayTeam, const League& league) const {

  Scores scores(DIM, DIM);

  for (int game = 0; game < nGames; ++game) {
    ScoreLine scoreLine = simulateScoreLine(homeTeam, awayTeam, league);
    ++scores.at(scoreLine.first, scoreLine.second);
  }

  return scores / nGames;
}

inline double Model::logLikelihoodSingleSeason(const Data& data, const Season& season) const {

  League league(season);

  EventsMap eventsMap = data.getEventsMap();
  MatchesMap matchesMap = data.getMatchesMap();

  map<Match, ScoreLine> matchScoresMap;

  double logLik = 0.0;

  for (const auto dateAndMatches : matchesMap) {

    Date date = dateAndMatches.first;

    vector<Match> matches = dateAndMatches.second;

    matchScoresMap.clear();

    //maps of team scores for each home/away team playing that day
    for (const Match match : matches) {
      matchScoresMap.insert(make_pair(match, make_pair(0, 0)));

      //each match begins at 0-0 i.e. a draw so each team gets 1 point
      league.currentMatch(match.getHomeTeam(), match.getAwayTeam(), 0, 0);
    }

    //go through each event integrating each teams rate + adding log of goal events
    double previousEventTime = 0.0;
    for (const Event event : eventsMap.at(date)) {

      double eventTime = event.getTime();

      //integral of rates for each team
      for (const Match match : matches) {
        Team matchHomeTeam = match.getHomeTeam();
        Team matchAwayTeam = match.getAwayTeam();

        ScoreLine matchScoreLine = matchScoresMap.at(match);

        int matchHomeScore = matchScoreLine.first;
        int matchAwayScore = matchScoreLine.second;

        DoublePair rateIntegrals = rateIntegral(previousEventTime, eventTime, matchHomeScore,
            matchAwayScore, matchHomeTeam, matchAwayTeam, league);

        logLik -= rateIntegrals.first;
        logLik -= rateIntegrals.second;
      }

      previousEventTime = eventTime;

      int eventNumber = event.getNumber();
      if (eventNumber != -1) {

        Team eventHomeTeam = event.getHomeTeam();
        Team eventAwayTeam = event.getAwayTeam();

        ScoreLine& eventScoreLine = matchScoresMap.at(Match(eventHomeTeam, eventAwayTeam));

        int& eventHomeScore = eventScoreLine.first;
        int& eventAwayScore = eventScoreLine.second;

        //log of rate at goal times
        if (eventNumber == 1) {
          logLik += rate(eventTime, eventHomeScore, eventAwayScore, eventHomeTeam, eventAwayTeam,
              league, true).first;
          ++eventHomeScore;
        } else if (eventNumber == 2) {
          logLik += rate(eventTime, eventHomeScore, eventAwayScore, eventHomeTeam, eventAwayTeam,
              league, true).second;
          ++eventAwayScore;
        }

        //update league positions after event
        league.currentMatch(eventHomeTeam, eventAwayTeam, eventHomeScore, eventAwayScore);
      } else if (eventTime == 1.0) {
        league.endMatches(matchScoresMap);
      }
    }
  }

  return logLik;
}

//temp
//inline double Model::logLikelihoodSingleSeason(const Data& data, const Season& season) const {
//
//	League league(season);
//
//	EventsMap eventsMap = data.getEventsMap();
//	MatchesMap matchesMap = data.getMatchesMap();
//
//	map<Match, ScoreLine> matchScoresMap;
//
//	double logLik = 0.0;
//
//	ofstream fileRates;
//	fileRates.open("rates.csv");
//
//	for (const auto dateAndMatches : matchesMap) {
//
//		Date date = dateAndMatches.first;
//
//		vector<Match> matches = dateAndMatches.second;
//
//		matchScoresMap.clear();
//
//		//temp to print out league before last 10 matches
//		for (const Match match : matches) {
//			Team matchHomeTeam = match.getHomeTeam();
//			Team matchAwayTeam = match.getAwayTeam();
//
//			if (matchHomeTeam == Team("Manchester City") && matchAwayTeam == Team("Queens Park Rangers")) {
//				std::cout << league << std::endl;
//			}
//		}
//
//		//maps of team scores for each home/away team playing that day
//		for (const Match match : matches) {
//			matchScoresMap.insert(make_pair(match, make_pair(0, 0)));
//
//			//each match begins at 0-0 i.e. a draw so each team gets 1 point
//			Team matchHomeTeam = match.getHomeTeam();
//			Team matchAwayTeam = match.getAwayTeam();
//
//			league.currentMatch(matchHomeTeam, matchAwayTeam, 0, 0);
//		}
//
//		//go through each event integrating each teams rate + adding log of goal events
//		double previousEventTime = 0.0;
//		for (const Event event : eventsMap.at(date)) {
//
//			double eventTime = event.getTime();
//
//			//integral of rates for each team
//			for (const Match match : matches) {
//				Team matchHomeTeam = match.getHomeTeam();
//				Team matchAwayTeam = match.getAwayTeam();
//
//				ScoreLine matchScoreLine = matchScoresMap.at(match);
//
//				int matchHomeScore = matchScoreLine.first;
//				int matchAwayScore = matchScoreLine.second;
//
//				if (matchHomeTeam == Team("Manchester City")
//						&& matchAwayTeam == Team("Queens Park Rangers")) {
//					for (double t = previousEventTime; t < eventTime; t += 0.001) {
//
//						DoublePair matchRates = rate(t, matchHomeScore, matchAwayScore, matchHomeTeam,
//								matchAwayTeam, league, false);
//						fileRates << t << "," << matchRates.first << "," << matchRates.second << std::endl;
//
//					}
//				}
//
//				DoublePair rateIntegrals = rateIntegral(previousEventTime, eventTime, matchHomeScore,
//						matchAwayScore, matchHomeTeam, matchAwayTeam, league);
//
//				logLik -= rateIntegrals.first;
//				logLik -= rateIntegrals.second;
//			}
//
//			previousEventTime = eventTime;
//
//			int eventNumber = event.getNumber();
//			if (eventNumber != -1) {
//
//				Team eventHomeTeam = event.getHomeTeam();
//				Team eventAwayTeam = event.getAwayTeam();
//
//				ScoreLine& eventScoreLine = matchScoresMap.at(Match(eventHomeTeam, eventAwayTeam));
//
//				int& eventHomeScore = eventScoreLine.first;
//				int& eventAwayScore = eventScoreLine.second;
//
//				//log of rate at goal times
//				if (eventNumber == 1) {
////					logLik += rate(eventTime, eventHomeScore, eventAwayScore, eventHomeTeam, eventAwayTeam,
////							league, true).first;
//					++eventHomeScore;
//				} else if (eventNumber == 2) {
////					logLik += rate(eventTime, eventHomeScore, eventAwayScore, eventHomeTeam, eventAwayTeam,
////							league, true).second;
//					++eventAwayScore;
//				}
//
//				//update league positions after event
//				league.currentMatch(eventHomeTeam, eventAwayTeam, eventHomeScore, eventAwayScore);
//			} else if (eventTime == 1.0) {
//				league.endMatches(matchScoresMap);
//			}
//		}
//	}
//	fileRates.close();
//
//	return logLik;
//}

inline double Model::logLikelihoodAllSeasons(const vector<Data>& allData) const {

  const int nSeasons = ALL_SEASONS.size();
  double logLik = 0.0;

  #pragma omp parallel for reduction(+:logLik)
  for (int i = 0; i < nSeasons; ++i) {
    logLik += logLikelihoodSingleSeason(allData[i], ALL_SEASONS[i]);
  }

  return logLik;
}

#endif /* MODEL_HPP_ */
