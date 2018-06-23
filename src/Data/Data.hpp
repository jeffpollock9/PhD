#ifndef DATA_HPP_
#define DATA_HPP_

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "Match.hpp"
#include "Event.hpp"
#include "../Defn/defn.hpp"

typedef map<Date, vector<Event> > EventsMap;
typedef map<Date, vector<Match> > MatchesMap;

//TODO: date should read column headers for correct indexes

class Data {
public:
	Data() {};
	Data(const string& matchFile, const string& eventFile);
	void readMatches(const string& fileName);
	void readEvents(const string& fileName);

	const EventsMap& getEventsMap() const {
		return eventsMap;
	}

	const MatchesMap& getMatchesMap() const {
		return matchesMap;
	}

private:
	MatchesMap matchesMap;
	EventsMap eventsMap;
};

inline Data::Data(const string& matchFile, const string& eventFile) {
	readMatches(matchFile);
	readEvents(eventFile);
}

inline void Data::readMatches(const string& fileName) {
	ifstream file(fileName.c_str());

	string line, field;
	vector<string> v;

	vector<string> date, homeTeam, awayTeam;

	getline(file, line); //ignore first line
	while (getline(file, line)) {
		v.clear();
		stringstream ss(line);

		while (getline(ss, field, ',')) {
			v.push_back(field);
		}

		date.push_back(v[2]);
		homeTeam.push_back(v[4]);
		awayTeam.push_back(v[5]);
	}

	int n = date.size();
	string thisDate;
	string previousDate;
	vector<Match> m;

	previousDate = date[0];

	for (int i = 0; i < n; ++i) {
		thisDate = date[i];
		if (thisDate == previousDate) {
			m.push_back(Match(homeTeam[i], awayTeam[i]));
		} else {
			matchesMap.insert(make_pair(previousDate, m));
			m.clear();
			m.push_back(Match(homeTeam[i], awayTeam[i]));
		}
		previousDate = thisDate;
	}
	if (m.size() > 0) {
		matchesMap.insert(make_pair(previousDate, m));
	}
}

inline void Data::readEvents(const string& fileName) {
	ifstream file(fileName.c_str());

	string line, field;
	vector<string> v;

	vector<string> date, homeTeam, awayTeam;
	IntegerVector event;
	DoubleVector time;

	getline(file, line); //ignore first line
	while (getline(file, line)) {
		v.clear();
		stringstream ss(line);

		while (getline(ss, field, ',')) {
			v.push_back(field);
		}

		date.push_back(v[0]);
		homeTeam.push_back(v[5]);
		awayTeam.push_back(v[6]);
		time.push_back(atoi(v[3].c_str()) / 90.0);
		event.push_back(atoi(v[2].c_str()));
	}

	int n = date.size();
	string thisDate;
	string previousDate;
	vector<Event> e;

	previousDate = date[0];

	for (int i = 0; i < n; ++i) {
		thisDate = date[i];
		if (thisDate == previousDate) {
			e.push_back(Event(homeTeam[i], awayTeam[i], event[i], time[i]));
		} else {
			eventsMap.insert(make_pair(previousDate, e));
			e.clear();
			e.push_back(Event(homeTeam[i], awayTeam[i], event[i], time[i]));
		}
		previousDate = thisDate;
	}
	if (e.size() > 0) {
		eventsMap.insert(make_pair(previousDate, e));
	}

}

#endif /* DATA_HPP_ */
