#ifndef DISTANCE_H_
#define DISTANCE_H_

#include "Vector.h"
#include <cmath>
#include <map>
#ifdef bet
	#include "../Data/BetData.h"
#else
	#ifdef getGamma
		#include "../Data/Data2010_2011.h"
	#else
		#include "../Data/Data2011_2012.h"
	#endif
#endif

struct D {
	int D1;
	int D2;
};

static std::map<IntegerVector, D> map;
static const unsigned int MAX_MAP_SIZE = 1000000;

inline int teamInPosition(const int position, const IntegerVector &Rorder) {

	for (int team = 0; team < 20; ++team) {
		if (position == Rorder[team]) {
			return team;
		}
	}

	throw std::runtime_error("team position not found");
}

inline void bubbleSwaps(int &D1, int &D2, DoubleVector &R) {

	D1 = 0;
	D2 = 0;

	DoubleVector Rsorted = R;
	std::sort(Rsorted.begin(), Rsorted.end());

	IntegerVector Rorder(20);
	for (int i = 0; i < 20; ++i) {
		for (int j = 0; j < 20; ++j) {
			if (R[i] == Rsorted[j]) {
				Rorder[i] = 20 - j;
				break;
			}
		}
	}

	if (map.count(Rorder)) {
		D1 = map.at(Rorder).D1;
		D2 = map.at(Rorder).D2;
		return;
	}

	int i, j, firstTeam, secondTeam;
	bool madeSwaps = true;

	while (madeSwaps) {

		madeSwaps = false;
		for (i = 1; i <= 19; ++i) {
			j = i + 1;

			firstTeam = teamInPosition(i, Rorder);
			secondTeam = teamInPosition(j, Rorder);

			if (LAST_POSITION_2011_2012[firstTeam] > LAST_POSITION_2011_2012[secondTeam]) {
				//swap two teams
				madeSwaps = true;
				++Rorder[firstTeam];
				--Rorder[secondTeam];
				if (LAST_POSITION_2011_2012[firstTeam] <= 17 && LAST_POSITION_2011_2012[secondTeam] <= 17) {
					++D1;
				} else {
					++D2;
				}
			}
		}
	}

	if (map.size() <= MAX_MAP_SIZE) {
		D d;
		d.D1 = D1;
		d.D2 = D2;
		map.insert(std::make_pair(Rorder, d));
	}
	return;
}

inline int sign(double x) {
	return x < 0 ? -1 : +1;
}

inline double kendallTau(DoubleVector x, DoubleVector y) {

	int numerator = 0;
	int n = x.size();

	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			numerator += sign(x[i] - x[j]) * sign(y[i] - y[j]);
		}
	}

	return numerator / (0.5 * n * (n - 1.0));
}

#endif /* DISTANCE_H_ */
