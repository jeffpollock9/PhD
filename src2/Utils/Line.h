
#ifndef LINE_H_
#define LINE_H_

#include<cmath>

class Line {
public:
	static int round(double x);

	static bool isWholeLine(double line);
	static bool isOneQuarterLine(double line);
	static bool isHalfLine(double line);
	static bool isThreeQuarterLine(double line);

};

inline int Line::round(double x) {
	return x > 0.0 ? std::floor(x + 0.5) : std::ceil(x - 0.5);
}

inline bool Line::isWholeLine(double line) {
	int x = std::abs(Line::round(line * 100));
	return x % 100 == 0;
}

inline bool Line::isOneQuarterLine(double line) {
	int x = std::abs(Line::round(line * 100));
	return x % 100 == 25;
}

inline bool Line::isHalfLine(double line) {
	int x = std::abs(Line::round(line * 100));
	return x % 100 == 50;
}

inline bool Line::isThreeQuarterLine(double line) {
	int x = std::abs(Line::round(line * 100));
	return x % 100 == 75;
}

#endif /* LINE_H_ */
