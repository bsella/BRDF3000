#ifndef PARAMETRISATIONTEST_H
#define PARAMETRISATIONTEST_H

#include "BaseTest.h"

class ParametrisationTest : public BaseTest{
public:
	ParametrisationTest();
private:
	static std::istringstream testCovariance(std::istream&);
};

#endif // PARAMETRISATIONTEST_H