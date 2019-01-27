#ifndef PARAMETRISATIONTEST_H
#define PARAMETRISATIONTEST_H

#include "BaseTest.h"

class ParametrisationTest : public BaseTest{
public:
    ParametrisationTest();
private:
	static std::istream& test1(std::istream&);
	static std::istream& test2(std::istream&);
	static std::istream& test3(std::istream&);
};

#endif // PARAMETRISATIONTEST_H