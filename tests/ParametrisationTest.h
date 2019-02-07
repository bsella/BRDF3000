#ifndef PARAMETRISATIONTEST_H
#define PARAMETRISATIONTEST_H

#include "BaseTest.h"

class ParametrisationTest : public BaseTest{
public:
	ParametrisationTest();
private:
	static std::istringstream testCovariance(std::istream&);
        static std::istringstream testComputeCovVector(std::istream&);
        static std::istringstream testReconstructWithoutMean(std::istream&);
        static std::istringstream testCenter(std::istream&);
        static std::istringstream testReconstructionError(std::istream&);
        static std::istringstream testReconstruct(std::istream&);
};

#endif // PARAMETRISATIONTEST_H
