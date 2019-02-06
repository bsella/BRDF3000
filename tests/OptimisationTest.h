#ifndef OPTIMISATIONTEST_H
#define OPTIMISATIONTEST_H

#include "BaseTest.h"


class OptimisationTest : public BaseTest{
public:
    OptimisationTest();

private:
    static std::istringstream shermanMorissonUpdate(std::istream& istr);
};

#endif // OPTIMISATIONTEST_H