#ifndef OPTIMISATIONTEST_H
#define OPTIMISATIONTEST_H

#include "BaseTest.h"
#include "Parametrisation/types.h"


class OptimisationTest : public BaseTest{
public:
    OptimisationTest();

private:
    using Scalar = long double;

    static std::istringstream testShermanMorissonUpdate(std::istream& istr);

    static std::istringstream testCost(std::istream& istr);

    static ChefDevr::Matrix<Scalar> testReadMatrix(std::istream &istr, unsigned int num_rows, unsigned int num_cols);


};

#endif // OPTIMISATIONTEST_H
