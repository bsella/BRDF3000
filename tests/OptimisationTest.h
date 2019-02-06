#ifndef OPTIMISATIONTEST_H
#define OPTIMISATIONTEST_H

#include "BaseTest.h"
#include "Parametrisation/types.h"


class OptimisationTest : public BaseTest{
public:
    OptimisationTest();

private:
    using Scalar = long double;

    static std::istringstream shermanMorissonUpdate(std::istream& istr);

    static std::istringstream cost(std::istream& istr);

    static ChefDevr::Matrix<Scalar> readMatrix(std::istream &istr, unsigned int num_rows, unsigned int num_cols);
};

#endif // OPTIMISATIONTEST_H