#include "OptimisationTest.h"

#include "Parametrisation/types.h"
#include "Optimisation/OptimisationSolver.h"



OptimisationTest::OptimisationTest(): BaseTest("Optimisation"){
    addTest(&shermanMorissonUpdate, "Sherman Morisson update 1", "../tests/data/shermanMorisson/shermanMorissonUpdateSet1",
            "../tests/data/shermanMorisson/shermanMorissonUpdateSet1_output");
    addTest(&shermanMorissonUpdate, "Sherman Morisson update 2", "../tests/data/shermanMorisson/shermanMorissonUpdateSet2",
            "../tests/data/shermanMorisson/shermanMorissonUpdateSet2_output");
    addTest(&shermanMorissonUpdate, "Sherman Morisson update 3", "../tests/data/shermanMorisson/shermanMorissonUpdateSet3",
            "../tests/data/shermanMorisson/shermanMorissonUpdateSet3_output");
}


std::istringstream OptimisationTest::shermanMorissonUpdate(std::istream& istr){
    using Scalar = long double;
    unsigned int num_rows, lv_num;
    Scalar detK, new_detK;
    ChefDevr::Matrix<Scalar> dummy{1, 1};
    ChefDevr::OptimisationSolver<Scalar> optimisation{0.1, dummy, 2};

    istr >> num_rows;
    istr >> lv_num;

    // Read K_minus_1
    ChefDevr::Matrix<Scalar> K_minus1{num_rows, num_rows};
    ChefDevr::Matrix<Scalar> new_K_minus1{num_rows, num_rows};

    for (unsigned int i = 0; i < num_rows; ++i) {
        for (unsigned int j = 0; j < num_rows; ++j) {
            istr >> K_minus1(i, j);
        }
    }

    istr >> detK;

    // Read diff_vector
    ChefDevr::Vector<Scalar> diff_vector{num_rows};
    for (unsigned int i = 0; i < num_rows; ++i) {
        istr >> diff_vector(i);
    }

    optimisation.shermanMorissonUpdate(K_minus1, new_K_minus1, detK, new_detK, lv_num, diff_vector);

    // Writing the output

    std::stringstream ret;
    ret.precision(100);

    for (unsigned int i = 0; i < num_rows; ++i) {
        for (unsigned int j = 0; j < num_rows; ++j) {
            ret << new_K_minus1(i, j) << ' ';
        }
        ret << std::endl;
    }

    ret << std::endl;
    ret << new_detK;

    return std::istringstream(ret.str());
}