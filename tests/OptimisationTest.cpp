#include "OptimisationTest.h"

#include "Parametrisation/types.h"
#include "Optimisation/OptimisationSolver.h"



OptimisationTest::OptimisationTest(): BaseTest("Optimisation"){
    addTest(&testShermanMorissonUpdate, "Sherman Morisson update 1", "../tests/data/shermanMorisson/shermanMorissonUpdateSet1",
            "../tests/data/shermanMorisson/shermanMorissonUpdateSet1_output");
    addTest(&testShermanMorissonUpdate, "Sherman Morisson update 2", "../tests/data/shermanMorisson/shermanMorissonUpdateSet2",
            "../tests/data/shermanMorisson/shermanMorissonUpdateSet2_output");
    addTest(&testCost, "Cost 1", "../tests/data/cost/costSet1",
            "../tests/data/cost/costSet1_output");
    addTest(&testCost, "Cost 2", "../tests/data/cost/costSet2",
            "../tests/data/cost/costSet2_output");
    addTest(&testCost, "Cost 3", "../tests/data/cost/costSet3",
            "../tests/data/cost/costSet3_output");
}


std::istringstream OptimisationTest::testShermanMorissonUpdate(std::istream& istr){
    unsigned int num_rows, lv_num;
    Scalar detK, new_detK;
    ChefDevr::Matrix<Scalar> dummy{1, 1};
    ChefDevr::OptimisationSolver<Scalar> optimisation{0.1, dummy, 2};

    istr >> num_rows;
    istr >> lv_num;

    const auto K_minus1 = readMatrix(istr, num_rows, num_rows);

    ChefDevr::Matrix<Scalar> new_K_minus1{num_rows, num_rows};

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


std::istringstream OptimisationTest::testCost(std::istream& istr) {
    unsigned int num_rows, d;
    Scalar detK, cost;

    istr >> num_rows;
    istr >> d;

    const auto Z = readMatrix(istr, num_rows, d);

    ChefDevr::OptimisationSolver<Scalar> optimisation{0.1, Z, 2};

    const auto K_minus1 = readMatrix(istr, num_rows, num_rows);

    istr >> detK;

    optimisation.cost(cost, K_minus1, detK);

    // Write Result
    std::stringstream ret;
    ret.precision(100);
    ret << cost;

    return std::istringstream(ret.str());
}


ChefDevr::Matrix<OptimisationTest::Scalar> OptimisationTest::testReadMatrix(std::istream &istr, unsigned int num_rows, unsigned int num_cols) {
    ChefDevr::Matrix<Scalar> matrix{num_rows, num_cols};

    for (unsigned int i = 0; i < num_rows; ++i) {
        for (unsigned int j = 0; j < num_cols; ++j) {
            istr >> matrix(i, j);
        }
    }

    return matrix;
}
