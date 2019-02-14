#include "OptimisationTest.h"

#include "Parametrisation/types.h"
#include "Optimisation/OptimisationSolver.h"
#include "BRDFReader/BRDFReader.h"

#include <iostream>
#include <string>
#include <fstream>



OptimisationTest::OptimisationTest(): BaseTest("Optimisation"){
    addTest(&testShermanMorissonUpdate, "Sherman Morisson update 1", "../tests/data/Optimisation/shermanMorisson/shermanMorissonUpdateSet1",
            "../tests/data/Optimisation/shermanMorisson/shermanMorissonUpdateSet1_output");
    addTest(&testShermanMorissonUpdate, "Sherman Morisson update 2", "../tests/data/Optimisation/shermanMorisson/shermanMorissonUpdateSet2",
            "../tests/data/Optimisation/shermanMorisson/shermanMorissonUpdateSet2_output");
    addTest(&testCost, "Cost 1", "../tests/data/Optimisation/cost/costSet1",
            "../tests/data/Optimisation/cost/costSet1_output");
    addTest(&testCost, "Cost 2", "../tests/data/Optimisation/cost/costSet2",
            "../tests/data/Optimisation/cost/costSet2_output");
    addTest(&testCost, "Cost 3", "../tests/data/Optimisation/cost/costSet3",
            "../tests/data/Optimisation/cost/costSet3_output");
}


std::istringstream OptimisationTest::testShermanMorissonUpdate(std::istream& istr){
    unsigned int num_rows, lv_num;
    Scalar detK, new_detK;
    ChefDevr::Matrix<Scalar> dummy{1, 1};
    ChefDevr::OptimisationSolver<Scalar> optimisation{1, 0.1, dummy, 2};

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
    const auto ZZt = Z * Z.transpose();

    ChefDevr::OptimisationSolver<Scalar> optimisation{d, 0.1, ZZt, 2};

    const auto K_minus1 = readMatrix(istr, num_rows, num_rows);

    istr >> detK;

    optimisation.cost(cost, K_minus1, detK);

    // Write Result
    std::stringstream ret;
    ret.precision(100);
    ret << cost;

    return std::istringstream(ret.str());
}


ChefDevr::Matrix<OptimisationTest::Scalar> OptimisationTest::readMatrix(std::istream &istr, unsigned int num_rows, unsigned int num_cols) {
    ChefDevr::Matrix<Scalar> matrix{num_rows, num_cols};

    for (unsigned int i = 0; i < num_rows; ++i) {
        for (unsigned int j = 0; j < num_cols; ++j) {
            istr >> matrix(i, j);
        }
    }

    return matrix;
}

std::istringstream OptimisationTest::testInitX(std::istream& istr) {
    ChefDevr::Matrix<Scalar> Z;
    uint latentDim;
    //load Z and latentDim

    ChefDevr::OptimisationSolver<Scalar> solver(1, Scalar(0), Z, latentDim);

    ChefDevr::Matrix<Scalar> ZZt;
    //load ZZt

    solver.initX(ZZt);

    std::stringstream ret;

    ret << solver.X;

    return std::istringstream(ret.str());
}
