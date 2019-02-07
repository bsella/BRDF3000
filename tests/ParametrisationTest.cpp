#include "ParametrisationTest.h"
#include "Parametrisation/Parametrisation.h"
#include "Parametrisation/types.h"
#include "BRDFReaderTest.h"

ParametrisationTest::ParametrisationTest(): BaseTest("Parametrisation") {
    addTest(&testCovariance, "Covariance1", "../tests/data/covTestSet1", "../tests/data/GT_covTestSet1");
    addTest(&testCovariance, "Covariance2", "../tests/data/covTestSet2", "../tests/data/GT_covTestSet2");
    addTest(&testCenter, "Center1", "../tests/data/centerTestSet1", "../tests/data/GT_centerTestSet1");
    addTest(&testCenter, "Center2", "../tests/data/centerTestSet2", "../tests/data/GT_centerTestSet2");
}

std::istringstream ParametrisationTest::testCovariance(std::istream& istr) {
    uint dim;
    istr >> dim;
    ChefDevr::Vector<double> A(dim), B(dim);
    for(uint i=0; i<dim; i++) {
        double tmp;
        istr >> tmp;
        A[i]=tmp;
    }
    for(uint i=0; i<dim; i++) {
        double tmp;
        istr >> tmp;
        B[i]=tmp;
    }
    return std::istringstream(std::to_string(ChefDevr::covariance(A,B)));
}

std::istringstream ParametrisationTest::testCenter(std::istream& istr) {
    uint w,h;
    istr >> w >> h;
    ChefDevr::Matrix<double> M(h,w);
    for(uint i=0; i<h; i++) {
        for(uint j=0; j<w; j++) {
            double tmp;
            istr >> tmp;
            M(i,j)=tmp;
        }
    }
    ChefDevr::RowVector<double> mean(h);
    ChefDevr::centerMat(M, mean);
    std::stringstream ret;
    ret << M;
    return std::istringstream(ret.str());
}
