#include "ParametrisationTest.h"
#include <Parametrisation/Parametrisation.h>
#include <Parametrisation/types.h>

ParametrisationTest::ParametrisationTest(): BaseTest("Parametrisation"){
	addTest(&testCovariance, "../tests/data/covTestSet1", "../tests/data/GT_covTestSet1");
}

std::istringstream ParametrisationTest::testCovariance(std::istream& istr){
	uint dim;
	istr >> dim;
	ChefDevr::Vector<double> A(10), B(10);
	for(uint i=0; i<dim; i++){
		double tmp;
		istr >> tmp;
		A[i]=tmp;
	}
	for(uint i=0; i<dim; i++){
		double tmp;
		istr >> tmp;
		B[i]=tmp;
	}
	return std::istringstream(std::to_string(ChefDevr::covariance(A,B)));
}