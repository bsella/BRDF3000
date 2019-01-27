#include "ParametrisationTest.h"

ParametrisationTest::ParametrisationTest(): BaseTest(){
	addTest(&test1, "../tests/data/test.txt", "../tests/data/test2.txt");
	addTest(&test1, "../tests/data/test2.txt", "../tests/data/test3.txt");
}

std::istream& ParametrisationTest::test1(std::istream& istr){
	return istr;
}