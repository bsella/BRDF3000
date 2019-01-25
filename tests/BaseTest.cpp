#include "BaseTest.h"

BaseTest::BaseTest(){}
BaseTest::~BaseTest(){}

bool BaseTest::doAllTests(std::ostream& out){
	bool successAll= true;
	uint testsSucceeded= 0, testIndex= 0;
	for(const auto& test : tests){
		std::string message;
		out << "Test " << testIndex << ": ";
		bool success;
		if(success = test(message)){
			testsSucceeded++;
			out << "\033[1;32m[OK";
		}
		else
			out << "\033[1;31m[KO";
		out << "] \033[0m"<< message << std::endl;
		testIndex++;
		successAll &= success;
	}
	if(successAll)
		out << "\033[1;32m[OK";
	else
		out << "\033[1;31m[KO";
	out << "] \033[0m" << testsSucceeded << " tests out of " << tests.size() << " have been executed successfuly" << std::endl;
	return successAll;
}