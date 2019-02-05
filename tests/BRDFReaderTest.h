#ifndef BRDFREADERTEST_H
#define BRDFREADERTEST_H

#include "BaseTest.h"

class BRDFReaderTest : public BaseTest{
public:
    BRDFReaderTest();

private:
    static std::istringstream createZ(std::istream& istr);
    static std::istringstream readBRDF(std::istream& istr);
};

#endif // BRDFREADERTEST_H