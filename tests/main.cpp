#include <iostream>

#include "BRDFReaderTest.h"
#include "OptimisationTest.h"
#include "ParametrisationTest.h"

int main(){
    ParametrisationTest pt;
    OptimisationTest ot;
    BRDFReaderTest bt;
    pt.doAllTests(std::cout);
    ot.doAllTests(std::cout);
    bt.doAllTests(std::cout);
    return 0;
}
