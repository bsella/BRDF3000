#include "BRDFReaderTest.h"
#include "../src/BRDFReader/BRDFReader.h"

BRDFReaderTest::BRDFReaderTest(): BaseTest("BRDFReader"){
    addTest(&readBRDF, "Read BRDF", "../tests/data/BRDFReader_data/inputTest.txt", "../tests/data/BRDFReader_data/brdf_output.txt");
}

std::istringstream BRDFReaderTest::readBRDF(std::istream& istr){
    ChefDevr::BRDFReader reader;
    std::string pathFile_brdf;
    unsigned int num_coefficients;

    istr >> pathFile_brdf;
    istr >> num_coefficients;

    ChefDevr::Vector<double> brdf = reader.read_brdf<double>(num_coefficients, pathFile_brdf.c_str());

    std::stringstream ret;
    ret.precision(20);
    ret << brdf.topRows<100>();

    return std::istringstream(ret.str());
}