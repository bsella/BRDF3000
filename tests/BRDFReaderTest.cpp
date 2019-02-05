#include "BRDFReaderTest.h"

#include "BRDFReader/BRDFReader.h"



BRDFReaderTest::BRDFReaderTest(): BaseTest("BRDFReader"){
    addTest(&readBRDF, "Read BRDF", "../tests/data/BRDFReader_data/inputBRDF.txt", "../tests/data/BRDFReader_data/brdf_output.txt");
    addTest(&createZ, "Create Z", "../tests/data/BRDFReader_data/inputSetBRDFs.txt", "../tests/data/BRDFReader_data/setBRDF_output.txt");
}

std::istringstream BRDFReaderTest::readBRDF(std::istream& istr){
    ChefDevr::BRDFReader reader;
    std::string pathFile_brdf;
    unsigned int num_coefficients;

    istr >> pathFile_brdf;
    istr >> num_coefficients;

    auto brdf = reader.read_brdf<double>(num_coefficients, pathFile_brdf.c_str());

    std::stringstream ret;
    ret.precision(20);
    ret << brdf.leftCols<100>();

    return std::istringstream(ret.str());
}

std::istringstream BRDFReaderTest::createZ(std::istream& istr){
    ChefDevr::BRDFReader reader;
    std::string pathBRDFs;

    istr >> pathBRDFs;

    const auto Z = reader.createZ<double>(pathBRDFs.c_str());

    std::stringstream ret;
    ret.precision(20);
    ret << Z.leftCols<100>();

    return std::istringstream(ret.str());
}