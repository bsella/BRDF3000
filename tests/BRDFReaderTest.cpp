#include "BRDFReaderTest.h"

#include "BRDFReader/BRDFReader.h"


using Scalar = long double;


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

    const ChefDevr::RowVector<Scalar> brdf = reader.read_brdf<Scalar>(num_coefficients, pathFile_brdf.c_str());

    std::stringstream ret;
    ret.precision(20);
    ret << brdf.leftCols<100>();

    return std::istringstream(ret.str());
}

std::istringstream BRDFReaderTest::createZ(std::istream& istr){
    ChefDevr::BRDFReader reader;
    std::string pathBRDFs;

    istr >> pathBRDFs;

    const ChefDevr::Matrix<Scalar> Z = reader.createZ<Scalar>(pathBRDFs.c_str());

    std::stringstream ret;
    ret.precision(20);
    ret << Z.leftCols<100>();

    return std::istringstream(ret.str());
}