#include "BRDFReaderTest.h"

#include "BRDFReader/BRDFReader.h"


using Scalar = long double;


BRDFReaderTest::BRDFReaderTest(): BaseTest("BRDFReader"){
    addTest(&readBRDF, "Read BRDF", "../tests/data/BRDFReader/inputBRDF.txt", "../tests/data/BRDFReader/brdf_output.txt");
    addTest(&createZ, "Create Z", "../tests/data/BRDFReader/inputSetBRDFs.txt", "../tests/data/BRDFReader/setBRDF_output.txt");
    addTest(&createZZt_centered, "Create ZZt centered", "../tests/data/BRDFReader/inputSetBRDFs.txt", "../tests/data/BRDFReader/ZZt_centered.txt");
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


std::istringstream BRDFReaderTest::createZZt_centered(std::istream& istr){
    ChefDevr::BRDFReader reader;
    std::string pathBRDFs;

    istr >> pathBRDFs;

    const ChefDevr::Matrix<Scalar> ZZt_centered = reader.createZZt_centered<Scalar>(pathBRDFs.c_str());

    std::stringstream ret;
    ret.precision(100);
    ret << ZZt_centered;

    return std::istringstream(ret.str());
}