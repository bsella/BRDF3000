/**
 * @file BRDFReader.hpp
 */

#include "Parametrisation/MERLReader.h"


namespace ChefDevr {
    using namespace Eigen;
    using namespace std;


    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZ(const char *fileDirectory) {
        extract_brdfFilePaths(fileDirectory);

        const auto num_brdfs = brdf_filePaths.size();
        Matrix<Scalar> Z{num_brdfs, MERLReader::num_coefficientsBRDF};

#pragma omp parallel for
        for (unsigned int i = 0; i < num_brdfs; ++i) {
            Z.row(i) = read_brdf<Scalar>(i);
        }

        return Z;
    }

    template <typename Scalar>
    Matrix<Scalar> BRDFReader::createZZt_centered(const char *fileDirectory, RowVector<Scalar> &meanBRDF) {
        extract_brdfFilePaths(fileDirectory);

        const auto num_brdfs = brdf_filePaths.size();
        Matrix<Scalar> ZZt_centered{num_brdfs, num_brdfs};
        meanBRDF.resize(MERLReader::num_coefficientsBRDF);

        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector<Scalar> brdf_first = read_brdf<Scalar>(i);

            ZZt_centered(i, i) = brdf_first.dot(brdf_first);
            meanBRDF += brdf_first;

#pragma omp parallel for
            for (unsigned int j = i + 1; j < num_brdfs; ++j) {
                const RowVector<Scalar> brdf_second = read_brdf<Scalar>(j);
                const Scalar coefficient = brdf_first.dot(brdf_second);
                ZZt_centered(i, j) = ZZt_centered(j, i) = coefficient;
            }
        }

        meanBRDF /= num_brdfs;

        RowVector<Scalar> brdf_brdfMean{num_brdfs};

#pragma omp parallel for
        for (unsigned int i = 0; i < num_brdfs; ++i) {
            const RowVector<Scalar> brdf = read_brdf<Scalar>(i);
            brdf_brdfMean(i) = brdf.dot(meanBRDF);
        }

        ZZt_centered.colwise() -= brdf_brdfMean.transpose();
        ZZt_centered.rowwise() -= brdf_brdfMean;
        ZZt_centered = ZZt_centered.array() + meanBRDF.dot(meanBRDF);

        return ZZt_centered;
    }

    template<typename Scalar>
    RowVector <Scalar> BRDFReader::read_brdf(unsigned int index_brdf) {
        const char* path(brdf_filePaths[index_brdf].c_str());
        return MERLReader::read_brdf<Scalar>(path);
    }

} // namespace ChefDevr
