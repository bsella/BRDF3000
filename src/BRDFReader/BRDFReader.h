#ifndef BRDFREADER_H
#define BRDFREADER_H

/**
 * @file BRDFReader.h
 */

#include <cstdio>
#include <vector>
#include "Parametrisation/types.h"
#include "../tests/BRDFReaderTest.h"


namespace ChefDevr {
    
    /**
    * @brief This class is used to read all the BRDF references and to sample them
    */
    class BRDFReader {
    public:

        BRDFReader() = default;
        ~BRDFReader() = default;

        /**
        * @brief Read all the BRDFs stored in a given directory
        * @param fileDirectory the path of the directory where all the BRDFs are stored
        * @return Non-centered Z BRDFs data matrix where each row represents a BRDF
        *
        * Initializes the list of BRDFs filePaths and filenames in the order in which they were read.
        */
        template<typename Scalar>
        Matrix<Scalar> createZ(const char *fileDirectory);

        /**
        * @brief Creates the centered ZZt matrix
        * @param[in] fileDirectory the path of the directory where all the BRDFs are stored
        * @param[out] meanBRDF The mean BRDF of the brdfs set
        * @return the centered ZZt matrix
        *
        * Initializes the list of BRDFs filePaths and filenames in the order in which they were read.
        *
        * As the set of BRDFs can be heavy, Z is not stored entirely inside the RAM.
        * Thus, a problem is not likely to occur if the RAM is too small compared to the set of BRDFs.
        */
        template <typename Scalar>
        Matrix<Scalar> createZZt_centered(const char *fileDirectory, RowVector<Scalar> &meanBRDF);

        /**
         * @return the list of BRDF filenames in the order in which they were read
         */
        inline const std::vector<std::string>& getBRDFFilenames () const {
            return brdf_filenames;
        }

        /**
         * @return the list of BRDF filePaths in the order in which they were read
         */
        inline const std::vector<std::string>& getBRDFFilePaths() const {
            return brdf_filePaths;
        }

        class BRDFReaderError : public std::runtime_error {
        public:
            explicit BRDFReaderError(const std::string& msg) :
                std::runtime_error(msg){}
        };

    private:

        /* ------------*/
        /* Attributes */
        /* ------------*/

        /**
         * @brief the list of BRDF filenames in the order in which they were read
         */
        std::vector<std::string> brdf_filenames;

        /**
         * @brief the list of BRDF filePaths in the order in which they were read
         */
        std::vector<std::string> brdf_filePaths;

        /**
         * @brief Initializes the list of BRDFs filenames and the list of BRDF filePaths in the order in which they are read.
         * @param fileDirectory the path of the directory where all the BRDFs are stored
         */
        void extract_brdfFilePaths(const char *fileDirectory);

        /**
         * @brief Read a BRDF from a file
         * @param index_brdf Index of the brdf to read
         * @return All the coefficients of a BRDF as a vector of scalars
         *
         * If the file is not found, returns an error
         */
        template<typename Scalar>
        RowVector <Scalar> read_brdf(unsigned int index_brdf);

        /* ------------*/
        /* Friends */
        /* ------------*/

        friend BRDFReaderTest;
    };
} // namespace ChefDevr

#include "BRDFReader.hpp"

#endif // BRDFREADER_H
