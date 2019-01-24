#ifndef BRDF_RECONSTRUCTED_H
#define BRDF_RECONSTRUCTED_H

#include "BRDFMeasuredMERL.h"

/**
 * @file BRDFReconstructed.h
 * A reconstructed BRDF, in the original space
 */


class BRDFReconstructed : public BRDFMeasuredMERL {
public:

    /**
     * @brief Save the BRDF in a file
     *
     * The file's format is the same as the MERL DataBase
     */
    void saveAsBRDF();

private:
};

#endif
