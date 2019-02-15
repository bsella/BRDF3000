#include "Albedo.h"
#include <cmath>

#include "Parametrisation/MERLReader.h"


/**
 * @brief Albedo.hpp
 */

namespace ChefDevr
{
    template <typename Scalar>
    void Albedo::computeAlbedo (
        const RowVector<Scalar>& brdf,
        double& r, double& g, double& b,
        unsigned int num_sampling)
    {
        computeAlbedoOpenMP(brdf, r, g, b, num_sampling);
    }

    template <typename Scalar>
    void Albedo::computeAlbedoOpenMP (
        const RowVector<Scalar>& brdf,
        double& r, double& g, double& b,
        const unsigned int num_sampling)
    {
        const double test(1/((2*M_PI)/3));
        auto& thnum = num_sampling;
        const unsigned int phnum = 4 * num_sampling;
        const double normalizeIrradiance(1./(thnum*thnum*phnum*phnum*M_PI));
        const double thstep = 0.5 * M_PI / num_sampling;
        const double phstep = 2.0 * M_PI / phnum;
        #pragma omp parallel reduction(+:r,g,b)
        {
            unsigned int thondex, thindex, phondex, phindex;
            double tho(0), thi(0), pho(0), phi(0);
            double costhi(0), coscos(0);
            double red, green, blue;
            #pragma omp for
            for (thindex=0; thindex < thnum; ++thindex)
            {
                thi += thstep;
                costhi = std::cos(thi);
                for (phindex=0; phindex < phnum; ++phindex)
                {
                    phi += phstep;
                    for (thondex=0; thondex < thnum; ++thondex)
                    {
                        tho += thstep;
                        coscos = costhi * std::cos(tho);
                        for (phondex=0; phondex < phnum; ++phondex)
                        {
                            pho += phstep;

                            MERLReader::lookup_brdf_val(brdf, thi, phi, tho, pho, red, green, blue);
                            r += red*coscos;
                            g += green*coscos;
                            b += blue*coscos;
                        }
                        pho = 0.0;
                    }
                    tho = 0.0;
                }
                phi = 0.0;
            }
        }
        r = std::min(r*normalizeIrradiance, 1.);
        g = std::min(g*normalizeIrradiance, 1.);
        b = std::min(b*normalizeIrradiance, 1.);
    }

    template <typename Scalar>
    void Albedo::computeAlbedoCuda (
        const RowVector<Scalar>& brdf,
        double& r, double& g, double& b,
        const unsigned int num_sampling)
    {
    }
} // namespace ChefDevr
