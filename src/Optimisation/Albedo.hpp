#include "Albedo.h" 
#include <cmath>

#include "../BRDFReader/BRDFReader.h"


/**
 * @brief Albedo.hpp
 */

namespace ChefDevr
{
    template <typename Scalar>
    void Albedo::computeAlbedo (const Vector<Scalar>& brdf, Color& albedo, unsigned int num_sampling)
    {
        computeAlbedoOpenMP(brdf, albedo, num_sampling);
    }

    template <typename Scalar>
    void Albedo::computeAlbedoOpenMP (const Vector<Scalar>& brdf, Color& albedo, unsigned int num_sampling)
    {
        const double one_over_pi(0.31830988618);
        double r(0), g(0), b(0);
        auto& thnum = num_sampling;
        const unsigned int phnum = 4 * num_sampling;
        const double thstep = 0.5 * M_PI / num_sampling;
        const double phstep = 2.0 * M_PI / phnum;
        #pragma omp parallel
        {
            unsigned int thondex, thindex, phondex, phindex;
            double tho, thi, pho, phi;
            #pragma omp for reduction(+:r,g,b)
            for (thindex=0; thindex < thnum; ++thindex)
            {
                thi = thstep*thindex;
                for (phindex=0; phindex < phnum; ++phindex)
                {
                    phi = phstep*phindex;
                    for (thondex=0; thondex < thnum; ++thondex)
                    {   
                        tho = thstep*thondex;
                        for (phondex=0; phondex < phnum; ++phondex)
                        {
                            pho = phstep*thindex;

                            double red, green, blue;
                            BRDFReader::lookup_brdf_val(brdf, thi, phi, tho, pho, red, green, blue);

                            r += red;
                            g += green;
                            b += blue;
                        }
                    }
                }
            }
        }
        albedo.r = r*one_over_pi;
        albedo.g = g*one_over_pi;
        albedo.b = b*one_over_pi;
    }

    template <typename Scalar>
    void Albedo::computeAlbedoCuda (const Vector<Scalar>& brdf, Color& albedo)
    {
    }
} // namespace ChefDevr
