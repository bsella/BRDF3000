#include "Albedo.h" 
#include <cmath>


/**
 * @brief Albedo.cpp
 */

namespace ChefDevr
{
    void Albedo::computeAlbedo (const ResampledBRDF& brdf, Color& albedo)
    {
        computeAlbedoOpenMP(brdf, albedo);
    }
    
    void Albedo::computeAlbedoOpenMP (const ResampledBRDF& brdf, Color& albedo)
    {
        const double one_over_pi(0.31830988618);
        double r(0), g(0), b(0);
        using thnum = brdf.thetaNum;
        using phnum = brdf.phiNum;
        using thstep = brdf.thetaStep;
        using phstep = brdf.phiStep;
        #pragma omp parallel
        {
            unsigned int thondex, thindex, phondex, phindex;
            double tho, thi, costhi, coscos;
            #pragma omp for reduction(+:r,g,b)
            for (thindex=0; thindex < thnum; ++thindex)
            {
                thi = thstep*thindex;
                costhi = std::cos(thi);
                for (phindex=0; phindex < phnum; ++phindex)
                {
                    for (thondex=0; thondex < thnum; ++thondex)
                    {   
                        tho = thstep*thondex;
                        coscos = std::cos(tho)*costhi;
                        for (phondex=0; phondex < phnum; ++phondex)
                        {
                            r += brdf.data[ thindex*phnum*thnum*phnum*3 + 
                                            phindex*thnum*phnum*3 +
                                            thondex*phnum*3 +
                                            phondex*3 ] * coscos;
                            
                            g += brdf.data[ thindex*phnum*thnum*phnum*3 + 
                                            phindex*thnum*phnum*3 +
                                            thondex*phnum*3 +
                                            phondex*3 + 1] * coscos;
                                            
                            b += brdf.data[ thindex*phnum*thnum*phnum*3 + 
                                            phindex*thnum*phnum*3 +
                                            thondex*phnum*3 +
                                            phondex*3 + 2] * coscos;
                        }
                    }
                }
            }
        }
        albedo.r = r*one_over_pi;
        albedo.g = g*one_over_pi;
        albedo.b = b*one_over_pi;
    }
    
    void Albedo::computeAlbedoCuda (const ResampledBRDF& brdf, Color& albedo)
    {
    }
} // namespace ChefDevr
