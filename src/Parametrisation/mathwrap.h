#include "types.h"

/**
 * @file mathwrap.h
 * This file provides wrapper for math functions to work for both boost and standard types
 */

namespace ChefDevr
{
    using float128 = boost::multiprecision::float128;
    
    
    template <typename Scalar>
    inline Scalar exp(Scalar x)
    {
        return std::exp(x);
    }
    
    template <>
    inline float128 exp<float128>(float128 x)
    {
        return boost::multiprecision::exp(x);
    }
    
    template <typename Scalar>
    inline Scalar log(Scalar x)
    {
        return std::log(x);
    }
    
    template <>
    inline float128 log<float128>(float128 x)
    {
        return boost::multiprecision::log(x);
    }
    
    template <typename Scalar>
    inline Scalar abs(Scalar x)
    {
        return std::abs(x);
    }
    
    template <>
    inline float128 abs<float128>(float128 x)
    {
        return boost::multiprecision::abs(x);
    }
    
} // namespace ChefDevr
