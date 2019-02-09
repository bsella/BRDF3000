#ifndef __MATHWRAP_H
#define __MATHWRAP_H

#include <boost/math/special_functions/fpclassify.hpp>
/**
 * @file mathwrap.h
 * This file provides wrapper for math functions to work for both boost and standard types
 */

using float128 = boost::multiprecision::float128;
namespace ChefDevr
{
    template <typename Scalar>
    inline Scalar exp(const Scalar& x)
    {
        return std::exp(x);
    }
    
    template <>
    inline float128 exp<float128>(const float128& x)
    {
        return boost::multiprecision::exp(x);
    }
    
    template <typename Scalar>
    inline Scalar log(const Scalar& x)
    {
        return std::log(x);
    }
    
    template <>
    inline float128 log<float128>(const float128& x)
    {
        return boost::multiprecision::log(x);
    }
    
    template <typename Scalar>
    inline Scalar abs(const Scalar& x)
    {
        return std::abs(x);
    }
    
    template <>
    inline float128 abs<float128>(const float128& x)
    {
        return boost::multiprecision::abs(x);
    }
} // namespace ChefDevr

inline bool isfinite(const float128& x) { return boost::math::isfinite(x); }

#endif // __MATHWRAP_H
