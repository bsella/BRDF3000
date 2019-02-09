#ifndef TYPES__H
#define TYPES__H

/**
 * @file types.h
 * This file defines types used in ChefDevr classes
 */


#ifdef DEBUG
#include <iostream>
#endif


#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/number.hpp>
#include <Eigen/Core>


namespace ChefDevr
{
    /** @tparam Scalar Type used for the scalar numbers stored */
    template<typename Scalar>
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    /** @tparam Scalar Type used for the scalar numbers stored */
    template<typename Scalar>
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    /** @tparam Scalar Type used for the scalar numbers stored */
    template<typename Scalar>
    using RowVector = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;
} // namespace ChefDevr

#endif //TYPES__H
