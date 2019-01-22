#ifndef TYPES__H
#define TYPES__H

#include <Eigen/Core>

namespace ChefDevr
{
    template<typename Scalar>
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename Scalar>
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    template<typename Scalar>
    using RowVector = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;
} // namespace ChefDevr

#endif //TYPES__H
