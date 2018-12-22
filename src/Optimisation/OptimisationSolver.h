#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include <Eigen/Core>

/// TODO: make template or macro to change floating point precision (float/double)

class OptimisationSolver {
public:
	explicit OptimisationSolver (const Eigen::MatrixXd& Z);
	/**
	 * @brief computeOptimisation
	 * Compute an optimised version of _X using the
	 * Hooke Jeeves method
	 * @param minStep
	 * Size of the step used on each iteration
	 * @return
	 * Optimised matrix
	 */
	Eigen::VectorXd computeOptimisation (float minStep) const;
private:
	const Eigen::MatrixXd& _X;
	static constexpr float reduceStep = .5f;

	float cost (uint d, Eigen::MatrixXd& K) const;

	/**
	 * @brief computeExplorationDisplacement
	 * Modify the X vector to find a better solution
	 * add and substract _reduceStep_ from each element
	 * and reevaluate the cost function to check for
	 * better solutions
	 * @return
	 * New X vector, with lower cost function
	 */
	Eigen::VectorXd computeExplorationDisplacement (const Eigen::VectorXd& X) const;

	/**
	 * @brief computeLearnDisplacement
	 * Apply a previously calculated displacement
	 * to each element of X
	 * @return
	 * New modified vector
	 */
	Eigen::VectorXd computeLearnDisplacement (const Eigen::VectorXd& X, const Eigen::VectorXd& Xpred) const;
};

#endif // OPTIMISATIONSOLVER_H
