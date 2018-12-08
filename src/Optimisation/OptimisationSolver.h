#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include <Eigen/Core>

/// TODO: make template or macro to change floating point precision (float/double)

class OptimisationSolver {
public:
	explicit OptimisationSolver (const Eigen::MatrixXd& Z);
	/**
	 * @brief computeOptimisation
	 * @param minStep
	 * Size of the step used on each iteration
	 * @return
	 * Optimised matrix
	 */
	Eigen::VectorXd computeOptimisation (float minStep) const;
private:
	const Eigen::MatrixXd& _Z;
	static constexpr float reduceStep = .5f;

	float cost (uint d, Eigen::MatrixXd& K) const;

	/**
	 * @brief computeExplorationDisplacement
	 * @param X
	 * @param L
	 * @return
	 */
	Eigen::VectorXd computeExplorationDisplacement (const Eigen::VectorXd& X, float& L) const;

	/**
	 * @brief computeLearnDisplacement
	 * @param X
	 * @param Xpred
	 * @param L
	 * @return
	 */
	Eigen::VectorXd computeLearnDisplacement (const Eigen::VectorXd& X, const Eigen::VectorXd& Xpred, float& L) const;
};

#endif // OPTIMISATIONSOLVER_H
