#ifndef OPTIMISATIONSOLVER_H
#define OPTIMISATIONSOLVER_H

#include <Eigen/Core>

template <typename Scalar>
class OptimisationSolver {
public:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> OptiMatrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> OptiVector;
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic> OptiRowVector;
    
    struct OptiResult
    {
        OptiMatrix K; // Mapping matrix
        OptiMatrix K_minus1; // Inverse mapping matrix
        OptiVector X; // Latent variables vector
    };
    
	OptimisationSolver() = delete;
	~OptimisationSolver() = delete;
    
    /**
	 * @brief  Compute an optimised mapping from BRDFs space to a latent space
     *         whose dimension is specified by dim
	 * @param  Z The BRDFs data matrix
     * @param  dim The dimension of latent space
     * @param  minStep Step beyond which convergence is considered being reached
	 * @return Structure containing optimized mapping, inverse mapping and latent variables
     * 
     * Uses Hook & Jeeves method to solve the optimisation
	 */
	static OptiResult computeOptimisation (
        const OptiMatrix& Z,
        const int dim,
        float minStep);
    
private:
	static constexpr float reduceStep = .5f;

    /**
     * @brief Computes the cost of the inverse mapping defined by K_minus1
     * @param detK Determinant of the forward mapping matrix K
     * @param K_minus1 Inverse mapping matrix
     * @param Z BRDFs data matrix
     * @param Zt Z transposed
     */
	static Scalar cost(
        const Scalar& detK,
        const OptiMatrix& K_minus1,
        const OptiMatrix& Z,
        const OptiMatrix& Zt);

	/**
	 * @brief  Modify the X vector to find a better solution
	 * @return New X vector, with lower cost function
	 * 
     * Adds and substracts _reduceStep_ from each element
	 * and reevaluate the cost function to check for better solutions
	 */
	static OptiResult computeExplorationDisplacement (const OptiResult& optiRes);

	/**
	 * @brief Apply a previously calculated displacement to each element of X
	 * @return New modified vector
	 */
	static OptiResult computeLearnDisplacement (const OptiResult& optiRes);
};

#endif // OPTIMISATIONSOLVER_H
