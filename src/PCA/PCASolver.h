#ifndef PCASOLVER_H
#define PCASOLVER_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

class PCASolver{
public:
	explicit PCASolver (const Eigen::MatrixXd& Z);

	/**
	 * @brief computePCA
	 * Generate a matix in latent space by applying the PCA method
	 * on the _Z matrix using the provided member functions
	 * @param latentDim
	 * Latent dimension (2 by default)
	 * @return
	 * A matrix "reduced" in latent dimension
	 */
	Eigen::MatrixXd computePCA (uint latentDim=2) const;
private:
	/**
	 * @brief _Z
	 * A matrix containing all the BRDFs (extremely large dimension)
	 */
	const Eigen::MatrixXd& _Z;
	/**
	 * @brief computeCovariance
	 * Takes one matrix and computes its corresponding Variance-
	 * Covariance matrix "Cov(M)"
	 * @return
	 * Variance-Covariance Matrix
	 */
	Eigen::MatrixXd computeCovariance (const Eigen::MatrixXd&) const;

	/**
	 * @brief computeEigenValues
	 * Takes one matrix and computes its eigenvalues
	 * @return
	 * Eigenvalues
	 */
	Eigen::MatrixXd computeEigenValues (const Eigen::MatrixXd&) const;

	/**
	 * @brief computeEigenVectors
	 * Takes one matrix and computes its eigenvectors
	 * @return
	 * Eigenvectors
	 */
	Eigen::MatrixXd computeEigenVectors (const Eigen::MatrixXd&) const;
};

#endif // PCASOLVER_H
