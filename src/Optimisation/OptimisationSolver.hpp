#include "OptimisationSolver.h"
#include "../Parametrisation/Parametrisation.h"

#include <cmath>


namespace ChefDevr
{
    template <typename Scalar>
    OptimisationSolver<Scalar>::OptimisationSolver(
        Scalar _minStep,
        Matrix<Scalar>& _Z,
        const unsigned int _latentDim) :
        
        minStep(_minStep),
        nb_data(_Z.cols()),
        Z(_Z),
        latentDim(_latentDim)
    {
        step = reduceStep;
        ZZt = Z*Z.transpose();
    }
    
    template <typename Scalar>
    typename OptimisationSolver<Scalar>::OptiResult OptimisationSolver<Scalar>::optimizeMapping ()
    {
        X = computePCA();
        // Compute K
        // Compute detK
        // Compute K_minus1
        // Compute cost
        centerMat<Scalar>(Z);
        costval = cost(K_minus1, detK);
        
        Vector<Scalar> new_X;
        
        // while ...
        do
        {
            X = exploratoryMove();
            while(false)
            {
                patternMove(new_X);
            }
        }while(false/*change this of course*/);
        
        return OptiResult();
    }
    
    template <typename Scalar>
    Scalar OptimisationSolver<Scalar>::cost(const Matrix<Scalar>& K_minus1, const Scalar& detK) const
    {
        Scalar trace(0);
        // Compute trace of K_minus1 * ZZt
        # pragma omp parallel for reduction(+:trace)
        for (unsigned int i = 0; i < ZZt.cols(); ++i){
            trace += K_minus1.row(i).dot(ZZt.col(i));
        }
        return Scalar(0.5) * nb_data * std::log(detK) + Scalar(0.5) * trace;
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::exploratoryMove ()
    {
        const auto& nbcoefs(X.rows());
        Scalar new_costval, new_detK;
        Vector<Scalar> new_X(X);
        Vector<Scalar> cov_vector(nb_data), new_cov_vector(nb_data);
        Matrix<Scalar> new_K_minus1(K_minus1.rows(), K_minus1.cols());
        Vector<Scalar> X_move(nbcoefs);
        unsigned int lv_num;
        
        for (unsigned int i(0); i < nbcoefs;++i)
        {
            lv_num = i/latentDim;
            X_move[i] = step;
            new_X[i] += step;
            computeCovarianceVector(new_cov_vector, new_X, lv_num);
            
            // Update K_minus1 and detK with Sherman-Morisson formula
            updateInverse(K_minus1, new_K_minus1, lv_num, new_cov_vector);
            new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
            // Update costval
            new_costval = cost(new_K_minus1, new_detK);
            
            if (new_costval > costval){
                X_move[i] = -step;
                new_X[i] -= Scalar(2)*step;
                new_cov_vector = computeCovVector(new_X, lv_num);
                
                // Update K_minus1 and detK with Sherman-Morisson formula
                updateInverse(K_minus1, new_K_minus1, lv_num, new_cov_vector);
                new_detK = updateDeterminant(new_K_minus1, lv_num, new_cov_vector);
                // Update cost
                new_costval = cost(new_K_minus1, new_detK);
                
                if (new_costval > costval){
                    X_move[i] = Scalar(0);
                    new_X[i] += step;
                }
                else{
                    costval = new_costval;
                    // cost has changed -> keep new_K_minus1 and new_detK
                    K_minus1 = new_K_minus1;
                    detK = new_detK;
                }
            }
            else{
                costval = new_costval;
                // cost has changed -> keep new_K_minus1 and new_detK
                K_minus1 = new_K_minus1;
                detK = new_detK;
            }
        }
        // X = X+X_move; <=>
        X = new_X;
        return X_move;
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::computeInverse (
        const Matrix<Scalar>& old_K_minus1,
        Matrix<Scalar>& new_K_minus1,
        unsigned int lv_num,
        Vector<Scalar>& cov_vector) const
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    Scalar OptimisationSolver<Scalar>::computeDeterminant (
        const Matrix<Scalar>& new_K_minus1,
        unsigned int lv_num,
        Vector<Scalar>& cov_vector) const
    {
        return Matrix<Scalar>();
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::patternMove (Vector<Scalar>& new_X, Matrix<Scalar>& new_K_minus1) const
    {
    }
    
    template <typename Scalar>
    Vector<Scalar> OptimisationSolver<Scalar>::computePCA ()
    {
        unsigned char i;
        auto B(Z.transpose()*Z);
        
        // Use EigenSolver to compute eigen values and vectors
        auto solver(Eigen::EigenSolver<Matrix<Scalar>>(B, true));
        solver.compute();
        auto eigenValues(solver.eigenvalues());
        auto eigenVectors(solver.eigenvectors());
        
        // Use std::sort to sort an indices vector in the same way eigen values should be ordered
        // initialize original index locations
        std::vector<size_t> idx(eigenValues.size());
        std::iota(idx.begin(), idx.end(), 0);
        // sort indexes based on comparing values in eigenValues
        std::sort(idx.begin(), idx.end(),
                [&eigenValues](size_t i1, size_t i2) {return eigenValues[i1] > eigenValues[i2];});
        
        // Build D the diagonal matrix with ordered eigen values in diagonal
        auto D(Matrix<Scalar>::Zero(latentDim, latentDim));
        for (i=0; i<latentDim; ++i)
        {
            // use sorted indices to retrieve eigen values in descending order
            D(i,i) = std::sqrt(eigenValues[idx[i]]);
        }
        
        // Build V the transposed matrix of ordered eigen vectors 
        Matrix<Scalar> V(latentDim, eigenVectors.cols());
        for (i=0; i<latentDim; ++i)
        {
            // use sorted indices to retrieve eigen vectors in descending order
            V.row(i) = eigenVectors.col(idx[i]);
        }
        
        // Latent variables for each BRDFs are now in column
        auto X(D*V);
        // Output X as column vector
        return X.reshaped();
    }
} // namespace ChefDevr
