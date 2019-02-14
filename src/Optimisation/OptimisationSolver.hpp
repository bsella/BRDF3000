#include "OptimisationSolver.h"
#include "Parametrisation/Parametrisation.h"

/**
 * @file OptimisationSolver.h
 */

#include <cmath>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <numeric>

namespace ChefDevr
{
    template <typename Scalar>
    OptimisationSolver<Scalar>::OptimisationSolver(
        const long num_BRDFCoefficients,
        Scalar _minStep,
        const Matrix<Scalar>& _ZZt,
        const unsigned int _latentDim) :
        minStep(_minStep),
        step(step0),
        nb_data(_ZZt.rows()),
        num_BRDFCoefficients{num_BRDFCoefficients},
        ZZt(_ZZt),
        latentDim(_latentDim),
        X(_ZZt.rows()*latentDim),
        X_move(_ZZt.rows()*latentDim),
        K_minus1(_ZZt.rows(), _ZZt.rows()){}
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::optimizeMapping ()
    {
        Vector<Scalar> new_X(latentDim*nb_data);
        Matrix<Scalar> new_K_minus1(nb_data, nb_data);
        Scalar new_costval, new_detK;
        
        // Init X
        initX(ZZt);
        
        // Compute K
        // (We use K_minus1 to store it because we don't need K anymore after)
        for (unsigned int i=0; i < nb_data; ++i)
        {
            
            computeCovVector<Scalar>(K_minus1.col(i).data(), X,
                             X.segment(latentDim*i, latentDim),
                             latentDim, nb_data);
        }
        
        // Compute detK
        detK = K_minus1.determinant();
        // Compute K_minus1 from K
        K_minus1 = K_minus1.inverse().eval();
        // Init cost
        cost(costval, K_minus1, detK);
        
        #ifdef DEBUG
        std::cout << "cost :" << costval << std::endl;
        std::cout << "detK :" << detK << std::endl;
        std::cout << "K_minus1 :" << std::endl << K_minus1 << std::endl;
        std::cout << "X :" << std::endl << X << std::endl << std::endl;
        #endif
        // Optimisation loop
        do
        {
            #ifdef DEBUG
            std::cout << " ================ Exploratory Move ================ " << std::endl;
            #endif
            if(exploratoryMove())
            {
                #ifdef DEBUG
                std::cout << "~~~ exploratory moved ~~~" << std::endl;
                std::cout << "cost :" << costval << std::endl;
                std::cout << "detK :" << detK << std::endl;
                std::cout << "K_minus1 :" << std::endl << K_minus1 << std::endl;
                std::cout << "X :" << std::endl << X << std::endl << std::endl;
                assert(X.minCoeff() > Scalar(-1) && X.maxCoeff() < Scalar(1));
                #endif
                do
                {   
                    if(patternMove(new_X, new_K_minus1, new_detK))
                    {
                        #ifdef DEBUG
                        std::cout << "~~~ pattern moved ~~~" << std::endl;
                        assert(X.minCoeff() > Scalar(-1) && X.maxCoeff() < Scalar(1));
                        #endif 
                        cost(new_costval, new_K_minus1, new_detK);
                        if (new_costval < costval)
                        {
                            #ifdef DEBUG
                            std::cout << "    effective ! :)" << std::endl;
                            #endif
                            costval = new_costval;
                            X.noalias() = new_X; 
                            K_minus1.noalias() = new_K_minus1;
                            detK = new_detK;
                            #ifdef DEBUG
                            std::cout << "cost :" << costval << std::endl;
                            std::cout << "detK :" << detK << std::endl;
                            std::cout << "K_minus1 :" << std::endl << K_minus1 << std::endl;
                            std::cout << "X :" << std::endl << X << std::endl << std::endl;
                            #endif
                        }
                        else
                        {
                            #ifdef DEBUG
                            std::cout << "____ineffective____" << std::endl;
                            #endif
                            break;
                        }
                    }
                    else
                    {
                        #ifdef DEBUG
                        std::cout << "~~~ pattern didn't move : out of bounds ~~~" << std::endl;
                        #endif 
                        break;
                    }
                }
                while(true);
            }
            step *= reduceStep;
        }while(step >= minStep);
        #ifdef DEBUG
        std::cout << " ================ OPTIMIZATION HAS CONVERGED ================ " << std::endl;
        std::cout << "cost :" << costval << std::endl;
        std::cout << "detK :" << detK << std::endl;
        std::cout << "K_minus1 :" << std::endl << K_minus1 << std::endl;
        std::cout << "X :" << std::endl << X << std::endl << std::endl;
        #endif
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::cost(Scalar& cost, const Matrix<Scalar>& K_minus1, const Scalar& detK) const
    {
        Scalar trace(0);
        // Compute trace of K_minus1 * ZZt
        # pragma omp parallel for reduction(+:trace)
        for (unsigned int i = 0; i < ZZt.cols(); ++i)
        {
            trace += K_minus1.row(i).dot(ZZt.col(i));
        }
        cost = Scalar(0.5) * (num_BRDFCoefficients * log(detK) + trace);
    }
    
    template <typename Scalar>
    bool OptimisationSolver<Scalar>::exploratoryMove ()
    {
        const auto& nbcoefs(X.rows());
        Scalar new_detK, new_costval(std::numeric_limits<Scalar>::infinity());
        Vector<Scalar> cov_vector(nb_data), diff_cov_vector(nb_data);
        Matrix<Scalar> new_K_minus1(K_minus1.rows(), K_minus1.cols());
        Vector<Scalar> X_move(nbcoefs);
        unsigned int lv_num;
        bool moved(false);
        
        for (unsigned int i(0); i < nbcoefs;++i)
        {
            new_costval = std::numeric_limits<Scalar>::infinity();
            lv_num = i/latentDim;
            computeCovVector<Scalar>(cov_vector.data(), X,
                             X.segment(latentDim*lv_num, latentDim),
                             latentDim, nb_data);
            
            X[i] += step;
            X_move[i] = step;
            if ( X[i] < Scalar(1)) // latent variable constraint
            {
                computeCovVector<Scalar>(diff_cov_vector.data(), X,
                                 X.segment(latentDim*lv_num,latentDim),
                                 latentDim, nb_data);
                diff_cov_vector = diff_cov_vector - cov_vector;
                
                // Update K_minus1 and detK with Sherman-Morisson formula
                shermanMorissonUpdate(K_minus1, new_K_minus1,
                                      detK, new_detK,
                                      lv_num, diff_cov_vector);
                // Update costval
                cost(new_costval, new_K_minus1, new_detK);
            }
            if (new_costval > costval)
            {
                X[i] -= Scalar(2)*step;
                X_move[i] = -step;
                if (X[i] > Scalar(-1)) // latent variable constraint
                {
                    computeCovVector<Scalar>(diff_cov_vector.data(), X,
                                     X.segment(latentDim*lv_num, latentDim),
                                     latentDim, nb_data);
                    diff_cov_vector = diff_cov_vector - cov_vector;
                    
                    // Update K_minus1 and detK with Sherman-Morisson formula
                    shermanMorissonUpdate(K_minus1, new_K_minus1, 
                                          detK, new_detK,
                                          lv_num, diff_cov_vector);
                    // Update cost
                    cost(new_costval, new_K_minus1, new_detK);
                }   
                if (new_costval > costval)
                {
                    X[i] += step;
                    X_move[i] = Scalar(0);
                }
                else
                {
                    costval = new_costval;
                    // cost has changed -> keep new_K_minus1 and new_detK
                    K_minus1.noalias() = new_K_minus1;
                    detK = new_detK;
                    moved = true;
                }
            }
            else
            {
                costval = new_costval;
                // cost has changed -> keep new_K_minus1 and new_detK
                K_minus1.noalias() = new_K_minus1;
                detK = new_detK;
                moved = true;
            }
        }
        return moved;
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::shermanMorissonUpdate (
        const Matrix<Scalar>& old_K_minus1,
        Matrix<Scalar>& new_K_minus1,
        const Scalar& old_detK,
        Scalar& new_detK,
        unsigned int lv_num,
        Vector<Scalar>& diff_cov_vector) const
    {
        #ifdef DEBUG
        assert(old_detK != 0.0);
        #endif
        
        Scalar centerCoeff(diff_cov_vector[lv_num]);
        
        // ===== One row modification =====
        Scalar dotp1(diff_cov_vector.dot(old_K_minus1.col(lv_num)) + Scalar(1));
        // Determinant update
        new_detK = dotp1 * old_detK;
        
        // Inverse update 
        new_K_minus1.noalias() =
            (old_K_minus1 - ( ((old_K_minus1.col(lv_num)*(diff_cov_vector.transpose()/dotp1))*old_K_minus1))
            ).eval();
        
        // ===== One column modification =====
        diff_cov_vector[lv_num] = Scalar(0);
        
        dotp1 = new_K_minus1.row(lv_num).dot(diff_cov_vector) + Scalar(1);
        
        // Determinant update
        new_detK = dotp1 * new_detK;
        
        // Inverse update
        new_K_minus1 = 
            (new_K_minus1 - ( new_K_minus1*(diff_cov_vector/dotp1)*new_K_minus1.row(lv_num))
            ).eval();
        
        diff_cov_vector[lv_num] = centerCoeff;
    }
    
    template <typename Scalar>
    bool OptimisationSolver<Scalar>::patternMove (Vector<Scalar>& new_X, Matrix<Scalar>& new_K_minus1, Scalar& new_detK) const
    {
        unsigned int i;
        new_X = X + X_move;
        if (new_X.minCoeff() > Scalar(-1) && new_X.maxCoeff() < Scalar(1))
        {
            // Compute new_K (in new_K_minus1 so we don't have to allocate more memory)
            for (i=0; i<nb_data; ++i){
                computeCovVector<Scalar>(new_K_minus1.col(i).data(), new_X,
                                 new_X.segment(i*latentDim, latentDim),
                                 latentDim, nb_data);
            }
            // Compute new_detK
            new_detK = new_K_minus1.determinant();
            // Compute new_K_minus1
            new_K_minus1 = new_K_minus1.inverse().eval();
            return true;
        }
        return false;
    }
    
    template <typename Scalar>
    void OptimisationSolver<Scalar>::initX (const Matrix<Scalar>& ZZt)
    {
        // Use EigenSolver to compute eigen vector/values decomposition
        Eigen::EigenSolver<Matrix<Scalar>> solver(ZZt, true);

        const Vector<Scalar>& eigenValues(solver.eigenvalues().real());
        const Matrix<Scalar>& eigenVectors(solver.eigenvectors().real());
         
        // Use std::sort to sort an indices vector in the same way eigen values should be ordered
        // initialize original index locations
        std::vector<size_t> idx(eigenValues.size());
        std::iota(idx.begin(), idx.end(), 0);
        // sort indexes based on comparing values in eigenValues
        std::function<size_t(size_t, size_t)> cmpIndices(
            [&eigenValues](size_t i1, size_t i2) {
                return eigenValues[i1] > eigenValues[i2];
            }
        );

        std::sort(idx.begin(), idx.end(),cmpIndices);
        
        // Build V the transposed matrix of ordered eigen vectors 
        Matrix<Scalar> V(latentDim, eigenVectors.cols());
        # pragma omp parallel for
        for (unsigned int i=0; i<latentDim; ++i)
        {
            // use sorted indices to retrieve eigen vectors in descending order
            V.row(i).noalias() = eigenVectors.col(idx[i]).transpose();
        }
        
        // X as column vector
        X  = Eigen::Map<Vector<Scalar>>(V.data(), latentDim*nb_data, 1);

        // normalize X
        X = X / (std::numeric_limits<Scalar>::epsilon() + std::max(abs(X.maxCoeff()), abs(X.minCoeff())));
    }
} // namespace ChefDevr
