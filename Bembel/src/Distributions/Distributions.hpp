
// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_DISTR_H_
#define BEMBEL_DISTR_H_


#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

/*
  We need a functor that can pretend it's const,
  but to be a good random number generator
  it needs mutable state.
*/
namespace Eigen {
namespace internal {
template<typename Scalar>
struct scalar_normal_dist_op
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

/*
  Draw nn samples from a size-dimensional normal distribution
  with a specified mean and covariance
*/

Eigen::MatrixXd generate_gaussian_samples(int nn, const Eigen::VectorXd mean,const Eigen::MatrixXd &covar, const int seed = 1)
{
  int size = mean.size(); // Dimensionality (rows)
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(seed); // Seed the rng


  Eigen::MatrixXd normTransform(size,size);

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise()
                           + mean;

  return samples;
}


Eigen::MatrixXd generate_gaussian_samples(int nn, const Eigen::VectorXd mean,const Eigen::SparseMatrix<double> &covar, const int seed = 1)
{
  int size = mean.size(); // Dimensionality (rows)
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(seed); // Seed the rng


  Eigen::MatrixXd normTransform(size,size);

  //Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > cholSolver;
  cholSolver.compute(covar);

  // We can only use the cholesky decomposition if
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise()
                           + mean;

  return samples;
}



Eigen::MatrixXd generate_gaussian_samples_sparse(int nn, const Eigen::VectorXd mean,const Eigen::SparseMatrix<double> &normTransform, const int seed = 1)
{
  int size = mean.size(); // Dimensionality (rows)
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(seed); // Seed the rng



  Eigen::MatrixXd samples = (normTransform
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise()
                           + mean;

  return samples;
}

Eigen::SparseMatrix<double> generate_sparse_norm_transform(const Eigen::SparseMatrix<double> &covar)
{
  int size = covar.rows(); // Dimensionality (rows)

  Eigen::SparseMatrix<double> normTransform(size,size);


  //Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > cholSolver;
  cholSolver.compute(covar);


  return cholSolver.matrixL();
}

Eigen::MatrixXd generate_gaussian_samples_sparse_covariance(int nn, const Eigen::VectorXd mean,const Eigen::SparseMatrix<double> &covar, const int seed = 1)
{
  Eigen::SparseMatrix<double> L = generate_sparse_norm_transform(covar);

  return generate_gaussian_samples_sparse(nn, mean,L, seed);
}

#endif
