

#ifndef RELIABLEDIAMONDNORM_UTILS_H
#define RELIABLEDIAMONDNORM_UTILS_H

#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>


#include <tomographer/tools/cxxutil.h>



// get the reduced state on A of rho_AB
template<typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
partial_trace_B(const Eigen::MatrixBase<Derived> & rho_AB, const int dimA)
{
  typedef Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixDynType;

  // dimB is deduced from the matrix size
  const int dimB = rho_AB.rows() / dimA;

  // consistency check: make sure that rho_AB is square and that the dimensions match
  tomographer_assert(rho_AB.rows() == rho_AB.cols());
  tomographer_assert(rho_AB.rows() == dimA*dimB);

  MatrixDynType rho_A = MatrixDynType::Zero(dimA, dimA);
  for (int i = 0; i < dimA; ++i) {
    for (int j = i + 1; j < dimA; ++j) {
      const typename MatrixDynType::Scalar x = rho_AB.block(dimB*i, dimB*j, dimB, dimB).trace();
      rho_A(i,j) = x;
      rho_A(j,i) = std::conj(x);
    }
    rho_A(i,i) = rho_AB.block(dimB*i, dimB*i, dimB, dimB).real().trace();
  }

  return rho_A;
}




// get the Choi matrix from the process matrix
template<typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
Choi_from_process_matrix(const Eigen::MatrixBase<Derived> & rho_XY, const int dimX, const int dimY)
{
  typedef typename Derived::Scalar ComplexScalar;
  typedef Eigen::Matrix<ComplexScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixDynType;

  // calculate rho_X -- reduced state on input system
  const MatrixDynType rho_X = partial_trace_B(rho_XY, dimX);
  
  Eigen::SelfAdjointEigenSolver<MatrixDynType> eigrhoX(rho_X);
  const MatrixDynType rhoXsqrtinvtimesEye =
    Eigen::kroneckerProduct(eigrhoX.operatorInverseSqrt(), MatrixDynType::Identity(dimY,dimY));

  return rhoXsqrtinvtimesEye * rho_XY.template selfadjointView<Eigen::Lower>() * rhoXsqrtinvtimesEye;//.adjoint();
}





#endif
