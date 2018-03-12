

#ifndef RELIABLEDIAMONDNORM_UTILS_H
#define RELIABLEDIAMONDNORM_UTILS_H

#include <cmath> // std::abs
#include <random>

#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>


#include <tomographer/tools/cxxutil.h>
#include <tomographer/tools/eigenutil.h>
#include <tomographer/tools/loggers.h>



typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> DynStride;




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






/** \brief A unitarily invariant random isometry
 *
 * Returns \a N randomly chosen orthogonal vectors packed together into columns of a
 * matrix \a V.  This makes \a V an isometry.
 *
 * The distribution is unitarily invariant, i.e. for any unitary matrix \a U, the outcome
 * \a U*V has the same probability as \a V.
 *
 * \param d the dimension of the space in which to generate random orthogonal vectors
 * \param N the number of orthogonal vectors to generate
 *
 * \param rng a std::random random number generator (such as std::mt19937)
 *
 * \param logger a reference to a logger (\ref pageLoggers) where we can log what we're
 *        doing.
 */
template<typename MatrixType, typename Rng, typename Logger>
inline MatrixType randomIsometry(int d, int N, Rng & rng, Logger & baselogger)
{
  auto logger = Tomographer::Logger::makeLocalLogger(TOMO_ORIGIN, baselogger);

  logger.longdebug("d=%d, N=%d", d, N);
  
  typedef typename MatrixType::Scalar Scalar;
  //typedef Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1> VectorType;

  // first, get a matrix of normally distributed random numbers
  MatrixType V(d,N);

  std::normal_distribution<> normdist(0.0, 1.0);
  V = Tomographer::Tools::denseRandom<MatrixType>(rng, normdist, V.rows(), V.cols());

  // perform Gram-Schmidt orthogonalization

  for (int j = 0; j < N; ++j) {

    auto v = V.col(j);

    for (int k = 0; k < j; ++k) {
      auto p = V.col(k).adjoint() * v;
      v -= p*V.col(k);
    }

    v /= v.norm();
  }

  // this should go into a test case
  tomographer_assert( (V.adjoint()*V - MatrixType::Identity(N,N)).norm()
                      < Eigen::NumTraits<Scalar>::dummy_precision() ) ;

  logger.longdebug([&](std::ostream& str) {
      str << "randomIsometry: got V = \n" << V << "\n"
  	  << "Check: V*V.adjoint() ==\n" << V*V.adjoint() << "\n"
  	  << "Check: V.adjoint()*V ==\n" << V.adjoint()*V;
    });

  return V;
}








template<typename VIsometryType,
         TOMOGRAPHER_ENABLED_IF_TMPL(!VIsometryType::IsRowMajor)>
inline Eigen::Map<const Eigen::Matrix<typename VIsometryType::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>,
                  0, DynStride>
remapIsometryToT(const VIsometryType & V, int dXY)
{
  typedef Eigen::Matrix<typename VIsometryType::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    MatrixType;
  
  // we have: V((jkl),i) |i>_A |jkl>_{BA'B'}
  //
  // we want: T((ij),(kl)) from V((jkl),i)
  //
  // (ij) = i*dJ+j
  // (jkl) = j*dK*dL+k*dL+l
  //
  // matrix is column-major, so:  (jkl),i == i*dJKL+(jkl) == (i*dJ+j)*dKL+(kl)
  //
  return Eigen::Map<const MatrixType,0,DynStride >( V.data(), dXY, dXY,
                                                    DynStride(1,dXY) );
}




template<typename MatrixType>
inline typename Eigen::NumTraits<typename MatrixType::Scalar>::Real
entanglement_fidelity(const MatrixType & E_Choi, const int dimX)
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic,
                        Eigen::internal::traits<MatrixType>::Options> SubMatrixType;

  // E is an unnormalized Choi matrix -> we need to normalize it to a
  // bipartite state corresponding to the process applied onto half a
  // normalized maximally entangled state
  //
  // Also, the overlap of A_XY with the unnormalized maximally entangled ket is
  //
  //     \sum_{i,j} A_XY(dX*i+i,dX*j+j) ,
  //
  // i.e., pick the submatrix with stride dX+1 starting at 0, and calculate
  // the sum of all its elements.

  Eigen::Map<const SubMatrixType,0,DynStride> submatrixentgl(
      E_Choi.data(), dimX, dimX,
      DynStride((dimX+1)*E_Choi.outerStride(), (dimX+1)*E_Choi.innerStride())
      );
  const Scalar overlapPhi =
    submatrixentgl.array().sum(); // DynStride defined in "utils.h"

  //std::cerr << "E_Choi = \n" << E_Choi << "\n";
  //std::cerr << "submatrixentgl = \n" << submatrixentgl << "\n";
  //std::cerr << "overlapPhi = " << overlapPhi << "\n";

  tomographer_assert(std::abs(overlapPhi.imag()) <= Eigen::NumTraits<Scalar>::dummy_precision());

  return overlapPhi.real() / RealScalar(dimX*dimX);
}





#endif
