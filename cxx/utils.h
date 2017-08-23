

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






template<typename IterCountIntType>
struct CtrlConvergedParams
{
  CtrlConvergedParams()
    : enabled(false),
      check_frequency_sweeps(-1),
      max_allowed_unknown(-1),
      max_allowed_unknown_notisolated(-1),
      max_allowed_not_converged(-1),
      max_add_run_iters(std::numeric_limits<double>::quiet_NaN())
  {
  }

  static inline CtrlConvergedParams fromPyDict(py::dict ctrl_converged_params,
                                               tpy::HistogramParams histogram_params)
  {
    CtrlConvergedParams c;
    c.enabled = ctrl_converged_params.attr("get")("enabled", true).cast<bool>();
    c.check_frequency_sweeps =
      ctrl_converged_params.attr("get")("check_frequency_sweeps", 1024).cast<int>();
    c.max_allowed_unknown =
      ctrl_converged_params.attr("get")("max_allowed_unknown",
                                        1+2*histogram_params.num_bins/100).cast<Eigen::Index>();
    c.max_allowed_unknown_notisolated =
      ctrl_converged_params.attr("get")("max_allowed_unknown_notisolated",
                                        1+histogram_params.num_bins/100).cast<Eigen::Index>();
    c.max_allowed_not_converged =
      ctrl_converged_params.attr("get")("max_allowed_not_converged",
                                        1+histogram_params.num_bins/200).cast<Eigen::Index>();
    c.max_add_run_iters = ctrl_converged_params.attr("get")("max_add_run_iters", 1.5).cast<double>();
    return c;
  }

  static inline CtrlConvergedParams fromPyDictWithGilAcq(py::dict ctrl_converged_params,
                                                         tpy::HistogramParams histogram_params)
  {
    py::gil_scoped_acquire gilacq;
    return fromPyDict(std::move(ctrl_converged_params), std::move(histogram_params));
  }

  bool enabled;
  IterCountIntType check_frequency_sweeps;
  Eigen::Index max_allowed_unknown;
  Eigen::Index max_allowed_unknown_notisolated;
  Eigen::Index max_allowed_not_converged;
  double max_add_run_iters;
};





#endif
