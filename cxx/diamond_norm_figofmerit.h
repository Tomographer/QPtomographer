
#ifndef DIAMOND_NORM_FIGOFMERIT_H
#define DIAMOND_NORM_FIGOFMERIT_H

#include <Eigen/Eigen>

#include <tomographer/tools/loggers.h>
#include <tomographer/tools/needownoperatornew.h>

#include "utils.h"


template<typename DMTypes_, typename DiamondNormSolverType_>
class DiamondNormToRefValueCalculator
  : public virtual Tomographer::Tools::NeedOwnOperatorNew<typename DMTypes_::MatrixType>::ProviderType
{
public:
  typedef DMTypes_ DMTypes;
  typedef typename DMTypes::MatrixType MatrixType;
  typedef typename DMTypes::MatrixTypeConstRef MatrixTypeConstRef;

  //  typedef DiamondNormSDPASolver<typename DMTypes::RealScalar, Tomographer::Logger::VacuumLogger>
  //  DiamondNormSolverType;
  typedef DiamondNormSolverType_ DiamondNormSolverType;

  typedef typename DiamondNormSolverType::RealScalarType ValueType;

  inline DiamondNormToRefValueCalculator(const DMTypes dmt, MatrixTypeConstRef E_ref_, int dimX,
                                         const typename DMTypes::RealScalar dnorm_epsilon)
    : E_ref(E_ref_),
      dnslv(dimX, dmt.dim()/dimX, dnorm_epsilon, Tomographer::Logger::vacuum_logger)
  {
    // make sure that dmt.dim() is indeed divisible by dimX
    tomographer_assert((int)dimX*dnslv.dimY() == (int)dmt.dim());
  }

  inline ValueType getValue(const MatrixTypeConstRef & T)
  {
    const int dimX = dnslv.dimX();
    const int dimY = dnslv.dimY();

    //const int dimXY = dnslv.dimXY();
    //typedef Eigen::Matrix<typename DMTypes::ComplexScalar, Eigen::Dynamic, Eigen::Dynamic> MatrixDynType;

    // TODO: act on & use only triangular part for taking the partial trace

    // calculate rho_{XY}
    const MatrixType rho_XY = T*T.adjoint();

    const MatrixType E = Choi_from_process_matrix(rho_XY, dimX, dimY);

    // finally, the difference between the two cpm's.
    MatrixType Delta = E - E_ref;

    // TODO: E might be trace-decreasing if rho_X doesn't have full support. However this
    // is generically not a problem because those states have zero measure --- ???

    return dnslv.calculate(Delta);
  }

private:
  const MatrixType E_ref;
  DiamondNormSolverType dnslv;
};




#endif
