
#ifndef DIAMOND_NORM_FIGOFMERIT_H
#define DIAMOND_NORM_FIGOFMERIT_H

#include <Eigen/Eigen>

#include <tomographer/tools/loggers.h>
#include <tomographer/tools/needownoperatornew.h>

#include "utils.h"


//
// ValueCalculator for "naive" bipartite state sampling version
//
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

    tomographer_assert( (Delta - Delta.adjoint()).norm()
                        < dimX*dimY*dimX*dimY*Eigen::NumTraits<typename DMTypes::RealScalar>::dummy_precision() );

    // TODO: E might be trace-decreasing if rho_X doesn't have full support. However this
    // is generically not a problem because those states have zero measure --- ???

    return dnslv.calculate(Delta);
  }

private:
  const MatrixType E_ref;
  DiamondNormSolverType dnslv;
};




//
// ValueCalculator for channel-space sampling version
//
template<typename ChannelTypes_, typename DiamondNormSolverType_>
class DiamondNormToRefChannelSpaceValueCalculator
  : public virtual Tomographer::Tools::NeedOwnOperatorNew<typename ChannelTypes_::MatrixType>::ProviderType
{
public:
  typedef ChannelTypes_ ChannelTypes;
  typedef typename ChannelTypes::MatrixType MatrixType;
  typedef typename ChannelTypes::MatrixTypeConstRef MatrixTypeConstRef;

  typedef DiamondNormSolverType_ DiamondNormSolverType;

  typedef typename DiamondNormSolverType::RealScalarType ValueType;

  inline DiamondNormToRefChannelSpaceValueCalculator(
      const ChannelTypes dmt, MatrixTypeConstRef E_ref_, int dimX,
      const double dnorm_epsilon,
      typename DiamondNormSolverType::BaseLoggerType & logger = Tomographer::Logger::vacuum_logger
      )
    : E_ref(E_ref_),
      dnslv(dimX, dmt.dim()/dimX, dnorm_epsilon, logger)
  {
    // make sure that dmt.dim() is indeed divisible by dimX
    tomographer_assert((int)dimX*dnslv.dimY() == (int)dmt.dim());
    // and that E_ref is Hermitian
    tomographer_assert( (E_ref - E_ref.adjoint()).norm()
                        < (dmt.dim()*dmt.dim()*
                           Eigen::NumTraits<typename ChannelTypes::RealScalar>::dummy_precision()) );
  }

  template<typename VIsometryType>
  inline ValueType getValue(const VIsometryType & Vpt)
  {
    tomographer_assert(Vpt.cols() == dnslv.dimX() && Vpt.rows() == dnslv.dimX()*dnslv.dimY()*dnslv.dimY());

    //const int dimX = dnslv.dimX();
    //const int dimY = dnslv.dimY();
    const int dimXY = dnslv.dimXY();

    MatrixType T(dimXY,dimXY);
    T = remapIsometryToT<MatrixType>(Vpt, dimXY);

    const auto E = T * T.adjoint();

    // finally, the difference between the two cpm's.
    MatrixType Delta(E_ref - E);

    // make sure Delta is hermitian(!)
    tomographer_assert( (Delta - Delta.adjoint()).norm()
                        < dimXY*dimXY*Eigen::NumTraits<typename ChannelTypes::RealScalar>::dummy_precision() );

    return dnslv.calculate(Delta);
  }

private:
  const MatrixType E_ref;
  DiamondNormSolverType dnslv;
};






#endif
