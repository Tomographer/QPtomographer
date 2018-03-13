
#ifndef ENTANGLEMENT_FIDELITY_FIGOFMERIT_H
#define ENTANGLEMENT_FIDELITY_FIGOFMERIT_H

#include <Eigen/Eigen>

#include <tomographer/tools/loggers.h>
#include <tomographer/tools/needownoperatornew.h>

#include "utils.h"


template<typename DMTypes_>
class EntglFidelityValueCalculator
{
public:
  typedef DMTypes_ DMTypes;
  typedef typename DMTypes::MatrixType MatrixType;
  typedef typename DMTypes::MatrixTypeConstRef MatrixTypeConstRef;

  typedef typename DMTypes::RealScalar ValueType;

  inline EntglFidelityValueCalculator(const DMTypes dmt, int dimX_)
    : dimX(dimX_)
  {
  }

  inline ValueType getValue(const MatrixTypeConstRef & T)
  {
    // TODO: act on & use only triangular part for taking the partial trace

    tomographer_assert(T.rows() == T.cols());
    const int dimXY = T.rows();

    // only makes sense to calculate entanglement fidelity for dimX == dimY
    tomographer_assert(dimX*dimX == dimXY) ;
    const int dimY = dimX;

    // calculate rho_{XY}
    const MatrixType rho_XY = T * T.adjoint();

    const MatrixType E = Choi_from_process_matrix(rho_XY, dimX, dimY);

    // calculate the overlap of the Choi state with the maximally entangled state.

    return entanglement_fidelity(E, dimX);
  }

private:
  const int dimX;
};







template<typename ChannelTypes_>
class EntglFidelityChannelSpaceValueCalculator
{
public:
  typedef ChannelTypes_ ChannelTypes;
  typedef typename ChannelTypes::MatrixType MatrixType;
  typedef typename ChannelTypes::MatrixTypeConstRef MatrixTypeConstRef;

  typedef typename ChannelTypes::RealScalar ValueType;

  inline EntglFidelityChannelSpaceValueCalculator(const ChannelTypes cht_)
    : cht(cht_)
  {
  }

  template<typename VIsometryType>
  inline ValueType getValue(const VIsometryType & Vpt) const
  {
    // TODO: act on & use only triangular part for taking the partial trace

    tomographer_assert(Vpt.cols() == cht.dimX());
    tomographer_assert(Vpt.rows() == cht.dimXY2());

    // only makes sense to calculate entanglement fidelity for dimX == dimY
    tomographer_assert(cht.dimX()*cht.dimX() == cht.dim()) ;

    MatrixType T(cht.dim(), cht.dim());
    T = remapIsometryToT<MatrixType>(Vpt, cht.dim());

    const MatrixType E( T * T.adjoint() );

    // calculate the overlap of the Choi state with the maximally entangled state.

    return entanglement_fidelity(E, cht.dimX());
  }

private:
  const ChannelTypes cht;
};




#endif
