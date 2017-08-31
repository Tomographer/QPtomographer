
#ifndef CHANNELSPACE_H
#define CHANNELSPACE_H

#include <random>

#include <unsupported/Eigen/MatrixFunctions> // matrix exponential

#include <tomographer/tomographer_version.h>
#include <tomographer/mhrw.h>
#include <tomographer/densedm/param_herm_x.h>
#include <tomographer/densedm/indepmeasllh.h>
#include <tomographer/tools/needownoperatornew.h>
#include <tomographer/tools/eigenutil.h> // denseRandom
//#include <tomographer/mathtools/random_unitary.h>


#include "utils.h"




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
  tomographer_assert( (V.adjoint()*V - MatrixType::Identity(N,N)).norm() < Eigen::NumTraits<Scalar>::dummy_precision() );

  logger.longdebug([&](std::ostream& str) {
      str << "randomIsometry: got V = \n" << V << "\n"
  	  << "Check: V*V.adjoint() ==\n" << V*V.adjoint() << "\n"
  	  << "Check: V.adjoint()*V ==\n" << V.adjoint()*V;
    });

  return V;
}











template<typename Scalar = double>
struct ChannelTypes : public Tomographer::DenseDM::DMTypes<Eigen::Dynamic,Scalar>
{
  typedef Tomographer::DenseDM::DMTypes<Eigen::Dynamic,Scalar> Base;

  using typename Base::RealScalar;
  using typename Base::ComplexScalar;

private:
  std::size_t dim_x;
  std::size_t dim_y;

public:
  ChannelTypes(std::size_t dim_x_, std::size_t dim_y_)
    : Base(dim_x_*dim_y_), dim_x(dim_x_), dim_y(dim_y_)
  {
  }

  inline std::size_t dimX() const { return dim_x; }
  inline std::size_t dimY() const { return dim_y; }
  inline std::size_t dimXY2() const { return Base::dim()*dim_y; }

};




typedef Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> DynStride;

template<typename VIsometryType,
         TOMOGRAPHER_ENABLED_IF_TMPL(!VIsometryType::IsRowMajor)>
inline Eigen::Map<const Eigen::Matrix<typename VIsometryType::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>,
                  0,DynStride>
remapIsometryToT(const VIsometryType & V, int dXY)
{
  typedef Eigen::Matrix<typename VIsometryType::Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixType;
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





enum JumpMode {
  RandHermExp,
  ElemRotations,
};


/** \brief Compiles with the MHWalker type interface (see Tomographer API documentation)
 *
 * \note the \a DMTypes class used by the DenseLLH must be a ChannelTypes type.
 */
template<typename DenseLLH_, typename Rng_, typename BaseLoggerType_>
class ChannelSpaceMHWalker
  : public Tomographer::Tools::NeedOwnOperatorNew<typename DenseLLH_::DMTypes::MatrixType>::ProviderType
{
public:
  typedef DenseLLH_ DenseLLH;
  typedef Rng_ Rng;
  typedef BaseLoggerType_ BaseLoggerType;

  typedef typename DenseLLH::DMTypes ChannelTypes;
#if TOMOGRAPHER_VERSION_MAJ >= 5
  typedef typename Tomographer::MHWalkerParamsStepSize<typename ChannelTypes::RealScalar> WalkerParams; // for Tomographer >=5
#else
  typedef typename ChannelTypes::RealScalar StepRealType; // for Tomographer <5
#endif
  typedef typename DenseLLH::LLHValueType FnValueType;

  static constexpr int UseFnSyntaxType = Tomographer::MHUseFnLogValue;
  
  static_assert(ChannelTypes::MatrixType::RowsAtCompileTime == Eigen::Dynamic &&
                ChannelTypes::MatrixType::RowsAtCompileTime == Eigen::Dynamic,
                "ChannelSpaceMHWalker: the code here requires the matrix type to be dynamic");
  static_assert((int)DenseLLH::LLHCalcType == (int)Tomographer::DenseDM::LLHCalcTypeX,
                "ChannelSpaceMHWalker: The DenseLLH object must have LLHCalcType==LLHCalcTypeX.");

  typedef typename ChannelTypes::MatrixType MatrixType;
  typedef typename ChannelTypes::RealScalar RealScalar;
  typedef typename ChannelTypes::ComplexScalar ComplexScalar;

  //! Type of the V isometry which will be the point type.
  typedef Eigen::Matrix<ComplexScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> VIsometryType;
  typedef VIsometryType PointType;

private:

  const DenseLLH & llh;
  Rng & rng;
  std::normal_distribution<RealScalar> normal_distr_rnd;
  const Tomographer::DenseDM::ParamX<ChannelTypes> px;

  Tomographer::Logger::LocalLogger<BaseLoggerType> locallogger;
  
  JumpMode mode;

public:
  ChannelSpaceMHWalker(const DenseLLH & llh_, Rng & rng_, BaseLoggerType & logger_, JumpMode mode_)
    : llh(llh_),
      rng(rng_),
      normal_distr_rnd(0.0, 1.0),
      px(llh.dmt),
      locallogger(TOMO_ORIGIN, logger_),
      mode(mode_)
  {
  }

  ChannelSpaceMHWalker(ChannelSpaceMHWalker && other) // move constructor
    : llh(std::move(other.llh)),
      rng(other.rng),
      normal_distr_rnd(0.0, 1.0),
      px(llh.dmt),
      locallogger(std::move(other.locallogger))
  {
  }
  
  inline PointType startPoint()
  {
    return randomIsometry<PointType>(llh.dmt.dimXY2(), llh.dmt.dimX(), rng, locallogger.parentLogger());
  }

  inline void init()
  {
    auto logger = locallogger.subLogger(TOMO_ORIGIN);
    logger.debug("Starting random walk");
  }

  inline void thermalizingDone()
  {
    auto logger = locallogger.subLogger(TOMO_ORIGIN);
    logger.debug("thermalizing done");
  }

  inline void done()
  {
    auto logger = locallogger.subLogger(TOMO_ORIGIN);
    logger.debug("Random walk done");
  }

  // for Tomographer 5: now argument is a struct with walker params
  template<typename MHWalkerParams>
  inline PointType jumpFn(const PointType & curV, MHWalkerParams && params)
  {
    return jumpFn(curV, params.step_size) ;
  }

  inline PointType jumpFn(const PointType & curV, RealScalar step_size) // not const -- updates the rng state
  {
    auto logger = locallogger.subLogger(TOMO_ORIGIN);

    PointType newV;

    switch (mode) {
    case RandHermExp:
      {
        logger.longdebug("Choosing random H, jump via e^{iH}") ;
        // Define choose a relatively small H, which will define the W as W=e^{iH}
        MatrixType A(llh.dmt.dimXY2(), llh.dmt.dimXY2());
        A = step_size*Tomographer::Tools::denseRandom<MatrixType>(rng, normal_distr_rnd,
                                                                  llh.dmt.dimXY2(), llh.dmt.dimXY2());
        MatrixType H = A + A.adjoint(); // make it Hermitian!
        newV = (ComplexScalar(0,1)*H).exp() * curV; // e^{iH}*V, w/ matrix exponential 
        break;
      }
    case ElemRotations:
      {
        logger.longdebug("Choosing random elementary Pauli-type rotation") ;

        newV = curV;

        for (int elem_rot_iter = 0; elem_rot_iter < 8*(int)llh.dmt.dimXY2(); ++elem_rot_iter) {
          logger.longdebug("iteration #%d", elem_rot_iter) ;

          // select two random indices, and randomly select whether we apply an elementary x, y or z rotation
          std::uniform_int_distribution<int> rnddist_i(0, llh.dmt.dimXY2()-1);
          std::uniform_int_distribution<int> rnddist_3(0, 2);
          int i = rnddist_i(rng);
          int j; do { j = rnddist_i(rng); } while (j == i); // choose also randomly, but different than i.
          if (i > j) { std::swap(i, j); } // ensure that i < j
          int xyz = rnddist_3(rng);
        
          logger.longdebug([&](std::ostream & stream) {
              stream << "dimXY2=" << llh.dmt.dimXY2() << " i=" << i << " j=" << j << " xyz=" << xyz;
            });

          RealScalar sina = step_size * normal_distr_rnd(rng); // * llh.dmt.dimXY2(); -- not with the several iterations
          if (sina < -1) { sina = -1; }
          if (sina >  1) { sina =  1; }
          RealScalar cosa = std::sqrt(1 - sina*sina);
          // apply transformation
          Eigen::Matrix<ComplexScalar,2,2> tr2d;
          // remember: e^{i\phi(\vec{n}\cdot\vec{\sigma})} = \cos\phi \Ident + i\sin\phi (\vec{n}\cdot\vec{\sigma})
          switch (xyz) {
          case 0: // X-type rotation
            tr2d(0,0) = cosa;                   tr2d(0,1) = ComplexScalar(0,sina);
            tr2d(1,0) = ComplexScalar(0,sina);  tr2d(1,1) = cosa;
            break;
          case 1: // Y-type rotation
            tr2d(0,0) = cosa;   tr2d(0,1) = sina;
            tr2d(1,0) = -sina;  tr2d(1,1) = cosa;
            break;
          case 2: // Z-type rotation
            tr2d(0,0) = ComplexScalar(cosa, sina);  tr2d(0,1) = 0;
            tr2d(1,0) = 0;                          tr2d(1,1) = ComplexScalar(cosa, -sina);
            break;
          default:
            tomographer_assert(false && "Invalid rotation type number sampled!");
          }
          logger.longdebug([&](std::ostream & stream) { stream << "cosa,sina=" << cosa << "," << sina
                                                               << "  tr2d = "  << std::setprecision(10) << tr2d; }) ;
          logger.longdebug([&](std::ostream & stream) { stream << "newV before update =\n" << newV; }) ;
          // do this:
          for (int k = 0; k < (int)llh.dmt.dimX(); ++k) {
            Eigen::Matrix<ComplexScalar,2,1> res;
            res(0) = tr2d(0,0) * newV(i,k) + tr2d(0,1) * newV(j,k) ;
            res(1) = tr2d(1,0) * newV(i,k) + tr2d(1,1) * newV(j,k) ;
            newV(i,k) = res(0);
            newV(j,k) = res(1);
          }
          //
          // getting out of hand complicated for what it's supposed to do:
          //
          // logger.longdebug([&](std::ostream & stream) { stream << "tr2d = " << tr2d; }) ;
          // typedef Eigen::Map<Eigen::Matrix<typename PointType::Scalar,2,Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic,1> >
          //   MappedType;
          // MappedType subm(newV.block(i, 0, 1, 1).data(), 2, llh.dmt.dimX(), Eigen::Stride<Eigen::Dynamic,1>(j-i));
          // logger.longdebug([&](std::ostream & stream) { stream << "subm = " << subm; }) ;
          // subm = tr2d * subm;
          // logger.longdebug([&](std::ostream & stream) { stream << "subm is now = " << subm; }) ;
          //
        }
        logger.longdebug([&](std::ostream & stream) { stream << "newV after =\n" << newV; });
        break;
      }
    default:
      {
        tomographer_assert(false && "Invalid jump function mode!") ;
      }
    };

    logger.longdebug([&](std::ostream & stream) {
        stream << "newV = ("<<newV.rows()<<","<<newV.cols()<<")\n"
               << newV << "\n"
               << "newV.adjoint()*newV = \n" << (newV.adjoint()*newV);
      });

    // make sure that all columns are indeed still orthogonal
    const double errok = Eigen::NumTraits<RealScalar>::dummy_precision() / 4; // "/4" --> get some leeway
    if ( ( newV.adjoint()*newV - MatrixType::Identity(llh.dmt.dimX(), llh.dmt.dimX()) ).norm() >= errok ) {
      MatrixType newVbad(newV.rows(), newV.cols()); newVbad = newV;
      // orthogonalize:
      for (int k = 0; k < (int)llh.dmt.dimX(); ++k) {
        if (k > 0) {
          auto projkdimspace = newV.block(0,0,newV.rows(),k) * newV.block(0,0,newV.rows(),k).adjoint();
          newV.col(k) -= projkdimspace * newV.col(k);
        }
        newV.col(k) /= newV.col(k).norm();
      }
      logger.debug([&](std::ostream & stream) {
          stream << "New V does not have orthogonal columns, orthogonalized: newV =\n"
                 << std::setprecision(10) << newVbad << "\n  to\n" << std::setprecision(10) << newV
                 << "\nNow (V*V' - Identity).norm() = "
                 << ( newV.adjoint()*newV - MatrixType::Identity(llh.dmt.dimX(), llh.dmt.dimX()) ).norm()
                 << "  errok = " << errok;
        });
    }

    return newV;
  }

  inline FnValueType fnLogVal(const PointType & V)
  {
    auto logger = locallogger.subLogger(TOMO_ORIGIN);
    logger.longdebug("fnLogVal()");
    logger.longdebug([&](std::ostream & stream) { stream << "V = \n" << V; });
    auto T = remapIsometryToT<MatrixType>(V, llh.dmt.dim());
    logger.longdebug([&](std::ostream & stream) { stream << "T = \n" << T; });
    auto Rho = T * T.adjoint();
    logger.longdebug([&](std::ostream & stream) { stream << "Rho = \n" << Rho; });
    // make sure we wrote the remapIsometryToT() function correctly--we must have partial_trace_B(Rho) = identity
    // CHECK:
    MatrixType RhoA = partial_trace_B(Rho, llh.dmt.dimX());
    logger.longdebug([&](std::ostream & stream) {
        stream << "Rho_A = \n" << RhoA << " ; should be identity";
      });
    tomographer_assert( (RhoA - MatrixType::Identity(llh.dmt.dimX(),llh.dmt.dimX())).norm()
                        <  Eigen::NumTraits<RealScalar>::dummy_precision() );
    //
    
    typename ChannelTypes::VectorParamType x = px.HermToX(Rho);
    // NOTE: we assume that the currect input state for each repetition is encoded into the list of POVM effects
    return llh.logLikelihoodX(x);
  }

};








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

  inline DiamondNormToRefChannelSpaceValueCalculator(const ChannelTypes dmt, MatrixTypeConstRef E_ref_, int dimX,
                                                     const double dnorm_epsilon,
                                                     typename DiamondNormSolverType::BaseLoggerType & logger
                                                     = Tomographer::Logger::vacuum_logger)
    : E_ref(E_ref_),
      dnslv(dimX, dmt.dim()/dimX, dnorm_epsilon, logger)
  {
    // make sure that dmt.dim() is indeed divisible by dimX
    tomographer_assert((int)dimX*dnslv.dimY() == (int)dmt.dim());
    // and that E_ref is Hermitian
    tomographer_assert( (E_ref - E_ref.adjoint()).norm()
                        < dmt.dim()*dmt.dim()*Eigen::NumTraits<typename ChannelTypes::RealScalar>::dummy_precision() );
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
