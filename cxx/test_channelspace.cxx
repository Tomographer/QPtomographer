
#include <cmath>

#include <string>
#include <iostream>
#include <random>
#include <algorithm>

#include <boost/math/constants/constants.hpp>

#define EIGEN_INITIALIZE_MATRICES_BY_NAN

#include <cmath>
#include <cstdlib>
#include <complex>
#include <type_traits>

#include <iostream>
#include <iomanip>

// define the exception class, but don't override eigen's eigen_assert() macro itself
#include <tomographer/tools/eigen_assert_exception.h>

#include <Eigen/Core>
#include <unsupported/Eigen/KroneckerProduct>

#define BOOST_TEST_MODULE test_channelspace
#include <boost/test/unit_test.hpp>

#include <tomographer/tools/boost_test_logger.h>
#include <tomographer/mathtools/pos_semidef_util.h>
#include <tomographer/densedm/dmtypes.h>

#include "channelspace.h"


// -----------------------------------------------------------------------------
// fixture(s)


struct MyStatsCollector {
  inline void init() { }
  inline void thermalizingDone() { }
  inline void done() { }
  template<typename PointType, typename MHRandomWalk>
  inline void processSample(int k, int n, const PointType & Vpt, double fnval, MHRandomWalk & rw)
  {
    BOOST_TEST_MESSAGE("New sample: V = \n" << Vpt) ;
  }
  template<typename... Args>
  inline void rawMove(Args&&...) const { }
};


template<int JumpMode_, int LogLevel = Tomographer::Logger::LONGDEBUG>
struct qubit_to_qutrit_fixture {
  
  typedef ChannelTypes<> ChannelTypes;
  typedef ChannelTypes::MatrixType MatrixType;
  typedef Tomographer::DenseDM::IndepMeasLLH<ChannelTypes> DenseLLH;

  ChannelTypes cht;
  DenseLLH llh;

  std::mt19937 rng;
  Tomographer::Logger::BoostTestLogger logger;

  // the channel MHWalker
  typedef ChannelSpaceMHWalker<DenseLLH, std::mt19937, Tomographer::Logger::BoostTestLogger>
    ChannelSpaceMHWalker;

  ChannelSpaceMHWalker mhwalker;

  static inline DenseLLH _make_llh(ChannelTypes cht)
  {
    DenseLLH llh(cht);

    MatrixType Zp(2,2);
    MatrixType Zm(2,2);
    MatrixType Xp(2,2);
    MatrixType Xm(2,2);

    MatrixType out(3,3);

    Zp << 1, 0,  0, 0;
    Zm << 0, 0,  0, 1;
    Xp << 1, 1,  1, 1;
    Xm << 1,-1, -1, 1;

    out << 1, 0, 0,  0, 0, 0,  0, 0, 0; // pure |0> state on qutrit

    llh.addMeasEffect(Eigen::kroneckerProduct(Zp, out), 50);
    llh.addMeasEffect(Eigen::kroneckerProduct(Zm, out), 20);
    llh.addMeasEffect(Eigen::kroneckerProduct(Xp, out), 10);
    llh.addMeasEffect(Eigen::kroneckerProduct(Xm, out), 90);

    return llh;
  }

  qubit_to_qutrit_fixture()
    : cht(2,3),
      llh(_make_llh(cht)),
      rng(0),
      logger(LogLevel),
      mhwalker(llh, rng, logger, (JumpMode)JumpMode_)
  {
  }

};

// because BOOST macros don't like commas in their fixture argument...
template<int JumpType_>
struct qubit_to_qutrit_fixture_DEBUG : public qubit_to_qutrit_fixture<JumpType_, Tomographer::Logger::DEBUG> { };


// -----------------------------------------------------------------------------
// test suites

BOOST_AUTO_TEST_SUITE(test_channelspacemhwalker);
// =============================================================================


BOOST_FIXTURE_TEST_CASE(simple_qubit_to_qutrit_eiH, qubit_to_qutrit_fixture<RandHermExp>)
{
  MyStatsCollector stats;
  typedef Tomographer::MHRandomWalk<std::mt19937, ChannelSpaceMHWalker, MyStatsCollector,
                                    Tomographer::Logger::BoostTestLogger> MHRandomWalk;
  MHRandomWalk rw(0.1, 10, 500, 1024, mhwalker, stats, rng, logger);
  rw.run();

  BOOST_CHECK( true ) ; // dummy test case -- for now, just check that the above compiles & runs
}
BOOST_FIXTURE_TEST_CASE(simple_qubit_to_qutrit_elr, qubit_to_qutrit_fixture<ElemRotations>)
{
  MyStatsCollector stats;
  typedef Tomographer::MHRandomWalk<std::mt19937, ChannelSpaceMHWalker, MyStatsCollector,
                                    Tomographer::Logger::BoostTestLogger> MHRandomWalk;
  MHRandomWalk rw(0.1, 10, 500, 1024, mhwalker, stats, rng, logger);
  rw.run();

  BOOST_CHECK( true ) ; // dummy test case -- for now, just check that the above compiles & runs
}

BOOST_FIXTURE_TEST_CASE(jump_fn_symmetric_eiH, qubit_to_qutrit_fixture_DEBUG<RandHermExp>)
{
  // start at a fixed starting point
  typedef ChannelSpaceMHWalker::PointType PointType;
  PointType V0 = PointType::Zero(llh.dmt.dimXY2(), llh.dmt.dimX());
  V0(0,0) = V0(1,1) = 1;

  // average over many steps starting from the same V0 and normalize --> should be V0
  // again

  const int N = 5000000;
  const double step_size = 0.1;

  PointType totalV = PointType::Zero(llh.dmt.dimXY2(), llh.dmt.dimX());

  for (int k = 0; k < N; ++k) {
    PointType newV = mhwalker.jumpFn(V0, step_size) ;
    totalV += newV;
  }
  BOOST_TEST_MESSAGE("total summed newV's = \n" << totalV) ;

  totalV /= (totalV.norm() / std::sqrt((double)llh.dmt.dimX()));

  BOOST_TEST_MESSAGE("normalized = \n" << totalV) ;

  double err = (totalV - V0).cwiseAbs().maxCoeff();
  double errok = 4*step_size/std::sqrt(N); // approx. 4 sigma
  BOOST_TEST_MESSAGE("err = " << err << ", errok = " << errok) ;
  BOOST_CHECK_LT(err, errok) ;
}
BOOST_FIXTURE_TEST_CASE(jump_fn_symmetric_elr, qubit_to_qutrit_fixture_DEBUG<ElemRotations>)
{
  // start at a fixed starting point
  typedef ChannelSpaceMHWalker::PointType PointType;
  PointType V0 = PointType::Zero(llh.dmt.dimXY2(), llh.dmt.dimX());
  V0(0,0) = V0(1,1) = 1;

  // average over many steps starting from the same V0 and normalize --> should be V0
  // again

  const int N = 5000000;
  const double step_size = 0.1;

  PointType totalV = PointType::Zero(llh.dmt.dimXY2(), llh.dmt.dimX());

  for (int k = 0; k < N; ++k) {
    PointType newV = mhwalker.jumpFn(V0, step_size) ;
    totalV += newV;
  }
  BOOST_TEST_MESSAGE("total summed newV's = \n" << totalV) ;

  totalV /= (totalV.norm() / std::sqrt((double)llh.dmt.dimX()));

  BOOST_TEST_MESSAGE("normalized = \n" << totalV) ;

  double err = (totalV - V0).cwiseAbs().maxCoeff();
  double errok = 4*step_size/std::sqrt(N); // approx. 4 sigma
  BOOST_TEST_MESSAGE("err = " << err << ", errok = " << errok) ;
  BOOST_CHECK_LT(err, errok) ;
}


// =============================================================================
BOOST_AUTO_TEST_SUITE_END() ;

