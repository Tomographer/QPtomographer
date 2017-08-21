
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

#define BOOST_TEST_MODULE test_diamond_norm
#include <boost/test/unit_test.hpp>

#include "testutils.h"

#include "diamond_norm_figofmerit.h"
#include "diamond_norm_sdpa.h"
#include "diamond_norm_scs.h"

#include <tomographer/tools/boost_test_logger.h>
#include <tomographer/mathtools/pos_semidef_util.h>
#include <tomographer/densedm/dmtypes.h>


// -----------------------------------------------------------------------------
// fixture(s)

template<int QuDim_>
struct identdefs
{
  static constexpr int QuDim = QuDim_;
  static constexpr int QuDim2 = QuDim*QuDim;


  typedef Eigen::Matrix<std::complex<double>,QuDim2,1> PureKetType;
  typedef Eigen::Matrix<std::complex<double>,QuDim2,QuDim2> ChoiMatrixType;
  typedef Eigen::Matrix<std::complex<double>,QuDim,QuDim> ReducedDMType;

  PureKetType Psi;
  ChoiMatrixType  Eident;
  ChoiMatrixType  Eclassident;

  identdefs()
    : Psi(PureKetType::Zero())
  {
    BOOST_TEST_MESSAGE("identdefs(), QuDim="<<QuDim<<", QuDim2="<<QuDim2);

    // initialize Psi
    for (std::size_t j = 0; j < QuDim; ++j) {
      Psi(j+QuDim*j) = 1;
    }

    // initialize Eident
    Eident = Psi * Psi.adjoint();

    validate_channel<double>(Eident, QuDim, QuDim, "Eident");

    // initialize Eclassident (classical identity channel == fully depolarizing)
    Eclassident = ChoiMatrixType::Zero();
    for (std::size_t j = 0; j < QuDim; ++j) {
      Eclassident(QuDim*j+j, QuDim*j+j) = 1;
    }

    validate_channel<double>(Eclassident, QuDim, QuDim, "Eclassident");

    BOOST_TEST_MESSAGE("identdefs() done");
  }
};
template<int QuDim_> constexpr int identdefs<QuDim_>::QuDim;
template<int QuDim_> constexpr int identdefs<QuDim_>::QuDim2;




// -----------------------------------------------------------------------------
// test suites

BOOST_AUTO_TEST_SUITE(test_diamond_norm)
// =============================================================================

BOOST_AUTO_TEST_SUITE(examples)

BOOST_FIXTURE_TEST_CASE(simple_qubit_sdpa, identdefs<2>)
{
  // This is now a test channel which is hopefully close to Eclassident.
  ChoiMatrixType E;
  E <<
    0.9, 0, 0, 0,
    0, 0.1, 0, 0,
    0,   0, 0, 0,
    0,   0, 0, 1;

  const ChoiMatrixType Delta = E - Eclassident;

  Tomographer::Logger::BoostTestLogger lg(Tomographer::Logger::LONGDEBUG);

  DiamondNormSDPASolver<double,Tomographer::Logger::BoostTestLogger,true> d(2, 2, lg);
  double value = d.calculate(Delta);

  BOOST_CHECK_CLOSE( value, 0.1, 1e-1/*percent*/ );
}
BOOST_FIXTURE_TEST_CASE(simple_qubit_scs, identdefs<2>)
{
  // This is now a test channel which is hopefully close to Eclassident.
  ChoiMatrixType E;
  E <<
    0.9, 0, 0, 0,
    0, 0.1, 0, 0,
    0,   0, 0, 0,
    0,   0, 0, 1;

  const ChoiMatrixType Delta = E - Eclassident;

  Tomographer::Logger::BoostTestLogger lg(Tomographer::Logger::LONGDEBUG);

  DiamondNormSCSSolver<double,Tomographer::Logger::BoostTestLogger,true> d(2, 2, lg);
  double value = d.calculate(Delta);

  BOOST_CHECK_CLOSE( value, 0.1, 1e-1/*percent*/ );
}

BOOST_FIXTURE_TEST_CASE(qutrit_case_sdpa, identdefs<3>)
{
  // This is now a test channel which is hopefully close to Eident.

  ChoiMatrixType E;
  E  << 
    0.9,0,  0,  0,  0.5,0,  0,  0,0.2,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,0.1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0.5,0,  0,  0,  1,  0,  0,  0,0.8,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
  0.2,  0,  0,  0,0.8,  0,  0,  0,  1 ;
  
  ChoiMatrixType Delta = E - Eident;

  BOOST_TEST_MESSAGE("\n\nAll good, continuing.\n"
                     << "Delta is = \n"
                     << Delta << "\n\n");

  Tomographer::Logger::BoostTestLogger lg(Tomographer::Logger::LONGDEBUG);

  DiamondNormSDPASolver<double,Tomographer::Logger::BoostTestLogger,true> d(3, 3, lg);
  double value = d.calculate(Delta);

  // calculated by solving the SDP with MATLAB/CVX
  BOOST_CHECK_CLOSE( value, 0.42666667, 1e-4/*percent*/ );
}
BOOST_FIXTURE_TEST_CASE(qutrit_case_scs, identdefs<3>)
{
  // This is now a test channel which is hopefully close to Eident.

  ChoiMatrixType E;
  E  << 
    0.9,0,  0,  0,  0.5,0,  0,  0,0.2,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,0.1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0.5,0,  0,  0,  1,  0,  0,  0,0.8,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
  0.2,  0,  0,  0,0.8,  0,  0,  0,  1 ;
  
  ChoiMatrixType Delta = E - Eident;

  BOOST_TEST_MESSAGE("\n\nAll good, continuing.\n"
                     << "Delta is = \n"
                     << Delta << "\n\n");

  Tomographer::Logger::BoostTestLogger lg(Tomographer::Logger::LONGDEBUG);

  DiamondNormSCSSolver<double,Tomographer::Logger::BoostTestLogger,true> d(3, 3, lg);
  double value = d.calculate(Delta);

  // calculated by solving the SDP with MATLAB/CVX
  BOOST_CHECK_CLOSE( value, 0.42666667, 1e-1/*percent*/ );
}

BOOST_AUTO_TEST_SUITE_END(); // examples

BOOST_FIXTURE_TEST_SUITE(diamond_norm_from_bistate, identdefs<2>)

BOOST_AUTO_TEST_CASE(dsolver_sdpa)
{
  Eigen::MatrixXcd true_E(4,4);
  true_E <<
    0.95,  0.,    0.,    0.95,
    0.,    0.05,  0.,    0.,
    0.,    0.,    0.,    0.,
    0.95,  0.,    0.,    1. ;

  DiamondNormSDPASolver<double, Tomographer::Logger::VacuumLogger> dnslv(2, 2, Tomographer::Logger::vacuum_logger);

  double val = dnslv.calculate(true_E - Eident);

  BOOST_TEST_MESSAGE("diamond-norm(true_E) = " << val);
  // value calculated using Python/QuTip
  BOOST_CHECK_CLOSE(val, 0.05, 1e-1/*percent*/) ;
}
BOOST_AUTO_TEST_CASE(dsolver_scs)
{
  Eigen::MatrixXcd true_E(4,4);
  true_E <<
    0.95,  0.,    0.,    0.95,
    0.,    0.05,  0.,    0.,
    0.,    0.,    0.,    0.,
    0.95,  0.,    0.,    1. ;

  DiamondNormSCSSolver<double, Tomographer::Logger::VacuumLogger> dnslv(2, 2, Tomographer::Logger::vacuum_logger);

  double val = dnslv.calculate(true_E - Eident);

  BOOST_TEST_MESSAGE("diamond-norm(true_E) = " << val);
  // value calculated using Python/QuTip
  BOOST_CHECK_CLOSE(val, 0.05, 1e-1/*percent*/) ;
}

BOOST_AUTO_TEST_CASE(example)
{
  Eigen::MatrixXcd rho(4,4);
  rho <<
    5.69952010e-01,   5.22989693e-03,   5.22989693e-03,   4.65354072e-01,
    5.22989693e-03,   3.00479897e-02,   4.79896928e-05,   4.77010307e-03,
    5.22989693e-03,   4.79896928e-05,   4.79896928e-05,   4.27010307e-03,
    4.65354072e-01,   4.77010307e-03,   4.27010307e-03,   3.99952010e-01;

  BOOST_TEST_MESSAGE("rho = \n" << rho );

  Eigen::MatrixXcd T = Tomographer::MathTools::safeOperatorSqrt<Eigen::MatrixXcd>(rho);
  BOOST_TEST_MESSAGE("T = \n" << T);

  Tomographer::DenseDM::DMTypes<4> dmt(4);

  {
    DiamondNormToRefValueCalculator<Tomographer::DenseDM::DMTypes<4>,
                                    DiamondNormSDPASolver<double, Tomographer::Logger::VacuumLogger> > valcalc(dmt, Eident, 2);

    double val = valcalc.getValue(T);
    // value calculated using Python/QuTip
    BOOST_CHECK_CLOSE(val, 0.05, 1e-1/*percent*/) ;
  }
  {
    DiamondNormToRefValueCalculator<Tomographer::DenseDM::DMTypes<4>,
                                    DiamondNormSCSSolver<double, Tomographer::Logger::VacuumLogger> > valcalc(dmt, Eident, 2);

    double val = valcalc.getValue(T);
    // value calculated using Python/QuTip
    BOOST_CHECK_CLOSE(val, 0.05, 1e-1/*percent*/) ;
  }
}

BOOST_AUTO_TEST_SUITE_END();


// =============================================================================
BOOST_AUTO_TEST_SUITE_END() ;

