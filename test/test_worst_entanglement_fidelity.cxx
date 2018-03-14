
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
#include <unsupported/Eigen/MatrixFunctions> // matrix sqrt()

#define BOOST_TEST_MODULE test_worst_entanglement_fidelity
#include <boost/test/unit_test.hpp>

#include "testutils.h"

#include "worstentglfidelity_figofmerit.h"
#include "channelspace.h"

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

  /*
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
  */
};
template<int QuDim_> constexpr int identdefs<QuDim_>::QuDim;
template<int QuDim_> constexpr int identdefs<QuDim_>::QuDim2;



// -----------------------------------------------------------------------------
// test suites

BOOST_AUTO_TEST_SUITE(test_worst_entanglement_fidelity)
// =============================================================================

BOOST_AUTO_TEST_SUITE(examples)

BOOST_FIXTURE_TEST_CASE(simple_qubit, identdefs<2>)
{
  // This is now a test channel which is hopefully close to Eclassident.
  ChoiMatrixType E;
  E <<
    0.9, 0, 0, 0,
    0, 0.1, 0, 0,
    0,   0, 0, 0,
    0,   0, 0, 1;

  Tomographer::Logger::BoostTestLogger logger(Tomographer::Logger::LONGDEBUG) ;
  typedef WorstEntglFidelitySCSSolver<SCS::scs_float, Tomographer::Logger::BoostTestLogger> Solver;
  Solver fidslv(2, Solver::factorizeChoiMatrix(E), 1e-8, logger);

  double fid = fidslv.calculate();

  // value calculated manually (Drafts&Calculations Vol VIII 3/13/2018 & MATLAB)
  BOOST_CHECK_CLOSE( fid, 1.-(1./1.9) /*0.4736842105263158*/, 1e-4/*percent*/ ) ;
}

BOOST_FIXTURE_TEST_CASE(qutrit_case, identdefs<3>)
{
  // This is now a test of another channel

  const auto J = std::complex<double>(0,1);

  ChoiMatrixType E;
  E  << 
    0.9,0,  0,  0,-0.5*J,0, 0,  0,0.2,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,0.1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
  0.5*J,0,  0,  0,  1,  0,  0,  0,0.8,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
  0.2,  0,  0,  0,0.8,  0,  0,  0,  1 ;
  
  Tomographer::Logger::BoostTestLogger logger(Tomographer::Logger::LONGDEBUG) ;
  typedef WorstEntglFidelitySCSSolver<SCS::scs_float, Tomographer::Logger::BoostTestLogger> Solver;
  Solver fidslv(3, Solver::factorizeChoiMatrix(E), 1e-7, logger);

  double fid = fidslv.calculate();

  // value calculated using MATLAB
  BOOST_CHECK_CLOSE( fid, 4.736842105263159e-01, 1e-4/*percent*/ );
}

BOOST_FIXTURE_TEST_CASE(qutrit_case_close, identdefs<3>)
{
  // This is now a test channel which is hopefully close to Eident.

  ChoiMatrixType E;
  E  << 
    0.9,0,  0,  0,  0.9,0,  0,  0,0.9,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,0.1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0.9,0,  0,  0,  0.9,0,  0,  0,0.9,
    0,  0,  0,  0,  0,0.1,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0.9,0,  0,  0,0.9,  0,  0,  0,  1 ;
  
  Tomographer::Logger::BoostTestLogger logger(Tomographer::Logger::LONGDEBUG) ;
  typedef WorstEntglFidelitySCSSolver<SCS::scs_float, Tomographer::Logger::BoostTestLogger> Solver;
  Solver fidslv(3, Solver::factorizeChoiMatrix(E), 1e-7, logger);

  double fid = fidslv.calculate();

  // value calculated using MATLAB
  BOOST_CHECK_CLOSE( fid, 8.999999999999997e-01, 1e-4/*percent*/ );
}


BOOST_FIXTURE_TEST_CASE(qutrit_case_ideal, identdefs<3>)
{
  // This is now an ideal test channel

  ChoiMatrixType E;
  E  << 
    1,  0,  0,  0,  1,  0,  0,  0,  1,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    1,  0,  0,  0,  1,  0,  0,  0,  1,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    1,  0,  0,  0,  1,  0,  0,  0,  1 ;
  
  Tomographer::Logger::BoostTestLogger logger(Tomographer::Logger::LONGDEBUG) ;
  typedef WorstEntglFidelitySCSSolver<SCS::scs_float, Tomographer::Logger::BoostTestLogger> Solver;
  Solver fidslv(3, Solver::factorizeChoiMatrix(E), 1e-7, logger);

  double fid = fidslv.calculate();

  // value calculated manually
  BOOST_CHECK_CLOSE( fid, 1.0, 1e-4/*percent*/ );
}

BOOST_AUTO_TEST_SUITE_END(); // examples

BOOST_FIXTURE_TEST_SUITE(worst_entanglement_fidelity_valcalc, identdefs<2>)

BOOST_FIXTURE_TEST_CASE(simple_qubit, identdefs<2>)
{
  // process matrix of the channel diag(0.9,0.1,0,1) applied to the state (0.8,0.2)
  ChoiMatrixType E;
  E <<
    0.72, 0, 0,   0,
    0, 0.08, 0,   0,
    0,    0, 0,   0,
    0,    0, 0, 0.2;

  ChoiMatrixType T = E.sqrt(); // a square root of E, to test the value calculator's getValue()

  typedef Tomographer::DenseDM::DMTypes<4> DMTypes;
  DMTypes dmt(4);

  WorstEntglFidelityValueCalculator<DMTypes> wefidcalc(dmt, 2, 1e-8);

  // value calculated manually
  BOOST_CHECK_CLOSE( wefidcalc.getValue(T), 1-1./1.9, 1e-4/*percent*/ );
}


BOOST_FIXTURE_TEST_CASE(simple_qubit_channelspace, identdefs<2>)
{
  // Choi matrix = diag(0.9,0.1,0,1)
  //
  // Here we need a Stinespring isometry:
  //
  //   |0>  -->  sqrt(0.9) * |0>_out |0>_ref + sqrt(0.1) * |1>_out |1>_ref
  //   |1>  -->  |1>_out |2>_ref
  //

  Eigen::Matrix<std::complex<double>,8,2> Vpt;
  Vpt.col(0)
    << std::sqrt(0.9),0,0,0, 0,std::sqrt(0.1),0,0 ;
  Vpt.col(1)
    << 0,0,0,0, 0,0,1.0,0 ;

  typedef ChannelTypes<double> MyChannelTypes;
  MyChannelTypes cht(2,2);

  WorstEntglFidelityChannelSpaceValueCalculator<MyChannelTypes> wefidcalc(cht, 1e-8);

  // value calculated manually
  BOOST_CHECK_CLOSE( wefidcalc.getValue(Vpt), 1.-1./1.9, 1e-4/*percent*/ );
}

BOOST_AUTO_TEST_SUITE_END() ;



// =============================================================================
BOOST_AUTO_TEST_SUITE_END() ;

