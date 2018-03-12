
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

#define BOOST_TEST_MODULE test_entanglement_fidelity
#include <boost/test/unit_test.hpp>

#include "testutils.h"

#include "entglfidelity_figofmerit.h"
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

BOOST_AUTO_TEST_SUITE(test_entanglement_fidelity)
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

  // value calculated manually
  BOOST_CHECK_CLOSE( entanglement_fidelity(E, 2), 0.475, 1e-4/*percent*/ );
}

BOOST_FIXTURE_TEST_CASE(qutrit_case, identdefs<3>)
{
  // This is now a test of another channel

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

  // value calculated manually
  BOOST_CHECK_CLOSE( entanglement_fidelity(E, 3), 5.9/9, 1e-4/*percent*/ );
}

BOOST_FIXTURE_TEST_CASE(qutrit_case_close, identdefs<3>)
{
  // This is now a test channel which is hopefully close to Eident.

  ChoiMatrixType E;
  E  << 
    0.9,0,  0,  0,  0.9,0,  0,  0,0.95,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,0.1,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0.9,0,  0,  0,  0.9,0,  0,  0,0.95,
    0,  0,  0,  0,  0,0.1,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,
    1,  0,  0,  0,0.9,  0,  0,  0,  1 ;
  
  // value calculated manually
  BOOST_CHECK_CLOSE( entanglement_fidelity(E, 3), 8.4/9, 1e-4/*percent*/ );
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
  
  // value calculated manually
  BOOST_CHECK_CLOSE( entanglement_fidelity(E, 3), 1.0, 1e-4/*percent*/ );
}

BOOST_AUTO_TEST_SUITE_END(); // examples

BOOST_FIXTURE_TEST_SUITE(entanglement_fidelity_valcalc, identdefs<2>)

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

  EntglFidelityValueCalculator<DMTypes> efidcalc(dmt,2);

  // value calculated manually
  BOOST_CHECK_CLOSE( efidcalc.getValue(T), 0.475, 1e-4/*percent*/ );
}


BOOST_FIXTURE_TEST_CASE(simple_qubit_channelspace, identdefs<2>)
{
  // Choi matrix = diag(0.9,0.1,0,1)
  //
  // Here we need a Stinespring isometry:
  //
  //   |0>  -->  sqrt(0.9) * |0>_out |0>_ref + sqrt(0.1) * |1>_out |0>_ref
  //   |1>  -->  |1>_out |1>_ref
  //
  // ref actually has four dimensions (that's what's needed in general) but we
  // don't need all of them for this particular channel.

  Eigen::Matrix<std::complex<double>,8,2> Vpt;
  Vpt.col(0)
    << std::sqrt(0.9),0,0,0, std::sqrt(0.1),0,0,0 ;
  Vpt.col(1)
    << 0,0,0,0, 0,1.0,0,0 ;

  typedef ChannelTypes<double> MyChannelTypes;
  MyChannelTypes cht(2,2);

  EntglFidelityChannelSpaceValueCalculator<MyChannelTypes> efidcalc(cht);

  // value calculated manually
  BOOST_CHECK_CLOSE( efidcalc.getValue(Vpt), 0.475, 1e-4/*percent*/ );
}

BOOST_AUTO_TEST_SUITE_END() ;


// =============================================================================
BOOST_AUTO_TEST_SUITE_END() ;

