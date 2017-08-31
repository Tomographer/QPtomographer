
#include <cmath>

#include <string>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <cstdlib>
#include <complex>
#include <type_traits>

#include <iostream>
#include <iomanip>

#define EIGEN_INITIALIZE_MATRICES_BY_NAN

// defines the exception class, but doesn't override eigen's eigen_assert() macro itself
#include <tomographer/tools/eigen_assert_exception.h>

#include <Eigen/Core>

#define BOOST_TEST_MODULE test_utils
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/output_test_stream.hpp>

#include "utils.h"

#include <tomographer/tools/boost_test_logger.h>
#include <tomographer/mathtools/pos_semidef_util.h>
#include <tomographer/densedm/dmtypes.h>


// -----------------------------------------------------------------------------
// fixture(s)


// -----------------------------------------------------------------------------
// test suites

BOOST_AUTO_TEST_SUITE(test_utils)
// =============================================================================


BOOST_AUTO_TEST_CASE(partial_trace)
{
  Eigen::MatrixXcd rho_AB(4,4);
  rho_AB <<
    0.8,  0, 0, 0.05,
    0,  0.1, 0, 0,
    0,    0, 0, 0,
    0.05, 0, 0, 0.1;
  Eigen::MatrixXcd rho_A = partial_trace_B(rho_AB, 2);

  BOOST_TEST_MESSAGE("rho_A = " << rho_A) ;

  Eigen::MatrixXcd true_rho_A(2,2);
  true_rho_A <<
    0.9, 0,
    0, 0.1;
  
  BOOST_CHECK((rho_A - true_rho_A).norm() < 1e-6) ;
}

BOOST_AUTO_TEST_CASE(partial_trace_4_2)
{
  Eigen::MatrixXcd E_AB(8,8);
  E_AB <<
    0.01, 0, 0, 0, 0,   0,   0, 0,
    0, 0.99, 0, 0, 0,   0,   0, 0,
    0,  0, 0.6, 0, 0,   0,   0, 0,
    0,  0, 0, 0.4, 0,   0,   0, 0,
    0,  0,   0, 0, 0.9, 0.3, 0, 0,
    0,  0,   0, 0, 0.3, 0.1, 0, 0,
    0,  0,   0, 0, 0,   0, 0.1, 0,
    0,  0,   0, 0, 0,   0,   0, 0.9 ;

  Eigen::MatrixXcd E_A = partial_trace_B(E_AB, 4);

  BOOST_TEST_MESSAGE("E_A = " << E_A) ;

  Eigen::MatrixXcd true_E_A = Eigen::Matrix4d::Identity();
  
  BOOST_CHECK((E_A - true_E_A).norm() < 1e-6) ;
}


BOOST_AUTO_TEST_CASE(t_Choi_from_process_matrix)
{
  Eigen::MatrixXcd rho(4,4);
  rho <<
    5.69952010e-01,   5.22989693e-03,   5.22989693e-03,   4.65354072e-01,
    5.22989693e-03,   3.00479897e-02,   4.79896928e-05,   4.77010307e-03,
    5.22989693e-03,   4.79896928e-05,   4.79896928e-05,   4.27010307e-03,
    4.65354072e-01,   4.77010307e-03,   4.27010307e-03,   3.99952010e-01;
  Eigen::MatrixXcd E = Choi_from_process_matrix(rho, 2, 2);

  BOOST_TEST_MESSAGE("E = \n" << E);

  Eigen::MatrixXcd true_E(4,4);
  true_E <<
    0.95,  0.,    0.,    0.95,
    0.,    0.05,  0.,    0.,
    0.,    0.,    0.,    0.,
    0.95,  0.,    0.,    1. ;

  BOOST_CHECK( (E - true_E).norm() < 1e-6 );

  BOOST_TEST_MESSAGE("DEBUG: rho_X = " << partial_trace_B(rho, 2));
}



// =============================================================================
BOOST_AUTO_TEST_SUITE_END()
