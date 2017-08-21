
#ifndef RELIABLEDIAMONDNORM_TESTUTILS_H
#define RELIABLEDIAMONDNORM_TESTUTILS_H

#include <Eigen/Eigen>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/output_test_stream.hpp>


template<typename RealScalarType = double>
inline
void validate_channel(const Eigen::Ref<
                          const Eigen::Matrix<std::complex<RealScalarType>, Eigen::Dynamic, Eigen::Dynamic>
                      > & E,
                      const std::size_t dim_x, const std::size_t dim_y,
                      const std::string debug_name)
{
  // just to make sure that E is a valid channel:
  BOOST_TEST_MESSAGE("Channel " << debug_name << " is = \n" << E << "\n");

  BOOST_CHECK(E.rows() == E.cols());
  BOOST_CHECK(dim_x * dim_y == (std::size_t)E.rows());

  typedef Eigen::Matrix<std::complex<RealScalarType>, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
  
  MatrixType E_ptA = MatrixType::Zero(dim_x,dim_x);
  for (std::size_t i = 0; i < dim_x; ++i) {
    for (std::size_t j = 0; j < dim_x; ++j) {
      for (std::size_t k = 0; k < dim_y; ++k) {
        E_ptA(i,j) += E(dim_y*i+k, dim_y*j+k);
      }
    }
  }
  BOOST_TEST_MESSAGE("tr_B("<<debug_name<<") is = \n" << E_ptA << "\n");
  BOOST_TEST_MESSAGE("tr("<<debug_name<<") is = " << E.real().trace() << "\n");
  auto eigvals = Eigen::SelfAdjointEigenSolver<MatrixType>(E).eigenvalues();
  BOOST_TEST_MESSAGE("Eigenvalues of "<<debug_name<<" are  = " << eigvals.transpose() << "\n");

  // make sure that everything is OK
  BOOST_CHECK( abs(E.real().trace() - (double)dim_x) < 1e-8 );
  BOOST_CHECK( (E_ptA - MatrixType::Identity(dim_x,dim_x)).norm() < 1e-8 );
  BOOST_CHECK( (eigvals.array()+1e-12 >= 0).all() );
}





#endif
