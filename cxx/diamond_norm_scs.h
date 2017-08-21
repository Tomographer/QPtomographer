
#ifndef DIAMOND_NORM_SCS_H
#define DIAMOND_NORM_SCS_H

#include <iostream>

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include <boost/math/constants/constants.hpp>

#include <tomographer/tools/cxxutil.h>

//
// NOTES: INSTRUCTIONS FOR COMPILING SCS: EDIT YOUR scs-X.Y.Z/scs.mk, AND FIND AND CHANGE
// THE FOLLOWING VARIABLES TO READ:
//
// ...
// CTRLC = 0
// ...
// USE_OPENMP = 0
// ...
// USE_LAPACK = 1
// ...
//


namespace SCS {
  extern "C" {
#   include <scs.h>
#   include <linSys.h>
#   include <linsys/amatrix.h>
  }
} // namespace SCS





namespace tomo_internal {
inline SCS::scs_int ltrilinindex(SCS::scs_int i, SCS::scs_int j, SCS::scs_int d)
{
  assert(i >= j && i < d);
  // Drafts & Calculations Vol VI ~60% 11/29/2016
  //return j*(d-j+1) + j*(j-1)/2 + i - j;
  return d*j - (j*j+j)/2 + i;
}
}



template<typename RealScalarType_, typename BaseLoggerType_ = Tomographer::Logger::VacuumLogger,
         bool CheckInputs_ = true>
class DiamondNormSCSSolver
{
public:
  typedef RealScalarType_ RealScalarType;
  typedef std::complex<RealScalarType> ComplexScalarType;

  static_assert( std::is_same<RealScalarType, SCS::scs_float>::value,
                 "DiamondNormSCSSolver only supports the scs_float scalar type SCS was compiled with.");

  typedef BaseLoggerType_ BaseLoggerType;

  static constexpr bool CheckInputs = CheckInputs_;

protected:

  const SCS::scs_int DimX;
  const SCS::scs_int DimY;
  const SCS::scs_int DimXY;

  const SCS::scs_float epsilon;

  Eigen::SparseMatrix<SCS::scs_float, Eigen::ColMajor, SCS::scs_int> A;
  Eigen::Matrix<SCS::scs_float, Eigen::Dynamic, 1> bVec;
  Eigen::Matrix<SCS::scs_float, Eigen::Dynamic, 1> cVec;
  
  SCS::Cone * cone;
  SCS::Data * data;

  SCS::Sol * sol;
  SCS::Info * info;

  SCS::Work * work;

  const std::vector<SCS::scs_int> sdpcones;

  Tomographer::Logger::LocalLogger<BaseLoggerType> _logger;

private:
  // helpers for the constructor
  
  inline static constexpr SCS::scs_int total_constraint_dof(SCS::scs_int DimX, SCS::scs_int /*DimY*/, SCS::scs_int DimXY) {
    // number of constraints degrees of freedom = dC1x + dC2x + dC3x = dC1x + dC2x * 2
    return DimX*(2*DimX+1) + DimXY*(2*DimXY+1) * 2;
  }

  inline void init_cone()
  {
    // zero cone {x | x = 0 } (dual to the free cone {x | x in R})
    cone->f = 0;
    // positive orthant {x | x >= 0}
    cone->l = 0;
    // second-order cone {(t,x) | ||x||_2 <= t}
    cone->qsize = 0;
    cone->q = NULL;
    // positive semidefinite cone { X | min(eig(X)) >= 0, X = X^T }
    cone->ssize = sdpcones.size();
    cone->s = const_cast<SCS::scs_int*>(sdpcones.data());
    // exponential cone {(x,y,z) | y e^(x/y) <= z, y>0 }
    cone->ep = 0;
    // dual exponential cone {(u,v,w) | −u e^(v/u) <= e w, u<0}
    cone->ed = 0;
    // power cone {(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}
    cone->psize = 0;
    // dual power cone {(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
    cone->p = NULL;
  }

  inline void init_data()
  {
    // store the A matrix in the data:
    data->m = A.rows(); // rows of A
    data->n = A.cols(); // columns of A

    data->A = (SCS::AMatrix*)std::calloc(1, sizeof(SCS::AMatrix));
    data->A->m = data->m;
    data->A->n = data->n;
    // use Eigen's support of Sparse matrices to get correct representation.
    data->A->x = A.valuePtr();
    data->A->i = A.innerIndexPtr();
    data->A->p = A.outerIndexPtr();

    data->b = bVec.data();

    data->c = cVec.data();

    // we should probably expose some API to tune these values -- these are hard-coded for now.
    data->stgs = (SCS::Settings*)std::malloc(sizeof(SCS::Settings));
    data->stgs->normalize = 1;
    data->stgs->scale = 5.0;
    data->stgs->rho_x = 1e-3;
    data->stgs->max_iters = 2500;
    data->stgs->eps = epsilon;
    data->stgs->alpha = 1.5;
    data->stgs->cg_rate = 2; // no effect as we use direct method and not indirect
    data->stgs->verbose = 0;
    data->stgs->warm_start = 0; // warm_start will be set after the first solution
  }

  inline void init_work()
  {
    work = SCS::scs_init(data, cone, info);
  }

public:
  DiamondNormSCSSolver(const DiamondNormSCSSolver & copy)
    : DimX(copy.DimX),
      DimY(copy.DimY),
      DimXY(copy.DimXY),
      epsilon(copy.epsilon),
      A(copy.A),
      bVec(copy.bVec), // make sure data is initialized, at least outside of the segment updated in calculate()
      cVec(copy.cVec),
      cone(NULL),
      data(NULL),
      sol(NULL),
      info(NULL),
      work(NULL),
      sdpcones(copy.sdpcones),
      _logger(copy._logger)
  {
    auto logger = _logger.subLogger(TOMO_ORIGIN);
    logger.debug("Constructing diamond norm SCS solver -- copy constructor");

    // allocate SCS structures
    cone = (SCS::Cone*)std::calloc(1, sizeof(SCS::Cone));
    data = (SCS::Data*)std::calloc(1, sizeof(SCS::Data));
    sol = (SCS::Sol*)std::calloc(1, sizeof(SCS::Sol));
    info = (SCS::Info*)std::calloc(1, sizeof(SCS::Info));

    // init the cone
    init_cone();

    // set up the SCS data for A (the Eigen sparse matrix A was copied from other object)
    init_data();

    // initializing working space for SCS
    init_work();
  }

  DiamondNormSCSSolver(DiamondNormSCSSolver && other)
    : DimX(other.DimX),
      DimY(other.DimY),
      DimXY(other.DimXY),
      epsilon(other.epsilon),
      A(std::move(other.A)),
      bVec(std::move(other.bVec)), // make sure data is initialized, at least outside of the segment updated in calculate()
      cVec(std::move(other.cVec)),
      cone(other.cone),
      data(other.data),
      sol(other.sol),
      info(other.info),
      work(other.work),
      sdpcones(std::move(other.sdpcones)),
      _logger(other._logger)
  {
    auto logger = _logger.subLogger(TOMO_ORIGIN);
    logger.debug("Constructing diamond norm SCS solver -- move constructor");

    // leave other in a well-defined, null state
    other.cone = NULL;
    other.data = NULL;
    other.sol = NULL;
    other.info = NULL;
    other.work = NULL;
  }

  
  DiamondNormSCSSolver(const SCS::scs_int sysDimX, const SCS::scs_int sysDimY, 
                       const SCS::scs_float epsilon_, BaseLoggerType & baselogger)
    : DimX(sysDimX),
      DimY(sysDimY),
      DimXY(sysDimX * sysDimY),
      epsilon(epsilon_),
      A(total_constraint_dof(DimX, DimY, DimXY), 1+DimXY*DimXY),
      bVec(total_constraint_dof(DimX, DimY, DimXY)),
      cVec(1+DimXY*DimXY),
      cone(NULL),
      data(NULL),
      sol(NULL),
      info(NULL),
      work(NULL),
      sdpcones{2*DimX, 2*DimXY, 2*DimXY},
      _logger(TOMO_ORIGIN, baselogger)
  {
    auto logger = _logger.subLogger(TOMO_ORIGIN);
    logger.debug("Constructing diamond norm SCS solver.");

    // allocate SCS structures
    cone = (SCS::Cone*)std::calloc(1, sizeof(SCS::Cone));
    data = (SCS::Data*)std::calloc(1, sizeof(SCS::Data));
    sol = (SCS::Sol*)std::calloc(1, sizeof(SCS::Sol));
    info = (SCS::Info*)std::calloc(1, sizeof(SCS::Info));
    
    // cones:
    // ------------------------------
    
    init_cone();

    // prepare the A matrix:
    // ------------------------------
    
    typedef Eigen::Triplet<SCS::scs_float> TT;
    std::vector<TT> AT; // prepare triplets for sparse A matrix

    // A has:
    //   **  number of rows = # of constraint degrees of freedom == the three SDP cones
    //       == dC1x + dC2x + dC3x rows  (dC1=2*DimX, dC2=2*DimXY, dC3=2*DimXY,  dC?x=dC?*(dC?+1)/2)
    //   **  number of columns = 1 + DimXY*DimXY = # of variable degrees of freedom
    //       == 1 + dzR + dzI  (dzR = dXY*(dXY+1)/2 , dzI = dXY*(dXY-1)/2)
    //

    // dimension of each constraint block (matrix dimension)
    const SCS::scs_int dC1 = 2*DimX;
    const SCS::scs_int dC2 = 2*DimXY;
    const SCS::scs_int dC3 = dC2;
    // # of degrees of freedom of each constraint block (== number of corresponding rows in A)
    const SCS::scs_int dC1x = dC1*(dC1+1)/2;
    const SCS::scs_int dC2x = dC2*(dC2+1)/2;
    const SCS::scs_int dC3x = dC2x;

    const SCS::scs_int con_offset_1 = 0;
    const SCS::scs_int con_offset_2 = con_offset_1+dC1x;
    const SCS::scs_int con_offset_3 = con_offset_2+dC2x;
  
    const SCS::scs_int var_alpha_offset = 0;
    const SCS::scs_int var_zRdiag_offset = var_alpha_offset + 1;
    const SCS::scs_int var_zRtri_offset = var_zRdiag_offset + DimXY;
    const SCS::scs_int var_zItri_offset = var_zRtri_offset + DimXY*(DimXY-1)/2;

    // A(0:dC1*,0) = A("block #1","alpha") = -vec(\Ident_{2dX})
    SCS::scs_int k = 0;
    for (SCS::scs_int j = 0; j < 2*DimX; ++j) {
      AT.push_back(TT(k, var_alpha_offset, SCS::scs_float(-1)));
      k += 2*DimX - j;
    }

    const SCS::scs_float SQRT2 = boost::math::constants::root_two<SCS::scs_float>();

    using tomo_internal::ltrilinindex;

    // A(0:dC1x,1:1+dzR) = A("block #1","zR"),
    // A(0:dC1x,1+dzR:1+dzR+dzI) = A("block #1","zI"),
    // A(dC1x:dC1x+dC2x,1:1+dzR) = A("block #2","zR"),
    // A(dC1x:dC1x+dC2x,1+dzR:1+dzR+dzI) = A("block #2","zI"),
    // A(dC1x+dC2x:dC1x+dC2x+dC3x,1:1+dzR) = A("block #3","zR"),
    // A(dC1x+dC2x:dC1x+dC2x+dC3x,1+dzR:1+dzR+dzI) = A("block #3","zI"),
    //
    SCS::scs_int trivarno = 0;
    for (SCS::scs_int ij = 0; ij < DimXY; ++ij) {
      const SCS::scs_int i = ij / DimY;
      const SCS::scs_int j = ij % DimY;
      for (SCS::scs_int ipjp = 0; ipjp < ij; ++ipjp) {
        const SCS::scs_int ip = ipjp / DimY;
        const SCS::scs_int jp = ipjp % DimY;
        // deal with z^R_{(ij),(i'j')} and z^I_{(ij),(i'j')}
        //std::cout << "Initializing i,j="<<i<<","<<j<<"  i',j'="<<ip<<","<<jp<<" (trivarno="<<trivarno<<") ... \n";
        // --------
        // BLOCK 1:
        //std::cout << "\tblock #1 ...\n";
        if (j == jp) {
          // zR:
          AT.push_back(TT(con_offset_1+ltrilinindex(i, ip, dC1), var_zRtri_offset+trivarno, SQRT2));
          AT.push_back(TT(con_offset_1+ltrilinindex(DimX+i, DimX+ip, dC1), var_zRtri_offset+trivarno, SQRT2));
          // zI:
          AT.push_back(TT(con_offset_1+ltrilinindex(DimX+i, ip, dC1), var_zItri_offset+trivarno, SQRT2));
          AT.push_back(TT(con_offset_1+ltrilinindex(DimX+ip, i, dC1), var_zItri_offset+trivarno, -SQRT2));
        }
        // BLOCK 2:
        //std::cout << "\tblock #2 ...\n";
        // zR:
        AT.push_back(TT(con_offset_2+ltrilinindex(ij, ipjp, dC2), var_zRtri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_2+ltrilinindex(DimXY+ij, DimXY+ipjp, dC2), var_zRtri_offset+trivarno, -SQRT2));
        // zI:
        AT.push_back(TT(con_offset_2+ltrilinindex(DimXY+ij, ipjp, dC2), var_zItri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_2+ltrilinindex(DimXY+ipjp, ij, dC2), var_zItri_offset+trivarno, SQRT2));
        // BLOCK 3 (same as block 2 for F_{z^R...} and F_{z^I...}:
        //std::cout << "\tblock #3 ...\n";
        // zR:
        AT.push_back(TT(con_offset_3+ltrilinindex(ij, ipjp, dC3), var_zRtri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_3+ltrilinindex(DimXY+ij, DimXY+ipjp, dC3), var_zRtri_offset+trivarno, -SQRT2));
        // zI:
        AT.push_back(TT(con_offset_3+ltrilinindex(DimXY+ij, ipjp, dC3), var_zItri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_3+ltrilinindex(DimXY+ipjp, ij, dC3), var_zItri_offset+trivarno, SQRT2));
        // --------
        ++trivarno;
      }
      // deal with z^R_{(ij),(ij)}:
      //std::cout << "Initializing diagonal i,j="<<i<<","<<j<<" = i',j'  ij=" << ij <<" ... \n";
      // BLOCK 1:
      AT.push_back(TT(con_offset_1+ltrilinindex(i, i, dC1), var_zRdiag_offset+ij, SCS::scs_float(1)));
      AT.push_back(TT(con_offset_1+ltrilinindex(DimX+i, DimX+i, dC1), var_zRdiag_offset+ij, SCS::scs_float(1)));
      // BLOCK 2:
      AT.push_back(TT(con_offset_2+ltrilinindex(ij, ij, dC2), var_zRdiag_offset+ij, SCS::scs_float(-1)));
      AT.push_back(TT(con_offset_2+ltrilinindex(DimXY+ij, DimXY+ij, dC2), var_zRdiag_offset+ij, SCS::scs_float(-1)));
      // BLOCK 3:
      AT.push_back(TT(con_offset_3+ltrilinindex(ij, ij, dC3), var_zRdiag_offset+ij, SCS::scs_float(-1)));
      AT.push_back(TT(con_offset_3+ltrilinindex(DimXY+ij, DimXY+ij, dC3), var_zRdiag_offset+ij, SCS::scs_float(-1)));
    }

    if (CheckInputs) {
      // perform some checks
      for (auto it = AT.begin(); it != AT.end(); ++it) {
        //std::cout << "A(" << it->row() << "," << it->col() << ") = " << it->value() << "\n";
        tomographer_assert(it->row() < dC1x+dC2x+dC3x);
        tomographer_assert(it->col() < 1+DimXY*DimXY);
      }
      //std::cout << "DimX=" << DimX << ", DimXY=" << DimXY
      //          << ", dC1=" << dC1 << ", dC1x=" << dC1x << ", dC2=" << dC2 << ", dC2x=" << dC2x << "\n";
      //std::cout << "dC1x+dC2x+dC3x=" << dC1x+dC2x+dC3x << ",  1+DimXY**2=" << 1+DimXY*DimXY << "\n";
    }

    // Finally, populate the A matrix.
    A.setFromTriplets(AT.begin(), AT.end());
    A.makeCompressed();


    // B vector

    // make sure the bVec is zeroed out, because the calcualte() function will only update
    // the relevant segement inside this vector.
    bVec.setZero();
    
    // C vector
    cVec.setZero();
    cVec(0) = SCS::scs_float(1); // minimize: alpha -- which is the first variable

    // finally, set up the SCS data structure
    init_data();

    // and the SCS working area
    init_work();
  }

  ~DiamondNormSCSSolver()
  {
    // free workspace variables
    if (work != NULL) {
      SCS::scs_finish(work);
    }
    // free other structures
    if (cone != NULL) {
      std::free(cone);
    }
    if (data != NULL) {
      if (data->A != NULL) {
        std::free(data->A);
      }
      if (data->stgs != NULL) {
        std::free(data->stgs);
      }
      std::free(data);
    }
    if (sol != NULL) {
      std::free(sol);
    }
    if (info != NULL) {
      std::free(info);
    }
  }
  
  inline SCS::scs_int dimX() const { return DimX; }
  inline SCS::scs_int dimY() const { return DimY; }
  inline SCS::scs_int dimXY() const { return DimXY; }

  RealScalarType calculate(const Eigen::Ref<const Eigen::Matrix<ComplexScalarType,
                                                                Eigen::Dynamic,Eigen::Dynamic> > & Delta)
  {
    auto logger = _logger.subLogger(TOMO_ORIGIN);
  
    logger.longdebug([&](std::ostream & stream) {
        stream << "Calculating the diamond norm of\n" << Delta;
      });
    
    if (CheckInputs) {
      tomographer_assert((SCS::scs_int)Delta.rows() == DimXY);
      tomographer_assert(Delta.cols() == Delta.rows());
      tomographer_assert( (Delta - Delta.adjoint()).norm() <
                          Delta.rows()*Delta.cols()*Eigen::NumTraits<RealScalarType>::dummy_precision());
    }

    const SCS::scs_int con_offset_2 = DimX*(2*DimX+1);
    const SCS::scs_int dC2x = DimXY*(2*DimXY+1);
    const SCS::scs_float SQRT2 = boost::math::constants::root_two<SCS::scs_float>();

    // copy the lower tri part of
    //
    //   - [ Delta_R, -Delta_I;
    //       Delta_I,  Delta_R ]
    //
    // column-wise onto bVec(con_offset_2:con_offset_2+dC2):
    //
    SCS::scs_int koff = con_offset_2;
    // offset for the second (lower right) DeltaR's columns in linear index
    SCS::scs_int koff2 = con_offset_2 + dC2x - DimXY*(DimXY+1)/2;
    for (SCS::scs_int j = 0; j < DimXY; ++j) {
      const SCS::scs_int col_below_len = DimXY-j-1;
      // DeltaR's j-th diagonal element, twice:
      bVec(koff2) = bVec(koff) = - Delta(j,j).real();
      ++koff;
      ++koff2;
      // DeltaR's lower part of j-th column, twice:
      bVec.segment(koff2, col_below_len) = bVec.segment(koff, col_below_len) =
        - Delta.real().block(j+1, j, col_below_len, 1) * SQRT2;
      koff += col_below_len;
      koff2 += col_below_len;
      // DeltaI's j-th column:
      bVec.segment(koff, DimXY) = - Delta.imag().col(j) * SQRT2;
      koff += DimXY;
    }

    //std::cout << "bVec = " << bVec << "\n";

    // because the SCS data structure directly points to the data in bVec, the SCS problem
    // is automatically updated in this way

    // Now solve!
    SCS::scs_int status = SCS::scs_solve(work, data, cone, sol, info);

    // Check that the SDP converged.
    if (status != SCS_SOLVED) {
      logger.warning("SCS could not solve the SDP: ended with status %s (code %d)", info->status, status);
    }

    SCS::scs_float dnorm = (info->pobj + info->dobj)/2;
    logger.longdebug("[%s] code %d solution is dnorm=%f", info->status, (int)status, (double)dnorm);

    // Use the current solution as a warm-start for future calls to calculate()
    if ( ! data->stgs->warm_start ) {
      data->stgs->warm_start = 1;
    }

    return dnorm;
  } // calculate()

}; // class DiamondNormSCSSolver

















#endif
