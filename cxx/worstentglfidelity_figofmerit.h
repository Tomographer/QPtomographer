
#ifndef WORST_ENTANGLEMENT_FIDELITY_FIGOFMERIT_H
#define WORST_ENTANGLEMENT_FIDELITY_FIGOFMERIT_H

#include <iostream>
#include <string>
#include <sstream>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

#include <tomographer/tools/loggers.h>
#include <tomographer/tools/needownoperatornew.h>
#include <tomographer/tools/cxxutil.h>

#include "qptomo_use_scs.h"
#include "utils.h"



template<typename Scalar>
struct my_to_string {
  std::string operator()(Scalar x) const
  {
    return std::to_string(x);
  }
};
template<typename Scalar>
struct my_to_string<std::complex<Scalar> > {
  std::string operator()(std::complex<Scalar> x) const
  {
    return std::to_string(x.real()) + " + " + std::to_string(x.imag()) + "*i";
  }
};



template<typename RealScalarType_, typename BaseLoggerType_ = Tomographer::Logger::VacuumLogger,
         bool CheckInputs_ = true>
class WorstEntglFidelitySCSSolver
{
public:
  typedef RealScalarType_ RealScalarType;
  typedef std::complex<RealScalarType> ComplexScalarType;

  typedef Eigen::Matrix<ComplexScalarType,Eigen::Dynamic,Eigen::Dynamic> MatrixType;

  static_assert( std::is_same<RealScalarType, scs_float>::value,
                 "WorstEntglFidelitySCSSolver only supports the scs_float scalar type SCS was compiled with.");

  typedef BaseLoggerType_ BaseLoggerType;

  static constexpr bool CheckInputs = CheckInputs_;

protected:

  const scs_int DimX;
  const scs_int DimXX;

  const scs_float epsilon;

  Eigen::SparseMatrix<scs_float, Eigen::ColMajor, scs_int> A;
  Eigen::Matrix<scs_float, Eigen::Dynamic, 1> bVec;
  Eigen::Matrix<scs_float, Eigen::Dynamic, 1> cVec;
  
  ScsCone * cone;
  ScsData * data;

  ScsSolution * sol;
  ScsInfo * info;

  ScsWork * work;

  const std::vector<scs_int> sdpcones;

  Tomographer::Logger::LocalLogger<BaseLoggerType> _logger;

private:
  // helpers for the constructor
 
  inline static constexpr scs_int total_constraint_dof(scs_int DimX,
                                                            scs_int DimXX) {
    // number of constraints degrees of freedom = dC0x + dC1x + dC2x
    return 1 + DimX*(2*DimX+1) +  (1+DimXX)*(2*(1+DimXX)+1) ;
  }
  
  inline void init_cone()
  {
    // zero cone {x | x = 0 } (dual to the free cone {x | x in R})
    cone->f = 1; // constraint: tr(rho) == 1
    // positive orthant {x | x >= 0}
    cone->l = 0;
    // second-order cone {(t,x) | ||x||_2 <= t}
    cone->qsize = 0;
    cone->q = NULL;
    // positive semidefinite cone { X | min(eig(X)) >= 0, X = X^T }
    cone->ssize = sdpcones.size();
    cone->s = const_cast<scs_int*>(sdpcones.data());
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

    data->A = (ScsMatrix*)std::calloc(1, sizeof(ScsMatrix));
    data->A->m = data->m;
    data->A->n = data->n;
    // use Eigen's support of Sparse matrices to get correct representation.
    data->A->x = A.valuePtr();
    data->A->i = A.innerIndexPtr();
    data->A->p = A.outerIndexPtr();

    data->b = bVec.data();

    data->c = cVec.data();

    // we should probably expose some API to tune these values -- these are
    // hard-coded for now.
    data->stgs = (ScsSettings*)std::malloc(sizeof(ScsSettings));
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
    work = scs_init(data, cone, info);
  }

public:
  WorstEntglFidelitySCSSolver(const WorstEntglFidelitySCSSolver & copy)
    : DimX(copy.DimX),
      DimXX(copy.DimXX),
      epsilon(copy.epsilon),
      A(copy.A),
      bVec(copy.bVec), // make sure data is initialized, at least outside of the
                       // segment updated in calculate()
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
    logger.debug("Constructing worst-case entanglement fidelity SCS solver -- copy constructor");

    // allocate SCS structures
    cone = (ScsCone*)std::calloc(1, sizeof(ScsCone));
    data = (ScsData*)std::calloc(1, sizeof(ScsData));
    sol = (ScsSolution*)std::calloc(1, sizeof(ScsSolution));
    info = (ScsInfo*)std::calloc(1, sizeof(ScsInfo));

    // init the cone
    init_cone();

    // set up the SCS data for A (the Eigen sparse matrix A was copied from other object)
    init_data();

    // initializing working space for SCS
    init_work();
  }

  WorstEntglFidelitySCSSolver(WorstEntglFidelitySCSSolver && other)
    : DimX(other.DimX),
      DimXX(other.DimXX),
      epsilon(other.epsilon),
      A(std::move(other.A)),
      bVec(std::move(other.bVec)), // make sure data is initialized, at least
                                   // outside of the segment updated in
                                   // calculate()
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
    logger.debug("Constructing worst-case entanglement fidelity SCS solver -- move constructor");

    // leave other in a well-defined, null state
    other.cone = NULL;
    other.data = NULL;
    other.sol = NULL;
    other.info = NULL;
    other.work = NULL;
  }

  
  WorstEntglFidelitySCSSolver(
      const scs_int sysDimX,
      const Eigen::Ref<const Eigen::Matrix<ComplexScalarType,Eigen::Dynamic,Eigen::Dynamic> > & M_YR,
      const scs_float epsilon_,
      BaseLoggerType & baselogger
      )
    : DimX(sysDimX),
      DimXX(sysDimX * sysDimX),
      epsilon(epsilon_),
      A(total_constraint_dof(DimX, DimXX), 1+DimX*DimX),
      bVec(total_constraint_dof(DimX, DimXX)),
      cVec(1+DimX*DimX),
      cone(NULL),
      data(NULL),
      sol(NULL),
      info(NULL),
      work(NULL),
      sdpcones{2*DimX, 2*(1+DimXX)},
      _logger(TOMO_ORIGIN, baselogger)
  {
    auto logger = _logger.subLogger(TOMO_ORIGIN);
    logger.debug("Constructing worst-case entanglement fidelity SCS solver.");


    // allocate SCS structures
    cone = (ScsCone*)std::calloc(1, sizeof(ScsCone));
    data = (ScsData*)std::calloc(1, sizeof(ScsData));
    sol = (ScsSolution*)std::calloc(1, sizeof(ScsSolution));
    info = (ScsInfo*)std::calloc(1, sizeof(ScsInfo));
    
    // util stuff
    using tomo_internal::ltrilinindex;
    const scs_float SQRT2 = boost::math::constants::root_two<scs_float>();
    const scs_float ONE = scs_float(1);


    // cones:
    // ------------------------------
    
    init_cone();

    
    //
    // We solve the following semidefinite program:
    //
    // Given: E_YR, factorized as E_YR = M_YR * M_YR'
    //
    // variables mu (real), rho (dX*dX complex hermitian)
    // 
    // minimize  mu
    //
    // constraint 0:  tr(rho) == 1
    //
    // constraint 1:  rho >= 0  (cplx. pos. semidef.)
    //
    // constraint 2:
    //    [              |             ]
    //    [      1       |  M vec(rho) ]
    //    [              |             ]    >=   0   (cplx. pos. semidef.)
    //    [--------------|-------------]
    //    [ vec(rho)'*M' |     mu      ]
    //
    //

    // Notes about the C++/SCS formulation & notation.
    //
    //  - rename "rho" -> "z" for variable names, it's simpler
    //

    // Variable degrees of freedom represent the following quantities:
    //
    // [ mu (1) |  zRdiag (dimX)  |   zRtri (dimX*(dimX-1)/2)   |   zItri (dimX*(dimX-1)/2)  ]
    //
    // zRtri and zItri are *strictly lower triangular elements* of rho(=z),
    // picked column-wise. (I.e., iterate linearly using  for (ip=0..d) { for (i=ip+1..d) { ... } })

    const scs_int dvmu = 1;
    const scs_int dvzRdiag = DimX;
    const scs_int dvzRtri = DimX*(DimX-1)/2;
    //const scs_int dvzItri = dvzRtri; // == DimX*(DimX-1)/2;

    const scs_int var_mu_offset = 0;
    const scs_int var_zRdiag_offset = var_mu_offset + dvmu;
    const scs_int var_zRtri_offset = var_zRdiag_offset + dvzRdiag;
    const scs_int var_zItri_offset = var_zRtri_offset + dvzRtri;

    // dimension of each semidefinite constraint block (matrix dimension)
    const scs_int dC1 = 2*DimX;
    const scs_int dC2 = 2*(1+DimXX);
    // # of degrees of freedom of each constraint block (== number of corresponding rows in A)
    const scs_int dC0x = 1;
    const scs_int dC1x = dC1*(dC1+1)/2;
    const scs_int dC2x = dC2*(dC2+1)/2;
    // offset for each constraint number
    const scs_int con_offset_0 = 0;
    const scs_int con_offset_1 = con_offset_0+dC0x;
    const scs_int con_offset_2 = con_offset_1+dC1x;

    // Objective:
    // ----------

    cVec.setZero();
    cVec(0) = ONE; // minimize mu


    // Prepare the A matrix and bVec:
    // ------------------------------
    
    typedef Eigen::Triplet<scs_float> TT;
    std::vector<TT> AT; // prepare triplets for sparse A matrix
    // total number of entries in the A matrix (just count the number of
    // push_back's below. Doesn't have to be exact, actually, an estimate would
    // be enough:
    AT.reserve(DimX + 2*dvzRdiag + 2*DimX*(DimX-1) + 2 + 8*DimX*DimX*DimXX );

    bVec.setZero();

    // A has:
    //   **  number of rows = # of constraint degrees of freedom == one zero cone + two SDP cones
    //       == 1 + dC1x + dC2x  (dC1=2*DimX,  dC2=2*(1+DimXX),  dC?x=dC?*(dC?+1)/2)
    //   **  number of columns = 1 + DimX*DimX = # of variable degrees of freedom
    //       == 1 + dzR + dzI  (dzR = DimX*(DimX+1)/2 , dzI = DimX*(DimX-1)/2)

    scs_int i, ip, jp, ij, k;
    scs_int trivarno;

    logger.longdebug([&](std::ostream & stream) {
        stream << "constraint 0 ... ";
      });

    // Constraint 0
    // ------------
    // A("blck #0", "zRdiag") = 1     [ for "tr(rho) == 1" ]
    for (k = 0; k < DimX; ++k) {
      AT.push_back(TT(0, var_zRdiag_offset+k, ONE));
    }

    bVec(0) = ONE;

    logger.longdebug([&](std::ostream & stream) {
        stream << "constraint 1 ... ";
      });

    // Constraint 1
    // ------------
    for (k = 0; k < dvzRdiag; ++k) {
      AT.push_back(TT(con_offset_1+ltrilinindex(k, k, dC1), var_zRdiag_offset+k, -ONE));
      AT.push_back(TT(con_offset_1+ltrilinindex(DimX+k, DimX+k, dC1), var_zRdiag_offset+k, -ONE));
    }
    trivarno = 0;
    for (ip = 0; ip < DimX; ++ip) {
      for (i = ip+1; i < DimX; ++i) {
        AT.push_back(TT(con_offset_1+ltrilinindex(i, ip, dC1), var_zRtri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_1+ltrilinindex(DimX+i, DimX+ip, dC1), var_zRtri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_1+ltrilinindex(DimX+i, ip, dC1), var_zItri_offset+trivarno, -SQRT2));
        AT.push_back(TT(con_offset_1+ltrilinindex(DimX+ip, i, dC1), var_zItri_offset+trivarno, SQRT2));
        ++trivarno;
      }
    }

    logger.longdebug([&](std::ostream & stream) {
        stream << "constraint 2 ... ";
      });

    // Constraint 2
    // ------------
    // the diagonal of this constraint is easy, place mu's at two places
    AT.push_back(TT(con_offset_2+ltrilinindex(DimXX, DimXX, dC2), var_mu_offset, -ONE));
    AT.push_back(TT(con_offset_2+ltrilinindex(1+2*DimXX, 1+2*DimXX, dC2), var_mu_offset, -ONE));
    // for the rest of the matrix we gotta work more
    // first place the diagonal elements rho(k,k) where they belong in A
    for (k = 0; k < dvzRdiag; ++k) {
      for (ij = 0; ij < DimXX; ++ij) {
        auto const valR = SQRT2*M_YR(DimX*k+k,ij).real();
        auto const valI = SQRT2*M_YR(DimX*k+k,ij).imag();
        // "-Re(<\rho|M^\dagger)" at (DimXX,0:DimXX) in PSD constraint matrix
        AT.push_back(TT(con_offset_2+ltrilinindex(DimXX, ij, dC2), var_zRdiag_offset+k, -valR));
        // "-Re(<\rho|M^\dagger)" at (1+2*DimXX,1+DimXX:1+2*DimXX) in PSD constraint matrix
        AT.push_back(TT(con_offset_2+ltrilinindex(1+2*DimXX, 1+DimXX+ij, dC2), var_zRdiag_offset+k, -valR));
        // "-Im(M|rho>) at (1+DimXX:1+2*DimXX,DimXX) in PSD constraint matrix
        AT.push_back(TT(con_offset_2+ltrilinindex(1+DimXX+ij, DimXX, dC2), var_zRdiag_offset+k, valI));
        // "-Im(<\rho|M^\dagger) = (Im(M|rho>))^T" at (1+2*DimXX,0:DimXX) in PSD constraint matrix
        AT.push_back(TT(con_offset_2+ltrilinindex(1+2*DimXX, ij, dC2), var_zRdiag_offset+k, -valI));
      }
    }
    // then place the remaining elements rho(i,j) chosen columnwise
    trivarno = 0;
    // iteration (ip,jp) over strict lower tri part of rho
    for (jp = 0; jp < DimX; ++jp) {
      for (ip = jp+1; ip < DimX; ++ip) {
        // iteration (ij) over corresponding row/column in PSD constraint where
        // we have to encode "Re(M vec(rho))" or "Im(M vec(rho))"
        for (ij = 0; ij < DimXX; ++ij) {
          const auto valR  = SQRT2*M_YR(DimX*ip+jp,ij).real(); // Re[M]_{(i'j'),(ij)}
          const auto valRt = SQRT2*M_YR(DimX*jp+ip,ij).real(); // Re[M]_{(j'i'),(ij)}
          const auto valI  = SQRT2*M_YR(DimX*ip+jp,ij).imag();
          const auto valIt = SQRT2*M_YR(DimX*jp+ip,ij).imag();
          // "-Re(<\rho|M^\dagger)" at (DimXX,0:DimXX) in PSD constraint matrix
          AT.push_back(TT(con_offset_2 + ltrilinindex(DimXX, ij, dC2), var_zRtri_offset+trivarno, -valR-valRt));
          AT.push_back(TT(con_offset_2 + ltrilinindex(DimXX, ij, dC2), var_zItri_offset+trivarno, -valI+valIt));
          // "-Re(<\rho|M^\dagger)" at (1+2*DimXX,1+DimXX:1+2*DimXX) in PSD constraint matrix
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+2*DimXX, 1+DimXX+ij, dC2), var_zRtri_offset+trivarno,
                          -valR-valRt));
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+2*DimXX, 1+DimXX+ij, dC2), var_zItri_offset+trivarno,
                          -valI+valIt));
          // "-Im(M|rho>) at (1+DimXX:1+2*DimXX,DimXX) in PSD constraint matrix
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+DimXX+ij, DimXX, dC2), var_zRtri_offset+trivarno,
                          valI+valIt));
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+DimXX+ij, DimXX, dC2), var_zItri_offset+trivarno,
                          -valR+valRt));
          // "-Im(<\rho|M^\dagger) = (Im(M|rho>))^T" at (1+2*DimXX,0:DimXX) in PSD constraint matrix
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+2*DimXX, ij, dC2), var_zRtri_offset+trivarno,
                          -valI-valIt));
          AT.push_back(TT(con_offset_2 + ltrilinindex(1+2*DimXX, ij, dC2), var_zItri_offset+trivarno,
                          valR-valRt));
        }
        ++trivarno;
      }
    }

    // constants: Identities on the diagonal except where there are mu's
    for (k = 0; k < DimXX; ++k) {
      bVec(con_offset_2+ltrilinindex(k, k, dC2)) = ONE;
      bVec(con_offset_2+ltrilinindex(1+DimXX+k, 1+DimXX+k, dC2)) = ONE;
    }
    
    if (CheckInputs) {
      // perform some checks
      for (auto it = AT.begin(); it != AT.end(); ++it) {
        //std::cout << "A(" << it->row() << "," << it->col() << ") = " << it->value() << "\n";
        tomographer_assert(it->row() < dC0x+dC1x+dC2x);
        tomographer_assert(it->col() < 1+DimX*DimX);
      }
    }

    logger.longdebug([&](std::ostream & stream) {
        stream << "DimX=" << DimX << ", DimXX=" << DimXX << ", dC0x=" << dC0x
               << ", dC1=" << dC1 << ", dC1x=" << dC1x << ", dC2=" << dC2 << ", dC2x=" << dC2x << "\n"
               << "dC0x+dC1x+dC2x=" << dC0x+dC1x+dC2x << ",  1+2*DimXX=" << 1+2*DimXX ;
      });

    // Finally populate the A matrix.
    A.setFromTriplets(AT.begin(), AT.end());
    A.makeCompressed();

    logger.longdebug([&](std::ostream & stream) {
        stream << "finished constructing the problem.  A =\n"
          // format A's rows into the cone structure of our problem
               << fmt_A_matrix()
               << "b = \n" << bVec << "\n"
               << "c = \n" << cVec;
      });
    
    // Finally^{\otimes 2} set up the SCS data structure
    init_data();

    // and the SCS working area
    init_work();
  }

  ~WorstEntglFidelitySCSSolver()
  {
    // free workspace variables
    if (work != NULL) {
      scs_finish(work);
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

  template<typename LoggerType>
  static inline MatrixType factorizeChoiMatrix(const Eigen::Ref<const MatrixType> & E_YR,
                                               LoggerType & logger)
  {
    // factorize E_YR into M*M'

    Eigen::LDLT<MatrixType> E_YR_ldlt(E_YR);
    MatrixType MM(E_YR.rows(), E_YR.cols()); // M = P^T * L * \sqrt{D}
    MM.setZero();
    MM = E_YR_ldlt.vectorD().array().sqrt().matrix().asDiagonal();
    MM = E_YR_ldlt.matrixL() * MM;
    MM = E_YR_ldlt.transpositionsP().transpose() * MM;

    // std::cerr << "DEBUG ACTIVATED USING SQRTM...\n";
    // MatrixType MM = MatrixType(E_YR).sqrt();

    logger.longdebug("WorstEntglFidelitySCSSolver::factorizeChoiMatrix()",
                     [&](std::ostream & stream) {
        stream << "Calculating the worst-case entanglement fidelity of E_YR =\n" << E_YR << "\n"
               << " ... factorized as M*M' with M =\n" << MM << "\n"
               // << " ... L = \n" << MatrixType(E_YR_ldlt.matrixL()) << "\n ... D = \n"
               // << E_YR_ldlt.vectorD() << "\n ... P = " << E_YR_ldlt.transpositionsP().indices() << "\n"
               << " ... M*M' = \n" << MM*MM.adjoint();
      });

    return MM;
  }

  static inline MatrixType factorizeChoiMatrix(const Eigen::Ref<const MatrixType> & E_YR)
  {
    return factorizeChoiMatrix(E_YR, Tomographer::Logger::vacuum_logger) ;
  }
  
  inline scs_int dimX() const { return DimX; }
  inline scs_int dimXX() const { return DimXX; }

  RealScalarType calculate()
  {
    // NOTE: unlike the diamond norm problem, here the process enters in the
    // main "A" matrix of the SCS problem.  This means that we can't just change
    // the b vector for each new channel and re-solve, we need to recompute the
    // full factorization again for each new channel.  Here the way this happens
    // is that this Solver instance is dedicated to a single channel.
    // Re-instantiate to calculate for a new channel.

    auto logger = _logger.subLogger(TOMO_ORIGIN);
  
    // Now solve!
    scs_int status = scs_solve(work, data, cone, sol, info);

    // Check that the SDP converged.
    if (status != SCS_SOLVED) {
      logger.warning("SCS could not solve the SDP: ended with status %s (code %d)", info->status, status);
    }

    scs_float fid = (info->pobj + info->dobj)/2;
    logger.longdebug("[%s] code %d solution is fid=%f", info->status, (int)status, (double)fid);

    return fid;
  } // calculate()


private:
  
  // debugging tool
  template<typename FormatFn>
  std::string fmt_A_matrix(FormatFn && fn)
  {
    //
    // Format the "A" matrix in a way which ease visual debugging
    //

    MatrixType Amat(A);

    std::vector<std::vector<std::string> > fmtall;
    fmtall.resize(Amat.rows());
    std::size_t w = 0;
    Eigen::Index i, j;
    for (i = 0; i < Amat.rows(); ++i) {
      fmtall[i].resize(Amat.cols());
      for (j = 0; j < Amat.cols(); ++j) {

        fmtall[i][j] = fn(Amat(i,j));

        if (w < fmtall[i][j].size()) {
          w = fmtall[i][j].size();
        }

      }
    }
    w += 3;

    fmtall.insert(fmtall.begin(), std::vector<std::string>());
    fmtall[0].resize(1+DimX*DimX);
    fmtall[0][0] = "mu  ";
    for (i = 0; i < DimX; ++i) {
      fmtall[0][1+i] = "rhoR("+std::to_string(i)+","+std::to_string(i)+")  ";
    }
    int varcount = 0;
    for (j = 0; j < DimX; ++j) {
      for (i = j+1; i < DimX; ++i) {
        fmtall[0][1+DimX+varcount] = "rhoR("+std::to_string(i)+","+std::to_string(j)+")  ";
        fmtall[0][1+DimX+DimX*(DimX-1)/2+varcount] = "rhoI("+std::to_string(i)+","+std::to_string(j)+")  ";
      }
    }

    std::vector<std::string> lines;
    lines.resize(fmtall.size());
    for (i = 0; i < (Eigen::Index)fmtall.size(); ++i) {
      lines[i] = std::string();
      for (j = 0; j < (Eigen::Index)fmtall[i].size(); ++j) {
        lines[i] += std::string(w - fmtall[i][j].size(), ' ') + fmtall[i][j];
      }
    }

    std::ostringstream stream;

    // dimension of each semidefinite constraint block (matrix dimension)
    const scs_int dC1 = 2*DimX;
    const scs_int dC2 = 2*(1+DimXX);
    // # of degrees of freedom of each constraint block (== number of corresponding rows in A)
    //const scs_int dC0x = 1;
    //const scs_int dC1x = dC1*(dC1+1)/2;
    //const scs_int dC2x = dC2*(dC2+1)/2;
    // offset for each constraint number
    //const scs_int con_offset_0 = 0;
    //const scs_int con_offset_1 = con_offset_0+dC0x;
    //const scs_int con_offset_2 = con_offset_1+dC1x;

    // the labels come first

    stream << std::string(lines[0].size()+7, '-') << "\n";
    stream << lines[0] << "\n";
    stream << std::string(lines[0].size()+7, '-') << "\n";

    // the zero cone first

    int offset = 1;

    stream << lines[offset] << "\n";
    stream << std::string(lines[offset].size()+7, '-') << "\n";

    ++offset;

    std::vector<int> conerealdims{ dC1, dC2 };

    for (int cone = 0; cone < 2; ++cone) {

      for (int k = conerealdims[cone]; k > 0; --k) {
        stream << lines[offset] << "  < #" << std::to_string(conerealdims[cone]-k) << "\n";
        for (int l = 1; l < k; ++l) {
          stream << lines[offset+l] << "\n";
        }
        offset += k;
      }
      stream << std::string(lines[0].size()+7, '-') << "\n";

    }

    return stream.str();
  }

  inline std::string fmt_A_matrix()
  {
    return fmt_A_matrix(my_to_string<typename MatrixType::Scalar>());
  }



}; // class WorstEntglFidelitySCSSolver



template<typename DMTypes_>
class WorstEntglFidelityValueCalculator
{
public:
  typedef DMTypes_ DMTypes;
  typedef typename DMTypes::MatrixType MatrixType;
  typedef typename DMTypes::MatrixTypeConstRef MatrixTypeConstRef;

  typedef typename DMTypes::RealScalar ValueType;

  inline WorstEntglFidelityValueCalculator(const DMTypes dmt, int dimX_, ValueType epsilon_)
    : dimX(dimX_), epsilon(epsilon_)
  {
  }

  inline ValueType getValue(const MatrixTypeConstRef & T) const
  {
    tomographer_assert(T.rows() == T.cols());
    const int dimXY = T.rows();

    // only makes sense to calculate worst-case entanglement fidelity for dimX == dimY
    tomographer_assert(dimX*dimX == dimXY) ;
    const int dimY = dimX;

    // calculate rho_{XY}
    const MatrixType rho_XY = T * T.adjoint();

    typedef Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      MatrixDynType;

    // calculate rho_X -- reduced state on input system
    const MatrixDynType rho_X = partial_trace_B(rho_XY, dimX);
  
    // need \sqrt{\rho_X}
    Eigen::SelfAdjointEigenSolver<MatrixDynType> eigrhoX(rho_X);
    const MatrixDynType rhoXsqrtinvtimesEye =
      Eigen::kroneckerProduct(eigrhoX.operatorInverseSqrt(), MatrixDynType::Identity(dimY,dimY));

    // rho^{-1/2}*T*T'*rho^{-1/2} == E, so we can choose a factorization of E as M = rho^{-1/2}*T
    const MatrixType M = rhoXsqrtinvtimesEye * T;

    //    const MatrixType E = Choi_from_process_matrix(rho_XY, dimX, dimY);


    // calculate the worst-case channel/entanglement fidelity by solving a
    // semidefinite program
    //
    // unlike with the diamond norm, we can't re-use the same matrix
    // factorization from one call to another because the SDP's main matrix
    // depends on the channel E itself. Hopefully this is no big deal
    typedef WorstEntglFidelitySCSSolver<scs_float>  MyWorstEntglFidSCSSolver;
    MyWorstEntglFidSCSSolver fidslv(dimX, M, epsilon, Tomographer::Logger::vacuum_logger);

    return fidslv.calculate();
  }

private:
  const int dimX;
  const ValueType epsilon;
};



template<typename ChannelTypes_>
class WorstEntglFidelityChannelSpaceValueCalculator
{
public:
  typedef ChannelTypes_ ChannelTypes;
  typedef typename ChannelTypes::MatrixType MatrixType;
  typedef typename ChannelTypes::MatrixTypeConstRef MatrixTypeConstRef;

  typedef typename ChannelTypes::RealScalar ValueType;

  inline WorstEntglFidelityChannelSpaceValueCalculator(const ChannelTypes cht_, ValueType epsilon_)
    : cht(cht_), epsilon(epsilon_)
  {
  }

  template<typename VIsometryType>
  inline ValueType getValue(const VIsometryType & Vpt) const
  {
    // TODO: act on & use only triangular part for taking the partial trace

    tomographer_assert(Vpt.cols() == cht.dimX());
    tomographer_assert(Vpt.rows() == cht.dimXY2());

    // only makes sense to calculate worst-case entanglement fidelity for dimX == dimY
    tomographer_assert(cht.dimX()*cht.dimX() == cht.dim()) ;

    MatrixType T(cht.dim(), cht.dim());
    T = remapIsometryToT<MatrixType>(Vpt, cht.dim());

    //    const MatrixType E( T * T.adjoint() ); -- not needed, we can specify
    //                                              the factorization directly!!

    // calculate the worst-case channel/entanglement fidelity by solving a
    // semidefinite program
    //
    // unlike with the diamond norm, we can't re-use the same matrix
    // factorization from one call to another because the SDP's main matrix
    // depends on the channel E itself. Hopefully this is no big deal
    WorstEntglFidelitySCSSolver<scs_float> fidslv(cht.dimX(), T, epsilon,
                                                  Tomographer::Logger::vacuum_logger);

    return fidslv.calculate();
  }

private:
  const ChannelTypes cht;
  const ValueType epsilon;
};




#endif
