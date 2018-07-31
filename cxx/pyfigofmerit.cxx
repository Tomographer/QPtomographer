
//
// QPtomographer.figofmerit -- Python module for calculating figures of merit
// like the diamond norm and the entanglement fidelity.
//

#include <tomographerpy/common.h>

#include <tomographer/tools/loggers.h>
#include <tomographer/densedm/dmtypes.h>
#include <tomographer/valuecalculator.h>
#include <tomographer/tomographer_version.h>

#include "diamond_norm_figofmerit.h"
#include "diamond_norm_scs.h"
#include "entglfidelity_figofmerit.h"
#include "worstentglfidelity_figofmerit.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;



static tpy::PyLogger * pylogger = nullptr;


scs_float figofmerit_diamond_norm_dist(
    const Eigen::Matrix<std::complex<scs_float>, Eigen::Dynamic, Eigen::Dynamic> & Delta_XY,
    const int dim_X,
    py::kwargs opt
    )
{
  tomographer_assert(dim_X > 0);

  const int dim_XY = Delta_XY.cols();
  const int dim_Y = dim_XY / dim_X;

  pylogger->longdebug("figofmerit_diamond_norm", [&](std::ostream & stream) {
      stream << "Delta_XY = \n" << Delta_XY
             << "\ndim_X, dim_Y, dim_XY = " << dim_X << ", " << dim_Y << ", " << dim_XY;
    });

  tomographer_assert(dim_X * dim_Y == dim_XY) ;

  scs_float epsilon = opt.attr("get")("epsilon"_s, 1e-6).cast<scs_float>();

  DiamondNormSCSSolver<scs_float, tpy::PyLogger> dnslv(dim_X, dim_Y, epsilon, *pylogger);

  scs_float val = dnslv.calculate(Delta_XY);

  pylogger->longdebug("figofmerit_diamond_norm", [&](std::ostream & stream) {
      stream << "  --> val = " << val;
    });

  return val;
}



tpy::RealScalar figofmerit_avg_entgl_fidelity2(
    const Eigen::Matrix<tpy::ComplexScalar, Eigen::Dynamic, Eigen::Dynamic> & E_XY,
    const int dim_X
    )
{
  tomographer_assert(dim_X > 0);

  const int dim_XY = E_XY.cols();
  const int dim_Y = dim_XY / dim_X;

  pylogger->longdebug("figofmerit_avg_entgl_fidelity2", [&](std::ostream & stream) {
      stream << "E_XY = \n" << E_XY
             << "\ndim_X, dim_Y, dim_XY = " << dim_X << ", " << dim_Y << ", " << dim_XY;
    });

  tomographer_assert(dim_X * dim_Y == dim_XY) ;

  if (dim_X != dim_Y) {
    throw py::value_error("average entanglement fidelity to identity channel: "
                          "input and output dimensions must be equal.");
  }

  return entanglement_fidelity(E_XY, dim_X);
}


scs_float figofmerit_worst_entgl_fidelity2(
    const Eigen::Matrix<std::complex<scs_float>, Eigen::Dynamic, Eigen::Dynamic> & E_XY,
    const int dim_X,
    py::kwargs opt)
{
  tomographer_assert(dim_X > 0);

  const int dim_XY = E_XY.cols();
  const int dim_Y = dim_XY / dim_X;

  pylogger->longdebug("figofmerit_worst_entgl_fidelity2", [&](std::ostream & stream) {
      stream << "E_XY = \n" << E_XY
             << "\ndim_X, dim_Y, dim_XY = " << dim_X << ", " << dim_Y << ", " << dim_XY;
    });

  tomographer_assert(dim_X * dim_Y == dim_XY) ;

  if (dim_X != dim_Y) {
    throw py::value_error("worst-case entanglement fidelity to identity channel: "
                          "input and output dimensions must be equal.");
  }

  scs_float epsilon = opt.attr("get")("epsilon"_s, 1e-6).cast<scs_float>();

  typedef WorstEntglFidelitySCSSolver<scs_float, tpy::PyLogger> Solver;
  Solver fidslv(dim_X, Solver::factorizeChoiMatrix(E_XY), epsilon, *pylogger);

  double fid = fidslv.calculate();

  pylogger->longdebug("figofmerit_worst_entgl_fidelity2", [&](std::ostream & stream) {
      stream << "  --> fid = " << fid;
    });

  return fid;
}



// tpy::RealScalar hello(
//     const Eigen::Matrix<tpy::ComplexScalar, Eigen::Dynamic, Eigen::Dynamic> & Delta_XY,
//     const int dim_X,
//     py::kwargs opt)
// {
//   const int dim_XY = Delta_XY.cols();
//   const int dim_Y = dim_XY / dim_X;
  
//   pylogger->longdebug("hello()", [&](std::ostream & stream) {
//       stream << "cols = " << Delta_XY.cols();
//     });

//   tomographer_assert(k == 1) ;

//   tpy::RealScalar x = opt.attr("get")("x", 0).cast<tpy::RealScalar>();
  
//   return k + Delta_XY.norm() + x;
// }


PYBIND11_MODULE(figofmerit, m)
{
  m.doc() = "Calculate various figures of merit including the diamond norm, "
    "the average entanglement fidelity and the worst-case entanglement fidelity";

  // make sure the tomographer module is loaded & initialized
  tpy::import_tomographer();

  pylogger = new tpy::PyLogger; // ownership will be taken over by Python/PyBind11

  pylogger->initPythonLogger("QPtomographer.figofmerit");
  auto logger = Tomographer::Logger::makeLocalLogger(TOMO_ORIGIN, *pylogger);

  logger.debug("INIT QPTOMOGRAPHER.FIGOFMERIT");

  // the main run call:
  // m.def(
  //     "hello", // function name
  //     &hello,
  //     "DOCSTRING",
  //     py::arg("m"),
  //     py::arg("k")
  //     );
  // the main run call:
  m.def(
      "diamond_norm_dist", // function name
      &figofmerit_diamond_norm_dist, // fn pointer
      // docstring:
      "Calculate half of the diamond norm of the given hermiticity-preserving map. [...]",
      // arguments:
      py::arg("Delta_XY"),
      py::arg("dim_X")
      );

  m.def(
      "avg_entgl_fidelity2", // fn name
      &figofmerit_avg_entgl_fidelity2, // fn pointer
      // docstring:
      "Calculate the (average) entanglement fidelity (squared) of E_XY to the identity channel. [...]",
      // arguments:
      py::arg("E_XY"),
      py::arg("dim_X")
      );

  m.def(
      "worst_entgl_fidelity2", // fn name
      &figofmerit_worst_entgl_fidelity2, // fn pointer
      // docstring:
      "Calculate the worst-case entanglement fidelity (squared) of E_XY to the identity channel. [...]",
      // arguments:
      py::arg("E_XY"),
      py::arg("dim_X")
      );

  logger.debug("defined function(s)") ;

  m.attr("cxxlogger") = pylogger; // ownership is transferred
  logger.debug("added cxxlogger") ;
}
