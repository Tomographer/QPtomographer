
//
// dnormtomo.bistates -- Python module for determining the distribution of the diamond
// norm for reliable quantum process tomography, "naive" version
//
// Requires the `tomographer` python package to be installed, version >= 5.0
//

#include <tomographerpy/common.h>

#include <tomographer/tools/loggers.h>
#include <tomographer/densedm/dmtypes.h>
#include <tomographer/densedm/indepmeasllh.h>
#include <tomographer/densedm/tspacellhwalker.h>
#include <tomographer/densedm/tspacefigofmerit.h>
#include <tomographer/mhrw.h>
#include <tomographer/mhrwtasks.h>
#include <tomographer/multiproc.h>
#include <tomographer/multiprocthreads.h>
#include <tomographer/tomographer_version.h>

#include <tomographer/mhrw_valuehist_tools.h>
#include <tomographer/mhrwstepsizecontroller.h>
#include <tomographer/mhrwvalueerrorbinsconvergedcontroller.h>

#include <tomographerpy/pyhistogram.h>
#include <tomographerpy/pymultiproc.h>
#include <tomographerpy/pymhrwtasks.h>
#include <tomographerpy/pymhrw.h>
#include <tomographerpy/pygil.h>

#include "diamond_norm_figofmerit.h"
#include "diamond_norm_scs.h"


#include "pybind11/pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;


static tpy::PyLogger * pylogger = nullptr;



// define an exception class for invalid inputs
TOMOGRAPHER_DEFINE_MSG_EXCEPTION(DNormBiStatesInvalidInputError, "Invalid Input: ") ;


//
// Data types for our quantum objects.  For the sake of the example, we just leave the
// size to be dynamic, that is, fixed at run time and not at compile time.
//
typedef Tomographer::DenseDM::DMTypes<Eigen::Dynamic, double> DMTypes;


//
// The class which will store our tomography data. Just define this as "DenseLLH" as a
// shorthand.
//
typedef Tomographer::DenseDM::IndepMeasLLH<DMTypes> DenseLLH;


//
// Our value calculator will simply be our DiamondNormToRefValueCalculator.
//
typedef DiamondNormToRefValueCalculator<
  DMTypes,
  DiamondNormSCSSolver<typename DMTypes::RealScalar>
  > ValueCalculator;


typedef std::mt19937 RngType;

typedef Tomographer::MHRWTasks::ValueHistogramTools::CDataBase<
  ValueCalculator, // our value calculator
  true, // use binning analysis
  Tomographer::MHWalkerParamsStepSize<tpy::RealType>, // MHWalkerParams
  RngType::result_type, // RngSeedType
  tpy::CountIntType, // IterCountIntType
  tpy::RealType, // CountRealType
  tpy::CountIntType // HistCountIntType
  >
  CDataBaseType;

//
// We need to define a class which adds the capacity of creating the "master"
// random walk object to the engine in
// Tomographer::MHRWTasks::ValueHistogramTasks, which take care of running the
// random walks etc. as needed.
//
struct OurCDataBiStates : public CDataBaseType
{
  OurCDataBiStates(
      const DenseLLH & llh_, // data from the the tomography experiment
      ValueCalculator valcalc, // the figure-of-merit calculator
      HistogramParams hist_params, // histogram parameters
      int binning_num_levels, // number of binning levels in the binning analysis
      tpy::MHRWParams mhrw_params, // parameters of the random walk
      RngType::result_type base_seed, // a random seed to initialize the random number generator
      py::dict ctrl_step_size_params_,
      py::dict ctrl_converged_params_)
    : CDataBaseType(
        valcalc, hist_params, binning_num_levels,
        MHRWParamsType(
            tpy::pyMHWalkerParamsFromPyObj<Tomographer::MHWalkerParamsStepSize<tpy::RealType> >(
                mhrw_params.mhwalker_params
                ),
            mhrw_params.n_sweep,
            mhrw_params.n_therm,
            mhrw_params.n_run
            ),
        base_seed
        ),
      llh(llh_),
      ctrl_step_size_params(ctrl_step_size_params_),
      ctrl_converged_params(ctrl_converged_params_)
  {
  }

  const DenseLLH llh;

  const py::dict ctrl_step_size_params;
  const py::dict ctrl_converged_params;

  // The result of a task run -- pass on to the default StatsResults type
  // provided by ValueHistogramTools.  Because we have several stats collectors
  // set, we need to pick out the result of our
  // "value-histogram-stats-collector", which is the first one in the
  // multiple-stats-collector object which we have create in
  // setupRandomWalkAndRun().  We thus pick out the first result in the tuple of
  // all stats collectors results, i.e. the one at index 0 (using std::get).
  struct MHRWStatsResultsType : public MHRWStatsResultsBaseType
  {
    template<typename... Types>
    MHRWStatsResultsType(std::tuple<ValueStatsCollectorResultType, Types...> && r)
      : MHRWStatsResultsBaseType(std::move(std::get<0>(r)))
    {
    }
  };

  //
  // This function is called automatically by the task manager/dispatcher via
  // MHRWTasks.  It should set up the random walk as required (at a minimum,
  // create a MHWalker instance and pass on the default value stats collector
  // from Tomographer::MHRWTasks::ValueHistogramTools), and run it.
  //
  // We should not forget to call run(), to actually run the random walk!
  //
  // Here, our example is slightly more involved as in "minimal_tomorun". In
  // addition, we'll include more stats collectors and set up some random walk
  // controllers.
  //
  template<typename Rng, typename LoggerType, typename ExecFn>
  inline void setupRandomWalkAndRun(Rng & rng, LoggerType & logger, ExecFn run) const
  {
    //
    // NOTE: This function can be called from multiple threads simultaneously as
    // the GIL (Global Interpreter Lock) is currently not held.  Don't call any
    // Python-related function without acquiring the Python GIL, see examples
    // below.  Also, Don't write to global variables.
    //
    // However, do not call logger methods while holding the GIL, as the logger
    // will itself try to acquire the GIL again causing a deadlock.
    //
    // A handy macro is TPY_EXPR_WITH_GIL( expression ), which evaluates the
    // given expression while holding the GIL and returns the result,
    // immediately re-releasing the GIL. Just surround any Python-related calls
    // with that macro.
    //

    auto val_stats_collector = createValueStatsCollector(logger);
    Tomographer::MHRWMovingAverageAcceptanceRatioStatsCollector<> movavg_accept_stats(
        TPY_EXPR_WITH_GIL( ctrl_step_size_params.attr("get")("num_samples", 8192).cast<int>() )
        );
    auto stats_collectors =
      Tomographer::mkMultipleMHRWStatsCollectors(val_stats_collector, movavg_accept_stats);

    auto therm_step_controller =
      Tomographer::mkMHRWStepSizeController<MHRWParamsType>(
          movavg_accept_stats,
          logger);

    std::tuple<int,Eigen::Index,Eigen::Index,Eigen::Index,double> stp =
      TPY_EXPR_WITH_GIL(
          std::make_tuple(
              ctrl_converged_params.attr("get")("check_frequency_sweeps", 1024).cast<int>(),
              ctrl_converged_params.attr("get")("max_allowed_unknown",
                                                1+2*histogram_params.num_bins/100).template cast<Eigen::Index>(),
              ctrl_converged_params.attr("get")("max_allowed_unknown_notisolated",
                                                1+histogram_params.num_bins/100).template cast<Eigen::Index>(),
              ctrl_converged_params.attr("get")("max_allowed_not_converged",
                                                1+histogram_params.num_bins/200).template cast<Eigen::Index>(),
              ctrl_converged_params.attr("get")("max_add_run_iters", 1.5).cast<double>()
              ) ) ;

    auto numsamples_controller =
      Tomographer::mkMHRWValueErrorBinsConvergedController(
          val_stats_collector,
          logger,
          std::get<0>(stp), // check_frequency_sweeps
          // max allowed: unknown, unk-not-isolated, not-converged
          std::get<1>(stp), std::get<2>(stp), std::get<3>(stp),
          // max add run iters
          std::get<4>(stp)
          );

    auto controllers =
      Tomographer::mkMHRWMultipleControllers(therm_step_controller, numsamples_controller);

    Tomographer::DenseDM::TSpace::LLHMHWalker<DenseLLH,Rng,LoggerType> mhwalker(
	llh.dmt.initMatrixType(),
	llh,
	rng,
	logger
	);
 
    run(mhwalker, stats_collectors, controllers);
  };

};



py::dict tomo_run_dnorm_bistates(
    const int dimX,
    const int dimY,
    const py::list& Emn,
    const Eigen::VectorXi& Nm,
    const tpy::HistogramParams& hist_params,
    const tpy::MHRWParams& mhrw_params,
    int binning_num_levels,
    const int num_repeats,
    py::object ref_channel_XY,
    py::object progress_fn,
    const int progress_interval_ms,
    const double dnorm_epsilon,
    py::dict ctrl_step_size_params,
    py::dict ctrl_converged_params
    )
{
  Tomographer::Logger::LocalLogger<tpy::PyLogger> logger(TOMO_ORIGIN, *pylogger);

  logger.debug("tomo_run_dnorm_bistates()");

  typedef DMTypes::MatrixType MatrixType;

  const int dimXY = dimX*dimY;

  DMTypes dmt(dimXY);

  DenseLLH llh(dmt);

  if (py::len(Emn) != (std::size_t)Nm.rows()) {
    throw DNormBiStatesInvalidInputError("Mismatch in number of measurements: len(Emn)="
                                         + std::to_string(py::len(Emn)) +
                                         " but Nm.rows()=" + std::to_string(Nm.rows()));
  }
  for (std::size_t k = 0; k < (std::size_t)Nm.rows(); ++k) {
    MatrixType POVMeffect = Emn[k].cast<MatrixType>();
    llh.addMeasEffect(POVMeffect, Nm(k), true);
  }

  logger.debug([&](std::ostream & ss) {
      ss << "\n\nExn: size="<<llh.Exn().size()<<"\n"
	 << llh.Exn() << "\n";
      ss << "\n\nNx: size="<<llh.Nx().size()<<"\n"
	 << llh.Nx() << "\n";
    });


  //
  // Prepare the figure of merit calculator:  diamond norm to the reference channel.
  //
  DMTypes::MatrixType mat_ref_channel_XY = DMTypes::MatrixType::Zero(dimXY, dimXY);

  if (ref_channel_XY.is_none()) {
    // if None is given, then use the unnormalized maximally entangled state = Choi matrix
    // of the identity channel.  (|Phi> = \sum_k |kk>)
    if (dimX != dimY) {
      // but this can only be done if dimX==dimY
      throw DNormBiStatesInvalidInputError("You must specify a reference channel ref_channel_XY if dimX != dimY");
    }
    for (int k = 0; k < dimX; ++k) {
      for (int k2 = 0; k2 < dimX; ++k2) {
        mat_ref_channel_XY(k + dimX*k, k2 + dimX*k2) = 1;
      }
    }
  } else {
    // extract an Eigen::Matrix from the python object
    mat_ref_channel_XY = ref_channel_XY.cast<DMTypes::MatrixType>();
  }

  logger.debug([&](std::ostream & stream) {
      stream << "Using the reference channel ref_channel_XY =\n" << mat_ref_channel_XY;
    });

  // ValueCalculator = DiamondNormToRefValueCalculator<DMTypes>
  ValueCalculator valcalc(dmt, mat_ref_channel_XY, dimX, dnorm_epsilon);


  // prepare the random walk tasks

  typedef Tomographer::MHRWTasks::MHRandomWalkTask<OurCDataBiStates, RngType>  OurMHRandomWalkTask;

  // seed for random number generator
  auto base_seed = std::chrono::system_clock::now().time_since_epoch().count();

  binning_num_levels = Tomographer::sanitizeBinningLevels(binning_num_levels, mhrw_params.n_run,
                                                          tpy::CountIntType(128), logger) ;

  OurCDataBiStates taskcdat(
      llh, valcalc, hist_params, binning_num_levels, mhrw_params, base_seed,
      ctrl_step_size_params, ctrl_converged_params
      );

  // worry about GIL.

  tpy::GilProtectedPyLogger logger_with_gil(logger.parentLogger(), false);

  Tomographer::MultiProc::CxxThreads::TaskDispatcher<OurMHRandomWalkTask,
                                                     OurCDataBiStates,
                                                     tpy::GilProtectedPyLogger>
    tasks(
      &taskcdat, // constant data
      logger_with_gil, // the main logger object -- automatically acquires the
                       // GIL for emitting messages
      num_repeats // num_runs
      );

  setTasksStatusReportPyCallback(tasks, progress_fn, progress_interval_ms,
                                 true /* GIL */);

  auto time_start = std::chrono::steady_clock::now();

  // release the GIL for multiprocessing
  {
    logger_with_gil.requireGilAcquisition(true);
    py::gil_scoped_release gil_release;

    // and run our tomo process

    try {

      tasks.run();

    } catch (tpy::PyFetchedException & pyerr) {

      // acquire GIL for PyErr_Restore()
      py::gil_scoped_acquire gil_acquire;

      pyerr.restorePyException();
      throw py::error_already_set();
      
    } catch (std::exception & e) {

      // another exception
      logger.debug("Exception: %s", e.what());
      throw; // error via pybind11

    }

  } // gil released scope
  logger_with_gil.requireGilAcquisition(false);

  auto time_end = std::chrono::steady_clock::now();

  logger.debug("Random walks done.");

  // delta-time, in seconds and fraction of seconds
  std::string elapsed_s = Tomographer::Tools::fmtDuration(time_end - time_start);

  py::dict res;

  // individual results from each task
  const auto & task_results = tasks.collectedTaskResults();
  // ... aggregated into a full averaged histogram
  auto aggregated_histogram = taskcdat.aggregateResultHistograms(task_results) ;

  res["final_histogram"] = tpy::HistogramWithErrorBars(aggregated_histogram.final_histogram);
  res["simple_final_histogram"] = tpy::HistogramWithErrorBars(aggregated_histogram.simple_final_histogram);
  res["elapsed_seconds"] = 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(
      time_end - time_start
      ).count();

  py::list runs_results;
  for (std::size_t k = 0; k < task_results.size(); ++k) {
    const auto & run_result = *task_results[k];
    runs_results.append(
        tpy::MHRandomWalkTaskResult(
            py::cast(tpy::ValueHistogramWithBinningMHRWStatsCollectorResult(run_result.stats_results)),
            tpy::MHRWParams(py::dict("step_size"_a=run_result.mhrw_params.mhwalker_params.step_size),
                            run_result.mhrw_params.n_sweep,
                            run_result.mhrw_params.n_therm,
                            run_result.mhrw_params.n_run),
            run_result.acceptance_ratio
            )
        );
  }
  res["runs_results"] = runs_results;

  // full final report
  std::string final_report;
  { std::ostringstream ss;
    Tomographer::MHRWTasks::ValueHistogramTools::printFinalReport(
        ss, // where to output
        taskcdat, // the cdata
        task_results, // the results
        aggregated_histogram // aggregated
        );
    final_report = ss.str();
    res["final_report"] = final_report;
  }

  // final report of runs only
  { std::ostringstream ss;
    Tomographer::MHRWTasks::ValueHistogramTools::printFinalReport(
        ss, // where to output
        taskcdat, // the cdata
        task_results, // the results
        aggregated_histogram, // aggregated
        0, // width -- use default
        false // don't print the histogram
        );
    res["final_report_runs"] = ss.str();
  }

  logger.debug([&](std::ostream & stream) {
      stream << final_report << "\n";
      stream << "Computation time: " <<  elapsed_s << "\n";
    });

  return res;
}





PYBIND11_PLUGIN(bistates)
{
  py::module m("bistates",
               "python interface to run random walk on bipartite states with diamond norm "
               "fig of merit on encoded channel");

  // make sure the tomographer module is loaded & initialized
  tpy::import_tomographer();

  pylogger = new tpy::PyLogger; // ownership will be taken over by Python/PyBind11

  pylogger->initPythonLogger("dnormtomo.bistates");
  auto logger = Tomographer::Logger::makeLocalLogger(TOMO_ORIGIN, *pylogger);

  logger.debug("INIT DNORMTOMO.BISTATES");

  // the main run call:
  m.def(
      "run", // function name
      &tomo_run_dnorm_bistates, // fn pointer
      py::arg("dimX"),
      py::arg("dimY"),
      py::arg("Emn") = py::list(),
      py::arg("Nm") = Eigen::VectorXi(),
      py::arg("hist_params") =  tpy::HistogramParams(),
      py::arg("mhrw_params") = tpy::MHRWParams(),
      py::arg("binning_num_levels") = -1,
      py::arg("num_repeats") = std::thread::hardware_concurrency(),
      py::arg("ref_channel_XY") = py::none(),
      py::arg("progress_fn") = py::none(),
      py::arg("progress_interval_ms") = (int)500,
      py::arg("dnorm_epsilon") = (double)1e-3,
      py::arg("ctrl_step_size_params") = py::dict(),
      py::arg("ctrl_converged_params") = py::dict(),
      "Docstring here"
      );

  logger.debug("defined the main function") ;

  m.attr("cxxlogger") = pylogger; // ownership is transferred

  logger.debug("added cxxlogger") ;

  m.attr("cxx_tomographer_version") = std::string(TOMOGRAPHER_VERSION);

  logger.debug([&](std::ostream & stream) { stream << "tomographer version is " << TOMOGRAPHER_VERSION ; }) ;


  logger.debug("DNormBiStatesInvalidInputError ...");

  py::register_exception<DNormBiStatesInvalidInputError>(m, "DNormBiStatesInvalidInputError");
    
  return m.ptr();
}
