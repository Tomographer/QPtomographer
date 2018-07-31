
//
// QPtomographer.channelspace -- Python module for determining the distribution
// of the diamond norm for reliable quantum process tomography
//
// Requires the `tomographer` python package to be installed, version >= 5.4
//

#include <tomographerpy/common.h>

#include <tomographer/tools/loggers.h>
#include <tomographer/densedm/dmtypes.h>
#include <tomographer/densedm/indepmeasllh.h>
#include <tomographer/densedm/tspacellhwalker.h>
#include <tomographer/densedm/tspacefigofmerit.h>
#include <tomographer/valuecalculator.h>
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
#include "entglfidelity_figofmerit.h"
#include "worstentglfidelity_figofmerit.h"

#include "channelspace.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

#include "pyutils.h"



static tpy::PyLogger * pylogger = nullptr;



// define an exception class for invalid inputs
TOMOGRAPHER_DEFINE_MSG_EXCEPTION(DNormChannelSpaceInvalidInputError, "Invalid Input: ") ;


//
// Data types for our quantum objects.  For the sake of the example, we just leave the
// size to be dynamic, that is, fixed at run time and not at compile time.
//
typedef ChannelTypes<double> MyChannelTypes;


//
// The class which will store our tomography data. Just define this as "DenseLLH" as a
// shorthand.
//
typedef Tomographer::DenseDM::IndepMeasLLH<MyChannelTypes> DenseLLH;



//
// The ValueCalculators -- use a MultiplexorValueCalculator to choose at
// run-time between a DiamondNormToRefChannelSpaceValueCalculator, a
// EntglFidelityChannelSpaceValueCalculator, a
// WorstEntglFidelityChannelSpaceValueCalculator, and a
// CallableChannelSpaceValueCalculators (any python callback)
//

class CallableChannelSpaceValueCalculator
{
public:
  typedef tpy::RealScalar ValueType;
  
  CallableChannelSpaceValueCalculator(py::object fn_) : fn(fn_) { }

  tpy::RealScalar getValue(const MyChannelTypes::VIsometryType & T) const
  {
    py::gil_scoped_acquire gil_acquire;
    qptomographer_checkPyException();
    tpy::RealScalar val = fn(py::cast(T)).cast<tpy::RealScalar>();
    qptomographer_checkPyException();
    return val;
  }

private:
  py::object fn;
};

typedef DiamondNormToRefChannelSpaceValueCalculator<
  MyChannelTypes,
  DiamondNormSCSSolver<MyChannelTypes::RealScalar>
  >
  DNormChannelSpaceValueCalculator;

typedef Tomographer::MultiplexorValueCalculator<
  tpy::RealScalar, // value type first, then:
  DNormChannelSpaceValueCalculator,
  EntglFidelityChannelSpaceValueCalculator<MyChannelTypes>,
  WorstEntglFidelityChannelSpaceValueCalculator<MyChannelTypes>,
  CallableChannelSpaceValueCalculator
  >
  MyValueCalculator;


typedef std::mt19937 RngType;

typedef Tomographer::MHRWTasks::ValueHistogramTools::CDataBase<
  MyValueCalculator, // our value calculator
  true, // use binning analysis
  Tomographer::MHWalkerParamsStepSize<tpy::RealScalar>, // MHWalkerParams
  RngType::result_type, // RngSeedType
  tpy::IterCountIntType, // IterCountIntType
  tpy::CountRealType, // CountRealType
  tpy::HistCountIntType // HistCountIntType
  >
  CDataBaseType;

//
// We need to define a class which adds the capacity of creating the "master"
// random walk object to the engine in
// Tomographer::MHRWTasks::ValueHistogramTools, which take care of running the
// random walks etc. as needed.
//
struct OurCDataChannelSpace : public CDataBaseType
{
  OurCDataChannelSpace(
      const DenseLLH & llh_, // data from the the tomography experiment
      MyValueCalculator valcalc, // the figure-of-merit calculator
      HistogramParams hist_params, // histogram parameters
      int binning_num_levels, // number of binning levels in the binning analysis
      tpy::MHRWParams mhrw_params, // parameters of the random walk
      RngType::result_type base_seed, // a random seed to initialize the random number generator
      py::dict ctrl_step_size_params_,
      py::dict ctrl_converged_params_)
    : CDataBaseType(
        valcalc, hist_params, binning_num_levels,
        MHRWParamsType(
            tpy::pyMHWalkerParamsFromPyObj<Tomographer::MHWalkerParamsStepSize<tpy::RealScalar> >(
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
      ctrl_converged_params(ctrl_converged_params_),
      chwalker_jump_mode(RandHermExp)
  {
  }

  const DenseLLH llh;

  const py::dict ctrl_step_size_params;
  const py::dict ctrl_converged_params;

  JumpMode chwalker_jump_mode;

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

    CtrlConvergedParams<int> stp =
      CtrlConvergedParams<int>::fromPyDictWithGilAcq(ctrl_converged_params, histogram_params);

    auto numsamples_controller =
      Tomographer::mkMHRWValueErrorBinsConvergedController(
          val_stats_collector,
          logger,
          stp.check_frequency_sweeps,
          // max allowed: unknown, unk-not-isolated, not-converged
          stp.max_allowed_unknown, stp.max_allowed_unknown_notisolated, stp.max_allowed_not_converged,
          // max add run iters
          stp.max_add_run_iters
          );

    auto controllers = Tomographer::mkMHRWMultipleControllers(therm_step_controller, numsamples_controller);

    ChannelSpaceMHWalker<DenseLLH,Rng,LoggerType> mhwalker(
	llh,
	rng,
	logger,
        chwalker_jump_mode
	);

    run(mhwalker, stats_collectors, controllers);
  };
};



py::dict tomo_run_dnorm_channels(py::kwargs kwargs)
{
  Tomographer::Logger::LocalLogger<tpy::PyLogger> logger(TOMO_ORIGIN, *pylogger);

  logger.debug("tomo_run_dnorm_channels()");

  logger.debug("level test: DEBUG");
  logger.longdebug("level test: LONGDEBUG");

  typedef MyChannelTypes::MatrixType MatrixType;

  auto pop_mandatory_kwarg = [&kwargs](std::string key) -> py::object {
    if (!kwargs.contains(py::cast(key))) {
      throw DNormChannelSpaceInvalidInputError("Missing required keyword argument `"+key+"'");
    }
    return kwargs.attr("pop")(py::cast(std::move(key)));
  };

  const int dimX = pop_mandatory_kwarg("dimX").cast<int>();
  const int dimY = pop_mandatory_kwarg("dimY").cast<int>();
  py::list Emn = pop_mandatory_kwarg("Emn").cast<py::list>();
  const Eigen::VectorXi Nm = pop_mandatory_kwarg("Nm").cast<Eigen::VectorXi>();
  const tpy::HistogramParams hist_params = pop_mandatory_kwarg("hist_params").cast<tpy::HistogramParams>();
  const tpy::MHRWParams mhrw_params = pop_mandatory_kwarg("mhrw_params").cast<tpy::MHRWParams>();
  py::object fig_of_merit = kwargs.attr("pop")("fig_of_merit"_s, py::none());
  py::object ref_channel_XY = kwargs.attr("pop")("ref_channel_XY"_s, py::none());
  const double sdp_epsilon = kwargs.attr("pop")("sdp_epsilon"_s, 1e-3).cast<double>();
  const int channel_walker_jump_mode = kwargs.attr("pop")("channel_walker_jump_mode"_s, (int)RandHermExp).cast<int>();
  int binning_num_levels = kwargs.attr("pop")("binning_num_levels"_s, -1).cast<int>();
  const int num_repeats = kwargs.attr("pop")("num_repeats"_s, std::thread::hardware_concurrency()).cast<int>();
  py::object progress_fn = kwargs.attr("pop")("progress_fn"_s, py::none());
  const int progress_interval_ms = kwargs.attr("pop")("progress_interval_ms"_s, 500).cast<int>();
  py::dict ctrl_step_size_params = kwargs.attr("pop")("ctrl_step_size_params"_s, py::dict()).cast<py::dict>();
  py::dict ctrl_converged_params = kwargs.attr("pop")("ctrl_converged_params"_s, py::dict()).cast<py::dict>();

  if (py::len(kwargs)) {
    throw DNormChannelSpaceInvalidInputError("Unknown extra arguments given: " +
                                             (", "_s.attr("join")(kwargs.attr("keys")())).cast<std::string>()) ;
  }

  MyChannelTypes dmt(dimX, dimY);

  const int dimXY = dmt.dim();

  DenseLLH llh(dmt);

  if ((std::size_t)py::len(Emn) != (std::size_t)Nm.rows()) {
    throw DNormChannelSpaceInvalidInputError("Mismatch in number of measurements: len(Emn)="
                                             + std::to_string(py::len(Emn)) +
                                             " but Nm.rows()=" + std::to_string(Nm.rows()));
  }
  for (std::size_t k = 0; k < (std::size_t)Nm.rows(); ++k) {
    MatrixType POVMeffect = Emn[k].cast<MatrixType>();
    llh.addMeasEffect(POVMeffect, Nm(k), true);
  }

  logger.debug([&](std::ostream & ss) {
      ss << "Measurement data is:\n\nExn: size="<<llh.Exn().size()<<"\n"
	 << llh.Exn() << "\n";
      ss << "\n\nNx: size="<<llh.Nx().size()<<"\n"
	 << llh.Nx() << "\n";
    });


  if (num_repeats < 1) {
    throw py::value_error("num_repeats must be >= 1") ;
  }

  logger.debug([&](std::ostream & stream) {
      stream << "About to consider figure of merit...";
    });

  //
  // Prepare the figure of merit calculator:  diamond norm to the reference channel.
  //
  MyChannelTypes::MatrixType mat_ref_channel_XY = MyChannelTypes::MatrixType::Zero(dimXY, dimXY);

  int fig_of_merit_multiplexor_idx = 0;

  if (!fig_of_merit.is_none() && py::hasattr(fig_of_merit, "__call__"_s)) {
    // custom callable. Consider this case first so we don't risk calling
    // '__eq__' on a callable object, which pybind11 doesn't like

    // got a callable as argument, the third possibility
    fig_of_merit_multiplexor_idx = 3;

    logger.debug([&](std::ostream & stream) {
        stream << "Using custom callable as figure of merit";
      });

  } else if (fig_of_merit.is_none() || fig_of_merit.attr("__eq__")("diamond-distance"_s).cast<bool>()) {

    // diamond norm distance, the default

    fig_of_merit_multiplexor_idx = 0;

    mat_ref_channel_XY =
      get_ref_channel<MyChannelTypes::MatrixType, DNormChannelSpaceInvalidInputError>(
          dimX, dimY, ref_channel_XY
          );

    logger.debug([&](std::ostream & stream) {
        stream << "Using diamond norm as figure of merit, with ref_channel_XY =\n" << mat_ref_channel_XY;
      });

  } else if (fig_of_merit.attr("__eq__")("entanglement-fidelity"_s).cast<bool>()) {

    // use entanglement fidelity, the second possibility
    fig_of_merit_multiplexor_idx = 1;

    logger.debug([&](std::ostream & stream) {
        stream << "Using entanglement fidelity as figure of merit";
      });

  } else if (fig_of_merit.attr("__eq__")("worst-entanglement-fidelity"_s).cast<bool>()) {

    // use entanglement fidelity, the second possibility
    fig_of_merit_multiplexor_idx = 2;

    logger.debug([&](std::ostream & stream) {
        stream << "Using entanglement fidelity as figure of merit";
      });

  } else {

    throw DNormChannelSpaceInvalidInputError(
            "Invalid value for `fig_of_merit' argument (see documentation)"
            );

  }

  MyValueCalculator valcalc(
      fig_of_merit_multiplexor_idx,
      // creator function for DNormChannelSpaceValueCalculator:
      [dmt,mat_ref_channel_XY,dimX,sdp_epsilon]() {
        // do NOT give a pylogger here, it's not thread-safe.
        return new DNormChannelSpaceValueCalculator(dmt, mat_ref_channel_XY, dimX, sdp_epsilon);
      },
      // creator function for EntanglmentFidelityChannelSpaceValueCalculator:
      [dmt]() {
        return new EntglFidelityChannelSpaceValueCalculator<MyChannelTypes>(dmt);
      },
      // creator function for EntanglmentFidelityChannelSpaceValueCalculator:
      [dmt,sdp_epsilon]() {
        return new WorstEntglFidelityChannelSpaceValueCalculator<MyChannelTypes>(dmt, sdp_epsilon);
      },
      // creator function for CallableChannelSpaceValueCalculator:
      [fig_of_merit]() {
        return new CallableChannelSpaceValueCalculator(fig_of_merit);
      }
      ) ;

  // // DNormValueCalculator is DiamondNormToRefValueCalculator<ChannelTypes>
  // DNormValueCalculator valcalc(dmt, mat_ref_channel_XY, 
  //                              dimX, dnorm_epsilon); // do NOT give a pylogger here, it's not thread-safe.

  // prepare the random walk tasks

  typedef Tomographer::MHRWTasks::MHRandomWalkTask<OurCDataChannelSpace, RngType>  OurMHRandomWalkTask;

  // seed for random number generator
  auto base_seed = std::chrono::system_clock::now().time_since_epoch().count();

  binning_num_levels = Tomographer::sanitizeBinningLevels(binning_num_levels, mhrw_params.n_run,
                                                          tpy::IterCountIntType(128), logger) ;

  OurCDataChannelSpace taskcdat(
      llh, valcalc, hist_params, binning_num_levels, mhrw_params, base_seed,
      ctrl_step_size_params, ctrl_converged_params
      );

  taskcdat.chwalker_jump_mode = JumpMode( channel_walker_jump_mode );

  // worry about GIL.

  tpy::GilProtectedPyLogger logger_with_gil(logger.parentLogger(), false);

  Tomographer::MultiProc::CxxThreads::TaskDispatcher<OurMHRandomWalkTask,
                                                     OurCDataChannelSpace,
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

      // acquire GIL for logger and PyErr_Restore()
      py::gil_scoped_acquire gil_acquire;

      logger.debug("Python exception set.");

      pyerr.restorePyException();
      throw py::error_already_set();
      
    } catch (std::exception & e) {

      // acquire GIL for logger
      py::gil_scoped_acquire gil_acquire;

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






PYBIND11_MODULE(channelspace, m)
{
  m.doc() = "Python interface to run channel space random walk with diamond norm fig of merit";

  // make sure the tomographer module is loaded & initialized
  tpy::import_tomographer();

  pylogger = new tpy::PyLogger; // ownership will be taken over by Python/PyBind11

  pylogger->initPythonLogger("QPtomographer.channelspace");
  auto logger = Tomographer::Logger::makeLocalLogger(TOMO_ORIGIN, *pylogger);

  logger.debug("INIT QPTOMOGRAPHER.CHANNELSPACE");

  logger.debug("initializing global objects and constants, cxxlogger...");

  m.attr("cxxlogger") = pylogger; // ownership is transferred
  logger.debug("cxxlogger done");

  m.attr("cxx_tomographer_version") = std::string(TOMOGRAPHER_VERSION);

  m.attr("RandHermExp") = py::cast((int)RandHermExp);
  m.attr("ElemRotations") = py::cast((int)ElemRotations);

  logger.longdebug("now to define run() ...");

  // the main run call:
  m.def(
      "run", // function name
      &tomo_run_dnorm_channels, // fn pointer
      // no positional arguments, all kwargs
      // docstring
      "Docstring here"
      );
  
  logger.debug("DNormChannelSpaceInvalidInputError ...");
  
  py::register_exception<DNormChannelSpaceInvalidInputError>(m, "DNormChannelSpaceInvalidInputError");
}
