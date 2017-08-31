

#ifndef RELIABLEDIAMONDNORM_PYUTILS_H
#define RELIABLEDIAMONDNORM_PYUTILS_H




template<typename IterCountIntType>
struct CtrlConvergedParams
{
  CtrlConvergedParams()
    : enabled(false),
      check_frequency_sweeps(-1),
      max_allowed_unknown(-1),
      max_allowed_unknown_notisolated(-1),
      max_allowed_not_converged(-1),
      max_add_run_iters(std::numeric_limits<double>::quiet_NaN())
  {
  }

  static inline CtrlConvergedParams fromPyDict(py::dict ctrl_converged_params,
                                               tpy::HistogramParams histogram_params)
  {
    CtrlConvergedParams c;
    c.enabled = ctrl_converged_params.attr("get")("enabled", true).cast<bool>();
    c.check_frequency_sweeps =
      ctrl_converged_params.attr("get")("check_frequency_sweeps", 1024).cast<int>();
    c.max_allowed_unknown =
      ctrl_converged_params.attr("get")("max_allowed_unknown",
                                        1+2*histogram_params.num_bins/100).cast<Eigen::Index>();
    c.max_allowed_unknown_notisolated =
      ctrl_converged_params.attr("get")("max_allowed_unknown_notisolated",
                                        1+histogram_params.num_bins/100).cast<Eigen::Index>();
    c.max_allowed_not_converged =
      ctrl_converged_params.attr("get")("max_allowed_not_converged",
                                        1+histogram_params.num_bins/200).cast<Eigen::Index>();
    c.max_add_run_iters = ctrl_converged_params.attr("get")("max_add_run_iters", 1.5).cast<double>();
    return c;
  }

  static inline CtrlConvergedParams fromPyDictWithGilAcq(py::dict ctrl_converged_params,
                                                         tpy::HistogramParams histogram_params)
  {
    py::gil_scoped_acquire gilacq;
    return fromPyDict(std::move(ctrl_converged_params), std::move(histogram_params));
  }

  bool enabled;
  IterCountIntType check_frequency_sweeps;
  Eigen::Index max_allowed_unknown;
  Eigen::Index max_allowed_unknown_notisolated;
  Eigen::Index max_allowed_not_converged;
  double max_add_run_iters;
};



#endif

