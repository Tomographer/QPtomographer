

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


template<typename MatrixType, typename ExceptionClass>
MatrixType get_ref_channel(const int dimX, const int dimY, py::object ref_channel_XY)
{
  // return the reference channel as specified by the python argument
  // ref_channel_XY; handle the case = None by defaulting to the identity
  // channel (only for dimX == dimY)

  if (ref_channel_XY.is_none()) {

    // if None is given, then use the unnormalized maximally entangled state = Choi matrix
    // of the identity channel.  (|Phi> = \sum_k |kk>)

    if (dimX != dimY) {
      // but this can only be done if dimX==dimY
      throw ExceptionClass(
          "You must specify a reference channel ref_channel_XY if dimX != dimY"
          );
    }

    MatrixType mat_ref_channel_XY(dimX*dimX,dimX*dimX);
    mat_ref_channel_XY.setZero();
    
    for (int k = 0; k < dimX; ++k) {
      for (int k2 = 0; k2 < dimX; ++k2) {
        mat_ref_channel_XY(k + dimX*k, k2 + dimX*k2) = 1;
      }
    }

    std::cerr << "mat_ref_channel_XY = \n" << mat_ref_channel_XY << "\n";

    return mat_ref_channel_XY;
  }

  // extract an Eigen::Matrix from the python object
  return ref_channel_XY.cast<MatrixType>();
}





// provide this function for compatibility with tomographer 5.4
// (tpy::checkPyException() was introduced in 5.5)
void qptomographer_checkPyException()
{
  if (PyErr_Occurred() != NULL || PyErr_CheckSignals() == -1) {
    fprintf(stderr, "DEBUG:: Python exception set, throwing PyFetchedException\n") ;
    throw tpy::PyFetchedException();
  }
}


#endif

