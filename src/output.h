#ifndef OUTPUT_H
#define OUTPUT_H

// ======== R via Rcpp ========
#ifdef USE_RCPP
  #include <Rcpp.h>
  #define COUT Rcpp::Rcout
  #define CERR Rcpp::Rcerr
  #define STOP(msg) Rcpp::stop(msg)
  #define WARN(msg) Rcpp::Rcerr << "Warning: " << msg << std::endl

// ======== Python via pybind11 ========
#elif defined(USE_PYBIND11)
  #include <iostream>
  #include <stdexcept>
  #include <pybind11/pybind11.h>
  namespace py = pybind11;
  #define COUT std::cout
  #define CERR std::cerr
  #define STOP(msg) throw std::runtime_error(msg)
  #define WARN(msg) do { \
    std::cerr << "Warning: " << msg << std::endl; \
    py::module::import("warnings").attr("warn")(msg); \
  } while(0)

// ======== Plain C++ ========
#else
  #include <iostream>
  #include <stdexcept>
  #define COUT std::cout
  #define CERR std::cerr
  #define STOP(msg) throw std::runtime_error(msg)
  #define WARN(msg) std::cerr << "Warning: " << msg << std::endl
#endif

#endif // OUTPUT_H
