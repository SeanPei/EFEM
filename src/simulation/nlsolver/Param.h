// Copyright (C) 2016-2020 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef EFEM_PARAM_H
#define EFEM_PARAM_H

#include "io/Config.h"
#include <Eigen/Core>
#include <memory>
#include <stdexcept> // std::invalid_argument

namespace EFEM {

///
/// Parameters to control the L-BFGS algorithm.
///
class NLSolverParam {
private:
    typedef Eigen::SparseMatrix<double> SparseMatrix;

public:
    ///
    /// Absolute tolerance for convergence test.
    /// This parameter determines the absolute accuracy \f$\epsilon_{abs}\f$
    /// with which the solution is to be found. A minimization terminates when
    /// \f$||g|| < \max\{\epsilon_{abs}, \epsilon_{rel}||x||\}\f$,
    /// where \f$||\cdot||\f$ denotes the Euclidean (L2) norm. The default value is
    /// \c 1e-5.
    ///
    double epsilon;
    ///
    /// Relative tolerance for convergence test.
    /// This parameter determines the relative accuracy \f$\epsilon_{rel}\f$
    /// with which the solution is to be found. A minimization terminates when
    /// \f$||g|| < \max\{\epsilon_{abs}, \epsilon_{rel}||x||\}\f$,
    /// where \f$||\cdot||\f$ denotes the Euclidean (L2) norm. The default value is
    /// \c 1e-5.
    ///
    double epsilon_rel;
    ///
    /// The maximum number of iterations.
    /// The process is terminated when the iteration count
    /// exceeds this parameter. Setting this parameter to zero continues an
    /// process until a convergence or error. The default value
    /// is \c 0.
    ///
    int max_iterations;
    ///
    /// The maximum number of trials for the line search.
    /// This parameter controls the number of function and residuals evaluations
    /// per iteration for the line search routine. The default value is \c 20.
    ///
    int max_linesearch;
    ///
    /// The minimum step length allowed in the line search.
    /// The default value is \c 1e-20. Usually this value does not need to be
    /// modified.
    ///
    double min_step;
    ///
    /// The maximum step length allowed in the line search.
    /// The default value is \c 1e+20. Usually this value does not need to be
    /// modified.
    ///
    double max_step;
    ///
    /// Add benchmark data or not
    /// default is false
    ///
    bool benchmark;
    ///
    /// For LBFGS
    /// The number of corrections to approximate the inverse Jacobi matrix.
    /// The L-BFGS routine stores the computation results of previous \ref m
    /// iterations to approximate the inverse Jacobi matrix of the current
    /// iteration. This parameter controls the size of the limited memories
    /// (corrections). The default value is \c 6. Values less than \c 3 are
    /// not recommended. Large values will result in excessive computing time.
    ///
    int m;
    double lns_fea;
    ///
    /// The simulation config.
    ///
    std::shared_ptr<Config> config;

public:
    ///
    /// Constructor for L-BFGS parameters.
    /// Default values for parameters will be set when the object is created.
    ///
    NLSolverParam()
    {
        epsilon = double(1e-5);
        epsilon_rel = double(1e-5);
        max_iterations = 0;
        max_linesearch = 20;
        min_step = double(1e-20);
        max_step = double(1e+20);
        benchmark = false;
        lns_fea = double(1.0);
        m = 6;
    }

    ///
    /// Checking the validity of L-BFGS parameters.
    /// An `std::invalid_argument` exception will be thrown if some parameter
    /// is invalid.
    ///
    inline void check_param() const
    {
        if (epsilon < 0)
            throw std::invalid_argument("'epsilon' must be non-negative");
        if (epsilon_rel < 0)
            throw std::invalid_argument("'epsilon_rel' must be non-negative");
        if (max_iterations < 0)
            throw std::invalid_argument("'max_iterations' must be non-negative");
        if (max_linesearch <= 0)
            throw std::invalid_argument("'max_linesearch' must be positive");
        if (min_step < 0)
            throw std::invalid_argument("'min_step' must be positive");
        if (max_step < min_step)
            throw std::invalid_argument("'max_step' must be greater than 'min_step'");
        // for LBFGS
        if (m <= 0)
            throw std::invalid_argument("'m' must be positive");
    }
};

}

#endif // PARAM_H
