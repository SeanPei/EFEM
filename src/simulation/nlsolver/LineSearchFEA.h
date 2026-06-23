// Copyright (C) 2016-2020 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef EFEM_LINE_SEARCH_FEA_H
#define EFEM_LINE_SEARCH_FEA_H

#include "Param.h"
#include <Eigen/Core>
#include <stdexcept> // std::runtime_error

namespace EFEM {

///
/// The backtracking line search algorithm for L-BFGS. Mainly for internal use.
///
class LineSearchFEA {
private:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

public:
    ///
    /// Line search by backtracking.
    ///
    /// \param f      A function object such that `f.residual(x, residual)` overwrites `residual` with
    ///               the residual.
    /// \param fx     In: The objective function value at the current point.
    ///               Out: The function value at the new point.
    /// \param x      Out: The new point moved to.
    /// \param residual   In: The current residual vector. Out: The residual at the
    ///               new point.
    /// \param step   In: The initial step length. Out: The calculated step length.
    /// \param drt    The current moving direction.
    /// \param xp     The current point.
    /// \param param  Parameters for the LBFGS algorithm
    ///
    template <typename Foo>
    static void LineSearch(Foo& f, Vector& x, Vector& residual,
        double& step,
        const Vector& drt, const Vector& xp,
        const NLSolverParam& param)
    {
        // Projection of residual on the search direction
        const double dg_init = residual.dot(drt);
        double curr_step = step;
        int iter;
        for (int iter = 0; iter < param.max_linesearch; iter++) {
            // x_{k+1} = x_k + step * d_k
            x.noalias() = xp + curr_step * drt;
            // Evaluate this candidate
            f.residual(x, residual);
            double fx = 0.5 * residual.dot(residual);
            double dg = residual.dot(drt);
            if (abs(dg) <= param.lns_fea * abs(dg_init)) {
                break;
            }

            if (iter >= param.max_linesearch)
                throw std::runtime_error("the line search routine reached the maximum number of iterations");

            if (step < param.min_step)
                throw std::runtime_error("the line search step became smaller than the minimum value allowed");

            if (step > param.max_step)
                throw std::runtime_error("the line search step became larger than the maximum value allowed");

            curr_step *= dg * dg_init > 0 ? 1.5 : 0.5;
        }
    }
};

} // namespace LBFGS

#endif // LINE_SEARCH_BACKTRACKING_H
