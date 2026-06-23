// Copyright (C) 2016-2020 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef EFEM_LBFGS_H
#define EFEM_LBFGS_H

#include "LBFGSMat.h"
#include "LineSearchFEA.h"
#include "Param.h"
#include <Eigen/Core>
#include <time.h>

namespace EFEM {

template <typename Foo>
class LBFGSSolver {
private:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Map<Vector> MapVec;
    typedef Eigen::SparseMatrix<double> SparseMatrix;
    typedef Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper> SparseMatrixSolver;

    const NLSolverParam m_param; // Parameters to control the LBFGS algorithm
    LBFGSMat m_bfgs; // Approximation to the Jacobi matrix
    Vector m_xp; // Old x
    Vector m_residual; // New residual
    Vector m_residualp; // Old residual
    Vector m_drt; // Moving direction
    Vector m_x_Initial;

    bool m_jacobi_available = false; // Jacobi available
    std::shared_ptr<SparseMatrixSolver> m_jacobi_solver_ptr; // Jacobi inverse solver

    std::vector<clock_t> m_time_benchmark;
    std::vector<double> m_f_benchmark;

    // Reset internal variables
    // n: dimension of the vector
    inline void reset(int n)
    {
        const int m = m_param.m;
        m_drt.resize(n);
        m_bfgs.reset(n, m);
        m_xp.resize(n);
        m_residual.resize(n);
        m_residualp.resize(n);
        m_time_benchmark.clear();
        m_f_benchmark.clear();
    }

public:
    ///
    /// Constructor for the L-BFGS solver.
    ///
    /// \param param An object of \ref Param to store parameters for the
    ///        algorithm
    ///
    LBFGSSolver(const NLSolverParam& param)
        : m_param(param)
    {
        m_param.check_param();
    }

    inline void updateApproxJacobi(Foo& f, Vector& x)
    {
        SparseMatrix jacobi;
        f.jacobi(x, jacobi);
        m_jacobi_available = false;
        m_jacobi_solver_ptr = std::make_shared<SparseMatrixSolver>();
        if (jacobi.nonZeros() > 0) {
            m_jacobi_solver_ptr->analyzePattern(jacobi);
            m_jacobi_solver_ptr->factorize(jacobi);
            m_jacobi_available = m_jacobi_solver_ptr->info() == Eigen::Success;
        }
        if (!m_jacobi_available) {
            Logger::warn("L_BFGS updateApproxJacobi: jacobi solve failed, info={}", m_jacobi_solver_ptr->info());
            m_jacobi_solver_ptr.reset();
        }
        m_bfgs.setApproxJacobiSolverPtr(m_jacobi_solver_ptr);
    }

    ///
    /// Minimizing a multivariate function using the L-BFGS algorithm.
    /// Exceptions will be thrown if error occurs.
    ///
    /// \param f  A function object such that `f(x, residual)` returns the
    ///           objective function value at `x`, and overwrites `residual` with
    ///           the residual.
    /// \param x  In: An initial guess. Out: The best point
    ///           found.
    /// \param fx Out: The objective function value at `x`.
    ///
    /// \return Number of iterations used.
    ///
    inline int solve(Foo& f, Vector& x)
    {
        using std::abs;

        // backup x
        m_x_Initial.noalias() = x;

        // Dimension of the vector
        const int n = x.size();
        reset(n);

        if (m_param.benchmark) {
            m_time_benchmark.push_back(clock());
        }
        // Evaluate function and compute residual
        f.residual(x, m_residual);
        if (m_param.benchmark) {
            m_f_benchmark.push_back(m_residual.norm());
        }
        double gnorm = m_residual.template lpNorm<Eigen::Infinity>();

        // Early exit if the initial x is already a minimizer
        if (gnorm <= m_param.epsilon || gnorm <= m_param.epsilon_rel * x.norm()) {
            return 1;
        }

        if (!m_jacobi_available && m_param.config->analysisConfig()->nl_solver == AnalysisConfig::NLSolver::LBFGS_AH) {
            updateApproxJacobi(f, x);
        }

        // Initial step size
        double step = double(1);
        // Initial direction
        if (m_jacobi_available) {
            m_drt.noalias() = -m_jacobi_solver_ptr->solve(m_residual);
        } else {
            m_drt.noalias() = -m_residual;
        }

        // Number of iterations used
        int k = 1;
        for (;;) {
            // Save the curent x and residual
            m_xp.noalias() = x;
            m_residualp.noalias() = m_residual;

            // Line search to update x and residual
            LineSearchFEA::LineSearch(f, x, m_residual, step, m_drt, m_xp, m_param);
            if (m_param.benchmark) {
                m_time_benchmark.push_back(clock());
                m_f_benchmark.push_back(m_residual.norm());
            }

            // New residual norm
            gnorm = m_residual.template lpNorm<Eigen::Infinity>();

            // Convergence test -- residual
            if (gnorm <= m_param.epsilon || gnorm <= m_param.epsilon_rel * x.template lpNorm<Eigen::Infinity>()) {
                return k;
            }
            // Maximum number of iterations
            if (m_param.max_iterations != 0 && k >= m_param.max_iterations) {
                return k;
            }

            // Update s and y
            // s_{k+1} = x_{k+1} - x_k
            // y_{k+1} = g_{k+1} - g_k
            m_bfgs.add_correction(x - m_xp, m_residual - m_residualp);
            // Recursive formula to compute d = -H * g
            m_bfgs.apply_Hv(m_residual, -double(1), m_drt);

            // Reset step = 1.0 as initial guess for the next line search
            step = double(1);
            k++;
        }

        return k;
    }

    inline void printBenchmark()
    {
        int length = m_f_benchmark.size();
        if (length == 0) {
            spdlog::info("no benchmark!");
            return;
        }
        double f_benchmark = m_f_benchmark[0] - m_f_benchmark[length - 1];
        for (int i = 0; i < length && i < 9000; i++) {
            spdlog::info("benchmark[{}]: f={}, time={}s", i,
                (m_f_benchmark[i] - m_f_benchmark[length - 1]) / f_benchmark,
                (double)(m_time_benchmark[i] - m_time_benchmark[0]) / CLOCKS_PER_SEC);
        }
    }
};

} // namespace LBFGS

#endif // LBFGS_H
