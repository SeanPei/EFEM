// Copyright (C) 2020 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef EFEM_BFGS_MAT_H
#define EFEM_BFGS_MAT_H

#include "Param.h"
#include "io/Logger.h"
#include "math/MathType.h"
#include <Eigen/Core>
#include <vector>

/// \cond

namespace EFEM {

//
// An *implicit* representation of the BFGS approximation to the Jacobi matrix B
//
// B = theta * I - W * M * W'
// H = inv(B)
//
// Reference:
// [1] D. C. Liu and J. Nocedal (1989). On the limited memory BFGS method for large scale optimization.
// [2] R. H. Byrd, P. Lu, and J. Nocedal (1995). A limited memory algorithm for bound constrained optimization.
//
class LBFGSMat {
private:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Ref<const Vector> RefConstVec;
    typedef std::vector<int> IndexSet;
    typedef Eigen::SparseMatrix<double> SparseMatrix;
    typedef Eigen::SimplicialLLT<SparseMatrix, Eigen::Upper> SparseMatrixSolver;

    int m_m; // Maximum number of correction vectors
    double m_theta; // theta * I is the initial approximation to the Jacobi matrix
    Matrix m_s; // History of the s vectors
    Matrix m_y; // History of the y vectors
    Vector m_ys; // History of the s'y values
    Vector m_alpha; // Temporary values used in computing H * v
    int m_ncorr; // Number of correction vectors in the history, m_ncorr <= m
    int m_ptr; // A Pointer to locate the most recent history, 1 <= m_ptr <= m
        // Details: s and y vectors are stored in cyclic order.
        //          For example, if the current s-vector is stored in m_s[, m-1],
        //          then in the next iteration m_s[, 0] will be overwritten.
        //          m_s[, m_ptr-1] points to the most recent history,
        //          and m_s[, m_ptr % m] points to the most distant one.
    bool m_jacobi_available = false; // Jacobi available
    std::shared_ptr<SparseMatrixSolver> m_jacobi_solver_ptr; // Jacobi inverse solver

public:
    // Constructor
    LBFGSMat() {}

    // Reset internal variables
    // n: dimension of the vector
    // m: maximum number of corrections to approximate the Jacobi matrix
    inline void reset(int n, int m)
    {
        m_m = m;
        m_theta = double(1);
        m_s.resize(n, m);
        m_y.resize(n, m);
        m_ys.resize(m);
        m_alpha.resize(m);
        m_ncorr = 0;
        m_ptr = m; // This makes sure that m_ptr % m == 0 in the first step
    }

    inline void setApproxJacobiSolverPtr(std::shared_ptr<SparseMatrixSolver> jacobi_solver_ptr)
    {
        m_jacobi_available = false;
        m_jacobi_solver_ptr = jacobi_solver_ptr;
        if (jacobi_solver_ptr) {
            m_jacobi_available = true;
        }
    }

    // Add correction vectors to the BFGS matrix
    inline void add_correction(const RefConstVec& s, const RefConstVec& y)
    {
        const int loc = m_ptr % m_m;

        m_s.col(loc).noalias() = s;
        m_y.col(loc).noalias() = y;

        // ys = y's = 1/rho
        const double ys = m_s.col(loc).dot(m_y.col(loc));
        m_ys[loc] = ys;

        m_theta = m_y.col(loc).squaredNorm() / ys;

        if (m_ncorr < m_m)
            m_ncorr++;

        m_ptr = loc + 1;
    }

    // Recursive formula to compute a * H * v, where a is a double, and v is [n x 1]
    // If m_jacobi_available = false, H0 = (1/theta) * I else H0 = m_jacobi^(-1) is the initial approximation to H
    // Algorithm 7.4 of Nocedal, J., & Wright, S. (2006). Numerical optimization.
    inline void apply_Hv(const Vector& v, const double& a, Vector& res)
    {
        res.resize(v.size());

        // L-BFGS two-loop recursion

        // Loop 1
        res.noalias() = a * v;
        int j = m_ptr % m_m;
        for (int i = 0; i < m_ncorr; i++) {
            j = (j + m_m - 1) % m_m;
            m_alpha[j] = m_s.col(j).dot(res) / m_ys[j];
            res.noalias() -= m_alpha[j] * m_y.col(j);
        }

        // Apply initial H0
        if (m_jacobi_available) {
            res = m_jacobi_solver_ptr->solve(res);
        } else {
            res /= m_theta;
        }

        // Loop 2
        for (int i = 0; i < m_ncorr; i++) {
            const double beta = m_y.col(j).dot(res) / m_ys[j];
            res.noalias() += (m_alpha[j] - beta) * m_s.col(j);
            j = (j + 1) % m_m;
        }
    }
};

} // namespace LBFGS

/// \endcond

#endif // BFGS_MAT_H
