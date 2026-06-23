#include "simulation/Integrator.h"
#include "io/Config.h"
#include "io/Logger.h"
#include "nlsolver/LBFGS.h"

EFEM::Integrator::Integrator(std::shared_ptr<Config> config, std::shared_ptr<Model> model, std::shared_ptr<Constraint> constraint)
{
    m_config = config;
    m_model = model;
    m_soft_objective.m_config = m_config;
    m_soft_objective.constraint = constraint;
    m_soft_objective.dimens = constraint->getDimens();
    NLSolverParam param;
    param.epsilon = m_config->analysisConfig()->residual_criterion;
    param.epsilon_rel = m_config->analysisConfig()->residual_criterion;
    param.max_iterations = m_config->analysisConfig()->max_solver_iters;
    param.max_linesearch = m_config->analysisConfig()->max_line_search;
    param.config = m_config;
    m_solver = std::make_shared<LBFGSSolver<SoftObjective>>(param);
}

bool EFEM::Integrator::integral()
{
    Eigen::VectorXd dx;
    dx.resize(m_soft_objective.dimens);
    dx.fill(0.0);
    int niter = m_solver->solve(m_soft_objective, dx);
    double dx_lp_norm = dx.template lpNorm<Eigen::Infinity>();
    double residual_lp_norm = m_soft_objective.nl_residual.template lpNorm<Eigen::Infinity>();
    if (!m_math.isReal(dx_lp_norm) || !m_math.isReal(residual_lp_norm)) {
        Logger::error("LBFGS solver error!");
        dx.setZero();
        return false;
    }

    bool converge = dx_lp_norm < m_config->analysisConfig()->dx_criterion;
    if (m_config->outputConfig()->show_log) {
        Logger::info("dx={}, residual={}, niter={}, r_evalutes={}, j_evalutes={}", dx_lp_norm, residual_lp_norm, niter,
            m_soft_objective.r_evalutes, m_soft_objective.j_evalutes);
    }
    return converge;
}

void EFEM::Integrator::printBenchmarch()
{
    m_solver->printBenchmark();
}

void EFEM::SoftObjective::residual(Eigen::VectorXd& dx, Eigen::VectorXd& residual)
{
    constraint->assignDx(dx);
    constraint->getResidual(residual);
    this->nl_residual = residual;
    this->r_evalutes++;
}

void EFEM::SoftObjective::jacobi(Eigen::VectorXd& dx, SparseMatrix& jacobi)
{
    constraint->assignDx(dx);
    buildJacobi(jacobi);
    this->j_evalutes++;
}

void EFEM::SoftObjective::buildJacobi(SparseMatrix& jacobi)
{
    jacobi.resize(dimens, dimens);
    std::vector<SparseMatrixTriplet> jacobi_triplets;
    constraint->getJacobi(jacobi_triplets);
    jacobi.setFromTriplets(jacobi_triplets.begin(), jacobi_triplets.end());
}