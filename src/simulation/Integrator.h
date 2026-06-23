#ifndef EFEM_SIMULATION_INTEGRATOR_H
#define EFEM_SIMULATION_INTEGRATOR_H

#include "fe/Constraint.h"
#include "io/Config.h"
#include "math/Math.h"
#include "model/Model.h"
#include "nlsolver/LBFGS.h"
#include <math/MathType.h>
#include <memory>
#include <vector>

namespace EFEM {

class SoftObjective {
public:
    int dimens;
    int r_evalutes = 0;
    int j_evalutes = 0;
    Eigen::VectorXd nl_residual;
    std::shared_ptr<Constraint> constraint;

    void residual(Eigen::VectorXd& dx, Eigen::VectorXd& residual);
    void jacobi(Eigen::VectorXd& dx, SparseMatrix& jacobi);

private:
    friend class Integrator;
    std::shared_ptr<Config> m_config;
    void buildJacobi(SparseMatrix& jacobi);
};

class Integrator {
public:
    Integrator(std::shared_ptr<Config> config, std::shared_ptr<Model> model, std::shared_ptr<Constraint> constraint);
    bool integral();
    int dimens() { return m_soft_objective.dimens; }
    void printBenchmarch();

private:
    Math m_math;
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
    SoftObjective m_soft_objective;
    std::shared_ptr<LBFGSSolver<SoftObjective>> m_solver;
};
}

#endif