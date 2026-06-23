#ifndef EFEM_CONSTRAINT_INERTIACONSTRAINT_H
#define EFEM_CONSTRAINT_INERTIACONSTRAINT_H

#include "fe/Constraint.h"

namespace EFEM {

class InertiaConstraint : public Constrainable {
public:
    InertiaConstraint(std::shared_ptr<Config> config, int dimens,
        std::vector<SparseMatrixTriplet>& consistent_inertia_matrix_triplets,
        std::vector<SparseMatrixTriplet>& lumped_inertia_matrix_triplets);
    void assignDx(Eigen::VectorXd dx);
    void conInitial();
    void conEvaluteResidual();
    void conGetResidual(Eigen::VectorXd& residual);
    void conEvaluteJacobi();
    void conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets);

private:
    Eigen::VectorXd m_con_residual;

    std::shared_ptr<Config> m_config;
    int m_dimens;
    std::vector<SparseMatrixTriplet> m_consistent_inertia_matrix_triplets;
    SparseMatrix m_consistent_inertia_matrix;
    std::vector<SparseMatrixTriplet> m_lumped_inertia_matrix_triplets;
    SparseMatrix m_lumped_inertia_matrix;
    Eigen::VectorXd m_dx;
    bool m_use_gravity;
    Eigen::VectorXd m_gravity;
    void buildGravity();
};
}

#endif