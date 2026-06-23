#include "InertiaConstraint.h"

EFEM::InertiaConstraint::InertiaConstraint(std::shared_ptr<Config> config, int dimens,
    std::vector<SparseMatrixTriplet>& inertia_consistent_matrix_triplets,
    std::vector<SparseMatrixTriplet>& lumped_inertia_matrix_triplets)
{
    constraint_type = ConstraintType::INERTIA;
    m_config = config;
    m_dimens = dimens;
    m_consistent_inertia_matrix_triplets = inertia_consistent_matrix_triplets;
    m_lumped_inertia_matrix_triplets = lumped_inertia_matrix_triplets;
}

void EFEM::InertiaConstraint::assignDx(Eigen::VectorXd dx)
{
    m_dx = dx;
}

void EFEM::InertiaConstraint::conInitial()
{
    m_dx.resize(m_dimens);
    m_dx.setZero();
    m_con_residual.resize(m_dimens);
    m_con_residual.setZero();

    m_consistent_inertia_matrix.resize(m_dimens, m_dimens);
    m_consistent_inertia_matrix.setFromTriplets(m_consistent_inertia_matrix_triplets.begin(), m_consistent_inertia_matrix_triplets.end());
    m_lumped_inertia_matrix.resize(m_dimens, m_dimens);
    m_lumped_inertia_matrix.setFromTriplets(m_lumped_inertia_matrix_triplets.begin(), m_lumped_inertia_matrix_triplets.end());

    buildGravity();
}

void EFEM::InertiaConstraint::conEvaluteResidual()
{
    m_con_residual.setZero();
    // gravity
    if (m_use_gravity) {
        Eigen::VectorXd gravity_residual;
        if (m_config->analysisConfig()->consistent_inertia) {
            gravity_residual = -m_consistent_inertia_matrix * m_gravity;
        } else {
            gravity_residual = -m_lumped_inertia_matrix * m_gravity;
        }
        m_con_residual += gravity_residual;
    }
}

void EFEM::InertiaConstraint::conGetResidual(Eigen::VectorXd& residual)
{
    residual += m_con_residual;
}

void EFEM::InertiaConstraint::conEvaluteJacobi()
{
}

void EFEM::InertiaConstraint::conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets)
{
}

void EFEM::InertiaConstraint::buildGravity()
{
    m_use_gravity = m_config->analysisConfig()->gravity > AnalysisConfig::Gravity::NON;
    Eigen::Vector3d gravity;
    if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::X) {
        gravity = Eigen::Vector3d(9.8, 0.0, 0);
    } else if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::NX) {
        gravity = Eigen::Vector3d(-9.8, 0.0, 0);
    } else if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::Y) {
        gravity = Eigen::Vector3d(0, 9.8, 0);
    } else if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::NY) {
        gravity = Eigen::Vector3d(0, -9.8, 0);
    } else if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::Z) {
        gravity = Eigen::Vector3d(0, 0, 9.8);
    } else if (m_use_gravity = m_config->analysisConfig()->gravity == AnalysisConfig::Gravity::NZ) {
        gravity = Eigen::Vector3d(0, 0, -9.8);
    }
    m_gravity.resize(m_dimens);
    m_gravity.setZero();
    for (int i = 0; i < m_dimens / 3; i++) {
        m_gravity.block<3, 1>(i * 3, 0) = gravity;
    }
}
