#ifndef EFEM_CONSTRAINT_EXTERNALCONSTRAINT_H
#define EFEM_CONSTRAINT_EXTERNALCONSTRAINT_H

#include "fe/Constraint.h"
#include "fe/element/IntegralElement.h"

namespace EFEM {

enum ExternalConstraintType {
    FIX,
    FORCE,
    PRESSURE
};

// 配置外部约束
class ExternalConstraintConfigurable {
public:
    ExternalConstraintType external_constraint_type;
    virtual void configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints) = 0;
};

// 通用约束
class ExternalConstraint : public Constrainable {
public:
    ExternalConstraint()
    {
        constraint_type = ConstraintType::EXTERNL;
    }
    ExternalConstraintType external_constraint_type;
};

// 通用约束
class ExternalConstraintConfig : public ConstraintConfigurable {
public:
    ExternalConstraintConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void configConstraint(std::vector<std::shared_ptr<Constrainable>>& constaints);

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
    std::vector<std::shared_ptr<ExternalConstraintConfigurable>> m_constraint_configs;
};

class ExternalFixConfig : public ExternalConstraintConfigurable {
public:
    ExternalFixConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints);

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
};

class ExternalForceConfig : public ExternalConstraintConfigurable {
public:
    ExternalForceConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints);

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
};

class ExternalForceConstraint : public ExternalConstraint {
public:
    ExternalForceConstraint(std::shared_ptr<Config> config, Eigen::Vector3d force);
    void addNode(std::shared_ptr<Node> node);
    void conInitial();
    void conUpdate();
    void conReset();
    void conEvaluteResidual();
    void conGetResidual(Eigen::VectorXd& residual);
    void conEvaluteJacobi();
    void conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets);

private:
    std::unordered_map<EFEM::EFEMId, Eigen::Vector3d> m_con_residuals;

    std::shared_ptr<Config> m_config;
    std::vector<std::shared_ptr<Node>> m_nodes;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_node_infos;
    Eigen::Vector3d m_force;
};

class ExternalPressureConfig : public ExternalConstraintConfigurable {
public:
    ExternalPressureConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints);

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
};

class ExternalPressureFace {
public:
    ExternalPressureFace(std::vector<std::shared_ptr<Node>> face);
    void getNodeDArea(int index, Eigen::Vector3d& d_area);
    double getNodeArea(int index);
    void getDAreaJacobi(std::vector<SparseMatrixTriplet> jacobi_triplet, double pressure);
    std::vector<std::shared_ptr<Node>>& getNodes()
    {
        return m_nodes;
    }
    std::vector<std::shared_ptr<ConstraintNodeInfo>>& getNodeInfos()
    {
        return m_node_infos;
    }
    std::vector<std::shared_ptr<Node>>& getTotalNodes()
    {
        return m_total_nodes;
    }
    std::vector<std::shared_ptr<ConstraintNodeInfo>>& getTotalNodeInfos()
    {
        return m_total_node_infos;
    }

private:
    void initial();
    void buildGaussPoint();

    std::vector<std::shared_ptr<Node>> m_nodes;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_node_infos;
    std::vector<std::shared_ptr<Node>> m_total_nodes;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_total_node_infos;
    std::shared_ptr<Element> m_element;
    std::shared_ptr<IntegralElement> m_integral_element;
    Eigen::MatrixXd m_gauss_local_coords;
    Eigen::VectorXd m_gauss_weight;
    Eigen::MatrixXd m_node_weight;
    std::vector<Eigen::Matrix3d> m_gauss_jacobi0;
    std::vector<Eigen::Matrix3d> m_gauss_jacobi0_inv;
    std::vector<Eigen::Matrix<double, 3, -1>> m_gauss_T;
    Eigen::Vector3d m_A0;
};

class ExternalPressureConstraint : public ExternalConstraint {
public:
    ExternalPressureConstraint(std::shared_ptr<Config> config, double value);
    ExternalPressureConstraint(std::shared_ptr<Config> config, Eigen::Vector3d pressure);
    void addPressureFace(std::shared_ptr<ExternalPressureFace> face);
    void conInitial();
    void conEvaluteResidual();
    void conGetResidual(Eigen::VectorXd& residual);
    void conEvaluteJacobi();
    void conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets);

private:
    std::unordered_map<EFEM::EFEMId, Eigen::Vector3d> m_con_residuals;
    std::vector<SparseMatrixTriplet> m_con_jacobi_triplets;

    std::shared_ptr<Config> m_config;
    bool m_direction;
    double m_pressure_value;
    Eigen::Vector3d m_pressure_dvalue;
    std::vector<std::shared_ptr<ExternalPressureFace>> m_pressure_faces;
};
}

#endif