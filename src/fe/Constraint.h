#ifndef EFEM_CONSTRAINT_CONSTRAINT_H
#define EFEM_CONSTRAINT_CONSTRAINT_H

#include "Eigen/Dense"
#include "core/Object.h"
#include "io/Config.h"
#include "math/MathType.h"
#include "model/Model.h"
#include <memory>
#include <vector>

namespace EFEM {

enum ConstraintType {
    ELEMENT,
    INERTIA,
    EXTERNL
};

class IntegralElement;
class InertiaConstraint;

class Constrainable : public Object {
public:
    ConstraintType constraint_type;
    virtual void conInitial() = 0;
    virtual void conEvaluteResidual() = 0;
    virtual void conGetResidual(Eigen::VectorXd& residual) = 0;
    virtual void conEvaluteJacobi() = 0;
    virtual void conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi) = 0;
};

class ConstraintConfigurable : public Object {
public:
    ConstraintType constraint_type;
    virtual void configConstraint(std::vector<std::shared_ptr<Constrainable>>& constaints)
        = 0;
};

class ConstraintNodeInfo : public Object {
public:
    int id = -1;
    bool fix = false;
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d dx = Eigen::Vector3d::Zero();
    double mass = 0;
    std::shared_ptr<Node> node;
};

class ConstraintElementInfo : public Object {
public:
    void buildIntegralElement();
    std::shared_ptr<IntegralElement> integral_element;
    std::shared_ptr<Element> element;
    std::vector<SparseMatrixTriplet> consistent_inertia_triplets;
    std::vector<SparseMatrixTriplet> lumped_inertia_triplets;
};

class Constraint {
public:
    Constraint(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void assignDx(Eigen::VectorXd& dx);
    void updateDx();
    void getResidual(Eigen::VectorXd& residual);
    void getJacobi(std::vector<SparseMatrixTriplet>& jacobi);
    int getDimens() { return m_dimens; };

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
    int m_dimens;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_infos;
    std::vector<std::shared_ptr<ConstraintConfigurable>> m_constraint_configs;
    std::vector<std::shared_ptr<Constrainable>> m_constraints;
    std::shared_ptr<InertiaConstraint> m_inertia_constraint;
};
}

#endif