#ifndef EFEM_CONSTRAINT_ELEMENTCONSTRAINT_H
#define EFEM_CONSTRAINT_ELEMENTCONSTRAINT_H

#include "fe/Constraint.h"
#include "fe/element/IntegralElement.h"
#include "material/Material.h"
#include <vector>

namespace EFEM {

class ElementConstraintConfig : public ConstraintConfigurable {
public:
    ElementConstraintConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model);
    void configConstraint(std::vector<std::shared_ptr<Constrainable>>& constaints);

private:
    std::shared_ptr<Config> m_config;
    std::shared_ptr<Model> m_model;
    std::unordered_map<int, std::shared_ptr<Material>> m_materials;

    void loadMaterials();
};

class ElementConstraint : public Constrainable {
public:
    ElementConstraint(std::shared_ptr<Config> config, std::shared_ptr<Element> element, std::shared_ptr<Material> material);
    void conInitial();
    void conEvaluteResidual();
    void conGetResidual(Eigen::VectorXd& residual);
    void conEvaluteJacobi();
    void conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets);

private:
    Eigen::Matrix<double, 3, -1> m_con_residuals;
    std::vector<SparseMatrixTriplet> m_con_jacobi_triplets;

    std::shared_ptr<Config> m_config;
    std::shared_ptr<Element> m_element;
    std::vector<std::shared_ptr<ConstraintNodeInfo>> m_infos;
    std::shared_ptr<Material> m_material;
    double m_volume = 0;
    std::shared_ptr<IntegralElement> m_integral_element;
};
}

#endif