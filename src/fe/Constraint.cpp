#include "fe/Constraint.h"
#include "fe/ElementConstraint.h"
#include "fe/ExternalConstraint.h"
#include "fe/InertiaConstraint.h"
#include "fe/element/IntegralElement.h"
#include "io/Logger.h"
#include <boost/algorithm/string.hpp>

void EFEM::ConstraintElementInfo::buildIntegralElement()
{
    if (element->nodes.size() == 4) { // tet_linear
        integral_element = std::make_shared<TetElement3D4>(element);
    } else if (element->nodes.size() == 8) { // hex_linear
        integral_element = std::make_shared<HexElement3D8>(element);
    } else if (element->nodes.size() == 10) { // tet_quadratic
        integral_element = std::make_shared<TetElement3D10>(element);
    } else if (element->nodes.size() == 20) { // hex_quadratic
        integral_element = std::make_shared<HexElement3D20>(element);
    }
}

EFEM::Constraint::Constraint(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    m_config = config;
    m_model = model;

    std::vector<std::shared_ptr<Body>> bodies = m_model->getBodies();
    // build ConstraintNodeInfo
    std::vector<std::shared_ptr<Node>> body_nodes;
    for (int i = 0; i < bodies.size(); i++) {
        std::vector<std::shared_ptr<Node>> b_nodes = bodies[i]->getNodes();
        body_nodes.insert(body_nodes.end(), b_nodes.begin(), b_nodes.end());
    }
    for (int i = 0; i < body_nodes.size(); i++) {
        std::shared_ptr<ConstraintNodeInfo> info = std::make_shared<ConstraintNodeInfo>();
        info->node = body_nodes[i];
        info->position = body_nodes[i]->position + body_nodes[i]->dx;
        body_nodes[i]->extra_info.insert({ "constraint_info", info });
    }
    // build ConstraintElementInfo
    std::vector<std::shared_ptr<Element>> body_elements;
    for (int i = 0; i < bodies.size(); i++) {
        std::vector<std::shared_ptr<Element>> b_elements = bodies[i]->getElements();
        body_elements.insert(body_elements.end(), b_elements.begin(), b_elements.end());
    }
    for (int i = 0; i < body_elements.size(); i++) {
        std::shared_ptr<ConstraintElementInfo> info = std::make_shared<ConstraintElementInfo>();
        info->element = body_elements[i];
        info->buildIntegralElement();
        body_elements[i]->extra_info.insert({ "constraint_info", info });
    }

    // build Constraints
    m_constraint_configs.push_back(std::make_shared<ElementConstraintConfig>(m_config, m_model));
    m_constraint_configs.push_back(std::make_shared<ExternalConstraintConfig>(m_config, m_model));
    for (int i = 0; i < m_constraint_configs.size(); i++) {
        m_constraint_configs[i]->configConstraint(m_constraints);
    }

    // update ConstraintNodeInfo
    int id = 0;
    for (int i = 0; i < body_nodes.size(); i++) {
        std::shared_ptr<ConstraintNodeInfo> info = std::static_pointer_cast<ConstraintNodeInfo>(body_nodes[i]->extra_info["constraint_info"]);
        if (!info->fix) {
            info->id = id;
            m_infos.push_back(info);
            id++;
        } else {
            info->id = -1;
        }
    }
    m_dimens = id * 3;

    // initial constaint
    for (int i = 0; i < m_constraints.size(); i++) {
        m_constraints[i]->conInitial();
    }

    // build inertia matrix
    std::vector<SparseMatrixTriplet> consistent_intertia_tripplets;
    std::vector<SparseMatrixTriplet> lumped_intertia_tripplets;
    for (int i = 0; i < body_elements.size(); i++) {
        std::shared_ptr<ConstraintElementInfo> info = std::static_pointer_cast<ConstraintElementInfo>(body_elements[i]->extra_info["constraint_info"]);
        consistent_intertia_tripplets.insert(consistent_intertia_tripplets.end(), info->consistent_inertia_triplets.begin(), info->consistent_inertia_triplets.end());
        lumped_intertia_tripplets.insert(lumped_intertia_tripplets.end(), info->lumped_inertia_triplets.begin(), info->lumped_inertia_triplets.end());
    }
    // build inertia constraint
    m_inertia_constraint = std::make_shared<InertiaConstraint>(m_config, m_dimens, consistent_intertia_tripplets, lumped_intertia_tripplets);
    m_inertia_constraint->conInitial();
    m_constraints.push_back(m_inertia_constraint);
}

void EFEM::Constraint::assignDx(Eigen::VectorXd& dx)
{
    m_inertia_constraint->assignDx(dx);
    for (int i = 0; i < m_infos.size(); i++) {
        auto info = m_infos[i];
        info->dx = dx.segment(info->id * 3, 3);
        info->position = info->node->position + info->node->dx + info->dx;
    }
}

void EFEM::Constraint::updateDx()
{
    for (int i = 0; i < m_infos.size(); i++) {
        auto info = m_infos[i];
        auto node = info->node;
        Eigen::Vector3d info_dx = info->dx;
        node->dx += info_dx;
        info->dx.setZero();
        info->position = node->position + node->dx;
    }
}

void EFEM::Constraint::getResidual(Eigen::VectorXd& residual)
{
    for (int i = 0; i < m_constraints.size(); i++) {
        m_constraints[i]->conEvaluteResidual();
    }

    residual.resize(m_dimens);
    residual.setZero();
    for (int i = 0; i < m_constraints.size(); i++) {
        m_constraints[i]->conGetResidual(residual);
    }
}

void EFEM::Constraint::getJacobi(std::vector<SparseMatrixTriplet>& jacobi)
{
    for (int i = 0; i < m_constraints.size(); i++) {
        m_constraints[i]->conEvaluteJacobi();
    }

    jacobi.clear();
    for (int i = 0; i < m_constraints.size(); i++) {
        m_constraints[i]->conGetJacobi(jacobi);
    }
}