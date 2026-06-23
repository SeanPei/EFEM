#include "fe/ElementConstraint.h"
#include "io/Logger.h"

EFEM::ElementConstraintConfig::ElementConstraintConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    constraint_type = ConstraintType::ELEMENT;
    m_config = config;
    m_model = model;
}

void EFEM::ElementConstraintConfig::configConstraint(std::vector<std::shared_ptr<Constrainable>>& constaints)
{
    loadMaterials();
    std::vector<std::shared_ptr<Element>> body_elements;
    std::vector<std::shared_ptr<Body>> bodies = m_model->getBodies();
    for (int i = 0; i < bodies.size(); i++) {
        std::vector<std::shared_ptr<Element>> b_element = bodies[i]->getElements();
        body_elements.insert(body_elements.end(), b_element.begin(), b_element.end());
    }
    for (int i = 0; i < body_elements.size(); i++) {
        std::shared_ptr<Element> element = body_elements[i];
        std::shared_ptr<Material> material = m_materials[element->material_id];
        if (!material) {
            Logger::error("Undefined material(id={}) in element(index={})", element->material_id, i);
            continue;
        }
        std::shared_ptr<ElementConstraint> element_constraint = std::make_shared<ElementConstraint>(m_config, element, material);
        constaints.push_back(element_constraint);
        double mass = element->volume() * material->density;
        if (mass < 4e-12) {
            mass = 4e-12;
        }
        for (auto node : element->nodes) {
            std::static_pointer_cast<ConstraintNodeInfo>(node->extra_info["constraint_info"])->mass += mass / element->nodes.size();
        }
    }
}

void EFEM::ElementConstraintConfig::loadMaterials()
{
    std::shared_ptr<Json> config_json = m_config->io()->configJson();
    if (!(config_json->find("material") != config_json->end() && config_json->find("material")->is_array())) {
        return;
    }
    std::vector<Json> material_jsons = config_json->find("material")->template get<std::vector<Json>>();

    for (auto mat_json : material_jsons) {
        if (mat_json["type"].is_string()) {
            std::string type = mat_json["type"].template get<std::string>();
            if (type == "neo_hookean") {
                auto material = std::make_shared<NeoHookean>();
                material->initFromJson(mat_json);
                m_materials.insert({ material->id, material });
            } else if (type == "mo_riv") {
                auto material = std::make_shared<MooneyRivlin>();
                material->initFromJson(mat_json);
                m_materials.insert({ material->id, material });
            } else if (type == "yeoh") {
                auto material = std::make_shared<Yeoh>();
                material->initFromJson(mat_json);
                m_materials.insert({ material->id, material });
            } else {
                Logger::info("unsupport material type: type={}", type);
            }
        }
    }
}

EFEM::ElementConstraint::ElementConstraint(std::shared_ptr<Config> config, std::shared_ptr<Element> element, std::shared_ptr<Material> material)
{
    constraint_type = ConstraintType::ELEMENT;
    m_element = element;
    m_config = config;
    m_material = material;
    m_volume = m_element->volume();
    for (int i = 0; i < m_element->nodes.size(); i++) {
        auto node = m_element->nodes[i];
        m_infos.push_back(std::static_pointer_cast<ConstraintNodeInfo>(node->extra_info["constraint_info"]));
    }

    auto element_info = std::static_pointer_cast<ConstraintElementInfo>(m_element->extra_info["constraint_info"]);
    m_integral_element = element_info->integral_element;
}

void EFEM::ElementConstraint::conInitial()
{
    auto element_info = std::static_pointer_cast<ConstraintElementInfo>(m_element->extra_info["constraint_info"]);

    Eigen::MatrixXd consistent_inertia_matrix;
    consistent_inertia_matrix.resize(m_infos.size() * 3, m_infos.size() * 3);
    consistent_inertia_matrix.setZero();
    int gp_num = m_integral_element->getGaussPointNum();
    for (int index = 0; index < gp_num; index++) {
        Eigen::VectorXd shape_funcs = m_integral_element->getGaussPointShapeFunc(index);
        Eigen::MatrixXd shape_matrix;
        shape_matrix.resize(3, 3 * shape_funcs.size());
        for (int i = 0; i < shape_funcs.size(); i++) {
            shape_matrix.block<3, 3>(0, i * 3) = shape_funcs[i] * Eigen::Matrix3d::Identity();
        }
        Eigen::MatrixXd gp_inertia_matrix = shape_matrix.transpose() * shape_matrix;
        consistent_inertia_matrix += m_integral_element->getGaussPointWeight(index) * m_material->density * gp_inertia_matrix;
    }
    consistent_inertia_matrix *= element_info->element->volume();

    Eigen::MatrixXd lumped_inertia_matrix;
    lumped_inertia_matrix.resize(m_infos.size() * 3, m_infos.size() * 3);
    lumped_inertia_matrix.setZero();
    double mass = element_info->element->volume() * m_material->density;
    if (mass < 4e-12) {
        mass = 4e-12;
    }
    mass = mass / element_info->element->nodes.size();
    for (int i = 0; i < m_infos.size(); i++) {
        lumped_inertia_matrix.block<3, 3>(i * 3, i * 3) = mass * Eigen::Matrix3d::Identity();
    }

    element_info->consistent_inertia_triplets.clear();
    element_info->lumped_inertia_triplets.clear();
    for (unsigned int i = 0; i < m_infos.size(); i++) {
        for (unsigned int j = 0; j < m_infos.size(); j++) {
            if (m_infos[i]->id >= 0 && m_infos[j]->id >= 0) {
                for (unsigned int row = 0; row < 3; row++) {
                    for (unsigned int col = 0; col < 3; col++) {
                        double consistent_temp = consistent_inertia_matrix(i * 3 + row, j * 3 + col);
                        if (std::abs(consistent_temp) > 0) {
                            element_info->consistent_inertia_triplets.push_back(SparseMatrixTriplet(m_infos[i]->id * 3 + row, m_infos[j]->id * 3 + col, consistent_temp));
                        }
                        double lmuped_temp = lumped_inertia_matrix(i * 3 + row, j * 3 + col);
                        if (std::abs(lmuped_temp) > 0) {
                            element_info->lumped_inertia_triplets.push_back(SparseMatrixTriplet(m_infos[i]->id * 3 + row, m_infos[j]->id * 3 + col, lmuped_temp));
                        }
                    }
                }
            }
        }
    }
}

void EFEM::ElementConstraint::conEvaluteResidual()
{
    int node_num = m_infos.size();
    m_con_residuals.resize(3, node_num);
    m_con_residuals.setZero();

    int gp_num = m_integral_element->getGaussPointNum();
    for (int index = 0; index < gp_num; index++) {
        Eigen::Matrix3d gp_F;
        double weight = this->m_volume * m_integral_element->getGaussPointWeight(index);
        gp_F = m_integral_element->getGaussPointF(index);
        m_con_residuals += weight * this->m_material->getFirstPiola(gp_F) * m_integral_element->getGaussPointT(index);
    }
}

void EFEM::ElementConstraint::conGetResidual(Eigen::VectorXd& residual)
{
    for (int i = 0; i < m_infos.size(); i++) {
        auto info = m_infos[i];
        if (info->id >= 0) {
            residual.segment(info->id * 3, 3) += m_con_residuals.block<3, 1>(0, i);
        }
    }
}

void EFEM::ElementConstraint::conEvaluteJacobi()
{
    int node_num = m_infos.size();
    Eigen::MatrixXd gp_H;
    gp_H.resize(node_num * 3, node_num * 3);
    gp_H.setZero();
    double ks = this->m_volume * m_material->elasticity_coeff;
    int gp_num = m_integral_element->getGaussPointNum();
    for (int index = 0; index < gp_num; index++) {
        Matrix3333 dPdF;
        dPdF.SetIdentity();
        double weight = ks * m_integral_element->getGaussPointWeight(index);
        Eigen::Matrix<double, 3, -1> gp_F_gradient = m_integral_element->getGaussPointT(index);
        Matrix3333 gp_dPdF = weight * dPdF;
        for (int k = 0; k < node_num; k++) {
            for (int i = 0; i < 3; i++) {
                Eigen::Matrix<double, 3, -1> gp_H_part;
                gp_H_part.resize(3, node_num);
                gp_H_part.setZero();
                for (int j = 0; j < 3; j++) {
                    gp_H_part += gp_F_gradient(j, k) * gp_dPdF(i, j) * gp_F_gradient;
                }
                for (int j = 0; j < node_num; j++) {
                    gp_H.block<3, 1>(j * 3, k * 3 + i) += gp_H_part.block<3, 1>(0, j);
                }
            }
        }
    }

    m_con_jacobi_triplets.clear();
    for (unsigned int i = 0; i < m_infos.size(); i++) {
        for (unsigned int j = 0; j < m_infos.size(); j++) {
            if (m_infos[i]->id >= 0 && m_infos[j]->id >= 0) {
                for (unsigned int row = 0; row < 3; row++) {
                    for (unsigned int col = 0; col < 3; col++) {
                        m_con_jacobi_triplets.push_back(SparseMatrixTriplet(m_infos[i]->id * 3 + row, m_infos[j]->id * 3 + col, gp_H(i * 3 + row, j * 3 + col)));
                    }
                }
            }
        }
    }
}

void EFEM::ElementConstraint::conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets)
{
    jacobi_triplets.insert(jacobi_triplets.end(), m_con_jacobi_triplets.begin(), m_con_jacobi_triplets.end());
}