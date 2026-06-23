#include "fe/ExternalConstraint.h"
#include "io/Logger.h"
#include "math/Math.h"
#include <boost/algorithm/string.hpp>

EFEM::ExternalConstraintConfig::ExternalConstraintConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    constraint_type = ConstraintType::ELEMENT;
    m_config = config;
    m_model = model;
}

void EFEM::ExternalConstraintConfig::configConstraint(std::vector<std::shared_ptr<Constrainable>>& constaints)
{
    std::shared_ptr<Json> config_json = m_config->io()->configJson();
    if (!(config_json->find("constraint") != config_json->end() && config_json->find("constraint")->is_array())) {
        return;
    }
    std::vector<Json> constraint_jsons = config_json->find("constraint")->template get<std::vector<Json>>();
    m_constraint_configs.clear();
    for (auto con_json : constraint_jsons) {
        if (con_json["type"].is_string()) {
            std::string type = con_json["type"].template get<std::string>();
            if (type == "fix") {
                std::shared_ptr<ExternalFixConfig> fix_config = std::make_shared<ExternalFixConfig>(m_config, m_model);
                fix_config->configConstraint(con_json, constaints);
                m_constraint_configs.push_back(fix_config);
            } else if (type == "force") {
                std::shared_ptr<ExternalForceConfig> force_config = std::make_shared<ExternalForceConfig>(m_config, m_model);
                force_config->configConstraint(con_json, constaints);
                m_constraint_configs.push_back(force_config);
            } else if (type == "pressure") {
                std::shared_ptr<ExternalPressureConfig> pressure_config = std::make_shared<ExternalPressureConfig>(m_config, m_model);
                pressure_config->configConstraint(con_json, constaints);
                m_constraint_configs.push_back(pressure_config);
            } else {
                Logger::info("unsupport constraint type: type={}", type);
            }
        }
    }
}

EFEM::ExternalFixConfig::ExternalFixConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    m_config = config;
    m_model = model;
    external_constraint_type = ExternalConstraintType::FIX;
}

void EFEM::ExternalFixConfig::configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints)
{
    std::vector<int> fixs;
    if (con_json["constraint_file"].is_string()) {
        std::string file = con_json["constraint_file"].template get<std::string>();
        if (file != "") {
            m_config->io()->readFile(file, [&](int index, const std::string line) {
                std::vector<std::string> s_c_list;
                boost::split(s_c_list, line, boost::is_any_of(", "), boost::token_compress_on);
                for (auto n_string : s_c_list) {
                    fixs.push_back(std::stoi(n_string) - 1);
                }
            });
        }
    }
    std::vector<int> model_ids;
    if (con_json["model_id"].is_array()) {
        model_ids = con_json["model_id"].template get<std::vector<int>>();
    } else if (con_json["model_id"].is_number()) {
        model_ids.push_back(con_json["model_id"].template get<int>());
    }
    if (model_ids.size() == 0) {
        model_ids.push_back(-1);
    }

    for (int i = 0; i < model_ids.size(); i++) {
        std::shared_ptr<Body> body = m_model->getBody(model_ids[i]);
        std::vector<std::shared_ptr<Node>> body_nodes = body->getNodes();
        for (int i = 0; i < fixs.size(); i++) {
            std::static_pointer_cast<ConstraintNodeInfo>(body_nodes[fixs[i]]->extra_info["constraint_info"])->fix = true;
        }
    }
}

EFEM::ExternalForceConfig::ExternalForceConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    m_config = config;
    m_model = model;
    external_constraint_type = ExternalConstraintType::FORCE;
}

void EFEM::ExternalForceConfig::configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints)
{
    Eigen::Vector3d force;
    if (con_json["value"].is_array()) {
        std::vector<double> force_vector = con_json["value"].template get<std::vector<double>>();
        force << force_vector[0], force_vector[1], force_vector[2];
    }
    std::vector<int> force_nodes;
    if (con_json["constraint_file"].is_string()) {
        std::string file = con_json["constraint_file"].template get<std::string>();
        if (file != "") {
            m_config->io()->readFile(file, [&](int index, const std::string line) {
                std::vector<std::string> s_c_list;
                boost::split(s_c_list, line, boost::is_any_of(", "), boost::token_compress_on);
                for (auto n_string : s_c_list) {
                    force_nodes.push_back(std::stoi(n_string) - 1);
                }
            });
        }
    }
    std::vector<int> model_ids;
    if (con_json["model_id"].is_array()) {
        model_ids = con_json["model_id"].template get<std::vector<int>>();
    } else if (con_json["model_id"].is_number()) {
        model_ids.push_back(con_json["model_id"].template get<int>());
    }
    if (model_ids.size() == 0) {
        model_ids.push_back(-1);
    }

    for (int i = 0; i < model_ids.size(); i++) {
        std::shared_ptr<Body> body = m_model->getBody(model_ids[i]);
        std::shared_ptr<ExternalForceConstraint> force_constraint = std::make_shared<ExternalForceConstraint>(m_config, force);
        std::vector<std::shared_ptr<Node>> body_nodes = body->getNodes();
        for (int i = 0; i < force_nodes.size(); i++) {
            force_constraint->addNode(body_nodes[force_nodes[i]]);
        }
        constaints.push_back(force_constraint);
    }
}

EFEM::ExternalForceConstraint::ExternalForceConstraint(std::shared_ptr<Config> config, Eigen::Vector3d force)
{
    m_config = config;
    m_force = force;
    external_constraint_type = ExternalConstraintType::FORCE;
}

void EFEM::ExternalForceConstraint::addNode(std::shared_ptr<Node> node)
{
    m_nodes.push_back(node);
    m_node_infos.push_back(std::static_pointer_cast<ConstraintNodeInfo>(node->extra_info["constraint_info"]));
}

void EFEM::ExternalForceConstraint::conInitial()
{
}

void EFEM::ExternalForceConstraint::conEvaluteResidual()
{
    m_con_residuals.clear();
    for (auto info : m_node_infos) {
        if (info->id >= 0) {
            if (m_con_residuals.find(info->id) == m_con_residuals.end()) {
                m_con_residuals[info->id] = Eigen::Vector3d::Zero();
            }
            m_con_residuals[info->id] += m_force;
        }
    }
}

void EFEM::ExternalForceConstraint::conGetResidual(Eigen::VectorXd& residual)
{
    for (auto pair : m_con_residuals) {
        residual.segment(pair.first * 3, 3) += -pair.second;
    }
}

void EFEM::ExternalForceConstraint::conEvaluteJacobi()
{
}

void EFEM::ExternalForceConstraint::conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets)
{
}

EFEM::ExternalPressureConfig::ExternalPressureConfig(std::shared_ptr<Config> config, std::shared_ptr<Model> model)
{
    m_config = config;
    m_model = model;
    external_constraint_type = ExternalConstraintType::PRESSURE;
}

void EFEM::ExternalPressureConfig::configConstraint(Json& con_json, std::vector<std::shared_ptr<Constrainable>>& constaints)
{
    std::vector<std::vector<int>> pressure_faces;
    double pressure_value = 0;
    Eigen::Vector3d pressure_dvalue = Eigen::Vector3d::Zero();
    if (con_json.find("value") != con_json.end() && con_json["value"].is_number()) {
        pressure_value = con_json.find("value")->template get<double>();
    } else if (con_json.find("value") != con_json.end() && con_json["value"].is_array()) {
        std::vector<double> p_vector = con_json["value"].template get<std::vector<double>>();
        pressure_dvalue << p_vector[0], p_vector[1], p_vector[2];
    }
    if (pressure_value <= 0 && pressure_dvalue.norm() <= 0) {
        return;
    }

    if (con_json["constraint_file"].is_string()) {
        std::string file = con_json["constraint_file"].template get<std::string>();
        if (file != "") {
            m_config->io()->readFile(file, [&](int index, const std::string line) {
                std::vector<std::string> s_c_list;
                boost::split(s_c_list, line, boost::is_any_of(", "), boost::token_compress_on);
                std::vector<int> face;
                for (int i = 0; i < s_c_list.size(); i++) {
                    int face_node_index = std::stoi(s_c_list[i]) - 1;
                    if (std::find(face.begin(), face.end(), face_node_index) == face.end()) {
                        face.push_back(face_node_index);
                    }
                }
                pressure_faces.push_back(face);
            });
        }
    }
    std::vector<int> model_ids;
    if (con_json["model_id"].is_array()) {
        model_ids = con_json["model_id"].template get<std::vector<int>>();
    } else if (con_json["model_id"].is_number()) {
        model_ids.push_back(con_json["model_id"].template get<int>());
    }
    if (model_ids.size() == 0) {
        model_ids.push_back(-1);
    }

    for (int i = 0; i < model_ids.size(); i++) {
        std::shared_ptr<Body> body = m_model->getBody(model_ids[i]);
        std::shared_ptr<ExternalPressureConstraint> pressure_constraint;
        if (pressure_value > 0) {
            pressure_constraint = std::make_shared<ExternalPressureConstraint>(m_config, pressure_value);
        } else {
            pressure_constraint = std::make_shared<ExternalPressureConstraint>(m_config, pressure_dvalue);
        }

        std::vector<std::shared_ptr<Node>> body_nodes = body->getNodes();
        for (int j = 0; j < pressure_faces.size(); j++) {
            std::vector<int> face = pressure_faces[j];
            std::vector<std::shared_ptr<Node>> face_nodes;
            for (int node_index : face) {
                face_nodes.push_back(body_nodes[node_index]);
            }
            std::shared_ptr<ExternalPressureFace> ex_press_face = std::make_shared<ExternalPressureFace>(face_nodes);
            pressure_constraint->addPressureFace(ex_press_face);
        }
        constaints.push_back(pressure_constraint);
    }
}

EFEM::ExternalPressureFace::ExternalPressureFace(std::vector<std::shared_ptr<Node>> face)
{
    m_nodes = face;
    initial();
}

void EFEM::ExternalPressureFace::initial()
{
    std::vector<std::shared_ptr<Element>> elements = m_nodes[0]->elements;
    for (int i = 0; i < m_nodes.size(); i++) {
        m_node_infos.push_back(std::static_pointer_cast<ConstraintNodeInfo>(m_nodes[i]->extra_info["constraint_info"]));
        if (i > 0) {
            for (std::vector<std::shared_ptr<Element>>::iterator it = elements.begin(); it != elements.end();) {
                auto e_find = std::find_if(m_nodes[i]->elements.begin(), m_nodes[i]->elements.end(), [&it](std::shared_ptr<Element> e) {
                    return e->id() == (*it)->id();
                });
                if (e_find == m_nodes[i]->elements.end()) {
                    it = elements.erase(it);
                } else {
                    it++;
                }
            }
        }
    }
    if (elements.size() == 1) {
        m_element = elements[0];
        auto element_info = std::static_pointer_cast<ConstraintElementInfo>(m_element->extra_info["constraint_info"]);
        m_integral_element = element_info->integral_element;
    } else {
        Logger::info("ExternalPressureFace error: multiple elements, size={}", elements.size());
    }

    m_total_nodes = m_element->nodes;
    for (auto t_node : m_total_nodes) {
        m_total_node_infos.push_back(std::static_pointer_cast<ConstraintNodeInfo>(t_node->extra_info["constraint_info"]));
    }

    buildGaussPoint();
}

void EFEM::ExternalPressureFace::buildGaussPoint()
{
    // linear tet
    if (m_nodes.size() == 3) {
        m_gauss_local_coords.resize(1, 3);
        m_gauss_weight.resize(1);
        m_gauss_weight(0) = 1;
        m_A0 = 0.5 * (m_nodes[2]->position - m_nodes[0]->position).cross(m_nodes[1]->position - m_nodes[0]->position);
        if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                return e->id() == this->m_element->nodes[0]->id();
            })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 0;
            m_gauss_local_coords(0, 1) = 1.0 / 3.0;
            m_gauss_local_coords(0, 2) = 1.0 / 3.0;
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[1]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 1.0 / 3.0;
            m_gauss_local_coords(0, 1) = 0;
            m_gauss_local_coords(0, 2) = 1.0 / 3.0;
            if ((this->m_element->nodes[1]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[2]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 1.0 / 3.0;
            m_gauss_local_coords(0, 1) = 1.0 / 3.0;
            m_gauss_local_coords(0, 2) = 0;
            if ((this->m_element->nodes[2]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[3]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 1.0 / 3.0;
            m_gauss_local_coords(0, 1) = 1.0 / 3.0;
            m_gauss_local_coords(0, 2) = 1.0 / 3.0;
            if ((this->m_element->nodes[3]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        }
    } else if (m_nodes.size() == 6) {
        m_gauss_local_coords.resize(3, 3);
        m_gauss_weight.resize(3);
        m_gauss_weight(0) = 1.0 / 3.0;
        m_gauss_weight(1) = 1.0 / 3.0;
        m_gauss_weight(2) = 1.0 / 3.0;
        m_A0 = 0.5 * (m_nodes[2]->position - m_nodes[0]->position).cross(m_nodes[1]->position - m_nodes[0]->position);
        if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                return e->id() == this->m_element->nodes[0]->id();
            })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 2.0 / 3.0;
            m_gauss_local_coords(0, 1) = 1.0 / 6.0;
            m_gauss_local_coords(0, 2) = 1.0 / 6.0;
            m_gauss_local_coords(1, 0) = 1.0 / 6.0;
            m_gauss_local_coords(1, 1) = 2.0 / 3.0;
            m_gauss_local_coords(1, 2) = 1.0 / 6.0;
            m_gauss_local_coords(2, 0) = 1.0 / 6.0;
            m_gauss_local_coords(2, 1) = 1.0 / 6.0;
            m_gauss_local_coords(2, 2) = 2.0 / 3.0;
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[1]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 0;
            m_gauss_local_coords(0, 1) = 2.0 / 3.0;
            m_gauss_local_coords(0, 2) = 1.0 / 6.0;
            m_gauss_local_coords(1, 0) = 0;
            m_gauss_local_coords(1, 1) = 1.0 / 6.0;
            m_gauss_local_coords(1, 2) = 2.0 / 3.0;
            m_gauss_local_coords(2, 0) = 0;
            m_gauss_local_coords(2, 1) = 1.0 / 6.0;
            m_gauss_local_coords(2, 2) = 1.0 / 6.0;
            if ((this->m_element->nodes[1]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[2]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 2.0 / 3.0;
            m_gauss_local_coords(0, 1) = 0;
            m_gauss_local_coords(0, 2) = 1.0 / 6.0;
            m_gauss_local_coords(1, 0) = 1.0 / 6.0;
            m_gauss_local_coords(1, 1) = 0;
            m_gauss_local_coords(1, 2) = 2.0 / 3.0;
            m_gauss_local_coords(2, 0) = 1.0 / 6.0;
            m_gauss_local_coords(2, 1) = 0;
            m_gauss_local_coords(2, 2) = 1.0 / 6.0;
            if ((this->m_element->nodes[2]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        } else if (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[3]->id();
                   })
            == m_nodes.end()) {
            m_gauss_local_coords(0, 0) = 2.0 / 3.0;
            m_gauss_local_coords(0, 1) = 1.0 / 6.0;
            m_gauss_local_coords(0, 2) = 0;
            m_gauss_local_coords(1, 0) = 1.0 / 6.0;
            m_gauss_local_coords(1, 1) = 2.0 / 3.0;
            m_gauss_local_coords(1, 2) = 0;
            m_gauss_local_coords(2, 0) = 1.0 / 6.0;
            m_gauss_local_coords(2, 1) = 1.0 / 6.0;
            m_gauss_local_coords(2, 2) = 0;
            if ((this->m_element->nodes[3]->position - m_nodes[0]->position).dot(m_A0) < 0) {
                m_A0 = -m_A0;
            }
        }
    } else if (m_nodes.size() == 4) {
        m_gauss_local_coords.resize(4, 3);
        m_gauss_weight.resize(4);
        m_gauss_weight(0) = 0.25;
        m_gauss_weight(1) = 0.25;
        m_gauss_weight(2) = 0.25;
        m_gauss_weight(3) = 0.25;
        Eigen::Vector3d A01 = 0.5 * (m_nodes[2]->position - m_nodes[0]->position).cross(m_nodes[1]->position - m_nodes[0]->position);
        Eigen::Vector3d A02 = 0.5 * (m_nodes[0]->position - m_nodes[3]->position).cross(m_nodes[2]->position - m_nodes[3]->position);
        if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                return e->id() == this->m_element->nodes[0]->id();
            }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[2]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 2) = 1.0;
            m_gauss_local_coords(1, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 2) = 1.0;
            m_gauss_local_coords(2, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 2) = 1.0;
            m_gauss_local_coords(3, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 2) = 1.0;
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[0]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[5]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 1) = 1.0;
            m_gauss_local_coords(0, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 1) = 1.0;
            m_gauss_local_coords(1, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 1) = 1.0;
            m_gauss_local_coords(2, 2) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 1) = 1.0;
            m_gauss_local_coords(3, 2) = std::sqrt(1.0 / 3.0);
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[0]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[7]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = 1.0;
            m_gauss_local_coords(0, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 0) = 1.0;
            m_gauss_local_coords(1, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 0) = 1.0;
            m_gauss_local_coords(2, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 2) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 0) = 1.0;
            m_gauss_local_coords(3, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 2) = std::sqrt(1.0 / 3.0);
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[3]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 1) = -1.0;
            m_gauss_local_coords(0, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 1) = -1.0;
            m_gauss_local_coords(1, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 1) = -1.0;
            m_gauss_local_coords(2, 2) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 1) = -1.0;
            m_gauss_local_coords(3, 2) = std::sqrt(1.0 / 3.0);
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[1]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -1.0;
            m_gauss_local_coords(0, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 0) = -1.0;
            m_gauss_local_coords(1, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 2) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 0) = -1.0;
            m_gauss_local_coords(2, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 2) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 0) = -1.0;
            m_gauss_local_coords(3, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 2) = std::sqrt(1.0 / 3.0);
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[4]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(0, 2) = -1.0;
            m_gauss_local_coords(1, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 1) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(1, 2) = -1.0;
            m_gauss_local_coords(2, 0) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(2, 2) = -1.0;
            m_gauss_local_coords(3, 0) = -std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 1) = std::sqrt(1.0 / 3.0);
            m_gauss_local_coords(3, 2) = -1.0;
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        }
        m_A0 = A01 + A02;
    } else if (m_nodes.size() == 8) {
        m_gauss_local_coords.resize(9, 3);
        m_gauss_weight.resize(9);
        m_gauss_weight(0) = 25.0 / 81.0;
        m_gauss_weight(1) = 40.0 / 81.0;
        m_gauss_weight(2) = 25.0 / 81.0;
        m_gauss_weight(3) = 40.0 / 81.0;
        m_gauss_weight(4) = 64.0 / 81.0;
        m_gauss_weight(5) = 40.0 / 81.0;
        m_gauss_weight(6) = 25.0 / 81.0;
        m_gauss_weight(7) = 40.0 / 81.0;
        m_gauss_weight(8) = 25.0 / 81.0;
        m_gauss_weight = 0.25 * m_gauss_weight;
        Eigen::Vector3d A01 = 0.5 * (m_nodes[2]->position - m_nodes[0]->position).cross(m_nodes[1]->position - m_nodes[0]->position);
        Eigen::Vector3d A02 = 0.5 * (m_nodes[0]->position - m_nodes[3]->position).cross(m_nodes[2]->position - m_nodes[3]->position);
        if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                return e->id() == this->m_element->nodes[0]->id();
            }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[2]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -1.0;
            m_gauss_local_coords(0, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 0) = -1.0;
            m_gauss_local_coords(1, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 2) = 0;
            m_gauss_local_coords(2, 0) = -1.0;
            m_gauss_local_coords(2, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 2) = std::sqrt(0.6);
            m_gauss_local_coords(3, 0) = -1.0;
            m_gauss_local_coords(3, 1) = 0;
            m_gauss_local_coords(3, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(4, 0) = -1.0;
            m_gauss_local_coords(4, 1) = 0;
            m_gauss_local_coords(4, 2) = 0;
            m_gauss_local_coords(5, 0) = -1.0;
            m_gauss_local_coords(5, 1) = 0;
            m_gauss_local_coords(5, 2) = std::sqrt(0.6);
            m_gauss_local_coords(6, 0) = -1.0;
            m_gauss_local_coords(6, 1) = std::sqrt(0.6);
            m_gauss_local_coords(6, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(7, 0) = -1.0;
            m_gauss_local_coords(7, 1) = std::sqrt(0.6);
            m_gauss_local_coords(7, 2) = 0;
            m_gauss_local_coords(8, 0) = -1.0;
            m_gauss_local_coords(8, 1) = std::sqrt(0.6);
            m_gauss_local_coords(8, 2) = std::sqrt(0.6);
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[0]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[5]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 2) = -1.0;
            m_gauss_local_coords(1, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 1) = 0;
            m_gauss_local_coords(1, 2) = -1.0;
            m_gauss_local_coords(2, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 1) = std::sqrt(0.6);
            m_gauss_local_coords(2, 2) = -1.0;
            m_gauss_local_coords(3, 0) = 0;
            m_gauss_local_coords(3, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(3, 2) = -1.0;
            m_gauss_local_coords(4, 0) = 0;
            m_gauss_local_coords(4, 1) = 0;
            m_gauss_local_coords(4, 2) = -1.0;
            m_gauss_local_coords(5, 0) = 0;
            m_gauss_local_coords(5, 1) = std::sqrt(0.6);
            m_gauss_local_coords(5, 2) = -1.0;
            m_gauss_local_coords(6, 0) = std::sqrt(0.6);
            m_gauss_local_coords(6, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(6, 2) = -1.0;
            m_gauss_local_coords(7, 0) = std::sqrt(0.6);
            m_gauss_local_coords(7, 1) = 0;
            m_gauss_local_coords(7, 2) = -1.0;
            m_gauss_local_coords(8, 0) = std::sqrt(0.6);
            m_gauss_local_coords(8, 1) = std::sqrt(0.6);
            m_gauss_local_coords(8, 2) = -1.0;
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[0]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[7]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 1) = 1.0;
            m_gauss_local_coords(0, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 1) = 1.0;
            m_gauss_local_coords(1, 2) = 0;
            m_gauss_local_coords(2, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 1) = 1.0;
            m_gauss_local_coords(2, 2) = std::sqrt(0.6);
            m_gauss_local_coords(3, 0) = 0;
            m_gauss_local_coords(3, 1) = 1.0;
            m_gauss_local_coords(3, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(4, 0) = 0;
            m_gauss_local_coords(4, 1) = 1.0;
            m_gauss_local_coords(4, 2) = 0;
            m_gauss_local_coords(5, 0) = 0;
            m_gauss_local_coords(5, 1) = 1.0;
            m_gauss_local_coords(5, 2) = std::sqrt(0.6);
            m_gauss_local_coords(6, 0) = std::sqrt(0.6);
            m_gauss_local_coords(6, 1) = 1.0;
            m_gauss_local_coords(6, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(7, 0) = std::sqrt(0.6);
            m_gauss_local_coords(7, 1) = 1.0;
            m_gauss_local_coords(7, 2) = 0;
            m_gauss_local_coords(8, 0) = std::sqrt(0.6);
            m_gauss_local_coords(8, 1) = 1.0;
            m_gauss_local_coords(8, 2) = std::sqrt(0.6);
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[0]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[3]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 2) = 1.0;
            m_gauss_local_coords(1, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 1) = 0;
            m_gauss_local_coords(1, 2) = 1.0;
            m_gauss_local_coords(2, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 1) = std::sqrt(0.6);
            m_gauss_local_coords(2, 2) = 1.0;
            m_gauss_local_coords(3, 0) = 0;
            m_gauss_local_coords(3, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(3, 2) = 1.0;
            m_gauss_local_coords(4, 0) = 0;
            m_gauss_local_coords(4, 1) = 0;
            m_gauss_local_coords(4, 2) = 1.0;
            m_gauss_local_coords(5, 0) = 0;
            m_gauss_local_coords(5, 1) = std::sqrt(0.6);
            m_gauss_local_coords(5, 2) = 1.0;
            m_gauss_local_coords(6, 0) = std::sqrt(0.6);
            m_gauss_local_coords(6, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(6, 2) = 1.0;
            m_gauss_local_coords(7, 0) = std::sqrt(0.6);
            m_gauss_local_coords(7, 1) = 0;
            m_gauss_local_coords(7, 2) = 1.0;
            m_gauss_local_coords(8, 0) = std::sqrt(0.6);
            m_gauss_local_coords(8, 1) = std::sqrt(0.6);
            m_gauss_local_coords(8, 2) = 1.0;
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[1]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 1) = -1.0;
            m_gauss_local_coords(0, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 1) = -1.0;
            m_gauss_local_coords(1, 2) = 0;
            m_gauss_local_coords(2, 0) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 1) = -1.0;
            m_gauss_local_coords(2, 2) = std::sqrt(0.6);
            m_gauss_local_coords(3, 0) = 0;
            m_gauss_local_coords(3, 1) = -1.0;
            m_gauss_local_coords(3, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(4, 0) = 0;
            m_gauss_local_coords(4, 1) = -1.0;
            m_gauss_local_coords(4, 2) = 0;
            m_gauss_local_coords(5, 0) = 0;
            m_gauss_local_coords(5, 1) = -1.0;
            m_gauss_local_coords(5, 2) = std::sqrt(0.6);
            m_gauss_local_coords(6, 0) = std::sqrt(0.6);
            m_gauss_local_coords(6, 1) = -1.0;
            m_gauss_local_coords(6, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(7, 0) = std::sqrt(0.6);
            m_gauss_local_coords(7, 1) = -1.0;
            m_gauss_local_coords(7, 2) = 0;
            m_gauss_local_coords(8, 0) = std::sqrt(0.6);
            m_gauss_local_coords(8, 1) = -1.0;
            m_gauss_local_coords(8, 2) = std::sqrt(0.6);
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        } else if ((std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                       return e->id() == this->m_element->nodes[6]->id();
                   }) == m_nodes.end())
            && (std::find_if(m_nodes.begin(), m_nodes.end(), [this](std::shared_ptr<Node> e) {
                   return e->id() == this->m_element->nodes[4]->id();
               }) == m_nodes.end())) {
            m_gauss_local_coords(0, 0) = 1.0;
            m_gauss_local_coords(0, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(0, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 0) = 1.0;
            m_gauss_local_coords(1, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(1, 2) = 0;
            m_gauss_local_coords(2, 0) = 1.0;
            m_gauss_local_coords(2, 1) = -std::sqrt(0.6);
            m_gauss_local_coords(2, 2) = std::sqrt(0.6);
            m_gauss_local_coords(3, 0) = 1.0;
            m_gauss_local_coords(3, 1) = 0;
            m_gauss_local_coords(3, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(4, 0) = 1.0;
            m_gauss_local_coords(4, 1) = 0;
            m_gauss_local_coords(4, 2) = 0;
            m_gauss_local_coords(5, 0) = 1.0;
            m_gauss_local_coords(5, 1) = 0;
            m_gauss_local_coords(5, 2) = std::sqrt(0.6);
            m_gauss_local_coords(6, 0) = 1.0;
            m_gauss_local_coords(6, 1) = std::sqrt(0.6);
            m_gauss_local_coords(6, 2) = -std::sqrt(0.6);
            m_gauss_local_coords(7, 0) = 1.0;
            m_gauss_local_coords(7, 1) = std::sqrt(0.6);
            m_gauss_local_coords(7, 2) = 0;
            m_gauss_local_coords(8, 0) = 1.0;
            m_gauss_local_coords(8, 1) = std::sqrt(0.6);
            m_gauss_local_coords(8, 2) = std::sqrt(0.6);
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A01) < 0) {
                A01 = -A01;
            }
            if ((this->m_element->nodes[6]->position - m_nodes[0]->position).dot(A02) < 0) {
                A02 = -A02;
            }
        }
        m_A0 = A01 + A02;
    }

    m_node_weight.resize(m_total_nodes.size(), m_gauss_local_coords.rows());
    for (int i = 0; i < m_gauss_local_coords.rows(); i++) {
        m_node_weight.col(i) = m_integral_element->getPointShapeFunc(m_gauss_local_coords.row(i)).transpose();
        m_gauss_jacobi0.push_back(m_integral_element->getPointJacobi(m_gauss_local_coords.row(i)));
        m_gauss_jacobi0_inv.push_back(m_gauss_jacobi0[i].inverse());

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> gp_gradient = m_integral_element->getPointShapeFuncLocalGradient(m_gauss_local_coords.row(i)).transpose();
        m_gauss_T.push_back(m_gauss_jacobi0_inv[i].transpose() * gp_gradient);
    }
}

void EFEM::ExternalPressureFace::getNodeDArea(int index, Eigen::Vector3d& d_area)
{
    d_area.setZero();
    for (int i = 0; i < m_gauss_local_coords.rows(); i++) {
        Eigen::Matrix3d F = m_integral_element->getPointJacobi(m_gauss_local_coords.row(i)) * m_gauss_jacobi0_inv[i];
        d_area += m_gauss_weight(i) * m_node_weight(index, i) * F.determinant() * F.inverse().transpose() * m_A0;
    }
}

double EFEM::ExternalPressureFace::getNodeArea(int index)
{
    double area = 0;
    for (int i = 0; i < m_gauss_local_coords.rows(); i++) {
        Eigen::Matrix3d F = m_integral_element->getPointJacobi(m_gauss_local_coords.row(i)) * m_gauss_jacobi0_inv[i];
        area += m_gauss_weight(i) * m_node_weight(index, i) * (F.determinant() * F.inverse().transpose() * m_A0).norm();
    }
    return area;
}

void EFEM::ExternalPressureFace::getDAreaJacobi(std::vector<SparseMatrixTriplet> jacobi_triplet, double pressure)
{
    int node_num = m_total_node_infos.size();
    Eigen::MatrixXd gp_H;
    gp_H.resize(node_num * 3, node_num * 3);
    gp_H.setZero();
    int gp_num = m_gauss_local_coords.rows();
    for (int index = 0; index < gp_num; index++) {
        Matrix3333 dPdF;
        Math math;
        EigenMatrix3 F = m_integral_element->getPointJacobi(m_gauss_local_coords.row(index)) * m_gauss_jacobi0_inv[index];
        double J = F.determinant();
        Matrix3333 IF;
        IF.SetIdentity();
        Eigen::Matrix3d FinvT = F.inverse().transpose();
        math.kroneckerProduct(dPdF, FinvT, FinvT);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dPdF(i, j) = dPdF(i, j) - FinvT * IF(i, j).transpose() * FinvT;
            }
        }
        dPdF = J * dPdF;

        double weight = m_gauss_weight[index];
        Eigen::Matrix<double, 3, -1> gp_F_gradient = m_gauss_T[index];
        Matrix3333 gp_dPdF = weight * dPdF;
        for (int k = 0; k < node_num; k++) {
            for (int i = 0; i < 3; i++) {
                Eigen::Matrix<double, 3, -1> gp_H_part;
                gp_H_part.resize(3, node_num);
                gp_H_part.setZero();
                for (int j = 0; j < 3; j++) {
                    gp_H_part += m_node_weight(k, index) * m_A0[j] * gp_dPdF(i, j) * gp_F_gradient;
                }
                for (int j = 0; j < node_num; j++) {
                    gp_H.block<3, 1>(j * 3, k * 3 + i) += gp_H_part.block<3, 1>(0, j);
                }
            }
        }
    }
    gp_H = pressure * gp_H;

    for (unsigned int i = 0; i < m_total_node_infos.size(); i++) {
        for (unsigned int j = 0; j < m_total_node_infos.size(); j++) {
            if (m_total_node_infos[i]->id >= 0 && m_total_node_infos[j]->id >= 0) {
                for (unsigned int row = 0; row < 3; row++) {
                    for (unsigned int col = 0; col < 3; col++) {
                        jacobi_triplet.push_back(SparseMatrixTriplet(m_total_node_infos[i]->id * 3 + row, m_total_node_infos[j]->id * 3 + col, gp_H(i * 3 + row, j * 3 + col)));
                    }
                }
            }
        }
    }
}

EFEM::ExternalPressureConstraint::ExternalPressureConstraint(std::shared_ptr<Config> config, double value)
{
    m_config = config;
    m_direction = false;
    m_pressure_value = value;
    external_constraint_type = ExternalConstraintType::PRESSURE;
}

EFEM::ExternalPressureConstraint::ExternalPressureConstraint(std::shared_ptr<Config> config, Eigen::Vector3d value)
{
    m_config = config;
    m_direction = true;
    m_pressure_dvalue = value;
    external_constraint_type = ExternalConstraintType::PRESSURE;
}

void EFEM::ExternalPressureConstraint::addPressureFace(std::shared_ptr<ExternalPressureFace> face)
{
    m_pressure_faces.push_back(face);
}

void EFEM::ExternalPressureConstraint::conInitial()
{
}

void EFEM::ExternalPressureConstraint::conEvaluteResidual()
{
    m_con_residuals.clear();

    for (auto face : m_pressure_faces) {
        std::vector<std::shared_ptr<ConstraintNodeInfo>> face_node_infos = face->getTotalNodeInfos();
        for (int i = 0; i < face_node_infos.size(); i++) {
            auto info = face_node_infos[i];
            Eigen::Vector3d norm;
            if (m_direction) {
                norm = m_pressure_dvalue;
                norm *= face->getNodeArea(i);
            } else {
                face->getNodeDArea(i, norm);
                norm *= m_pressure_value;
            }
            if (info->id >= 0) {
                if (m_con_residuals.find(info->id) == m_con_residuals.end()) {
                    m_con_residuals[info->id] = Eigen::Vector3d::Zero();
                }
                m_con_residuals[info->id] += norm;
            }
        }
    }
}

void EFEM::ExternalPressureConstraint::conGetResidual(Eigen::VectorXd& residual)
{
    for (auto pair : m_con_residuals) {
        residual.segment(pair.first * 3, 3) += -pair.second;
    }
}

void EFEM::ExternalPressureConstraint::conEvaluteJacobi()
{
}

void EFEM::ExternalPressureConstraint::conGetJacobi(std::vector<SparseMatrixTriplet>& jacobi_triplets)
{
}