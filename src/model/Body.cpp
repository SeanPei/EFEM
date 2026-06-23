#include "model/Body.h"
#include "io/Logger.h"
#include "math/Math.h"
#include <algorithm>

EFEM::Body::Body(const int id, const std::shared_ptr<Config> config, const std::vector<std::pair<std::string, std::vector<double>>>& poses, const std::string& mesh_file)
{
    this->id = id;
    clock_t loadStartTime, loadEndTime;
    loadStartTime = clock();
    Logger::info("Body Start: body load");
    Math math;
    m_nodes.clear();
    m_elements.clear();
    m_render_nodes.clear();
    config->io()->loadMesh(
        mesh_file,
        [&](int index, const Eigen::Vector3d& node_pos) {
            buildNode(index, math.poseTransform(poses, node_pos));
        },
        [&](int index, const std::vector<int>& node_list) {
            buildElement(index, node_list);
        });
    for (int i = 0; i < m_elements.size(); i++) {
        m_elements[i]->init();
    }
    loadEndTime = clock();
    Logger::info("End: Body load, use {} s", (double)(loadEndTime - loadStartTime) / CLOCKS_PER_SEC);
}

double EFEM::Body::volume()
{
    double volume = 0;
    for (int i = 0; i < m_elements.size(); i++) {
        volume = volume + m_elements[i]->volume();
    }
    return volume;
}

void EFEM::Body::buildNode(int index, const Eigen::Vector3d& node_pos)
{
    auto node = std::make_shared<Node>();
    node->index = index;
    node->position = node_pos;
    node->dx.setZero();
    m_nodes.push_back(node);
}

void EFEM::Body::buildElement(int index, const std::vector<int>& raw_node_list)
{
    auto element = std::make_shared<Element>();
    if (raw_node_list.size() <= 5) { // tet_linear
        for (int i = 0; i < 4; i++) {
            int node_index = raw_node_list[i] - 1;
            auto node = m_nodes[node_index];
            element->nodes.push_back(node);
            element->node_indexs.push_back(node->index);
            node->elements.push_back(element);
        }
        element->material_id = raw_node_list.size() > 4 ? raw_node_list[4] : -1;
    } else if (raw_node_list.size() <= 9) { // hex_linear
        for (int i = 0; i < 8; i++) {
            int node_index = raw_node_list[i] - 1;
            auto node = m_nodes[node_index];
            element->nodes.push_back(node);
            element->node_indexs.push_back(node->index);
            node->elements.push_back(element);
        }
        element->material_id = raw_node_list.size() > 8 ? raw_node_list[8] : -1;
    } else if (raw_node_list.size() <= 11) { // tet_quadratic
        for (int i = 0; i < 10; i++) {
            int node_index = raw_node_list[i] - 1;
            auto node = m_nodes[node_index];
            element->nodes.push_back(node);
            element->node_indexs.push_back(node->index);
            node->elements.push_back(element);
        }
        element->material_id = raw_node_list.size() > 10 ? raw_node_list[10] : -1;
    } else if (raw_node_list.size() <= 21) { // hex_quadratic
        for (int i = 0; i < 20; i++) {
            int node_index = raw_node_list[i] - 1;
            auto node = m_nodes[node_index];
            element->nodes.push_back(node);
            element->node_indexs.push_back(node->index);
            node->elements.push_back(element);
        }
        element->material_id = raw_node_list.size() > 20 ? raw_node_list[20] : -1;
    }
    m_elements.push_back(element);
}
